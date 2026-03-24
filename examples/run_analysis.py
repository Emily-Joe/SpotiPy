"""
SpotiPy: run_analysis.py
========================
Full CLV analysis pipeline for a solar active region.

Usage:
python run_analysis.py
(params.txt must be present in the same directory)
"""

import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
import astropy.units as u

from spotipy.downloading import download_series, load_file_list, ask_yn
from spotipy.aligning import align_images
from spotipy.limbdarkening_removal import remove_limb_darkening, get_header_geometry
from spotipy.segmentation import get_masks, get_aia_masks
from spotipy.tracking import track_spots, refine_centering, strip
from spotipy.gui_tools import get_manual_coordinates

# =====================================================================
# CONFIGURATION
# =====================================================================

def load_config(path="params.txt"):
    if not os.path.exists(path):
        print(f"[ERROR] Config file not found: {path}")
        sys.exit(1)
    config = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            k, v = k.strip(), v.strip()
            try: config[k] = int(v)
            except ValueError:
                try: config[k] = float(v)
                except ValueError: config[k] = v
    return config

def update_params_file(path, x_val, y_val):
    with open(path) as f:
        lines = f.readlines()
    with open(path, "w") as f:
        for line in lines:
            if line.startswith("X_START_ARCSEC"):
                f.write(f"X_START_ARCSEC = {x_val:.2f}\n")
            elif line.startswith("Y_START_ARCSEC"):
                f.write(f"Y_START_ARCSEC = {y_val:.2f}\n")
            else:
                f.write(line)

# =====================================================================
# FITS HELPERS
# =====================================================================

def read_fits(path):
    try:
        with fits.open(path) as hdul:
            if len(hdul) > 1 and hdul[1].data is not None:
                return hdul[1].data, hdul[1].header
            return hdul[0].data, hdul[0].header
    except Exception as e:
        print(f" [WARN] Could not read {os.path.basename(path)}: {e}")
        return None, None

def find_closest_file(target_time, file_list):
    best_path, best_dt = None, 9e99
    for path in file_list:
        _, hdr = read_fits(path)
        if hdr is None or "DATE-OBS" not in hdr:
            continue
        dt = abs((Time(hdr["DATE-OBS"]) - target_time).to_value("s"))
        if dt < best_dt:
            best_dt, best_path = dt, path
    return best_path

def build_mu_map(shape, header, x0, y0):
    (solar_cx, solar_cy), r_pix = get_header_geometry(header)[:2]
    if r_pix is None:
        return np.full(shape, np.nan)
    H, W = shape
    yy, xx = np.indices((H, W))
    X = (x0 + xx) - solar_cx
    Y = (y0 + yy) - solar_cy
    with np.errstate(invalid="ignore"):
        mu = np.sqrt(np.clip(1.0 - (np.sqrt(X**2 + Y**2) / r_pix)**2, 0.0, 1.0))
    return mu

def get_heliographic_coords(img, hdr, x0, y0, frame_size):
    """Extract heliographic lon/lat of crop center using SunPy."""
    try:
        import sunpy.map
        hdr_copy = hdr.copy()
        hdr_copy['cunit1'] = 'arcsec'
        hdr_copy['cunit2'] = 'arcsec'
        smap = sunpy.map.Map(img, hdr_copy)
        cx_pix = x0 + frame_size / 2
        cy_pix = y0 + frame_size / 2
        coord = smap.pixel_to_world(cx_pix * u.pix, cy_pix * u.pix)
        tx = coord.Tx.to(u.arcsec).value
        ty = coord.Ty.to(u.arcsec).value
        lon = coord.heliographic_stonyhurst.lon.to(u.deg).value
        lat = coord.heliographic_stonyhurst.lat.to(u.deg).value
        return tx, ty, lon, lat
    except Exception:
        return 0.0, 0.0, 0.0, 0.0

# =====================================================================
# OVERLAY PNG HELPERS
# =====================================================================

def save_hmi_overlay(crop_norm, umbra_m, penum_m, out_path, crop_extent=None):
    disp = np.clip(np.nan_to_num(crop_norm), 0, 2.0) / 2.0
    rgb = np.dstack([disp, disp, disp])
    rgb[penum_m & ~umbra_m] = [1.0, 1.0, 0.0]
    rgb[umbra_m] = [1.0, 0.0, 1.0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    kw = dict(origin='lower', extent=crop_extent) if crop_extent else dict(origin='lower')
    axes[0].imshow(disp, cmap="gray", vmin=0, vmax=1, **kw)
    axes[0].set_title("HMI Raw")
    axes[1].imshow(rgb, **kw)
    axes[1].set_title("HMI Mask")
    if crop_extent:
        for ax in axes:
            ax.set_xlabel("x [arcsec]")
            ax.set_ylabel("y [arcsec]")
    else:
        for ax in axes: ax.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def save_aia_overlay(crop_aia, spot_m, qs_m, plage_m, net_m, out_path, crop_extent=None, ondisk_m=None):
    disp = np.nan_to_num(crop_aia.astype(float))
    valid = disp[disp > 0.1]
    lo, hi = (np.percentile(valid, [1, 99]) if valid.size > 10 else (0, 1))
    if hi <= lo: hi = lo + 1e-6
    base = np.clip((disp - lo) / (hi - lo), 0, 1)
    rgb = np.dstack([base, base, base])
    ondisk = ondisk_m if ondisk_m is not None else np.ones(crop_aia.shape, dtype=bool)
    if qs_m is not None: rgb[qs_m.astype(bool) & ondisk] = [0.0, 1.0, 0.0]
    if net_m is not None: rgb[net_m.astype(bool) & ondisk] = [0.0, 0.0, 1.0]
    if plage_m is not None: rgb[plage_m.astype(bool) & ondisk] = [1.0, 0.0, 0.0]
    if spot_m is not None: rgb[spot_m.astype(bool) & ondisk] = [1.0, 0.0, 1.0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    kw = dict(origin='lower', extent=crop_extent) if crop_extent else dict(origin='lower')
    axes[0].imshow(base, cmap="gray", **kw)
    axes[0].set_title("AIA Raw")
    axes[1].imshow(rgb, **kw)
    axes[1].set_title("AIA Segmentation")
    if crop_extent:
        for ax in axes:
            ax.set_xlabel("x [arcsec]")
            ax.set_ylabel("y [arcsec]")
    else:
        for ax in axes: ax.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

# =====================================================================
# CLV FIT & PLOT HELPERS
# =====================================================================

def fit_clv(mu, intensity, poly_order=3, min_mu=0.15):
    mask = np.isfinite(mu) & np.isfinite(intensity) & (mu >= min_mu)
    if np.sum(mask) < poly_order + 1:
        return None, None
    coeffs_raw = np.polyfit(mu[mask], intensity[mask], poly_order)
    val_at_1 = np.polyval(coeffs_raw, 1.0)
    coeffs_norm = coeffs_raw / val_at_1 if val_at_1 != 0 else coeffs_raw
    return coeffs_raw, coeffs_norm

def _plot_clv(data_dict, colors, title, out_path, poly_order, min_mu):
    fig, ax = plt.subplots(figsize=(7, 5))
    xs = np.linspace(min_mu, 1.0, 300)
    for cat, col in colors.items():
        ml, il = data_dict[cat]
        if not ml: continue
        mu_arr = np.concatenate(ml)
        i_arr = np.concatenate(il)
        ax.scatter(mu_arr, i_arr, s=1, alpha=0.03, color=col, label=cat)
        cr, _ = fit_clv(mu_arr, i_arr, poly_order, min_mu)
        if cr is not None:
            ax.plot(xs, np.polyval(cr, xs), color=col, linewidth=2.5)
    ax.set_xlabel(r"$\mu = \cos\,\theta$")
    ax.set_ylabel("Intensity / Value")
    ax.set_title(title)
    ax.set_xlim(0, 1)
    try:
        leg = ax.legend(loc="best", frameon=True, framealpha=0.9)
        for lh in leg.legend_handles: lh.set_alpha(1.0)
    except AttributeError:
        ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

def save_hmi_clv(data_dict, obs_key, out_path, poly_order, min_mu):
    _plot_clv(data_dict, {"Umbra": "magenta", "Penumbra": "gold", "Full": "cyan"},
              f"HMI CLV — {obs_key}", out_path, poly_order, min_mu)

def save_aia_clv(data_dict, obs_key, out_path, poly_order, min_mu):
    _plot_clv(data_dict, {"Quiet": "gray", "Plage": "red", "Network": "steelblue"},
              f"AIA Masks / HMI {obs_key} CLV", out_path, poly_order, min_mu)

# =====================================================================
# SAVE / LOAD RAW RESULTS
# =====================================================================

def save_results(data_dict, out_dir, poly_order, min_mu):
    os.makedirs(out_dir, exist_ok=True)
    for cat, (mu_list, i_list) in data_dict.items():
        if not mu_list: continue
        mu_all = np.concatenate(mu_list)
        i_all = np.concatenate(i_list)
        np.savetxt(os.path.join(out_dir, f"{cat}_raw.txt"),
                   np.column_stack([mu_all, i_all]), header="mu intensity")
        cr, cn = fit_clv(mu_all, i_all, poly_order, min_mu)
        if cr is not None:
            np.savetxt(os.path.join(out_dir, f"{cat}_fit_coeffs_raw.txt"), cr[None, :])
            np.savetxt(os.path.join(out_dir, f"{cat}_fit_coeffs_norm.txt"), cn[None, :])

def load_results(data_dict, out_dir, categories):
    for cat in categories:
        f = os.path.join(out_dir, f"{cat}_raw.txt")
        if os.path.exists(f):
            d = np.loadtxt(f).reshape(-1, 2)
            data_dict[cat][0].append(d[:, 0])
            data_dict[cat][1].append(d[:, 1])

# =====================================================================
# MAIN
# =====================================================================

def main():
    # --- 0. Config ---
    p = load_config("params.txt")
    noaa        = p["NOAA_NUMBER"]
    email       = p["EMAIL"]
    start_date  = p["START_DATE"]
    days        = int(p["DAYS"])
    cadence_h   = int(p["CADENCE"])
    FRAME_SIZE  = int(p.get("FRAME_SIZE", 400))
    half        = FRAME_SIZE // 2
    MIN_MU      = float(p.get("MIN_MU_FOR_FIT", 0.15))
    HMI_POLY    = int(p.get("HMI_POLY_ORDER", 2))
    AIA_POLY    = int(p.get("AIA_POLY_ORDER", 5))

    print(f"\n{'='*55}")
    print(f" SpotiPy Pipeline | NOAA {noaa}")
    print(f"{'='*55}\n")

    # --- 1. Directories ---
    base_dir = f"data_NOAA_{noaa}"
    ROOT     = f"results_NOAA_{noaa}"
    #hmi_keys = ["Ic", "M", "V", "Ld", "Lw"]
    hmi_keys = ["Ic"]


    series_to_download = {
        "HMI_NoLD": ("hmi.Ic_noLimbDark_720s", None),
        "Ic":       ("hmi.Ic_720s",             None),
        "M":        ("hmi.M_720s",              None),
        "V":        ("hmi.V_720s",              None),
        "Ld":       ("hmi.Ld_720s",             None),
        "Lw":       ("hmi.Lw_720s",             None),
        "AIA":      ("aia.lev1_uv_24s",         1700),
    }

    data_dirs = {name: os.path.join(base_dir, name) for name in series_to_download}
    dirs = {
        "aia_aligned": os.path.join(ROOT, "AIA_Aligned"),
        "aia_nolbd":   os.path.join(ROOT, "AIA_NoLimbDark"),
        "mask_hmi":    os.path.join(ROOT, "Masks_HMI"),
        "mask_aia":    os.path.join(ROOT, "Masks_AIA"),
        "clv_hmi":     os.path.join(ROOT, "CLV_Plots_HMI"),
        "clv_aia":     os.path.join(ROOT, "CLV_Plots_AIA"),
    }

    hmi_res_dirs = {key: os.path.join(ROOT, "Results_HMI", key) for key in hmi_keys}
    aia_res_dirs = {key: os.path.join(ROOT, "Results_AIA", key) for key in hmi_keys}

    for d in list(dirs.values()) + list(hmi_res_dirs.values()) + list(aia_res_dirs.values()):
        os.makedirs(d, exist_ok=True)

    # --- 2. Download / load ---
    file_map = {}
    if ask_yn("Download data?"):
        for name, (series, wl) in series_to_download.items():
            file_map[name] = download_series(
                series=series, start_date=start_date, days=days,
                cadence_h=cadence_h, email=email,
                out_dir=data_dirs[name], wavelength_angstrom=wl,
            )
    else:
        print("Loading existing files from saved lists...")
        for name, (series, _) in series_to_download.items():
            safe = series.replace(".", "_")
            lpath = os.path.join(data_dirs[name], f"{safe}_files.txt")
            file_map[name] = load_file_list(lpath, data_dirs[name])
            print(f"  {name}: {len(file_map[name])} files")

    if not file_map.get("HMI_NoLD"):
        print("[ERROR] No HMI NoLimbDark data found. Exiting.")
        sys.exit(1)

    # --- 3. GUI coordinate selection (MUST happen before switching backend) ---
    if ask_yn("\nManually select starting coordinates?"):
        result = get_manual_coordinates(file_map["HMI_NoLD"][0])
        if result:
            x_new, y_new = result
            update_params_file("params.txt", x_new, y_new)
            p["X_START_ARCSEC"] = x_new
            p["Y_START_ARCSEC"] = y_new
            print(f"  Coordinates updated: x={x_new:.2f}'', y={y_new:.2f}''")

    x_start = float(p.get("X_START_ARCSEC", 0))
    y_start = float(p.get("Y_START_ARCSEC", 0))

    # Switch to non-interactive Agg backend — GUI is done
    matplotlib.use("Agg")
    plt.switch_backend("Agg")

    # --- 4. Differential rotation tracking ---
    print("\nCalculating differential rotation tracks...")
    nold_files = [os.path.basename(f) for f in file_map["HMI_NoLD"]]
    tracks = track_spots(nold_files, data_dirs["HMI_NoLD"], x_start, y_start)
    print(f"  Tracks computed for {len(tracks)} timesteps.")

    # --- Strip visualization ---
    print("Generating time-summed strip...")
    strip(
        series           = nold_files,
        directory        = data_dirs["HMI_NoLD"],
        tracks           = tracks,
        strip_height_arcsec = 400,
        frame_size       = FRAME_SIZE,
        overlay          = True,
        animate          = False,
        save_path        = ROOT,
    )

    # --- 5. Skip-if-exists check ---
    skip_fits = False
    example = os.path.join(hmi_res_dirs["Ic"], "Umbra_raw.txt")
    if os.path.exists(example):
        if not ask_yn("\nFound existing results. Recalculate from FITS?"):
            skip_fits = True
            print("  Skipping FITS loop — loading existing results.")

    # --- 6. Data storage ---
    hmi_cats = ["Umbra", "Penumbra", "Full"]
    aia_cats = ["Quiet", "Plage", "Network"]
    DATA_HMI = {key: {cat: [[], []] for cat in hmi_cats} for key in hmi_keys}
    DATA_AIA = {key: {cat: [[], []] for cat in aia_cats} for key in hmi_keys}

    # Containers for .npz output
    all_hmi_context    = []
    all_hmi_seg        = []
    all_aia_context    = []
    all_aia_seg        = []
    all_times          = []
    all_x_arcsec       = []
    all_y_arcsec       = []
    all_lon            = []
    all_lat            = []

    # --- 7. Main loop ---
    if not skip_fits:
        save_masks = ask_yn("\nSave mask overlay PNGs?")
        print("\nStarting main processing loop...\n")
        n_steps = len(file_map["HMI_NoLD"])

        for i, nold_path in enumerate(file_map["HMI_NoLD"]):
            print(f"  [{i+1}/{n_steps}] {os.path.basename(nold_path)}")

            d_nold, hdr_nold = read_fits(nold_path)
            if d_nold is None:
                print("  [SKIP] Could not read HMI NoLD."); continue

            t_ref = Time(hdr_nold["DATE-OBS"])
            (solar_cx, solar_cy), r_pix = get_header_geometry(hdr_nold)
            if r_pix is None:
                print("  [SKIP] Missing solar geometry."); continue

            base_name = os.path.basename(nold_path).replace(".fits", "")

            # Align AIA
            best_aia = find_closest_file(t_ref, file_map["AIA"])
            if best_aia is None:
                print("  [SKIP] No matching AIA file."); continue
            aligned_aia_path = os.path.join(dirs["aia_aligned"], f"{base_name}_AIA_aligned.fits")
            if not align_images(best_aia, nold_path, aligned_aia_path):
                print("  [SKIP] AIA alignment failed."); continue
            d_aia, _ = read_fits(aligned_aia_path)

            # AIA limb darkening removal
            flat_aia_path = os.path.join(dirs["aia_nolbd"], f"{base_name}_AIA_nolbd.fits")
            if not os.path.exists(flat_aia_path):
                flat_aia = remove_limb_darkening(d_aia, center=(solar_cx, solar_cy), radius_pix=r_pix)
                fits.writeto(flat_aia_path, flat_aia, hdr_nold, overwrite=True)
            else:
                flat_aia, _ = read_fits(flat_aia_path)

            # Crop
            x_arc, y_arc = tracks[i]
            cx = int(round(solar_cx - (x_arc / hdr_nold["CDELT1"])))
            cy = int(round(solar_cy - (y_arc / hdr_nold["CDELT2"])))
            H, W = d_nold.shape
            rx0 = int(np.clip(cx - half, 0, W - FRAME_SIZE))
            ry0 = int(np.clip(cy - half, 0, H - FRAME_SIZE))
            cx, cy = refine_centering(d_nold[ry0:ry0+FRAME_SIZE, rx0:rx0+FRAME_SIZE], rx0, ry0, FRAME_SIZE)
            x0 = int(np.clip(cx, 0, W - FRAME_SIZE))
            y0 = int(np.clip(cy, 0, H - FRAME_SIZE))

            # Arcsec extent for overlay axes
            _dx = hdr_nold["CDELT1"]
            _dy = hdr_nold["CDELT2"]
            x_left  = (x0 - solar_cx) * _dx
            x_right = (x0 + FRAME_SIZE - solar_cx) * _dx
            y_bot   = (y0 - solar_cy) * _dy
            y_top   = (y0 + FRAME_SIZE - solar_cy) * _dy
            crop_extent = [x_right, x_left, y_top, y_bot]

            crop_nold = np.rot90(d_nold[y0:y0+FRAME_SIZE, x0:x0+FRAME_SIZE], 2)
            crop_aia  = np.rot90(flat_aia[y0:y0+FRAME_SIZE, x0:x0+FRAME_SIZE], 2)
            nold_med  = np.nanmedian(crop_nold)
            crop_norm = crop_nold / nold_med if nold_med > 0 else crop_nold

            # mu map
            mu_map  = np.rot90(build_mu_map((FRAME_SIZE, FRAME_SIZE), hdr_nold, x0, y0), 2)
            ondisk_m = mu_map > 0

            # HMI segmentation
            hmi_masks = get_masks(
                crop_norm,
                disk_mask      = ondisk_m,
                umbra_range    = (p.get("UMBRA_MIN", 10),    p.get("UMBRA_MAX", 55)),
                penumbra_range = (p.get("PENUMBRA_MIN", 75), p.get("PENUMBRA_MAX", 120)),
            )
            umbra_m = hmi_masks["umbra"]
            penum_m = hmi_masks["penumbra"]
            full_m  = hmi_masks["spot"]
            all_dark = hmi_masks['all_dark']

            # Save HMI FITS + overlay PNG
            fits.writeto(os.path.join(dirs["mask_hmi"], f"{base_name}_spotmask.fits"),
                         full_m.astype(np.uint8), hdr_nold, overwrite=True)
            if save_masks:
                save_hmi_overlay(crop_norm, umbra_m, penum_m,
                                 os.path.join(dirs["mask_hmi"], f"{base_name}_overlay.png"),
                                 crop_extent=crop_extent)

            # AIA segmentation
            aia_masks = get_aia_masks(
                crop_aia,
                spot_mask=full_m,
                    all_dark_mask=all_dark,
                plage_excess_pct = float(p.get("PLAGE_EXCESS_PCT", 20.0)),
                qs_tol_pct       = float(p.get("QUIET_SUN_TOL_PCT", 15.0)),
                min_area         = int(p.get("MIN_PLAGE_AREA", 450)),
            )
            qs_m    = aia_masks["qs"]
            plage_m = aia_masks["plage"]
            net_m   = aia_masks["network"]

            # Save AIA FITS + overlay PNG
            seg_map = np.zeros_like(crop_aia, dtype=np.uint8)
            seg_map[full_m]  = 1
            seg_map[plage_m] = 2
            seg_map[net_m]   = 3
            seg_map[qs_m]    = 4
            fits.writeto(
                os.path.join(dirs["mask_aia"],
                             os.path.basename(best_aia).replace(".fits", "_segmentation.fits")),
                seg_map, hdr_nold, overwrite=True,
            )
            if save_masks:
                save_aia_overlay(crop_aia, full_m, qs_m, plage_m, net_m,
                                 os.path.join(dirs["mask_aia"],
                                              os.path.basename(best_aia).replace(".fits", "_overlay.png")),
                                 crop_extent=crop_extent, ondisk_m=ondisk_m)

            # Collect frames for .npz
            tx, ty, lon, lat = get_heliographic_coords(d_nold, hdr_nold, x0, y0, FRAME_SIZE)
            all_hmi_context.append(crop_norm.astype(np.float32))
            all_hmi_seg.append(seg_map)  # reuse combined HMI+AIA seg
            all_aia_context.append(crop_aia.astype(np.float32))
            all_aia_seg.append(seg_map)
            all_times.append(t_ref.iso)
            all_x_arcsec.append(tx)
            all_y_arcsec.append(ty)
            all_lon.append(lon)
            all_lat.append(lat)

            # Extract for all 5 HMI observables
            for obs_key in hmi_keys:
                obs_path = find_closest_file(t_ref, file_map[obs_key])
                if obs_path is None: continue
                d_obs, _ = read_fits(obs_path)
                if d_obs is None: continue
                crop_obs = np.rot90(d_obs[y0:y0+FRAME_SIZE, x0:x0+FRAME_SIZE], 2)

                for cat, mask in [("Umbra", umbra_m), ("Penumbra", penum_m), ("Full", full_m)]:
                    if mask is None or not np.any(mask): continue
                    mu_v = mu_map[mask]; i_v = crop_obs[mask]
                    ok = np.isfinite(mu_v) & np.isfinite(i_v)
                    DATA_HMI[obs_key][cat][0].append(mu_v[ok])
                    DATA_HMI[obs_key][cat][1].append(i_v[ok])

                for cat, mask in [("Quiet", qs_m), ("Plage", plage_m), ("Network", net_m)]:
                    if mask is None or not np.any(mask): continue
                    mu_v = mu_map[mask]; i_v = crop_obs[mask]
                    ok = np.isfinite(mu_v) & np.isfinite(i_v)
                    DATA_AIA[obs_key][cat][0].append(mu_v[ok])
                    DATA_AIA[obs_key][cat][1].append(i_v[ok])

        # --- 8. Save results ---
        print("\nSaving results...")
        for obs_key in hmi_keys:
            save_results(DATA_HMI[obs_key], hmi_res_dirs[obs_key], HMI_POLY, MIN_MU)
            save_results(DATA_AIA[obs_key], aia_res_dirs[obs_key], AIA_POLY, MIN_MU)
        print("  Done.")

        # --- Save .npz ---
        if all_times:
            npz_path = os.path.join(ROOT, f"spotipy_NOAA_{noaa}.npz")
            print(f"\nSaving {npz_path}...")
            np.savez_compressed(
                npz_path,
                hmi_context    = np.array(all_hmi_context),
                hmi_seg        = np.array(all_hmi_seg),
                aia_context    = np.array(all_aia_context),
                aia_seg        = np.array(all_aia_seg),
                times          = np.array(all_times),
                x_arcsec       = np.array(all_x_arcsec),
                y_arcsec       = np.array(all_y_arcsec),
                longitude      = np.array(all_lon),
                latitude       = np.array(all_lat),
            )
            print(f"  Saved {len(all_times)} frames → {npz_path}")

    else:
        # --- Load existing results ---
        print("\nLoading existing results...")
        for obs_key in hmi_keys:
            load_results(DATA_HMI[obs_key], hmi_res_dirs[obs_key], hmi_cats)
            load_results(DATA_AIA[obs_key], aia_res_dirs[obs_key], aia_cats)

    # --- 9. Generate CLV plots ---
    print("\nGenerating CLV plots...")
    for obs_key in hmi_keys:
        save_hmi_clv(DATA_HMI[obs_key], obs_key,
                     os.path.join(dirs["clv_hmi"], f"HMI_CLV_{obs_key}.png"), HMI_POLY, MIN_MU)
        save_aia_clv(DATA_AIA[obs_key], obs_key,
                     os.path.join(dirs["clv_aia"], f"AIA_CLV_{obs_key}.png"), AIA_POLY, MIN_MU)
        print(f"  {obs_key}: HMI + AIA CLV plots saved.")

    print(f"\n{'='*55}")
    print(f" Pipeline complete! Results saved to: {ROOT}/")
    print(f"{'='*55}\n")

if __name__ == "__main__":
    main()
