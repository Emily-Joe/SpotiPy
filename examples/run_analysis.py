"""
SpotiPy Pipeline: Full Implementation
-----------------------------------
1. Downloads HMI and AIA data
2. Prompts for manual GUI clicker to update coordinates.
3. Aligns AIA to HMI and saves intermediate FITS.
4. Removes AIA Limb Darkening and saves intermediate FITS.
5. Generates Masks (HMI: Umbra/Penumbra | AIA: Quiet/Plage/Network).
6. Extracts mu & intensity arrays, saves to .txt, and plots scatter+CLV fits.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time, TimeDelta

# --- Import SpotiPy Modules ---
from spotipy.downloading import download_data, ask_yn
from spotipy import (
    get_masks,
    get_aia_masks,
    remove_limb_darkening,
    track_spots,
    refine_centering,
    align_images,
    get_manual_coordinates
)

def load_config(path):
    config = {}
    if not os.path.exists(path):
        print(f"Error: Config file '{path}' not found.")
        sys.exit(1)
    with open(path, 'r') as f:
        for line in f:
            if line.strip().startswith('#') or '=' not in line: continue
            key, val = line.split('=', 1)
            try:
                config[key.strip()] = float(val.strip()) if '.' in val else int(val.strip())
            except ValueError:
                config[key.strip()] = val.strip()
    return config

def update_params_file(path, x_val, y_val):
    """Rewrites params.txt with the newly clicked coordinates."""
    with open(path, 'r') as f:
        lines = f.readlines()
    with open(path, 'w') as f:
        for line in lines:
            if line.startswith('X_START_ARCSEC'):
                f.write(f"X_START_ARCSEC={x_val:.2f}\n")
            elif line.startswith('Y_START_ARCSEC'):
                f.write(f"Y_START_ARCSEC={y_val:.2f}\n")
            else:
                f.write(line)

def calculate_mu_grid(shape, center, radius):
    y, x = np.indices(shape)
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    with np.errstate(invalid='ignore'):
        mu = np.sqrt(1.0 - (r / radius)**2)
    mu[np.isnan(mu)] = 0.0
    return mu

def collect_I_mu(img, mu, mask):
    v = img[mask]
    m = mu[mask]
    valid = np.isfinite(v) & np.isfinite(m)
    return m[valid], v[valid]

def plot_clv_with_scatter(mu_arrays, I_arrays, labels, colors, title, outpng):
    fig, ax = plt.subplots(figsize=(7,5))
    xs = np.linspace(0.15, 1.0, 100)

    for mu, I, lab, col in zip(mu_arrays, I_arrays, labels, colors):
        if mu is not None and mu.size > 0:
            ax.scatter(mu, I, s=1, alpha=0.03, c=col, label=f"{lab} pts")

            mask = (mu >= 0.15)
            if np.sum(mask) > 3:
                coeffs = np.polyfit(mu[mask], I[mask], 2)
                ys = np.polyval(coeffs, xs)
                ax.plot(xs, ys, color=col, linewidth=3, linestyle='-', label=f"{lab} fit")

    ax.set_xlabel(r"$\mu = \cos\,\theta$")
    ax.set_ylabel("Intensity / Value")
    ax.set_title(title)
    ax.set_xlim(0, 1)

    leg = ax.legend(loc='best', frameon=True)
    for lh in leg.legendHandles:
        lh.set_alpha(1.0)

    fig.tight_layout()
    fig.savefig(outpng, dpi=200)
    plt.close(fig)

def main():
    p = load_config("params.txt")
    print(f"--- Pipeline: NOAA {p['NOAA_NUMBER']} ---")

    observables = {
        'Ic': 'hmi.Ic_45s', 'M': 'hmi.M_45s', 'V': 'hmi.V_45s',
        'Ld': 'hmi.Ld_45s', 'Lw': 'hmi.Lw_45s', 'AIA': 'aia.lev1_uv_24s'
    }

    ROOT = f"results_NOAA_{p['NOAA_NUMBER']}"
    dirs = {
        "aia_aligned": os.path.join(ROOT, "AIA_Aligned"),
        "aia_nolbd": os.path.join(ROOT, "AIA_NoLimbDark"),
        "clv_plots": os.path.join(ROOT, "CLV_Plots"),
        "text_data": os.path.join(ROOT, "Text_Data")
    }
    for d in dirs.values(): os.makedirs(d, exist_ok=True)

    base_data_dir = f"data_NOAA_{p['NOAA_NUMBER']}"
    file_map = {}

    if ask_yn("\nDo you want to search and download new data?"):
        file_map = download_data(
            observables=observables,
            start_date=p['START_DATE'],
            days=p['DAYS'],
            cadence=p['CADENCE'],
            email=p['EMAIL'],
            base_dir=base_data_dir
        )
    else:
        print("\nSkipping download. Loading existing files...")
        for obs_name in observables.keys():
            dir_path = os.path.join(base_data_dir, obs_name)
            if os.path.exists(dir_path):
                files = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith('.fits')]
                files.sort()
                file_map[obs_name] = files
                print(f"Found {len(files)} files for {obs_name}.")
            else:
                print(f"Warning: Directory {dir_path} not found.")
                file_map[obs_name] = []

    # Check if we have required files before proceeding
    if not file_map.get('Ic'):
        print("\nError: No continuum (Ic) data found to start tracking. Exiting.")
        sys.exit(1)

    # --- GUI CLICKER INTEGRATION ---
    if ask_yn("\nDo you want to manually select the sunspot coordinates?"):
        first_ic = file_map['Ic'][0]
        x_new, y_new = get_manual_coordinates(first_ic)
        print(f"Selected Coordinates: X={x_new:.2f}, Y={y_new:.2f}")
        update_params_file("params.txt", x_new, y_new)
        p['X_START_ARCSEC'] = x_new
        p['Y_START_ARCSEC'] = y_new

    print(f"\n--- Calculating Tracks ---")
    ref_files = file_map['Ic']
    t0 = Time(p['START_DATE'])
    dt = TimeDelta(float(p['CADENCE']) * 3600, format='sec')
    times = [t0 + i*dt for i in range(len(ref_files))]

    coords = track_spots(p['START_DATE'],
                         (p.get('X_START_ARCSEC',0), p.get('Y_START_ARCSEC',0)),
                         p['DAYS'], times)

    DATA = {key: {'Umbra': [[],[]], 'Penumbra': [[],[]], 'Quiet': [[],[]], 'Plage': [[],[]], 'Network': [[],[]]}
            for key in observables if key != 'AIA'}

    sz = p.get('FRAME_SIZE', 400) // 2

    print(f"\n--- Processing {len(ref_files)} Timesteps ---")
    for i in range(len(ref_files)):
        missing = any(i >= len(file_map[key]) for key in observables)
        if missing:
            print(f"Step {i+1}: Missing data. Skipping.")
            continue

        base_name = os.path.basename(file_map['Ic'][i]).replace('.fits', '')
        print(f"Step {i+1}: {base_name}...", end="\r")

        with fits.open(file_map['Ic'][i]) as h:
            d_Ic = h[1].data if len(h)>1 else h[0].data
            header_Ic = h[1].header if len(h)>1 else h[0].header

        with fits.open(file_map['AIA'][i]) as h:
            d_AIA = h[1].data if len(h)>1 else h[0].data
            header_AIA = h[1].header if len(h)>1 else h[0].header

        aligned_aia_path = os.path.join(dirs["aia_aligned"], f"{base_name}_AIA_aligned.fits")
        if not os.path.exists(aligned_aia_path):
            aligned_aia = align_images(d_AIA, header_AIA, header_Ic)
            fits.writeto(aligned_aia_path, aligned_aia, header_Ic, overwrite=True)
        else:
            aligned_aia = fits.getdata(aligned_aia_path)

        flat_aia_path = os.path.join(dirs["aia_nolbd"], f"{base_name}_AIA_nolbd.fits")
        if not os.path.exists(flat_aia_path):
            flat_aia = remove_limb_darkening(aligned_aia, center=(2048,2048), radius_pix=1900)
            fits.writeto(flat_aia_path, flat_aia, header_Ic, overwrite=True)
        else:
            flat_aia = fits.getdata(flat_aia_path)

        tx, ty = coords[i] if i < len(coords) else (0,0)
        rough_cx = int(2048 + (tx / 0.6))
        rough_cy = int(2048 + (ty / 0.6))

        rough_crop = d_Ic[max(0, rough_cy-sz):min(4096, rough_cy+sz),
                          max(0, rough_cx-sz):min(4096, rough_cx+sz)]
        cx, cy = refine_centering(rough_crop, rough_cx, rough_cy, p.get('FRAME_SIZE', 400))

        def get_crop(data):
            return data[max(0, cy-sz):min(4096, cy+sz), max(0, cx-sz):min(4096, cx+sz)]

        crop_Ic = get_crop(d_Ic)
        crop_AIA = get_crop(flat_aia)

        hmi_masks = get_masks(crop_Ic,
                              umbra_range=(p['UMBRA_MIN'], p['UMBRA_MAX']),
                              penumbra_range=(p['PENUMBRA_MIN'], p['PENUMBRA_MAX']))
        aia_masks = get_aia_masks(crop_AIA)

        mu_map = calculate_mu_grid(crop_Ic.shape, center=(sz - (cx-2048), sz - (cy-2048)), radius=1900)

        for obs_key in ['Ic', 'M', 'V', 'Ld', 'Lw']:
            with fits.open(file_map[obs_key][i]) as h:
                d_obs = h[1].data if len(h)>1 else h[0].data

            crop_obs = get_crop(d_obs)
            if obs_key in ['M', 'V']:
                crop_obs = np.abs(crop_obs)

            for cat, mask in [('Umbra', hmi_masks['umbra']), ('Penumbra', hmi_masks['penumbra'])]:
                m, v = collect_I_mu(crop_obs, mu_map, mask)
                DATA[obs_key][cat][0].append(m)
                DATA[obs_key][cat][1].append(v)

            aia_regions = [('Quiet', aia_masks.get('quiet')),
                           ('Plage', aia_masks.get('plage')),
                           ('Network', aia_masks.get('network'))]

            for cat, mask in aia_regions:
                if mask is not None:
                    m, v = collect_I_mu(crop_obs, mu_map, mask)
                    DATA[obs_key][cat][0].append(m)
                    DATA[obs_key][cat][1].append(v)

    print("\nGenerating scatter plots and text files...")
    for obs_key in ['Ic', 'M', 'V', 'Ld', 'Lw']:
        obs_dir = os.path.join(dirs["text_data"], obs_key)
        os.makedirs(obs_dir, exist_ok=True)

        plot_mu, plot_I, plot_lbl, plot_col = [], [], [], []
        regions = [('Umbra', 'magenta'), ('Penumbra', 'yellow'), ('Quiet', 'gray'), ('Plage', 'orange'), ('Network', 'red')]

        for cat, color in regions:
            m_list = DATA[obs_key][cat][0]
            i_list = DATA[obs_key][cat][1]

            if m_list:
                mu_all = np.concatenate(m_list)
                I_all = np.concatenate(i_list)

                np.savetxt(os.path.join(obs_dir, f"lbd_{cat}_raw.txt"), np.column_stack([mu_all, I_all]))

                plot_mu.append(mu_all)
                plot_I.append(I_all)
                plot_lbl.append(cat)
                plot_col.append(color)

        plot_clv_with_scatter(
            plot_mu, plot_I, plot_lbl, plot_col,
            title=f"CLV ({obs_key} - Abs)",
            outpng=os.path.join(dirs["clv_plots"], f"CLV_{obs_key}.png")
        )

    print(f"\n Analysis Complete. Files saved in {ROOT}/")

if __name__ == "__main__":
    main()
