"""
SpotiPy Multi-Observable Pipeline (Full 5-Parameter Version)
------------------------------------------------------------
Workflow:
1. Download Intensity (Ic), Magnetogram (M), Doppler (V), Line Depth (Ld), Line Width (Lw).
2. Segment using Intensity (Create masks for Umbra/Penumbra).
3. Apply these masks to ALL 5 observables.
4. Generate CLV Candle Plots for all 5 physics.
"""

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time, TimeDelta

# --- Import SpotiPy Modules ---
import spotipy
from spotipy.segmentation import get_masks
from spotipy.processing import remove_limb_darkening
from spotipy.tracking import track_spots

def load_config(path):
    """Parses params.txt."""
    config = {}
    if not os.path.exists(path):
        print(f"❌ Error: Config file '{path}' not found.")
        sys.exit(1)
    with open(path, 'r') as f:
        for line in f:
            if line.strip().startswith('#') or '=' not in line: continue
            key, val = line.split('=', 1)
            key, val = key.strip(), val.strip()
            try:
                if '.' in val: config[key] = float(val)
                else: config[key] = int(val)
            except ValueError:
                config[key] = val
    return config

def calculate_mu_grid(shape, center, radius):
    """Generates a grid of mu (cosine theta) values."""
    y, x = np.indices(shape)
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    with np.errstate(invalid='ignore'):
        mu = np.sqrt(1.0 - (r / radius)**2)
    mu[np.isnan(mu)] = 0.0
    return mu

def save_candle_plot(mu, data, mask, label, unit, color, filename):
    """Generates a CLV candle plot for a specific physical quantity."""
    # Filter: valid pixels inside the mask
    # We use np.abs(data) for things like Magnetograms to show magnitude
    valid = mask.ravel() & (mu.ravel() > 0.1) & np.isfinite(data.ravel())

    if np.sum(valid) < 50: return # Skip if not enough data

    x = mu.ravel()[valid]
    y = data.ravel()[valid]

    # Binning (0.1 to 1.0)
    bins = np.linspace(0.1, 1.0, 10)
    centers = 0.5 * (bins[1:] + bins[:-1])
    binned_data = [y[(x >= b) & (x < bins[i+1])] for i, b in enumerate(bins[:-1])]

    # Plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.boxplot(binned_data, positions=centers, widths=0.05,
               patch_artist=True,
               boxprops=dict(facecolor=color, alpha=0.5),
               showfliers=False, medianprops=dict(color='black'))

    ax.set_title(f"CLV: {label}")
    ax.set_xlabel(r"$\mu = \cos(\theta)$")
    ax.set_ylabel(f"{label} [{unit}]")
    ax.set_xlim(0, 1.1)
    ax.grid(True, alpha=0.3)

    plt.savefig(filename, dpi=100)
    plt.close(fig)

def main():
    p = load_config("params.txt")
    print(f"--- 🌞 5-Observable Pipeline: NOAA {p['NOAA_NUMBER']} ---")

    # 1. Define ALL 5 Observables (Name: JSOC Series)
    observables = {
        'Ic': 'hmi.Ic_45s',   # Continuum Intensity
        'M':  'hmi.M_45s',    # Magnetogram
        'V':  'hmi.V_45s',    # Dopplergram
        'Ld': 'hmi.Ld_45s',   # Line Depth
        'Lw': 'hmi.Lw_45s'    # Line Width
    }

    file_map = {} # Stores lists of files: {'Ic': [...], 'M': [...], ...}

    # 2. Download Data loop
    try:
        for obs_name, series in observables.items():
            print(f"\n⬇️  Downloading {obs_name} ({series})...")
            files = spotipy.download_data(
                start_time=p['START_DATE'],
                duration_days=p['DAYS'],
                instrument=series,
                cadence=f"{p['CADENCE']}h",
                email=p['EMAIL'],
                output_dir=f"data_NOAA_{p['NOAA_NUMBER']}/{obs_name}"
            )
            files.sort()
            file_map[obs_name] = files

    except Exception as e:
        print(f"❌ Download Error: {e}")
        sys.exit(1)

    # 3. Track Spot
    print(f"\n--- Calculating Tracks ---")
    # Generate time steps based on the first observable (Ic)
    ref_files = file_map['Ic']
    t0 = Time(p['START_DATE'])
    dt = TimeDelta(float(p['CADENCE']) * 3600, format='sec')
    times = [t0 + i*dt for i in range(len(ref_files))]

    coords = track_spots(p['START_DATE'],
                         (p.get('X_START_ARCSEC',0), p.get('Y_START_ARCSEC',0)),
                         p['DAYS'], times)

    # 4. Process Loop
    res_dir = f"results_NOAA_{p['NOAA_NUMBER']}"
    if not os.path.exists(res_dir): os.makedirs(res_dir)

    print(f"\n--- Processing {len(ref_files)} Timesteps ---")

    for i in range(len(ref_files)):
        # Ensure we have files for all 5 observables at this step
        current_files = {}
        missing = False
        for key in observables:
            if i < len(file_map[key]):
                current_files[key] = file_map[key][i]
            else:
                missing = True

        if missing:
            print(f"   Step {i+1}: Missing data for some observables. Skipping.")
            continue

        base_name = os.path.basename(current_files['Ic']).replace('.fits', '')
        print(f"   Step {i+1}: {base_name}...", end="\r")

        # --- A. Load & Align Intensity (Reference) ---
        with fits.open(current_files['Ic']) as h:
            d_Ic = h[1].data if len(h)>1 else h[0].data

        # Calculate Crop Center from Tracking
        tx, ty = coords[i] if i < len(coords) else (0,0)
        cx = int(2048 + (tx / 0.6))
        cy = int(2048 + (ty / 0.6))
        sz = p.get('FRAME_SIZE', 400) // 2

        # Helper to crop
        def get_crop(data):
            # Bounds check
            y1, y2 = max(0, cy-sz), min(4096, cy+sz)
            x1, x2 = max(0, cx-sz), min(4096, cx+sz)
            return data[y1:y2, x1:x2]

        # --- B. Load All Other Physics ---
        with fits.open(current_files['M']) as h:  d_M  = h[1].data if len(h)>1 else h[0].data
        with fits.open(current_files['V']) as h:  d_V  = h[1].data if len(h)>1 else h[0].data
        with fits.open(current_files['Ld']) as h: d_Ld = h[1].data if len(h)>1 else h[0].data
        with fits.open(current_files['Lw']) as h: d_Lw = h[1].data if len(h)>1 else h[0].data

        # --- C. Segmentation (Using Intensity) ---
        d_Ic_clean = remove_limb_darkening(d_Ic, center=(2048,2048), radius_pix=1900)
        crop_Ic = get_crop(d_Ic_clean)

        # Generate Masks
        masks = get_masks(crop_Ic,
                          umbra_range=(p['UMBRA_MIN'], p['UMBRA_MAX']),
                          penumbra_range=(p['PENUMBRA_MIN'], p['PENUMBRA_MAX']))

        # --- D. Generate 5 Plots ---
        # 1. Mu Grid for the crop
        mu_map = calculate_mu_grid(crop_Ic.shape,
                                   center=(sz - (cx-2048), sz - (cy-2048)),
                                   radius=1900)

        out_path = os.path.join(res_dir, f"{base_name}")

        # 1. Intensity (I_c)
        save_candle_plot(mu_map, crop_Ic, masks['umbra'],
                         "Intensity", "Norm", "red", f"{out_path}_CLV_Ic.png")

        # 2. Magnetogram (M) - Absolute Value
        save_candle_plot(mu_map, np.abs(get_crop(d_M)), masks['umbra'],
                         "Magnetic Field", "Gauss", "blue", f"{out_path}_CLV_Mag.png")

        # 3. Dopplergram (V)
        save_candle_plot(mu_map, get_crop(d_V), masks['umbra'],
                         "Velocity", "m/s", "green", f"{out_path}_CLV_Vel.png")

        # 4. Line Depth (Ld)
        save_candle_plot(mu_map, get_crop(d_Ld), masks['umbra'],
                         "Line Depth", "Depth", "purple", f"{out_path}_CLV_Ld.png")

        # 5. Line Width (Lw)
        save_candle_plot(mu_map, get_crop(d_Lw), masks['umbra'],
                         "Line Width", "Width", "orange", f"{out_path}_CLV_Lw.png")

    print(f"\n✅ Analysis Complete. 5-Parameter Plots saved in {res_dir}/")

if __name__ == "__main__":
    main()
