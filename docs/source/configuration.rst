Configuration (params.txt)
==========================

The research pipeline is controlled via ``params.txt``. This text-based configuration ensures that analysis parameters are documented and reproducible.

Target and Query Parameters
---------------------------
* **NOAA_NUMBER**: Active Region identifier (e.g., ``13034``).
* **START_DATE**: ISO 8601 start time (e.g., ``2019-04-07T18:17:00``).
* **DAYS**: Duration of the observation window.
* **CADENCE**: Time step between downloaded frames (hours).
* **EMAIL**: JSOC-registered email address.

Trajectory and Extraction
-------------------------
* **X_START_ARCSEC / Y_START_ARCSEC**: Initial helioprojective coordinates of the target feature.
* **FRAME_SIZE**: The pixel dimensions of the square extraction window.

Physical Segmentation Thresholds
--------------------------------
HMI segmentation utilizes an 8-bit scale where the Quiet Sun median is approximately 128.

* **UMBRA_MIN / MAX**: Thresholds for umbral pixels (~8% to 43% of Quiet Sun level).
* **PENUMBRA_MIN / MAX**: Thresholds for penumbral pixels (~59% to 94% of Quiet Sun level).
* **PLAGE_EXCESS_PCT**: Percentage excess above the Quiet Sun median (typically ``20.0``).
* **QUIET_SUN_TOL_PCT**: Tolerance for network classification (typically ``15.0``).
* **MIN_PLAGE_AREA**: Minimum contiguous pixel area required for plage classification to filter transient UV noise.
