Usage Guide
===========

SpotiPy supports two primary research workflows: an automated end-to-end analysis pipeline and a modular library for custom solar physics applications.

Pipeline Mode
-------------
The ``run_analysis.py`` script executes the full scientific pipeline: **Download → Align → Limb Darkening Removal → Segmentation → Tracking**.

1. **Configuration:** The pipeline is driven by the ``params.txt`` file (see :doc:`configuration`).
2. **Execution:** Run ``python run_analysis.py``.
3. **Interactive Refinement:** The pipeline will prompt for manual coordinate selection using ``gui_tools`` to define the initial feature position.
4. **Active Tracking:** The system utilizes a differential rotation model to predict feature movement, supplemented by the ``refine_centering`` function which "locks on" to the darkest pixels of the sunspot umbra to account for local proper motions.

**Scientific Output:**
Results are saved to ``results_NOAA_{noaa}/``, including Center-to-Limb Variation (CLV) scatter plots and a consolidated ``.npz`` archive containing all cropped frames, masks, and heliographic coordinates.

Library Mode
------------
Individual components can be integrated into existing workflows. For example, to apply the SpotiPy segmentation logic to an external HMI continuum array:

.. code-block:: python

   from spotipy.segmentation import get_masks

   # Generate umbra and penumbra masks using 8-bit mapped thresholds
   masks = get_masks(
       hmi_data,
       umbra_range=(10, 55),
       penumbra_range=(75, 120)
   )

For detailed information on all available functions, see the :doc:`autoapi/index` reference.
