Installation
============

SpotiPy requires Python 3.9 or later. It is recommended to install the package in an editable mode so that any local modifications to the pipeline are immediately reflected in your environment.

Prerequisites: JSOC Registration
--------------------------------
SpotiPy's downloading module interacts directly with the Joint Science Operations Center (JSOC) to retrieve SDO/HMI and AIA data. To authorize these downloads, you **must** register your email address with JSOC.

1. Visit the `JSOC Email Registration page <http://jsoc.stanford.edu/ajax/register_email.html>`_.
2. Register the email address you intend to use.
3. Save this email address in your ``params.txt`` file (see :doc:`configuration`).

Installing SpotiPy
------------------
Clone the repository from GitHub and install it using ``pip``:

.. code-block:: bash

   git clone https://github.com/Emily-Joe/SpotiPy.git
   cd SpotiPy
   pip install -e .

The ``-e`` flag installs the package in "editable" mode. This will automatically install all required scientific dependencies, including ``numpy``, ``astropy``, ``sunpy``, and ``scipy``.
