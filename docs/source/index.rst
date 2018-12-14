.. pyscses documentation master file, created by
   sphinx-quickstart on Fri Mar  9 11:37:43 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyscses - PYthon Space Charge Site Explicit Solver
====================================================

``pyscses`` is a Python module that implements a site-explicit, one-dimensional Poisson-Boltzmann solver, used for modelling ionic space charge properties in solid materials. Space charge properties such as electrostatic potential, charge density and charge carrier distributions over the space charge region can be calculated using the Poisson-Boltzmann equation from the input of defect segregation energies and atomically resolved charge carrier positions. The grain boundary resistivity and activation energy can be calculated by extending the model using the calculated charge carrier distributions. ``pyscses`` also accounts for different approximations typically assumed when space charge formation is considered. These approximations include site explicit vs. continuum modelling, Mott-Schokkty (single mobile defect species) and Gouy-Chapman (all defect species mobile) conditions, and whether the charge of the non-defective species should be considered. Full mathematical derivations, definitions and example code can be found on github in the `userguide <https://github.com/georgiewellock/PYSCSES/blob/master/userguides/notebooks/userguide.ipynb>`.

API documentation can be found `here <https://gwpb.readthedocs.io/en/latest/>`

Source code is available as a git repository at `https://github.com/georgiewellock/PYSCSES/tree/master/pyscses <https://github.com/georgiewellock/PYSCSES/tree/master/pyscses>`

Scientific context
=================
In polycrystalline solid materials, grain boundaries and interfaces exist separating different crystalline domains. The structural distortion at these interfaces causes segregation of charge carriers to, or away from the grain boundary core. Due to this, the grain boundary core carries a net charge which causes the depletion or accumulation of charge carriers in the regions adjacent, known as space charge regions. Due to the variation on charge carrier concentrations, the ionic conductivity of the material can be strongly affected by the presence of grain boundaries.

Tests
=====

Jupyter notebooks that can be used to check the output of the calculations can be found in

::
    PYSCSES/tests/test_notebooks

The test notebooks can be found on github `here <https://github.com/georgiewellock/PYSCSES/tree/master/tests/test_notebooks>`.

.. toctree::
   :caption: API documentation
   :maxdepth: 3

   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
