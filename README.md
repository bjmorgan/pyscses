# `pyscses` - PYthon Space-Charge Site-Explicit Solver

<img src="https://github.com/bjmorgan/pyscses/blob/master/logo.png" width="200px"/>

[![PyPI version](https://badge.fury.io/py/pyscses.svg)](https://badge.fury.io/py/pyscses)
[![DOI](https://zenodo.org/badge/90385184.svg)](https://zenodo.org/badge/latestdoi/90385184)
[![Documentation Status](https://readthedocs.org/projects/pyscses/badge/?version=latest)](https://pyscses.readthedocs.io/en/latest/?badge=latest)

`pyscses` is a Python module that implements a site-explicit, one-dimensional Poisson-Boltzmann solver, used for modelling ionic space charge properties in solid materials. Space charge properties such as electrostatic potential, charge density and charge carrier distributions over the space charge region can be calculated using the Poisson-Boltzmann equation from the input of defect segregation energies and atomically resolved charge carrier positions. The grain boundary resistivity and activation energy can be calculated by extending the model using the calculated charge carrier distributions. `pyscses` also accounts for different approximations typically assumed when space charge formation is considered. These approximations include site explicit vs. continuum modelling, Mott-Schottky (single mobile defect species) and Gouy-Chapman (all defect species mobile) conditions, and whether the charge of the non-defective species should be considered. Full mathematical derivations, definitions and example code can be found in the [userguide](https://github.com/bjmorgan/pyscses/blob/master/userguides/notebooks/userguide.ipynb).

API documentation can be found [here](https://pyscses.readthedocs.io/en/latest/)

Source code is available as a git repository at [https://github.com/bjmorgan/pyscses/tree/master/pyscses](https://github.com/bjmorgan/pyscses/tree/master/pyscses)
  
## Installation

The simplest way to install `pyscses` is to use `pip` to install from [PyPI](https://pypi.org/project/pyscses/)
```
pip install pyscses
```

Alternatively, you can download the latest release from [GitHub](https://github.com/bjmorgan/pyscses/releases), and install directly:
```
cd pyscses
pip install -e .
```
which installs an editable (-e) version of pyscses in your userspace.

Or clone the latest version from [GitHub](https://github.com/bjmorgan/pyscses/releases) with
```
git clone git@github.com:bjmorgan/pyscses.git
```
and install the same way.
```
cd pyscses
pip install -e .
```
## Tests

Jupyter notebooks that can be used to check the output of the calculations can be found in
```
pyscses/tests/test_notebooks
```
The test notebooks can be found on github [here](https://github.com/bjmorgan/pyscses/tree/master/tests/test_notebooks).

There are four directories with varying conditions. In each there is a Jupyter notebook which can be run. The input for the calculations is stored in the `input_data` directory and the output for the calculations will be stored in the `generated_data` directory and can be compared to the verified data in the `expected_outputs` directory. A list of the input parameters used in the notebooks is reiterated in each of the four test folders in the `input_parameters` file. 

## Documentation
Once installed, the `pyscses` code is imported into, and run using a [Jupyter notebook](http://jupyter-notebook.readthedocs.io/en/latest/#).
An overview of the capabilities of `pyscses`, along with example code for running the code and varying the simulation conditions can be found in
```
pyscses/userguides/userguide.ipynb
```
or the Jupyter notebook can be found on github [here](https://github.com/bjmorgan/pyscses/blob/master/userguides/notebooks/userguide.ipynb) .

API documentation is available [here](https://pyscses.readthedocs.io/en/latest/) .

## Scientific context

In polycrystalline solid materials, grain boundaries and interfaces exist separating different crystalline domains. The structural distortion at these interfaces causes segregation of charge carriers to, or away from the grain boundary core. Due to this, the grain boundary core carries a net charge which causes the depletion or accumulation of charge carriers in the regions adjacent, known as space charge regions. Due to the variation on charge carrier concentrations, the ionic conductivity of the material can be strongly affected by the presence of grain boundaries.

## Citing `pyscses`

This code can be cited as:

Wellock, Georgina L., & Morgan, Benjamin J. (2019). *pyscses - a PYthon Space-Charge Site-Explicit Solver* Zenodo. http://doi.org/10.5281/zenodo.2536867

### BibTeX

```
@misc{wellock_georgina_l_2019_2536901,
  author       = {Wellock, Georgina L. and
                  Morgan, Benjamin J.},
  title        = {{pyscses - a PYthon Space-Charge Site-Explicit 
                   Solver}},
  month        = jan,
  year         = 2019,
  doi          = {10.5281/zenodo.2536901},
  url          = {https://doi.org/10.5281/zenodo.2536867}
}
```

