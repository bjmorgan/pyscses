[![Documentation Status](https://readthedocs.org/projects/gwpb/badge/?version=latest)](https://gwpb.readthedocs.io/en/latest/?badge=latest)

# `pyscses` - PYthon Space Charge Site Explicit Solver

<img src="logo.png" width="100"/>

`pyscses` is a Python module that implements a one-dimensional Poisson-Boltzmann solver, used for modelling space charge properties in solid materials. 
In polycrystalline solid materials, grain boundaries and interfaces exist separating different crystalline domains. The structural distortion at these interfaces causes segregation of charge carriers to, or away from the grain boundary core. Due to this, the grain boundary core carries a net charge which causes the depletion or accumulation of charge carriers in the regions adjacent, known as space charge regions. Due to the variation on charge carrier concentrations, the ionic conductivity of the material can be strongly affected by the presence of grain boundaries.

Taking input of defect segregation energies and atomically resolved charge carrier positions, `pyscses` can calculate space charge properties including the electrostatic potential, charge density and charge carrier distributions over the space charge region using the Poisson-Boltzmann equation. Using the calculated charge carrier distributions, the grain boundary resistivity and activation energy can be calculated. `pyscses` also accounts for different approximations typically assumed when space charge formation is considered. Full mathematical derivations, definitions and example code can be found in the [userguide](https://github.com/georgiewellock/PYSCSES/blob/master/userguides/notebooks/userguide.ipynb).

API documentation can be found [here](https://gwpb.readthedocs.io/en/latest/)

Source code is available as a git repository at [https://github.com/georgiewellock/PYSCSES/tree/master/pyscses](https://github.com/georgiewellock/PYSCSES/tree/master/pyscses)
  
## Installation

```
pip install pyscses
```

Or download the latest release from [GitHub](https://github.com/georgiewellock/PYSCSES/releases), and install
```
cd PYSCSES
python setup.py install
```

Or clone the latest development version
```
git clone git@github.com:georgiewellock/PYSCSES.git
```
and install the same way.
```
cd PYSCSES
python setup.py install 
```
## Tests

Some limited unit testing is available and can be run manually. First enter the directory with the test files
```
cd PYSCSES/tests/unit_tests
```
and then the tests can be run using
```
python -m unittest discover
```

Jupyter notebooks that can be used to check the output of the calculations can be found in
```
PYSCSES/tests/test_notebooks
```
The test notebooks can be found on github [here](https://github.com/georgiewellock/PYSCSES/tree/master/tests/test_notebooks).

There are four directories with varying conditions. In each there is a Jupyter notebook which can be run. The input for the calculations is stored in the `input_data` directory and the output for the calculations will be stored in the `generated_data` directory and can be compared to the verified data in the `expected_outputs` directory. A list of the input parameters used in the notebooks is reiterated in each of the four test folders in the `input_parameters` file. 

## Documentation
Once installed, the `pyscses` code is imported into, and run using a [Jupyter notebook](http://jupyter-notebook.readthedocs.io/en/latest/#).
An overview of the capabilities of `pyscses`, along with example code for running the code and varying the simulation conditions can be found in
```
PYSCSES/userguides/userguide.ipynb
```
or the Jupyter notebook can be found on github [here](https://github.com/georgiewellock/PYSCSES/blob/master/userguides/notebooks/userguide.ipynb) .

API documentation is available [here](https://gwpb.readthedocs.io/en/latest/) .
## Citing `pyscses`
