# `pyscses` - PYthon Space-Charge Site-Explicit Solver

<img src="https://github.com/bjmorgan/pyscses/blob/master/logo.png" width="200px"/>

[![status](http://joss.theoj.org/papers/803ed6dd19f453819bdd3ed9ceadf3b3/status.svg)](http://joss.theoj.org/papers/803ed6dd19f453819bdd3ed9ceadf3b3)
[![PyPI version](https://badge.fury.io/py/pyscses.svg)](https://badge.fury.io/py/pyscses)
[![DOI](https://zenodo.org/badge/90385184.svg)](https://zenodo.org/badge/latestdoi/90385184)
[![Documentation Status](https://readthedocs.org/projects/pyscses/badge/?version=latest)](https://pyscses.readthedocs.io/en/latest/?badge=latest)

`pyscses` is a Python package for modelling ionic space-charges in solid electrolytes. Its primary use is to calculate equilibrium distributions of point-charge atomic defects within one-dimensional &ldquo;Poisson-Boltzmann&rdquo;-like mean-field models. These calculations take as inputs a set of defect site positions, within a specific crystal structure, and the associated defect segregation energies. `pyscses` can also be used to calculate ionic transport properties (space-charge resistivities and activation energies) for these equilibrium defect distributions.

One approach to modelling space-charge formation in solid electrolytes is to consider defects as ideal point-charges embedded in a continuum dielectric, and to calculate equilibrium defect distributions by solving mean-field &ldquo;Poisson-Boltzmann&rdquo;-like equations. While numerical solutions to the 1D Poisson-Boltzmann equation are relatively simple to implement, published results are typically obtined using private closed-source codes, making it difficult to reproduce results or to test the effect of different approximations included in specific models. ``pyscses`` provides an open-source Python package for modelling space-charge formation in solid electrolytes, within a 1D Poisson-Boltzmann-like formalism.We are currently using ``pyscses`` in our own research into space-charge formation in solid electrolytes for fuel cells and lithium-ion batteries, and hope that this open-source resource will support reproducible research practices in future studies in this area.

``pyscses`` implements a range of numerical models within the general scheme of solving a modified Poisson-Boltzmann equation on a 1D grid:
- Continuum (regular grid) and site-explicit (irregular grid) models.
- Periodic and Dirichlet boundary conditions.
- &ldquo;Mott-Schottky&rdquo; and &ldquo;Gouy-Champman&rdquo; condtions. These are implemented by setting the mobilities of different defect species. In the case of Mott-Schottky conditions, all but one defect species have a mobility of zero.
- Inclusion of &ldquo;lattice&rdquo;site charges to account for non-defective species in the crystal structure.

Properties that can be calculated include:
- Defect mole fractions.
- Charge density.
- Electrostatic potential.
- Parallel and perpendicular grain boundary resistivities \[[1](#HwangEtAl_JElectroceram1999)\].
- Grain boundary activation energies \[[2](#Kim_PhysChemChemPhys2016)\].

Full mathematical derivations, definitions and example code can be found in the [userguide](https://github.com/bjmorgan/pyscses/blob/master/userguides/notebooks/userguide.ipynb). 

A more detailed overview of the code and its capabilities, and of the scientific context of modelling space-charge regions in solids, are given in the [JOSS paper](http://joss.theoj.org/papers/803ed6dd19f453819bdd3ed9ceadf3b3).

## Documentation
### Userguides
Introductory userguides are contained in the [userguide](https://github.com/bjmorgan/pyscses/blob/master/userguides/README.md). Each guide is presented as a [Jupyter notebook](http://jupyter-notebook.readthedocs.io/en/latest/#). The userguides cover the theory behind the Poisson-Boltzmann solver, how to set up a Jupyter notebook to run a calculation and examples of running the calculation under different approximations. These examples include site-explicit versus continuum models, Mott-Schottky (single mobile defect species) and Gouy-Chapman (all defect species mobile) conditions,  and running the solver on real data.

These userguides can also be found in the directory:
```
pyscses/userguides/userguide.ipynb
```
For online viewing of these userguides, we recommend using [nbviewer](https://nbviewer.jupyter.org). The links below open nbviewer versions of each userguide notebook.
- [Theory](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Theory.ipynb)
- [Setting up the notebook](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Setting_up.ipynb)
- [Running the calculation](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Running.ipynb)
- [Example 1 - continuum vs. site explicit and boundary conditions](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_1_BC.ipynb)
- [Example 2 - Mott-Schottky vs. Gouy-Chapman conditions](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_2_MSGC.ipynb)
- [Example 3 - Calculating the grain boundary resistivity and activation energy](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_3_Res.ipynb)
- [Example 4 - Comparison with &ldquo;experimental&rdquo; values](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_4_MSapp.ipynb)
- [Example 5 - Using the pyscses with real data](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_5_real_data.ipynb)
- [Example 6 - Storing the output data](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_6_store_data.ipynb)

### API Documentation
API documentation can be found [here](https://pyscses.readthedocs.io/en/latest/)

## Installation
Source code is available as a git repository at [https://github.com/bjmorgan/pyscses/tree/master/pyscses](https://github.com/bjmorgan/pyscses/tree/master/pyscses)
  
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

The directory `tests/test_notebooks` contains a set of Jupyter notebooks with specified outputs, that can be run to test code functionality. The test notebooks can be found on GitHub [here](https://github.com/bjmorgan/pyscses/tree/master/tests/test_notebooks).

There are four directories with varying conditions. In each there is a Jupyter notebook which can be run. The input for the calculations is stored in the `input_data` directory and the output for the calculations will be stored in the `generated_data` directory and can be compared to the verified data in the `expected_outputs` directory. A list of the input parameters used in the notebooks is reiterated in each of the four test folders in the `input_parameters` file. 

### Unit tests

Limited unit tests are contained in the top `tests` directory. These can be run using
```
cd tests
python -m unittest discover
```
Automated unit testing of the latest commit happens [here](https://travis-ci.org/bjmorgan/pyscses/).

## Contributing

### Bugs reports and feature requests

If you think you have found a bug, please report it on the [Issue Tracker](https://github.com/bjmorgan/pyscses/issues). This is also the place to propose ideas for new features or ask questions about the design of pyscses. Poor documentation is considered a bug, but please be as specific as possible when asking for improvements.

### Code contributions

We welcome your help in improving and extending the package with your own contributions. This is managed through GitHub pull requests; for external contributions we prefer the "fork and pull" workflow, while core developers use branches in the main repository:

- First open an [Issue](https://github.com/bjmorgan/pyscses/issues) to discuss the proposed contribution. This discussion might include how the changes fit pyscses' scope and a general technical approach.
- Make your own project fork and implement the changes there. Please keep your code style compliant with PEP8.
- Open a [pull request](https://github.com/bjmorgan/pyscses/pulls) to merge the changes into the main project. A more detailed discussion can take place there before the changes are accepted.

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

# References
1. <a name="HwangEtAl_JElectroceram1999" />[Hwang, J.-H., McLachlan, D. S., & Mason, T. O. (1999). Brick layer model analysis of nanoscale-to-microscale cerium dioxide. *J. Electroceram.*, **3** (1), 7–16.](https://dx.doi.org/10.1023/A:1009998114205)
1. <a name="Kim_PhysChemChemPhys2016" />[Kim, S., & Maier, J. (2002). On the conductivity mechanism of nanocrystalline ceria. *J. Electrochem. Soc.*, **149**(10), J73–J83.](https://dx.doi.org/10.1149/1.1507597)

