[![Documentation Status](https://readthedocs.org/projects/gwpb/badge/?version=latest)](https://gwpb.readthedocs.io/en/latest/?badge=latest)

# `pyscses` - PYthon Space Charge Site Explicit Solver

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

Some limited unit testing is available and can be run manually
`cd tests/unit_tests`
using
`python -m unittest discover`.

Jupyter notebooks that can be used to check the output of the calculations can be found in
`pyscses/tests/test_notebooks`
github link [here](https://github.com/georgiewellock/PYSCSES/tree/master/tests/test_notebooks).
There are four directories with varying conditions. In each there is a Jupyter notebook which can be run. The input for the calculations is stored in the `input_data` directory and the output for the calculations will be stored in the `generated_data` directory and can be compared to the verified data in the `expected_outputs` directory. A list of the input parameters used in the notebooks is reiterated in each of the four test folders in the `input_parameters` file. 

## Documentation
Once installed, the `pyscses` code is imported into, and run using a [Jupyter notebook](http://jupyter-notebook.readthedocs.io/en/latest/#).
An overview of the capabilities of `pyscses`, along with example code for running the code and varying the simulation conditions can be found [here](https://github.com/georgiewellock/PYSCSES/blob/master/userguides/notebooks/userguide.ipynb) .

API documentation is available [here](https://gwpb.readthedocs.io/en/latest/) .
## Citing `pyscses`
