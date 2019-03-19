# Test Notebooks

These test notebooks run through typical calculations that can be performed to calculate space charge properties using classically calculated data for a Gd-doped ceria system. To verify the results from running the calculations, the outputs from running each notebook will be stored in a `generated_outputs` directory and can be compared to the controlled outputs previously calculated in an `expected_outputs` directory. The notebooks run the calculations using a variety of typical calculation assumptions, including site explicit and continuum modelling, Mott-Schottky and Gouy-Chapman approximations and the inclusion of a site charge term. 

Running these notebooks requires additional Python dependencies. These can be installed using
```
cd tests/test_notebooks
pip install -r requirements.txt
```

## Notebooks

- [Test 1 - Site explicit modelling under Mott-Schottky conditions.](test_1/test_notebook_1.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/tests/test_notebooks/test_1/test_notebook_1.ipynb))
- [Test 2 - Continuum modelling under Mott-Schottky conditions.](test_2/test_notebook_2.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/tests/test_notebooks/test_2/test_notebook_2.ipynb))
- [Test 3 - Site explicit modelling under Gouy-Chapman conditions.](test_3/test_notebook_3.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/tests/test_notebooks/test_3/test_notebook_3.ipynb))
- [Test 4 - Site explicit modelling under Mott-Schottky conditions including the additional site charge.](test_4/test_notebook_4.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/tests/test_notebooks/test_4/test_notebook_4.ipynb))

## Running times

Some of these notebooks contain relatively intensive calculations. If you rerun these from scratch, in some cases please be prepared to wait for the calculation to reach completion.  
Typical running times for each notebook (4 GHz i7 iMac):

- Test 1: 2 minutes
- Test 2: 6.5 minutes
- Test 3: 27 minutes
- Test 4: 3 minutes

