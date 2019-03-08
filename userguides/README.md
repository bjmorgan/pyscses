# Userguides

These userguides cover the theory behind the Poisson-Boltzmann solver and go through setting up the Jupyter notebooks in order to calculate the space charge properties. Following this, the different approximations are explored and examples given to show how the calculated properties vary concluding with running the solver on some real crystalline data. 

Running these userguides requires additional Python dependencies. These can be installed using
```
cd userguides
pip install -r requirements.txt
```

## Notebooks

Each notebook can be viewed using [nbviewer](https://nbviewer.jupyter.org) by following each [nbviewer]() link.
- [Theory](notebooks/Theory.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Theory.ipynb))
- [Setting up the notebook](notebooks/Setting_up.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Setting_up.ipynb))
- [Running the calculation](notebooks/Running.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Running.ipynb))
- [Example 1 - continuum vs. site explicit and boundary conditions](notebooks/Ex_1_BC.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_1_BC.ipynb))
- [Example 2 - Mott-Schottky vs. Gouy-Chapman conditions](notebooks/Ex_2_MSGC.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_2_MSGC.ipynb))
- [Example 3 - Calculating the grain boundary resistivity and activation energy](notebooks/Ex_3_Res.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_3_Res.ipynb))
- [Example 4 - Comparison with &ldquo;experimental&rdquo; values](notebooks/Ex_4_MSapp.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_4_MSapp.ipynb))
- [Example 5 - Using the pyscses with real data](notebooks/Ex_5_real_data.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_5_real_data.ipynb))
- [Example 6 - Storing the output data](notebooks/Ex_6_store_data.ipynb) ([nbviewer](https://nbviewer.jupyter.org/github/bjmorgan/pyscses/blob/master/userguides/notebooks/Ex_6_store_data.ipynb))

## Running times

Some of these notebooks contain relatively intensive calculations. If you rerun these from scratch, in some cases please be prepared to wait for the calculation to reach completion.  
Typical running times for each notebook (4 GHz i7 iMac):
- Example 1: 81 minutes
- Example 2: 32 minutes
- Example 3: 3 minutes
- Example 4: 3 minutes
- Example 5: 2 minutes
- Example 6: 2 minutes

## Notes
Calculations run using Gouy-Chapman conditions take significantly longer to run than those with Mott-Schottky conditions. This is due to the Poisson-Boltzmann solver being run repeatedly in order to minimise the difference between the input and output mole fractions.

Calculations run using a continuum approximation take longer with increasing numbers of grid points. The default in the userguides is 1000 grid points.
