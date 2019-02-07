# Userguides

These userguides cover the theory behind the Poisson-Boltzmann solver and go through setting up the Jupyter notebooks in order to calculate the space charge properties. Following this, the different approximations are explored and examples given to show how the calculated properties vary concluding with running the solver on some real crystalline data. 

To run these userguides, there are additional dependencies. These can be installed using

```
pip install -r requirements.txt
```

from inside the `userguides` directory.

**Notebooks**

[Theory](notebooks/Theory.ipynb)

[Setting up the notebook](notebooks/Setting_up.ipynb)

[Running the calculation](notebooks/Running.ipynb)

[Example 1 - continuum vs. site explicit and boundary conditions](notebooks/Ex_1_BC.ipynb)

[Example 2 - Mott-Schottky vs. Gouy-Chapman conditions](notebooks/Ex_2_MSGC.ipynb)

[Example 3 - Calculating the grain boundary resistivity and activation energy](notebooks/Ex_3_Res.ipynb)

[Example 4 - Comparison with 'experimental' values](notebooks/Ex_4_MSapp.ipynb)

[Example 5 - Using the pyscses with real data](notebooks/Ex_5_real_data.ipynb)

[Example 6 - Storing the output data](notebooks/Ex_6_store_data.ipynb)

**Run times**

Expected run times for each notebook ( Using an iMac with a 4 Ghz i7 processor):

Example 1: 81 minutes

Example 2: 32 minutes

Example 3: 3 minutes

Example 4: 3 minutes

Example 5: 2 minutes

Example 6: 2 minutes

**Note** - Calculations run using Gouy-Chapman conditions take significantly longer to run. This is due to the Poisson-Boltzmann solver being run repeatedly in order to minimise the difference between the input and output mole fractions.

**Note** - Calculations run using a continuum approximation take longer depending on the number of points. Default in the userguides = 1000.
