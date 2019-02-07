`pyscses` is a Python module that implements a site-explicit, one-dimensional Poisson-Boltzmann solver, used for modelling ionic space charge properties in solid materials.

Space charge properties such as electrostatic potential, charge density and charge carrier distributions over the space charge region can be calculated using the Poisson-Boltzmann equation from the input of defect segregation energies and atomically resolved charge carrier positions. The grain boundary resistivity and activation energy can be calculated by extending the model using the calculated charge carrier distributions. `pyscses` also accounts for different approximations typically assumed when space charge formation is considered. These approximations include site explicit vs. continuum modelling, Mott-Schokkty (single mobile defect species) and Gouy-Chapman (all defect species mobile) conditions.

These userguides cover the theory behind the Poisson-Boltzmann solver and go through setting up the Jupyter notebooks in order to calculate the space charge properties. Following this, the different approximations are explored and examples given to show how the calculated properties vary concluding with running the solver on some real crystalline data. 

To run these userguides, there are additional dependencies. These can be installed using

```
pip install -r requirements.txt
```

from inside the `userguides` directory.

Expected run times for each notebook ( Using an iMac with a 4 Ghz i7 processor):
Example 1: 81 minutes
Example 2: 32 minutes
Example 3: 3 minutes
Example 4: 3 minutes
Example 5: 2 minutes
Example 6: 2 minutes

Note - Calculations run using Gouy-Chapman conditions take significantly longer to run. This is due to the Poisson-Boltzmann solver being run repeatedly in order to minimise the difference between the input and output mole fractions.

Note - Calculations run using a continuum approximation take longer depending on the number of points. Default in the userguides = 1000.
