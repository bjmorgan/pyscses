---
title: 'pyscses: python space charge site explicit solver'
tags:
 - Python
 - poisson-boltzmann
 - space charge formation
 - solid electrolytes
authors:
 - name: Georgina L. Wellock
 - email: g.l.wellock@bath.ac.uk
 - email: georgiewellock@gmail.com
 - affiliation: "1"
 - name: Benjamin J. Morgan
 - email: b.j.morgan@bath.ac.uk
 - orcid: 0000-0002-3056-8233
 - affiliation: "1"
affiliations:
 - name: Department of Chemistry, University of Bath, Claverton Down, UK, BA2 7AY
   index: 1
date: 15 December 2018
bibliography: paper.bib
---

# Summary

``pyscses`` is an open-source python package for numerical modelling of ionic space-charge regions in crystalline solids. Ionic space-charges are microscopic regions with a local excess or deficiency of charged point defects, giving these regions a net electrostatic charge. Space-charge regions form at crystallographic interfaces, such as grain boundaries or heterointerfaces between two different materials. The formation of space-charges is due to segregation of charged defects to, or from, these interfaces, and the associated redistribution of defects in adjacent regions of the crystal, due to defect&endash;defect electrostatic interactions. 

The segregation of defects to, or from, space-charge regions can produce local defect concentrations that strongly deviate from the average &ldquo;bulk&rdquo; values in a material. This can have significant consequences for key material properties. For example, in solid electrolytes, where ionic defects act as mobile charge carriers, space-charge formation at grain boundaries is associated with large changes in ionic conductivity both parallel and perpendicular to the grain boundaries, due to phenomena such as enhanced defect concentrations in grain bounday cores, or reduced defect concentrations in adjacent &ldquo;depletion&rquod; space-charge regions. In the case of solid electrolytes, understanding space-charge formation at grain boundaries and interfaces is a key challenge in developing a theoretical description of the role of crystalline microstructure on macroscopic ion transport (ionic conductivities).

``pyscses`` considers simple one-dimensional models of crystallographic interfaces, and calculates equilibrium defect distributions by solving a modified Poisson-Boltzmann equation. The driving force for defect segregation to or from the interface is described by defect segregation energies, which are defined as

TODO: GEORGIE, can you add the equation for the segregation energies?

These defect segregation energies are typically obtained from standard first-principles or classical atomistic calculations, and are then used as inputs for calculations performed with ``pyscses``. ``pyscses`` can solve conventional ``continuum`` models, using a regular 1D grid, and ``site-explicit`` models, using an irregular 1D grid, where grid points correspond to explicit crystallographic site positions for the material being modelled. Each site can have an associated segregation energy for one or more defect species, which avoids the need to define arbitrary &ldquo;core&rdquo; and &ldquo;bulk&rdquo; regions in the model.

# The numerical model
The default model implemented in ``pyscses`` calculates the equilibrium distribution of point charge defect species distributed on a 1D grid, with a defect electrochemical potential `$\mu^o_{i,x}$` of
$$
\mu^o_{i,x} + RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}} \right) + z_i F \Phi_x.
$$
The `$RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}}$` term gives the chemical potential for a non-interacting lattice-gas with site exclusion (maximum site occupancy of 1).

Finding the equilibrium distribution of defects is equivalent to solving a modified Poisson-Boltzmann equation, which can be derived by requiring the electrochemical potentials at each site to be equal to a set of reference bulk electrochemical potentials.
$$
\mu^o_{i,x} + RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}} \right) + z_i F \Phi_x = \mu^o_{i,\infty} + RT\ln \left(\frac{c_{i,\infty}}{1- c_{i,\infty}}\right) + z_i F \Phi_{\infty},
$$
the defect concentrations can be described,
$$
c_i = \frac{c_\infty \exp\left(\frac{-z_i\Phi_x + \mu_i}{kT}\right)}{1+ c_{\infty} \left(\exp\left( \frac{-z_i\Phi_x + \mu_i}{kT} \right) -1 \right) } . 
$$
The charge density is proportional to the concentration given by
$$
\rho = \sum_i c_i z_i F.
$$
The electrostatic potential can be calculated using Poisson's equation,
$$
\nabla^2 \Phi = \frac{ -\rho } { \epsilon \epsilon_0 }.
$$
``pyscses`` solves this as a second order partial differential equation using a second order finite difference approximation for each site.

TODO: Add the various model types you can solve that have not yet been discussed, e.g. Dirichlet vs. Periodic boundary conditions. Ability to selevtively fix specific defects, allowing modelling of &ldquo;Mott-Schottky&rdquo; and &ldquo;Gouy-Chapman&rdquo; conditions. Anything else I have forgotten?

# Typical workflow

TODO: GEORGIE, can you briefly describe here how you run a standard GB calculation. I would include figures from your poster showing the segregation energies arranged on the 1D grid, and some example outputs, e.g. space-charge potential, defect concentrations, and charge densities.

# Calculated properties

In addition to allowing calculation of defect, charge, and potential across an interface, ``pyscses`` can calculate interface (grain boundary) resistivities and activation energies. 

## Grain boundary resistivities

Perpendicular grain boundary resistivities are calculated by modelling the system as a set of resistors in series:
$$
r_{GB}^\perp = \frac{c_{i,\infty}\mu_{i,\infty}z_i}{\sum_{i} c_{i,x}\mu_{i,x}z_i}.
$$
Parallel grain boundary resistivities are calculated by modelling the system as a set of resistors in parallel:
$$
r_{GB}^\parallel = \frac{1}{\frac{ c_{i,x}\mu{i,x}z_i}{c_{i,\infty}mu_{i,\infty}z_i}}
$$
TODO: is the parallel equation missing a sum?
 
The temperature dependance of these resistivity ratios is then used to calculate the grain boundary activation energy by applying an Arrhenius equation.

## Activation energies
TODO: say something about these

# Approximations and limitations
``pyscses`` implements a modification of the Poisson-Boltzmann equation, which assumes that defects are non-interacting except for point-charge electrostatics and site exclusion.

TODO:
- resistivitity and activation energy calculations assume fixed mobilities for defects

# Acknowledgements

TODO: G. L. W. Acknowledgements (EPSRC grant code?)
B. J. M. acknowledges support from the Royal Society (Grant No. UF130329).

# References
