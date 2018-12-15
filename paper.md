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

``pyscses`` is a python module that implements a site-explicit, one dimensional Poisson-Boltzmann solver that can be used to model ionic space charge formation in solid materials. 
Space charge analyses typically treat space charge formation using continuum modelling [@ cite wahnstrom maier norby ]. ``pyscses`` takes explicitly calculated defect segregation energies and atomically resolved defect positions to calculate space charge properties. These properties include the electrostatic potential, charge density and defect distributions over the space charge region. The use of site explicit modelling allows the calculation of grain boundary resistivities and activation energies, which take into account the explicit grain boundary structure.

## Space charge formation
In polycrystalline solid materials, grain boundaries and interfaces exist separating different crystalline domains. The structural distortion means that defect concentrations and mobilities may deviate from their bulk values in localised regions of space, such as those close to grain boundaries. Defects typically segregate to, or away from the core of the grain boundary. This results in a grain boundary core which carries a net charge and an accumulation or depletion of defects in the regions adjacent to the grain boundary core, known as the space charge regions. Due to the relationship between conductivity and concentration, $ \sigma_i = c_i \mu_i z_i $, variation in defect concentration in these regions can strongly affect the the ionic conductivity of the material.  

## ``pyscses``
``pyscses`` implements a one dimensional Poisson-Boltzmann solver to model ionic space charge formation. The approach combines Boltzmann statistics with Poisson's equation to give a self-consistent method for calculating the electrostatic potential, charge density and defect distributions across the space charge region at equilibrium. This approach is available in ``pyscses``, however due to the site-explicit nature of crystalline materials, the default method in ``pyscses`` uses Fermi-Dirac statistics to give a more realistic expression of defect distributions in boundary layers. 
By equating the electrochemical potentials between the bulk and boundary layer at equilibrium, 

$$
\mu^o_{i,x} + RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}} \right) + z_i F \Phi_x = \mu^o_{i,\infty} + RT\ln \left(\frac{c_{i,\infty}}{1- c_{i,\infty}}\right) + z_i F \Phi_{\infty},
$$

the defect concentrations can be described,

$$
c_i = \frac{c_\infty \exp\left(\frac{-z_i\Phi_x + \mu_i}{kT}\right)}{1+ c_{\infty} \left(\exp\left( \frac{-z_i\Phi_x + \mu_i}{kT} \right) -1 \right) } . 
$$

The charge density is proportional to the concentration given by

$$
\rho = \sum_i c_i z_i F,
$$

and then the electrostatic potential can be calculated as given in Poisson's equation,

$$
\nabla^2 \Phi = \frac{ -\rho } { \epsilon \epsilon_0 }.
$$

``pyscses`` solves this as a second order partial differential equation using a second order finite difference approximation for each site,

$$
-\rho = \frac{ -2( \Delta x_2 \Phi_{x-1} - ( \Delta x_1 + \Delta x_2 ) \Phi_0 + \Delta x_1 \Phi_{x+1} ) } { \Delta x_1 \Delta x_2 ( \Delta x_1 + \Delta x_2 )  },
$$

combined with taking the sites as a set of linear equations,

$$
A = \Delta x_2 \Phi_{i-1} - ( \Delta x_1 + \Delta x_2 ) \Phi_0 + \Delta x_1 \Phi_{i+1},
$$

$$
p = \frac{2} {\Delta x_1 \Delta x_2 ( \Delta x_1 + \Delta x_2 )},
$$

$$
b = - \rho,
$$

$$
L =( A^T \cdot p )^T,
$$

$$
\Phi = L^{-1}\vec{b},
$$

which is solved using matrix inversion implemented with ``scipy.sparse``.

``pyscses`` extends the model by calculating the grain boundary resistivity and activation energies from the calculated defect distributions. 

Taking each site as a resistor in series, the perpendicular grain boundary resistivity is calculated by taking the ratio of the sum of the resistivity in the space charge region to the resistivity in the bulk,

$$
\r_{GB}^{+} = \frac{c_{i,\infty}\mu_{i,\infty}z_i}{\sum{i} c_{i,x}\mu_{i,x}z_i}.
$$

Taking each site as a resistor in parallel the parallel grain boundary is calculated by taking the inverse of the ratio of the sum of the conductivity in the space charge region to the conductivity in the bulk,
$$
\r_{GB}^{=} = \frac{1}{\frac{ c_{i,x}\mu{i,x}z_i}{c_{i,\infty}mu_{i,\infty}z_i}}
$$
 
The temperature dependance of these resistivity ratios is then used to calculate the grain boundary activation energy by applying an Arrhenius equation.

## Approximations
``pyscses`` accounts for a number of different approximations typically assumed when space charge formation is considered. These approximations include applying site-explicit modelling, where each site is taken at the atomically resolved defect position; or continuum modelling where the explicitly calculated defect segregation energies are interpolated onto a regular grid, and it is assumed that defects exist at all points. Another approximation considered is applying  Mott-Schottky conditions, where only one defect species is considered mobile in the system or Gouy-Chapman conditions, where all defect species are considered mobile in the system. ``pyscses`` can also consider the effect of the charge of the non-defective species being included when calculating the charge density over the space charge region and the boundary conditions of the calculation can be defined as either periodic or Dirichlet when the calculation is set up. 

## Limitations
While ``pyscses`` executes a computationally inexpensive and conceptually simple method to calculating space charge properties the model involves some assumptions about the defect behaviour. These assumptions include that the system is in the dilute limit, including only electrostatic interactions and that defects are point charges with no mass. 

## Acknowledgements


