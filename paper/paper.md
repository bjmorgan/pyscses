---
title: 'pyscses: A python space-charge site-explicit solver'
tags:
 - Python
 - poisson-boltzmann
 - space-charges
 - solid electrolytes
authors:
 - name: Georgina L. Wellock
   email: g.l.wellock@bath.ac.uk
   email: georgiewellock@gmail.com
   orcid: 0000-0002-1068-5889
   affiliation: "1"
 - name: Benjamin J. Morgan
   email: b.j.morgan@bath.ac.uk
   orcid: 0000-0002-3056-8233
   affiliation: "1, 2"
affiliations:
 - name: Department of Chemistry, University of Bath, Claverton Down, UK, BA2 7AY
   index: 1
 - name: The Faraday Institution, Quad One, Harwell Science and Innovation Campus, Didcot, UK
   index: 2
date: 15 December 2018
bibliography: paper.bib
---

# Summary

``pyscses`` is an open-source python package for numerical modelling of ionic space-charge regions in crystalline solids. Ionic space-charges are microscopic regions with a local excess or deficiency of charged point defects, giving these regions a net electrostatic charge[@Frenkel_KinTheLiq1946]. Space-charge regions form at crystallographic interfaces, such as grain boundaries or heterointerfaces between two different materials. The formation of space-charges is due to segregation of charged defects to, or from, these interfaces, and the associated redistribution of defects in adjacent regions of the crystal, due to defect&ndash;defect electrostatic interactions [@Maier_BerBunsenges1984; @Maier_JPhysChemSol1985; @Maier_SolStatIonics2003; @ChiangEtAl_ApplPhysLett1996; @KimAndMaier_JElectrochemSoc2002]. 

The segregation of defects to, or from, space-charge regions can produce local defect concentrations that strongly deviate from the average &ldquo;bulk&rdquo; values in a material[@FleigAndMaier_SolidStateIonics1996]. This can have significant consequences for key material properties. For example, in solid electrolytes, where ionic defects act as mobile charge carriers, space-charge formation at grain boundaries is associated with large changes in ionic conductivity both parallel and perpendicular to the grain boundaries, due to phenomena such as enhanced defect concentrations[@JamnikAndMaier_SolidStateIonics1999] in grain bounday cores, or reduced defect concentrations in adjacent &ldquo;depletion&rdquo; space-charge regions. In the case of solid electrolytes, understanding space-charge formation at grain boundaries and interfaces is a key challenge in developing a theoretical description of the role of crystalline microstructure on macroscopic ion transport (ionic conductivities)[@KimEtAl_PhysChemChemPhys2003].

``pyscses`` considers simple one-dimensional models of crystallographic interfaces, and calculates equilibrium defect distributions by solving a modified Poisson-Boltzmann equationi[@Maier_ProgSolStatChem1995; @DeSouzaEtAl_SolidStateIonics2011]. The driving force for defect segregation to or from the interface is described by defect segregation energies, which are defined as

$$
\Delta E_\mathrm{seg}^{i,x} = E_\mathrm{f}^{i,x} - E_\mathrm{f}^{i, \infty}
$$

where $E_\mathrm{f}^{i,x} = E_\mathrm{tot}^{i,\mathrm{lattice}} - E_\mathrm{tot}^\mathrm{lattice}$.

# The numerical model
The default model implemented in ``pyscses`` calculates the equilibrium distribution of point charge defect species distributed on a 1D grid [@HelgeeEtAl_FuelCells2012; @PolfusEtAl_SolStatIonics2016; @DeSouzaEtAl_SolidStateIonics2011; @NymanEtAl_ApplPhysLet2012], with a defect electrochemical potential $\mu^o_{i,x}$ of
$$
\mu^o_{i,x} + RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}} \right) + z_i F \Phi_x.
$$
The $RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}}\right)$ term gives the chemical potential for a non-interacting lattice-gas with site exclusion (maximum site occupancy of 1).

Finding the equilibrium distribution of defects is equivalent to solving a modified Poisson-Boltzmann equation, which can be derived by requiring the electrochemical potentials at each site to be equal to a set of reference bulk electrochemical potentials [@Maier_IonicsTextbook2005].
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

``pyscses`` can enforce a number of different approximations commonly assumed when space charge formation in considered, including:

- Continuum and site-explicit modelling.

- periodic and Dirichlet boundary conditions.

- &ldquo;Mott-Schottky&rdquo; and &ldquo;Gouy-Champman&rdquo; condtions (by selecting specific defect species as &ldquo;fixed&rdquo;, i.e. immobile).

- Inclusion of site charge for non-defective species.

# Typical workflow

To run a standard calculation on a grain boundary, the calculated defect segregation energies and explicit atomic position are projected onto a one-dimensional grid. If the calculation is run using a continuum approximation, the defect segregation energies and atomic positions are interpolated onto a regular grid.

![A crystal structure is generated and the explicit site positions are projected to create a one-dimensional grid. The segregation energy of each defect can then be associated with the appropriate grid point.](Figures/seg_energies_joss.pdf)

The Poisson-Boltzmann equation is solved self-consistently, giving the electrostatic potential, charge density and defect mole fractions across the space charge region.

![Example output for Gd-doped ceria. Comparison between the output from the Poisson-Boltzmann solver using continuum and site-explicit modelling](Figures/continuum_vs_se_joss_MS.pdf)

## Calculated properties

In addition to allowing calculation of defect, charge, and potential across an interface, ``pyscses`` can calculate interface (grain boundary) resistivities [@HwangEtAl_JElectroceram1999] and activation energies. 

## Grain boundary resistivities

Perpendicular grain boundary resistivities, $r_\mathrm{GB}^\perp$, are calculated by modelling the system as a set of resistors in series:
$$
r_\mathrm{GB}^\perp = \frac{c_{i,\infty}\mu_{i,\infty}z_i}{\sum_{i} c_{i,x}\mu_{i,x}z_i}.
$$
Parallel grain boundary resistivities, $r_\mathrm{GB}^\parallel$, are calculated by modelling the system as a set of resistors in parallel:
$$
r_\mathrm{GB}^\parallel = \frac{1}{\frac{\sum_{i} c_{i,x}\mu{i,x}z_i}{c_{i,\infty}mu_{i,\infty}z_i}}
$$
 
## Activation energies
``pyscses`` can be used to calculate the grain boundary contribution to a single defect activation energy by running the Poisson-Boltzmann and resistivity ratio calculations at a range of different temperatures and applying an Arrhenius equation. 

$$
\ln{r_\mathrm{GB}} = \frac{-E_\mathrm{act}}{kT}
$$

# Approximations and limitations
- ``pyscses`` implements a modification of the Poisson-Boltzmann equation, which assumes that defects are non-interacting except for point-charge electrostatics and site exclusion.

- The resistivitity and activation energy calculations implemented in ``pyscses`` assume fixed mobilities for defects.

# Acknowledgements

G. L. W. acknowledges support from the EPSRC Doctoral Training Partnership at the University of Bath (EPSRC grant code?).
B. J. M. acknowledges support from the Royal Society (Grant No. UF130329).

# References
