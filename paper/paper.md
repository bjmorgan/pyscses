---
title: 'pyscses: a PYthon Space-Charge Site-Explicit Solver'
tags:
 - Python
 - poisson-boltzmann
 - space-charge
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
date: 8 January 2019
bibliography: paper.bib
---

# Summary

``pyscses`` is an open-source Python package for modelling ionic space-charge regions in crystalline solids. Ionic space-charges are microscopic regions with a local excess or deficiency of charged point defects, giving these regions a net charge [@Frenkel_KinTheLiq1946]. Space-charge regions can form at crystallographic interfaces, such as grain boundaries or heterointerfaces between two different materials. The formation of space-charges is due to segregation of charged defects to, or from, these interfaces, and the associated redistribution of defects in adjacent regions of the crystal, due to defect&ndash;defect electrostatic interactions [@Maier_BerBunsenges1984; @Maier_JPhysChemSol1985; @Maier_SolStatIonics2003; @ChiangEtAl_ApplPhysLett1996; @KimAndMaier_JElectrochemSoc2002]. 

The segregation of defects to, or from, space-charge regions can produce local defect concentrations that strongly deviate from the average &ldquo;bulk&rdquo; values in a material [@FleigAndMaier_SolidStateIonics1996]. This can have significant consequences for key material properties. For example, in solid electrolytes, where ionic defects act as mobile charge carriers, space-charge formation at grain boundaries is associated with large changes in ionic conductivity, both parallel and perpendicular to the grain boundaries, due to effects such as enhanced defect concentrations in the grain boundary core [@JamnikAndMaier_SolidStateIonics1999], or reduced defect concentrations in adjacent &ldquo;depletion&rdquo; space-charge regions. In the case of solid electrolytes, understanding space-charge formation at grain boundaries and interfaces is a key challenge in developing a theoretical description of the role of crystalline microstructure on macroscopic ion transport (ionic conductivities) [@KimEtAl_PhysChemChemPhys2003].

``pyscses`` considers simple one-dimensional models of crystallographic interfaces, and calculates equilibrium defect distributions by solving a modified Poisson-Boltzmann equation [@Maier_ProgSolStatChem1995; @DeSouzaEtAl_SolidStateIonics2011]. The driving force for point-defect segregation to or from the interface is described by sets of defect segregation energies, $\Delta E_\mathrm{seg}^{i,x}$, defined as
$$
\Delta E_\mathrm{seg}^{i,x} = E_\mathrm{f}^{i,x} - E_\mathrm{f}^{i, \infty}
$$
i.e. the difference in defect formation energy for defect species $i$ at site $x$ compared to a reference site in the &ldquo;bulk&rdquo; of the crystal.

# The numerical model
The default model implemented in ``pyscses`` calculates the equilibrium distribution of point-charge defect species on a 1D grid [@HelgeeEtAl_FuelCells2012; @PolfusEtAl_SolStatIonics2016; @DeSouzaEtAl_SolidStateIonics2011; @NymanEtAl_ApplPhysLet2012]. The electrochemical potential, $\mu^o_{i,x}$, of defect species $i$ at site $x$, with site occupancy $c_{i,x}$ is 
$$
\mu^o_{i,x} + RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}} \right) + z_i F \Phi_x.
$$
The $RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}}\right)$ term gives the chemical potential for a non-interacting lattice-gas with site exclusion (maximum site occupancy of 1). $\Phi_x$ is the electrosatic potential at site $x$, and $F$ is the Faraday constant.

Finding the equilibrium distribution of defects described by this electrochemical potential function is equivalent to solving a modified Poisson-Boltzmann equation, which can be derived by requiring the electrochemical potentials at each site to be equal to a set of reference bulk electrochemical potentials [@Maier_IonicsTextbook2005].
$$
\mu^o_{i,x} + RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}} \right) + z_i F \Phi_x = \mu^o_{i,\infty} + RT\ln \left(\frac{c_{i,\infty}}{1- c_{i,\infty}}\right) + z_i F \Phi_{\infty}.
$$
This equilibrium condition can be rearranged to give the defect concentration (site occupancy) at site $x$:
$$
c_i = \frac{c_\infty \exp\left(\frac{-z_i\Phi_x + \mu_i}{kT}\right)}{1+ c_{\infty} \left(\exp\left( \frac{-z_i\Phi_x + \mu_i}{kT} \right) -1 \right) } . 
$$
The charge density at each site is proportional to a sum over defect concentrations given by
$$
\rho = \frac{1}{V_\mathrm{site}}\sum_i c_i z_i F.
$$
The electrostatic potential at all sites can then be calculated using Poisson's equation,
$$
\nabla^2 \Phi = \frac{ -\rho } { \epsilon \epsilon_0 }.
$$
The equilibrium defect distribution corresponds to the self-consistent solution to these three previous equations.
``pyscses`` solves this set of equations as a second order partial differential equation using a second order finite difference approximation for each site.

Within the general framework of solving this modified 1D Poisson-Boltzmann equation,``pyscses`` can follow a number of different modelling approaches:

- Continuum (regular grid) and site-explicit (irregular grid) models.
- Periodic and Dirichlet boundary conditions.
- &ldquo;Mott-Schottky&rdquo; and &ldquo;Gouy-Champman&rdquo; condtions. These are implemented by setting the mobilities of different defect species. In the case of Mott-Schottky conditions, all but one defect species have a mobility of zero.
- Inclusion of &ldquo;lattice&rdquo;site charges to account for non-defective species in the crystal structure.

# Typical workflow

The necessary input to model space-charge formation at a grain boundary is a set of defect site positions and segregation energies, projected onto a one-dimensional grid (see Figure 1). For calculations using a &ldquo;continuum&rdquo; (regular) grid, the defect segregation energies and atomic positions are interpolated onto a regular grid.

![(Top) An example crystal structure for a grain boundary in CeO<sub>2</sub>. The $x$ coordinate of each potential defect site (orange spheres) is used to construct a one-dimensional &ldquo;site-explicit&rdquo; grid. Defect segregation energies calculated using e.g. atomistic modelling methods are used to assign segregation energies to every grid point (bottom).](Figures/seg_energies_joss.pdf)

`pyscses` uses these input data to solve the self-consistent modified Poisson-Boltzmann equation. The calculated outputs include the equilibrium electrostatic potential, charge density, and defect mole fractions (site occupancies) across the space charge region (Figure 2).

![Example outputs (electrostatic potentials, charge densities, and site occupancies) for a grain boundary in Gd-doped CeO<sub>2</sub>. The left and right pairs of panels show equivalent results calculated using continuum and site-explicit models.](Figures/continuum_vs_se_joss_MS.pdf)

## Calculated properties

In addition to allowing calculation of defect, charge, and potential across an interface, ``pyscses`` can calculate interface (grain boundary) resistivities [@HwangEtAl_JElectroceram1999] and activation energies. 

### Grain boundary resistivities

Perpendicular grain boundary resistivities, $r_\mathrm{GB}^\perp$, are calculated by modelling the system as a set of resistors in series:
$$
r_\mathrm{GB}^\perp = \frac{c_{i,\infty}\mu_{i,\infty}z_i}{\sum_{i} c_{i,x}\mu_{i,x}z_i}.
$$
Parallel grain boundary resistivities, $r_\mathrm{GB}^\parallel$, are calculated by modelling the system as a set of resistors in parallel:
$$
r_\mathrm{GB}^\parallel = \frac{1}{\frac{\sum_{i} c_{i,x}\mu{i,x}z_i}{c_{i,\infty}mu_{i,\infty}z_i}}
$$
 
### Activation energies
``pyscses`` can be used to calculate the grain boundary contribution to an ionic conductivity activation energy. This procedure consists of solving a space-charge model at a series of temperatures, and numerically differentiating the logrithmic grain-boundary resistance [@Kim_PhysChemChemPhys2016]:

$$
E_\mathrm{act}=\frac{ \mathrm{d} \ln r_\mathrm{GB} } { \mathrm{d} \frac{1}{kT}}.
$$

# Approximations and limitations
- ``pyscses`` implements a modified Poisson-Boltzmann equation, which assumes that defects are non-interacting except for point-charge electrostatics and site exclusion.
- The resistivitity and activation energy calculations implemented in ``pyscses`` assume mobilities of each defect species are independent of the local crystal structure.

# Acknowledgements

G. L. W. acknowledges support from the EPSRC Doctoral Training Partnership at the University of Bath (EPSRC grant code?).
B. J. M. acknowledges support from the Royal Society (Grant No. UF130329).

# References
