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
Solid electrolytes are often crystalline materials with high ionic conductivities, often used in solid state electrochemical devices. Most crystalline materials are characterised by defects, imperfections in the crystal that violate the perfect symmetries. Many solid electrolytes are polycrystalline with grain boundaries separating individual locally ordered regions with different orientations. In ionic crystals, thermodynamic arguments show that charged point defects will segregate to, or away from, interfaces and grain boundaries, and that there will be an associated redistribution of point defects in the regions adjacent. This redistribution leads to the regions having a net charge - known as space charge regions [@Maier_BerBunsenges1984; @Maier_JPhysChemSol1985; @Maier_SolStatIonics2003; @ChiangEtAl_ApplPhysLett1996; @KimAndMaier_JElectrochemSoc2002]. The segregation of defects can produce local defect concentrations that strongly deviate from the average &ldquo;bulk&rdquo; values in a material [@FleigAndMaier_SolidStateIonics1996], which can result in large variations in ionic conductivities[@JamnikAndMaier_SolidStateIonics1999]. It is therefore important to have a good understanding of defect behaviour and space charge formation to improve ionic conductivities in solid electrolytes for use in electrochemical devices. ``pyscses`` is an open-source Python package for modelling ionic space-charge regions in crystalline solids. While the theory of space charge formation is well documented, there is currently no resources available for calculating space charge properties in solid electrolytes using explicitly defined defect positions and defect segregation energies. 

# The numerical model
``pyscses`` considers simple one-dimensional models of crystallographic interfaces, and calculates equilibrium defect distributions by solving a modified Poisson-Boltzmann equation [@Maier_ProgSolStatChem1995; @DeSouzaEtAl_SolidStateIonics2011]. The driving force for point-defect segregation to or from the interface is described by sets of defect segregation energies, $\Delta E_\mathrm{seg}^{i,x}$, defined as
$$
\Delta E_\mathrm{seg}^{i,x} = E_\mathrm{f}^{i,x} - E_\mathrm{f}^{i, \infty}
$$
i.e. the difference in defect formation energy for defect species $i$ at site $x$ compared to a reference site in the &ldquo;bulk&rdquo; of the crystal.

The default model implemented in ``pyscses`` calculates the equilibrium distribution of point-charge defect species on a 1D grid [@HelgeeEtAl_FuelCells2012; @PolfusEtAl_SolStatIonics2016; @DeSouzaEtAl_SolidStateIonics2011; @NymanEtAl_ApplPhysLet2012]. Finding the equilibrium distribution of defects is described by solving a modified Poisson-Boltzmann equation, which can be derived by requiring the electrochemical potentials at each site to be equal to a set of reference bulk electrochemical potentials [@Maier_IonicsTextbook2005].
$$
\mu^o_{i,x} + RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}} \right) + z_i F \Phi_x = \mu^o_{i,\infty} + RT\ln \left(\frac{c_{i,\infty}}{1- c_{i,\infty}}\right) + z_i F \Phi_{\infty}.
$$

Within the general framework of solving this modified 1D Poisson-Boltzmann equation, ``pyscses`` can follow a number of different modelling approaches:

- Continuum (regular grid) and site-explicit (irregular grid) models.
- Periodic and Dirichlet boundary conditions.
- &ldquo;Mott-Schottky&rdquo; and &ldquo;Gouy-Champman&rdquo; condtions. These are implemented by setting the mobilities of different defect species. In the case of Mott-Schottky conditions, all but one defect species have a mobility of zero.
- Inclusion of &ldquo;lattice&rdquo;site charges to account for non-defective species in the crystal structure.

To calculate a number of space charge properties:

- Defect mole fractions.
- Charge density.
- Electrostatic potential.
- Parallel and perpendicular grain boundary resistivities [@HwangEtAl_JElectroceram1999].
- Grain boundary activation energies [@Kim_PhysChemChemPhys2016].

<!---# Typical workflow

The necessary input to model space-charge formation at a grain boundary is a set of defect site positions and segregation energies, projected onto a one-dimensional grid (see Figure 1). For calculations using a &ldquo;continuum&rdquo; (regular) grid, the defect segregation energies and atomic positions are interpolated onto a regular grid.

![(Top) An example crystal structure for a grain boundary in CeO<sub>2</sub>. The $x$ coordinate of each potential defect site (orange spheres) is used to construct a one-dimensional &ldquo;site-explicit&rdquo; grid. Defect segregation energies calculated using e.g. atomistic modelling methods are used to assign segregation energies to every grid point (bottom).](Figures/seg_energies_joss.pdf)

`pyscses` uses these input data to solve the self-consistent modified Poisson-Boltzmann equation. The calculated outputs include the equilibrium electrostatic potential, charge density, and defect mole fractions (site occupancies) across the space charge region (Figure 2).

![Example outputs (electrostatic potentials, charge densities, and site occupancies) for a grain boundary in Gd-doped CeO<sub>2</sub>. The left and right pairs of panels show equivalent results calculated using continuum and site-explicit models.](Figures/continuum_vs_se_joss_MS.pdf)
-->
# Approximations and limitations
- ``pyscses`` implements a modified Poisson-Boltzmann equation, which assumes that defects are non-interacting except for point-charge electrostatics and site exclusion.
- The resistivitity and activation energy calculations implemented in ``pyscses`` assume mobilities of each defect species are independent of the local crystal structure.

# Acknowledgements

B. J. M. acknowledges support from the Royal Society (Grant No. UF130329).

# References
