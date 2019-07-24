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
`pyscses` is a Python package for modelling ionic space-charges in solid electrolytes. Its primary use is to calculate equilibrium distributions of point-charge atomic defects within one-dimensional &ldquo;Poisson-Boltzmann&rdquo;-like mean-field models. These calculations take as inputs a set of defect site positions, within a specific crystal structure, and the associated defect segregation energies. `pyscses` can also be used to calculate ionic transport properties (space-charge resistivities and activation energies) for these equilibrium defect distributions.

One approach to modelling space-charge formation in solid electrolytes is to consider defects as ideal point-charges embedded in a continuum dielectric, and to calculate equilibrium defect distributions by solving mean-field &ldquo;Poisson-Boltzmann&rdquo;-like equations [@Franceschetti_SolStatIonics1981; @GuoAndWaser_ProgMaterSci2006; @NymanEtAl_ApplPhysLett2012; @LindmanEtAl_SolStatIonics2013a; @PolfusEtAl_SolStatIonics2016; @HelgeeEtAl_FuelCells2013]. While numerical solutions to the 1D Poisson-Boltzmann equation are relatively simple to implement, published results are typically obtained using private closed-source codes, making it difficult to reproduce results or to test the effect of different approximations included in specific models. ``pyscses`` provides an open-source Python package for modelling space-charge formation in solid electrolytes, within a 1D Poisson-Boltzmann-like formalism.We are currently using ``pyscses`` in our own research into space-charge formation in solid electrolytes for fuel cells and lithium-ion batteries, and hope that this open-source resource will support reproducible research practices in future studies in this area [@SandveEtAl_PLoSComputBiol2013].

# Scientific Context
Crystalline solids consist of periodic arrangements of atoms that vibrate about fixed lattice sites. In most solids, atoms rarely move between lattice sites, and long ranged atomic diffusion is slow. Solid electrolytes are notable because they contain atoms that can rapidly move between lattice sites, giving high diffusion coefficients. This unusual property makes them useful in solid-state electrochemical devices such as fuel cells and solid-state batteries [@MahatoEtAl_ProgMaterSci2015; @ZhengEtAl_JPowerSources2018]. Many practical solid electrolytes are polycrystalline: they contain multiple crystalline domains, with varied orientations, that meet at grain boundaries. 

All (poly)crystalline materials contain defects&mdash;these are structural imperfections in the crystal that break local crystal symmetries, such as vacancies or interstitials. Thermodynamic arguments predict that defects will spontaneously segregate towards, or away from, interfaces and grain boundaries. In ionic crystals, such as solid electrolytes, these defects carry electric charge, and their segregation produces a build up of charge at grain boundaries. This, in turn, causes a redistribution of charged defects in adjacent crystalline regions, to form locally charged &ldquo;space-charge&rdquo; regions [@Maier_BerBunsenges1984; @Maier_JPhysChemSol1985; @Maier_SolStatIonics2003; @ChiangEtAl_ApplPhysLett1996; @KimAndMaier_JElectrochemSoc2002]. In space-charge regions, defect concentrations may differ significantly from average &ldquo;bulk&rdquo; values in a single crystal [@FleigAndMaier_SolidStateIonics1996]. Ionic conductivities depend on local defect concentrations, and space-charge regions with reduced or enhanced defect numbers can therefore significantly affect the ionic conductivities of solid electrolytes [@JamnikAndMaier_SolidStateIonics1999]. Understanding the defect segregation and space-charge formation in solid electrolytes is therefore key to optimising the properties of these materials for use in electrochemical devices.

# Numerical Model
``pyscses`` considers simple one-dimensional models of crystallographic interfaces, and calculates equilibrium defect distributions by solving a modified Poisson-Boltzmann equation [@Maier_ProgSolStatChem1995; @DeSouzaEtAl_SolidStateIonics2011], which can be derived by considering the condition that at equilibrium the electrochemical potential for a given defect species is constant [@Maier_IonicsTextbook2005]:
$$
\mu^o_{i,x} + RT\ln \left( \frac{c_{i,x}}{1-c_{i,x}} \right) + z_i F \Phi_x = \mu^o_{i,\infty} + RT\ln \left(\frac{c_{i,\infty}}{1- c_{i,\infty}}\right) + z_i F \Phi_{\infty}.
$$
The thermodynamic driving force for point-defect segregation to or from the interface is described by defect segregation energies, $\left\{\Delta E_\mathrm{seg}^{i,x}\right\}$
$$
\Delta E_\mathrm{seg}^{i,x} = E_\mathrm{f}^{i,x} - E_\mathrm{f}^{i, \infty}
$$
i.e., $\Delta E_\mathrm{seg}^{i,x}$ is the difference in defect formation energy for defect species $i$ at site $x$ compared to a reference site in the &ldquo;bulk&rdquo; of the crystal.

Within the general framework of solving this modified 1D Poisson-Boltzmann equation, ``pyscses`` implements a range of numerical models:

- Continuum (regular grid) and site-explicit (irregular grid) models.
- Periodic and Dirichlet boundary conditions.
- &ldquo;Mott-Schottky&rdquo; and &ldquo;Gouy-Chapman&rdquo; conditions. These are implemented by setting the mobilities of different defect species. In the case of Mott-Schottky conditions, all but one defect species have mobilities of zero.
- Inclusion of &ldquo;lattice-site&rdquo; charges to account for non-defective species in the crystal structure.

Properties that can be calculated include:

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
## Approximations and Limitations
- The modified Poisson-Boltzmann model implemented in ``pyscses`` assumes that defects only interact via point-charge electrostatics and volume exclusion. 
- The resistivitity and activation energy calculations implemented in ``pyscses`` assume that defect mobilities are independent of the local crystal structure.

# Acknowledgements

B. J. M. acknowledges support from the Royal Society (Grant No. UF130329), and from the Faraday Institution (FIRG003).

# References
