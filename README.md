# cfd-ense
This repository currently contains a sequential (rhoSimpleEnseFoam) and a block-coupled (blockCoupledSenseFoam) steady-state OpenFOAM / foam-extend solver for Sambasivam's [^4] set of extended Navier-Stokes equations (SENSE) as well as a few testcases.

## Background

Flows of ideal gases characterized by a certain range of Knudsennumbers ($Kn > 0.01$) cannot be accurately predicted by the classical Navier-Stokes equations (CNSE). This expresses itself in a deviation of the calculated from the experimentally determined mass flow rates, with higher values being observed for the latter. First approaches in the recent decades to take the missing physics into account for by extending the CNSE by new terms were pursued by Brenner [^5], which were then followed by those of several other scientists [^13]. One group of them derived a model that was found to be more easily alignable with fundamental physical principles.  Durst, Gomes and Sambasivam [^1] termed their newly derived model as "extended Navier-Stokes equations" ("ENSE"). While the terminologies and initial approaches appear to be different, from a mathematical point of view Brenner's and "ENSE" models are very similar: At several places an additional diffusive flux term gets introduced into the Navier-Stokes system. In both cases this flux contains a density gradient, which becomes significant under certain conditions, e.g. large Knudsen numbers, is otherwise commonly neglectable and not considered in the classical equations. However, the Brenner-flux appears to be more general and should be applicable to liquids aswell [^13]. The "ENSE"-flux on the other hand is more specifically related to the kinetic theory of (ideal) gases and also contains an additional temperature-gradient-term, also known as "Soret-diffusion". Therefore, the main focus for this repository lies on the "ENSE" model now or - to be more precise - on the variant derived by Sambasivam in his Ph.D. thesis (SENSE).

## Validity
The SENSE has already been validated for a few use-cases so far, e.g. certain laminar flows with Knudsen numbers of up to $O(1)$. First simulations of turbulent flows and shockwaves are promising [^4]. By including models for an effective mean free path it seems to be able to predict the transitional and free molecular flow regime ($Kn > 10$) [^7]. Using the $G_m$ - plot and the exact analytical model/method as proposed by Dongari et al. [^7] the curvatures of the calculated graph were found to fit indeed very well with those of the experimental curve of Ewart et al. [^12] for mean Knudsen numbers of up to 80. However, there was overall an obvious offset visible. This offset is also clearly observable in almost every single one of the over 30 comparisons with experiments carried out by Shen [^8].

## Solvers

A very brief introduction into the SENSE model and OpenFOAM programming with respect to the rhoSimpleEnseFoam solver is given in [^10]. The derivation of the SENSE model itself was done by Sambasivam in his Ph.D. thesis [^4].

### rhoSimpleEnseFoam

#### Features

- Pressure-based, steady-state and sequential solver for of the SENSE
- If the isothermal switch in the SIMPLE subdictionary is set to "yes", solving the energy equation will be omitted
- Naming of velocity fields is consistend with paper [^10] where this work is presented
- Based on rhoSimpleFoam in [OpenFOAM v2006 / v2212 ](https://www.openfoam.com)

### blockCoupledSenseFoam

#### Features

- Pressure-based, steady-state and block-coupled solver for the SENSE
- Solves a block-coupled compressible U-p-e system. Achieves implicit coupling between the flow variables - e.g. the temperature dependent diffusion term in the continuity equation gets also implicitly discretized
- Based on [foam-extend 5.0](https://sourceforge.net/projects/foam-extend/)

More background for solvers and testcases is given in the subfolders readmes.

\
\
(c) Johannes Schwarz, 2023

[^1]: Durst F, Gomes J, Sambasivam R. Thermofluiddynamics: Do We Solve the
Right Kind of Equations? In: Turbulence heat and mass transfer 5. Proceed-
ings of the international symposium on turbulence heat and mass transfer.
2006, http://dx.doi.org/10.1615/ICHMT.2006.TurbulHeatMassTransf.10

[^2]: Chakraborty Suman, Durst Franz. Derivations of extended Navier-Stokes
equations from upscaled molecular transport considerations for compress-
ible ideal gas flows: Towards extended constitutive forms. Phys Fluids
2007;19. http://dx.doi.org/10.1063/1.2759531

[^3]: Durst F, Filimonov D, Sambasivam R. Treatments of micro-channel flows
revisited: Continuum versus rarified gas considerations. J Inst Eng (India):
Ser C 2020;101. http://dx.doi.org/10.1007/s40032-020-00586-3

[^4]: Sambasivam Rajamani. Extended Navier-Stokes Equations:
Derivations and Applications to Fluid Flow Problems (Ph.D. thesis), Friedrich-
Alexander-Universität Erlangen-Nürnberg; 2012

[^5]: Brenner Howard. Navier–Stokes revisited. Physica A 2004;349:60–132. http://dx.doi.org/10.1016/j.physa.2004.10.034.

[^6]: Dongari, Nishanth & Sharma IITK, Ashutosh & Durst, F.. (2009). Pressure-driven diffusive gas flows in micro-channels: From the Knudsen to the continuum regimes. Microfluidics and Nanofluidics. 6. 679-692. https://doi.org/10.1007/S10404-008-0344-Y.

[^7]: Dongari Nishanth, Durst F, Chakraborty Suman. Predicting microscale
gas flows and rarefaction effects through extended Navier–Stokes–Fourier
equations from phoretic transport considerations. Microfluid Nanofluid
2010;9:831–46. http://dx.doi.org/10.1007/s10404-010-0604-5.

[^8]: [Shen Di. Viscous Compressible Flow Through a Micro-Conduit: Slip-Like Flow Rate with No-Slip Boundary Condition (Ph.D. thesis), Arizona State University; 2019](https://hdl.handle.net/2286/R.I.54955)

[^9]: Maxwell J. VII. On stresses in rarified gases arising from inequalities of temperature.
Philosophical Transactions of the Royal Society of London
1879;170:231-256
https://doi.org/10.1098/rstl.1879.0067

[^10]: Johannes Schwarz, Kristjan Axelsson, Daniel Anheuer, Martin Richter, Johanna Adam, Martin Heinrich, Rüdiger Schwarze. An OpenFOAM solver for the extended Navier–Stokes equations. SoftwareX. Volume 22. 2023. 101378. ISSN 2352-7110. https://doi.org/10.1016/j.softx.2023.101378.

[^11]: S Kokou Dadzie, Howard Brenner. Predicting enhanced mass flow rates in gas microchannels using nonkinetic models. Physical Review E. 2012. http://dx.doi.org/10.1103/PhysRevE.86.036318

[^12]: Timothée Ewart, Pierre Perrier, Irina A. Graur, J. Gilbert Méolans. Mass flow rate measurements in a microchannel, from hydrodynamic to near free molecular regimes. Journal of Fluid Mechanics. 2007, http://dx.doi.org/10.1017/S0022112007006374

[^13]: Howard Brenner. Beyond Navier-Stokes. International Journal of Engineering Science. 2012. http://dx.doi.org/10.1016/j.ijengsci.2012.01.006