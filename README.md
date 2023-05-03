# cfd-ense
OpenFOAM solver for investigating the extended Navier-Stokes equations (ENSE) as proposed by Durst, Sambasivam, Chakraborty, Dongari, Brenner et al. [^1][^2][^3][^4][^5][^6][^7].

## Background

Flows of ideal gases characterized by a certain range of Knudsennumbers ($Kn > \sim 0.01$) can not be accurately predicted by the classical Navier-Stokes equations (CNSE). This expresses itself in a deviation of the calculated from the experimentally determined mass flow rates, where the latter yields higher results. First attempts to take the missing physics into account for by extending the CNSE by new terms seem to be credited to Brenner [^5], which seem to have inspired others to take similar approaches. One of them was found to be more easily alignable with fundamental physical principles. Durst, Sambasivam et al. [^1][^2][^3][^4] initially named their newly derived model the extended Navier-Stokes equations. However, in the end these approaches are very similar and as Dongari et al. [^6] had proven in case of isothermal flows Brenners [^5] and Durst et al. [^1] models are identical. 


## Core principle
Considering the following stationary, isothermal NSE:

**Continuity equation:**

$$\nabla \cdot \left( \rho \bm{u} \right) = 0$$

**Momentum conservation equations:**

$$\nabla \cdot \left(\rho \bm{u}^C \bm{u}^C\right)=- \nabla p - \nabla \cdot \bm{\underline{\tau}}$$

* The continuity equation needs to be fulfilled for a total mass flux $\rho \bm{u}$
* The following types of mass fluxes with different root causes are identified:
    * Diffusion $\rho \bm{u}^D$: Considering the kinetic theory of gases, the statistical motion of particles (internal energy) along with the presence of density (and/or temperature) gradients will result in a macroscopic mass flux (not considered in the classical equations)
    * Convection $\rho \bm{u}^C$: Fluid gets driven/accelerated by external forces/mechanical work in the macroscopic environment
* Convective/transient acceleration terms are formulated in terms of the convective velocity $\bm{u}^C$ because macroscopic momentum is not a property of the diffusive velocity $\bm{u}^D$ 
* Transport of momentum by mass diffusion gets considered when deriving the shear stress tensor $\bm{\underline{\tau}}$
* The superposition of $\bm{u}^C$ and $\bm{u}^D$ is assumed : $\bm{u} = \bm{u}^C + \bm{u}^D$

## Validity
The ENSE has been validated for a limited range of use cases so far, e.g. isothermal, laminar flows with Knudsen numbers of up to ~1. The superposition assumption should be questioned & challenged. First attempts to solve turbulent flows and shockwaves seem to be promising [^4]. By including models for an effective mean free path it seems to be able to predict Knudsen diffusion in the free molecular regime $Kn > \sim 10$ [^7]. However, a comparision with over 30 experiments by Shen [^8] pointed out that there still seems to be something missing as the calculated massflow has some kind of offset starting in the slip flow regime. The quantity of this offset varies from experiment to experiment. As for now (05/2023) it is assumed that this is attributed to reflection of particles at the relatively smooth walls of microchannels as proposed by Maxwell [^9], hence the no-slip assumption for $\bm{u}^C$ might be insufficient here.


## rhoSimpleEnseFoam

Introduced in [^10].

### Features

- Pressure-based, steady-state and sequential solver for of the ENSE
- If the isothermal switch in the SIMPLE subdictionary is set to "yes", solving the energy equation will be omitted
- Naming of velocity fields is consistend with paper [^10] where this work is presented
- Based on rhoSimpleFoam in [OpenFOAM v2006 / v2212 ](https://www.openfoam.com)

### Compilation and usage

A guide for installation of OpenFOAM v2212 on Linux can be found [here](https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled/debian)

Installation on Windows seems to be most advantageously achieved when using WSL along with Ubuntu.

After sourcing the OpenFOAM v2006 / v2212 environment in a shell, clone the repository and compile the solver:
```bash
git clone https://github.com/imfd-stroemungsmechanik/cfd-ense.git
cd cfd-ense
wmake ./applications/rhoSimpleEnseFoam
```
 
To run the microchannel test case, create the mesh and start the solver:
```bash
cd ./run/microchannel
blockMesh
rhoSimpleEnseFoam
```
The most accurate reproduction of the residuals/underrelaxation-factors presented in the paper [^10] can be achieved by running it on one computational core and using version v2006. Background for solvers and testcases is given in the subfolders readmes.
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
