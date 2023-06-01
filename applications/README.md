# OpenFOAM solver

## rhoSimpleEnseFoam

Introduced in [^10]
Maxwell-slip and Smoluchowski-jump boundary conditions are adapted from rhoCentralFoam and could be used in addition, although the former comes without curvature terms so far

### Features

- Pressure-based, steady-state and sequential solver for of the (S)ENSE
- If the isothermal switch in the SIMPLE subdictionary is set to "yes", solving the energy equation will be omitted
- Naming of velocity fields is consistend with paper where this work is presented
- Based on rhoSimpleFoam in [OpenFOAM v2006 / v2212](https://www.openfoam.com)

### Compilation and usage on Linux

A guide for installation of OpenFOAM v2212 on Linux can be found [here](https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled/debian)

Installation on Windows seems to be most advantageously achieved when using WSL along with Ubuntu.

After sourcing the OpenFOAM v2006 / v2212 environment in a shell, clone the repository and compile the solver:
```bash
git clone https://github.com/imfd-stroemungsmechanik/cfd-ense.git
cd cfd-ense
wmake ./applications/rhoSimpleEnseFoam
```

### blockCoupledSenseFoam

#### Features

- Pressure-based, steady-state and block-coupled solver for the (S)ENSE
- Solves a block-coupled compressible U-p-e system. Achieves implicit coupling between the flow variables - e.g. the temperature dependent diffusion term in the continuity equation gets also implicitly discretized
- Utilizes most basic coupling strategies
- Succesfully tested for both BFS verification cases and a validation case (will be uploaded soon). Microchannel not yet available
- Based on [foam-extend 5.0](https://sourceforge.net/projects/foam-extend/)

\
\
(c) Johannes Schwarz, Martin Heinrich 2023

[^10]: Johannes Schwarz, Kristjan Axelsson, Daniel Anheuer, Martin Richter, Johanna Adam, Martin Heinrich, Rüdiger Schwarze. An OpenFOAM solver for the extended Navier–Stokes equations. SoftwareX. Volume 22. 2023. 101378. ISSN 2352-7110. https://doi.org/10.1016/j.softx.2023.101378.