# cfd-ense

## rhoSimpleEnseFoam

### Features

- Pressure-based, steady-state and sequential solver for the compressible flow of gases
- Solves the extended Navier-Stokes equations (ENSE) by [Durst et al.](https://doi.org/10.1007/s40032-020-00586-3) to include mass diffusion for the simulation of flows with high Knudsen numbers
- Maxwell-slip and Smoluchowski-jump boundary conditions are adapted from rhoCentralFoam and could be used in addition, although the former comes without curvature terms so far
- If the isothermal switch in the SIMPLE subdictionary is set to "yes", solving the energy equation will be omitted. 
- Solves for the convective velocity U. Diffusive velocity Ud and total velocity Ut are updated explicitly in each outer iteration
- Based on rhoSimpleFoam in [OpenFOAM v2006](https://www.openfoam.com)
- Developed by Johannes Schwarz, minor updates by Martin Heinrich

### Compilation and usage on Linux

A guide for compilation and installation of OpenFOAM v2006 on Linux can be found [here](https://develop.openfoam.com/Development/openfoam/-/blob/OpenFOAM-v2006/doc/Build.md)

After sourcing the OpenFOAM v2006 environment in a shell, clone the repository and compile the solver:
```bash
git clone https://github.com/imfd-stroemungsmechanik/cfd-ense.git
cd cfd-ense
wmake ./applications/rhoSimpleEnseFoam
```
 
To run the microchannel test case on one CPU, create the mesh and start the solver:
```bash
cd ./run/microchannel
blockMesh
rhoSimpleEnseFoam
```
