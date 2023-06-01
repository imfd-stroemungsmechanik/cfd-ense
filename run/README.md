# Testcases

## Verification

Numerical/analytical ENSE results available for verification purposes: [^4]

* Microchannel - Page 50 ff.
* BackwardStep_Kn0.018 - Page 101 ff.
* BackwardStep_Kn1 - Page 134 ff.

## Validation

Experimental / DSMC ... data available for validating the model itself


[^4]: Sambasivam Rajamani. Extended Navier-Stokes Equations:
Derivations and Applications to Fluid Flow Problems (Ph.D. thesis), Friedrich-
Alexander-Universität Erlangen-Nürnberg; 2012

### Instructions

To run e.g. the microchannel test case, change into the corresponding directory and execute the Allrun script. Otherwise manually create the mesh and start the solver. E.g. for rhoSimpleEnseFoam:
```bash
cd ./verification/ofseq/microchannel
blockMesh
rhoSimpleEnseFoam
```
Optionally, one can also decompose the domain after mesh creation and execute the solver on multiple CPUs, e.g. on 8 as follows (The decomposeParDict in the systems folder needs to be configured accordingly):

```bash
decomposePar
mpirun -np 8 rhoSimpleEnseFoam -parallel
```

The most accurate reproduction of the residuals/underrelaxation-factors presented in the paper can be achieved by running it on one CPU and using version v2006. There are several post processing scripts available in OpenFOAM, e.g. an evaluation of the mass flow rate on certain boundaries, see for example the controlDict in the microchannel testcase. Otherwise ParaView is a very powerful tool that is commonly used for visualization/post processing of OpenFOAM output.
