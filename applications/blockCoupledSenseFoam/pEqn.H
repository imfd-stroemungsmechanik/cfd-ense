// Pressure equation

fvScalarMatrix pEqn
(
    fvm::div(phiPsi, p)
);

// Calculate div U coupling.
BlockLduSystem<vector, scalar> UInp(fvm::UDiv(fvc::interpolate(rho), U));
UpEEqn.insertBlockCoupling(3, 0, UInp, false);
pEqn -= fvc::div(phi);


// ENSE:

fvScalarMatrix epEqn
(
    fvm::laplacian(nmuByp, p)
);

pEqn += epEqn;

// Temperature dependent diffusion - e coupling
fvScalarMatrix eInp(constrainmdT(fvm::laplacian(muBy2e, thermo.e())));
UpEEqn.insertEquationCoupling(3, 4, eInp);


pEqn.setReference(pRefCell, pRefValue);
UpEEqn.insertEquation(3, pEqn);
