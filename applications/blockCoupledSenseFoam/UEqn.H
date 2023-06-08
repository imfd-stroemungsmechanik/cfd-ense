// Momentum equation

fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + turbulence->divDevRhoReff()
  //- fvm::laplacian(thermo.mu(), U)
  //- fvc::div(thermo.mu()*dev2(T(fvc::grad(U))))
);

// Calculate grad p coupling matrix. VV, 9/June/2014
BlockLduSystem<vector, vector> pInU(fvm::grad(p));

surfaceScalarField U0f = fvc::interpolate(U.component(0));
surfaceScalarField U1f = fvc::interpolate(U.component(1));
surfaceScalarField U2f = fvc::interpolate(U.component(2));

// div(phi, U) - p coupling
fvScalarMatrix pInU0(fvm::div(phiPsi*U0f, p));
fvScalarMatrix pInU1(fvm::div(phiPsi*U1f, p));
fvScalarMatrix pInU2(fvm::div(phiPsi*U2f, p));
UEqn -= fvc::div(phi, U);


// ENSE terms:

UEqn +=
(
    fvm::div(phimd, U)
  + fvc::div(dev2(md*U))
);

// div(phimd, U) - p and e coupling
pInU0 += fvm::laplacian(nmuByp*U0f, p);
pInU1 += fvm::laplacian(nmuByp*U1f, p);
pInU2 += fvm::laplacian(nmuByp*U2f, p);
fvScalarMatrix EInU0(constrainmdT(fvm::laplacian(muBy2e*U0f, thermo.e())));
fvScalarMatrix EInU1(constrainmdT(fvm::laplacian(muBy2e*U1f, thermo.e())));
fvScalarMatrix EInU2(constrainmdT(fvm::laplacian(muBy2e*U2f, thermo.e())));
UpEEqn.insertEquationCoupling(0, 4, EInU0);
UpEEqn.insertEquationCoupling(1, 4, EInU1);
UpEEqn.insertEquationCoupling(2, 4, EInU2);
UEqn -= fvc::div(phimd, U);


// Under-relax momentum. Note this will destroy the H and A
UEqn.relax();

UpEEqn.insertEquation(0, UEqn);
#   include "addBlockCoupledBC.H"
UpEEqn.insertBlockCoupling(0, 3, pInU, true);
UpEEqn.insertEquationCoupling(0, 3, pInU0);
UpEEqn.insertEquationCoupling(1, 3, pInU1);
UpEEqn.insertEquationCoupling(2, 3, pInU2);
