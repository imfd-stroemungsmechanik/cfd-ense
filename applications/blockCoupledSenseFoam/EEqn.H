// Energy equation

{
    volScalarField& he = thermo.e();

    const volTensorField tauT = -thermo.mu()*(fvc::grad(U) + dev2(T(fvc::grad(U)))) + U*md + dev2(md*U);

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      - fvm::laplacian(turbulence->alphaEff(), he)
      //+ fvc::div(p*U) - Introduced by the coupling below
      + fvc::div(tauT & U)
    );

    const surfaceScalarField ef = fvc::interpolate(thermo.e());

    // div(phi, he) - U and p coupling
    {
        fvScalarMatrix pInE(fvm::div(ef*phiPsi, p));
        BlockLduSystem<vector, scalar> UInE(fvm::UDiv(rhof*ef, U));
        UpEEqn.insertEquationCoupling(4, 3, pInE);
        UpEEqn.insertBlockCoupling(4, 0, UInE, false);
        EEqn -= 2*fvc::div(phi, he);
    }

    // div(phi, k) - U and p coupling
    {
        const surfaceScalarField Kf = fvc::interpolate(K);
        BlockLduSystem<vector, scalar> UInE(fvm::UDiv(rhof*Kf, U));
        fvScalarMatrix pInE(fvm::div(Kf*phiPsi,p));
        UpEEqn.insertEquationCoupling(4, 3, pInE);
        UpEEqn.insertBlockCoupling(4, 0, UInE, false);
        EEqn -= fvc::div(phi, K);
    }

    // div(p*U) - U and p coupling
    {
        BlockLduSystem<vector, scalar> UInE(fvm::UDiv(fvc::interpolate(p), U));
        fvScalarMatrix pInE(fvm::div(phi/rhof, p));
        UpEEqn.insertBlockCoupling(4, 0, UInE, false);
        UpEEqn.insertEquationCoupling(4, 3, pInE);
        EEqn -= fvc::div(p*U);
    }


    // ENSE:

    EEqn +=
    (
        fvm::div(phimdp, he) + constrainmdT(fvm::laplacian(0.5*thermo.mu(), he))
      - K * fvc::div(phimd)
    );

    // div(phimdp, he) - p coupling
    {
        fvScalarMatrix pInE(fvm::laplacian(nmuByp*ef, p));
        UpEEqn.insertEquationCoupling(4, 3, pInE);
        EEqn -= fvc::div(phimdp, he);
    }


    EEqn.relax();
    UpEEqn.insertEquation(4, EEqn);
}
