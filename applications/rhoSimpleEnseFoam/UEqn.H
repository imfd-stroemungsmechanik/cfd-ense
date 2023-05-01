    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(Uc);   

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, Uc)
      + MRF.DDt(rho, Uc)
      + turbulence->divDevRhoReff(Uc)

      + fvm::div(phiMd, Uc)
      + fvc::div(dev2(md*Uc))
     ==
        fvOptions(rho, Uc)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(Uc);
