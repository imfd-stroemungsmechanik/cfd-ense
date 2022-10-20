    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);   

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)

      + fvm::div(phiMd, U)
      + fvc::div((md*U)-2*sph(md*U))
     ==
        fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
