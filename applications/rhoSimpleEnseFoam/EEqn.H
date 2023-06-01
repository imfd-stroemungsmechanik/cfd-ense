{
    // Must be internal energy e
    volScalarField& he = thermo.he();

    const volScalarField K("K", 0.5*magSqr(U));
    const volTensorField tauT = -thermo.mu()*(fvc::grad(U) + dev2(T(fvc::grad(U)))) + U*md + dev2(md*U);

    fvScalarMatrix divmdThe(constrainmdT(fvm::laplacian(thermo.mu()*0.5, he)));

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + fvc::div(phi, K)
      - fvm::laplacian(turbulence->alphaEff(), he)
      //- fvm::laplacian(thermo.mu(), he)
      //- fvm::laplacian(thermo.kappa()/thermo.Cv(), he)

      + fvm::div(phimdp, he) + divmdThe
      - K * fvc::div(phimd)
      + fvc::div(p*U)
      + fvc::div(tauT & U)
     ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    phimdT = divmdThe.flux() * fvc::interpolate(1/(thermo.T()*thermo.Cv()));
}
