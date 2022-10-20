{
    
    //surfaceScalarField phiT = phi + phiMd;
    volScalarField& he = thermo.he(); 
    volScalarField K("K", 0.5*magSqr(U)); 

    volTensorField tauT = - thermo.mu()*(fvc::grad(U)+dev2(T(fvc::grad(U)))) + md*U + U*md - (2 * sph(md*U));

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + fvc::div(phi, K)        
      - fvm::laplacian(turbulence->alphaEff(), he)
      + fvm::div(phiMd,he)
      //+ fvc::div(md*he)
      - K * fvc::div(md)
      + fvc::div(p*U)
      + fvc::div(tauT & U)
     ==
        fvOptions(rho, he)
    );

//    fvScalarMatrix EEqn
//    (
//        fvm::div(phiT, he)
//      + (
//            he.name() == "e"
//          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
//          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
//        )
//      - fvm::laplacian(turbulence->alphaEff(), he)
//     ==
//        fvOptions(rho, he)
//    );



    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
