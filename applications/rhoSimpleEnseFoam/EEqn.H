{   
    // Must be internal energy e
    volScalarField& he = thermo.he(); 
    
    const volScalarField K("K", 0.5*magSqr(Uc));  
    const volTensorField tauT = -thermo.mu()*(fvc::grad(Uc)+dev2(T(fvc::grad(Uc))))+Uc*md+dev2(md*Uc); 
    
    // Assemble div(mdT*e) and constrain face normal component of mdT to zero on wall boundaries 
    fvScalarMatrix divMdTHe(fvm::laplacian(thermo.mu()*0.5, he));
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchi)
    {
        if(patches[patchi].type() == "wall")
        {
            const label size = he.boundaryField()[patchi].size();    
            divMdTHe.internalCoeffs()[patchi] = Field<scalar>(size, pTraits<scalar>::zero);
            divMdTHe.boundaryCoeffs()[patchi] = Field<scalar>(size, pTraits<scalar>::zero);       
        }
    }    

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + fvc::div(phi, K)
      - fvm::laplacian(turbulence->alphaEff(), he)
      //- fvm::laplacian(thermo.mu(), he)
      //- fvm::laplacian(thermo.kappa()/thermo.Cv(), he)
      
      + fvm::div(phiMdp, he) + divMdTHe
      - K * fvc::div(phiMd)
      + fvc::div(p*Uc)
      + fvc::div(tauT & Uc)
     ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    phiMdT = divMdTHe.flux() * fvc::interpolate(1/(thermo.T()*thermo.Cv()));
}
