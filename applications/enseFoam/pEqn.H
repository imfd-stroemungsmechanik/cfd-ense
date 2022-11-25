volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
tUEqn.clear();

bool closedVolume = false;

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

// Assemble second part of continuity equation from new ENSE terms
fvScalarMatrix epEqn
(
    fvm::laplacian(-thermo.mu()/p, p)
  + fvc::div(phiMdT)
);

if (simple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );

    phiHbyA -= fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

    while (simple.correctNonOrthogonal())
    {
        // Assemble CNSE part of pressure equation
        fvScalarMatrix cpEqn
        (
            fvc::div(phiHbyA)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf, p)          
          ==
            fvOptions(psi, p, rho.name())
        );

        // Combine both pressure equation assemblies
        fvScalarMatrix pEqn
        (
            cpEqn
          + epEqn
        );

        // Relax the pressure equation to maintain diagonal dominance
        pEqn.relax();

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();
        
        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + cpEqn.flux();
            phiMdp = epEqn.flux();
        }
    }
}
else
{
    closedVolume = adjustPhi(phiHbyA, U, p);    

    while (simple.correctNonOrthogonal())
    {
        // Assemble CNSE part of pressure equation
        fvScalarMatrix cpEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)         
          ==
            fvOptions(psi, p, rho.name())
        );

        // Combine both pressure equation assemblies
        fvScalarMatrix pEqn
        (
            cpEqn
          + epEqn
        );        

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();
        
        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + cpEqn.flux();
            phiMdp = epEqn.flux();
        }
    }
}

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

bool pLimited = pressureControl.limit(p);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

if (pLimited || closedVolume)
{
    p.correctBoundaryConditions();
}

rho = thermo.rho();

if (!simple.transonic())
{
    rho.relax();
}
