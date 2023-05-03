rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), Uc, p));
tUEqn.clear();

bool closedVolume = false;

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

volScalarField rhorAtU("rhorAtU", rho*rAtU);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, Uc, phiHbyA, rhorAtU, MRF);

volScalarField nMuByp = -thermo.mu()/p;

if (simple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );

    phiHbyA +=
        fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf()
      - fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

    HbyA -= (rAU - rAtU)*fvc::grad(p);

    while (simple.correctNonOrthogonal())
    {
        // Assemble CNSE part of pressure equation
        fvScalarMatrix cpEqn
        (
            fvc::div(phiHbyA)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAtU, p)
         ==
            fvOptions(psi, p, rho.name())
        );  

        // Assemble second part of continuity equation from new ENSE terms
        fvScalarMatrix epEqn
        (
            fvm::laplacian(nMuByp, p)
          + fvc::div(phiMdT)
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
    closedVolume = adjustPhi(phiHbyA, Uc, p);

    phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU)*fvc::grad(p);

    while (simple.correctNonOrthogonal())
    {
        // Assemble CNSE part of pressure equation
        fvScalarMatrix cpEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rhorAtU, p)
          ==
            fvOptions(psi, p, rho.name())
        );

        // Assemble second part of continuity equation from new ENSE terms
        fvScalarMatrix epEqn
        (
            fvm::laplacian(nMuByp, p)
          + fvc::div(phiMdT)
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

Uc = HbyA - rAtU*fvc::grad(p);

Uc.correctBoundaryConditions();
fvOptions.correct(Uc);

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
