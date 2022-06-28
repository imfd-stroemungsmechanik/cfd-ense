Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, simple.dict());

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvOptions.H"

volVectorField mdT
(
    IOobject
    (
        "mdT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.mu() * 1/(2*thermo.T()) * fvc::grad(thermo.T())
);
#include "correctMdTBC.H"

volVectorField mdp
(
    IOobject
    (
        "mdp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -thermo.mu() * (1/p * fvc::grad(p))
);

volVectorField md
(
    IOobject
    (
        "md",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),   
    mdT + mdp//dimensionedVector("md", dimensionSet(1,-2,-1,0,0,0,0), Foam::vector(0,0,0))
);

volVectorField Ud
(
    IOobject
    (
        "Ud",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    md/rho
);

volVectorField Ut
(
    IOobject
    (
        "Ut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U + Ud
);

surfaceScalarField phiMd = fvc::interpolate(md) & mesh.Sf();

volScalarField massRes
(
    IOobject
    (
        "massRes",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(phi) + fvc::div(phiMd)
);
