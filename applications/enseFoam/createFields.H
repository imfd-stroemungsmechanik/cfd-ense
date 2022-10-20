Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "e");

mesh.setFluxRequired("e");

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

// Convective velocity
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

// Soret diffusion 
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
#include "constrainMdT.H"


// Fick diffusion
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


// Entire self-diffusion
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
    mdT + mdp
);


// Diffusive velocity
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

// Total velocity 
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

// Diffusive face-flux fields
surfaceScalarField phiMdT = fvc::flux(mdT);
surfaceScalarField phiMdp = fvc::flux(mdp);
surfaceScalarField phiMd = phiMdT + phiMdp;

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
    fvc::div(phi + phiMd)
);
