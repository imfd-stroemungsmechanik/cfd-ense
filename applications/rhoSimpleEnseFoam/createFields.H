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
        "Uc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(rho*U) & mesh.Sf()
);

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

bool isothermal = simple.dict().getOrDefault("isothermal", false);

// Soret's diffusion term
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
    mesh,
    dimensionedVector("mdT", dimensionSet(1,-2,-1,0,0,0,0), Foam::vector(0,0,0))
);

if(!isothermal)
{
    mdT = constrainmdT(thermo.mu()*1/(2*thermo.T())*fvc::grad(thermo.T()));
}


// Fick's diffusion term
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
    -thermo.mu()*(1/p*fvc::grad(p))
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
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U + Ud
);

// Diffusive face-flux fields
surfaceScalarField phimdT = fvc::flux(mdT);
surfaceScalarField phimdp = fvc::flux(mdp);
surfaceScalarField phimd = phimdT + phimdp;

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
    fvc::div(phi + phimd)
);
