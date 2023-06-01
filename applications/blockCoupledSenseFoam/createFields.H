Info<< "Reading thermophysical properties\n" << endl;

autoPtr<basicPsiThermo> pThermo
(
    basicPsiThermo::New(mesh)
);
basicPsiThermo& thermo = pThermo();
//thermo.validate(args.executable(), "e");

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();

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

    Info << "Reading field U\n" << endl;
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


#   include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    // Block vector field for velocity (first three entries), pressure (fourth
    // entry) and energy (fifth entry)
    Info << "Creating field UpE\n" << endl;
    volVector5Field UpE
    (
        IOobject
        (
            "UpE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector5("zero", dimless, vector5::zero)
    );

// Temperature dependent diffusion term
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
    constrainmdT(thermo.mu()*1/(2*thermo.T())*fvc::grad(thermo.T()))
);

// Pressure dependent diffusion term
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
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U + Ud
);

surfaceScalarField phimdp = fvc::interpolate(mdp) & mesh.Sf();
surfaceScalarField phimdT = fvc::interpolate(mdT) & mesh.Sf();
surfaceScalarField phimd  = phimdp + phimdT;

//Info<< "Creating field rAU\n" << endl;

//(
//    IOobject
//    (
//        "rAU",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//    mesh,
//    runTime.deltaT()
//);

mesh.schemesDict().setFluxRequired(thermo.e().name());
mesh.schemesDict().setFluxRequired(p.name());
mesh.schemesDict().setFluxRequired("pcorr");
