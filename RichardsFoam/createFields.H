    Info<< "Reading transportProperties\n" << endl;

// reading of the physical constant associated with the considered problem.
// Localisation of the data within the considered case: constant/transportProperties  
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar K
    (
        transportProperties.lookup("K")
    );
 
    dimensionedScalar alpha
    (
        transportProperties.lookup("alpha")
    );

    dimensionedScalar thetas
    (
        transportProperties.lookup("thetas")
    );

    dimensionedScalar thetar
    (
        transportProperties.lookup("thetar")
    );

    dimensionedScalar n
    (
        transportProperties.lookup("n")
    );
   
    dimensionedScalar C
    (
        transportProperties.lookup("C")
    );

    dimensionedScalar S
    (
        transportProperties.lookup("S")
    );

// declaration of the variable and results fields

// Water velocity field [m/s]
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

// Water saturation field [-]
    Info<< "Reading field theta\n" << endl;
    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

// Water pressure field [m] - field of resolution.
    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

// Field of residuals for the Picard loop [m].
    Info<< "Reading field err\n" << endl;
    volScalarField err
    (
        IOobject
        (
            "err",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

// Dimensionned unit scalar field [m].
    Info<< "Reading field vuz\n" << endl;
    volVectorField vuz
    (
        IOobject
        (
            "vuz",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

// Dimensionless unit vertical upward vector field.
    Info<< "Reading field usf\n" << endl;
    volScalarField usf
    (
        IOobject
        (
            "usf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "createPhi.H"

