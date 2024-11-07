#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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


Info<< "Reading field UMean\n" << endl;
volVectorField UMean
(
    IOobject
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field p_rghMean\n" << endl;
volScalarField p_rghMean
(
    IOobject
    (
        "p_rghMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

//////////////////////////////////////////////////////////////////////////////// 
volScalarField TKE ///////////
(
    IOobject
    (
        "TKE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "TKE", 
           rho.dimensions()*dimVelocity*dimVelocity,
           0.0
      )
);

volScalarField pPrime /////////
(
    IOobject
    (
        "pPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "pPrime", 
           p.dimensions(),
           scalar(0.0)
      )    
);

volVectorField TMF ////////////////////
(
    IOobject
    (
        "TMF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedVector
      (
           "TMF", 
           U.dimensions()*rho.dimensions(),
           vector (0,0,0)
      )
);

//~ 

volSymmTensorField devSGSstress //////////////
(
    IOobject
    (
        "devSGSstress",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "devSGSstress", 
           p.dimensions(),
           symmTensor::zero
      )
);

volSymmTensorField SGSstress //////////////
(
    IOobject
    (
        "SGSstress",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "SGSstress", 
           p.dimensions(),
           symmTensor::zero
      )
);

volSymmTensorField stressTensor //////////////
(
    IOobject
    (
        "stressTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "stressTensor", 
           p.dimensions(),
           symmTensor::zero
      )
);

volSymmTensorField ReyStr ///////////////
(
    IOobject
    (
        "ReyStr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "ReyStr", 
           U.dimensions()*U.dimensions()*rho.dimensions(),
           symmTensor::zero
      )
);

volVectorField turbDiff /////////////
(
    IOobject
    (
        "turbDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "turbDiff", 
           U.dimensions()*U.dimensions()*U.dimensions()*rho.dimensions(),
           vector (0,0,0)
      )
);

volVectorField pressDiff //////////////
(
    IOobject
    (
        "pressDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "pressDiff", 
           p.dimensions()*U.dimensions(),
           vector (0,0,0)
      )
);

volScalarField pressDila //////////
(
    IOobject
    (
        "pressDila",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "pressDila", 
           p.dimensions()*U.dimensions()/dimLength,
           0.0
      )
);

volScalarField pDiffDill //////////
(
    IOobject
    (
        "pDiffDill",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "pDiffDill", 
           p.dimensions()*U.dimensions()/dimLength,
           0.0
      )
);

volSymmTensorField strainTensor ///////////////
(
    IOobject
    (
        "strainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "strainTensor", 
           dimVelocity/dimLength,
           symmTensor::zero
      )
);

volVectorField viscDiff ///////////////
(
    IOobject
    (
        "viscDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "viscDiff", 
           rho.dimensions()*dimVelocity*dimVelocity*dimVelocity,
           vector (0,0,0)
      )
);

volScalarField viscDiss //////////////
(
    IOobject
    (
        "viscDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "viscDiss", 
           rho.dimensions()*dimVelocity*dimVelocity*dimVelocity/dimLength,
           0.0
      )
);

volVectorField SGSDiff /////////////
(
    IOobject
    (
        "SGSDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "SGSDiff", 
           rho.dimensions()*dimVelocity*dimVelocity*dimVelocity,
           vector (0,0,0)
      )
);

volVectorField SGSDif2 /////////////
(
    IOobject
    (
        "SGSDif2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "SGSDif2", 
           rho.dimensions()*dimVelocity*dimVelocity*dimVelocity,
           vector (0,0,0)
      )
);

volScalarField SGSDiss ///////////
(
    IOobject
    (
        "SGSDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "SGSDiss", 
           rho.dimensions()*dimVelocity*dimVelocity*dimVelocity/dimLength,
           0.0
      )
);

volVectorField UPrime ///////////
(
    IOobject
    (
        "UPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedVector
      (
           "UPrime", 
           dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);

volScalarField nut ///////////
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "nut", 
           dimVelocity*dimLength,
           0.0
      )
);

volScalarField nu_lam ///////////
(
    IOobject
    (
        "nu_lam",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "nu_lam", 
           dimVelocity*dimLength,
           0.0
      )
);

Info<< "Reading field Reynolds Stress tensor\n" << endl; // It is written in the time folders, calculated by TurbulenceFields FO

volSymmTensorField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
//////////////////////////////////////////////////////////////////////////////// 

#include "createMRF.H"
#include "createFvOptions.H"
