/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    myInterFoam

Group
    grpMultiphaseSolvers

Description
    While interFoam is a Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach, with optional mesh motion
    and mesh topology changes including adaptive re-meshing, interFoamPlus provide additional
    utilities for calculating of resolved TKE budjets.
    
    The aformentioned terms are acquired when LES approach (explicit here, SGS-related terms are
    skipped for iLES) is chosen for the turbulence modeling. The additional lines ignore the surface
    tension, and use the unconditioned temporal averaging with Reynolds decomposition framework for
    analysing the incompressible highly variable density turbilence (IHVDT) flows.In other words, it is
    a Reynolds-Favre-LES approach.
    
    Please note that some terms are calculated when the simulation is finished through "interFoamPost".
    These terms are: (i) TKE productions and (ii) TKE advection.
    
    In order to run this solver properly, the following fields also (in addition to default ones) need
    to be available in the time folders:
    UMean and p_rghMean. The function objects (FO) can provide them. Thus, one may run interFoam for the transition
    part of simulation and then interFoamPlus afterwards. 
    
	The Favre velocity Utilda used in the TMF/TKE diagnostics is computed inside the solver at each time-step.
	We form the instantaneous resolved momentum rhoU, apply a volume-weighted 1-ring spatial filter
	(cell + immediate face-neighbours) to both rhoU and rho and compute Utilda = rhoU/rho. The filter is implemented
	as a volume-weighted average to respect non-uniform cell volumes;
	
References:

1- Saeedi, M. & Wang, B-C., 2015. "Large-eddy simulation of turbulent flow and dispersion over a
   matrix of wall-mounted cubes." Phys. Fluids, 27, 115104.
 
2) Hendrikson, K. & Yue, D.K.-P., 2019. "Wake behind a three-dimensional dry transom stern.
   Part 2. Analysis and modelling of incompressible highly variable density turbulence."
   J. Fluid Mech., 875, 884–913.
   
3) Asgari, E. & Tadjfar, M., 2017. “Assessment of Four Inﬂow Conditions on Large-Eddy
   Simulation of a Gently Curved Backward-Facing Step.” J. Turbul., 18(1), 61–86.    

Date: 20 January 2023

M. Reza Rastan, https://scholar.google.com/citations?user=L5XYxF0AAAAJ&hl=en&oi=ao
Postdoctoral Research Assistant
Department of Engineering Mechanics
School of Naval Architecture, Ocean and Civil Engineering
Shanghai Jiao Tong University, Shanghai, China
rastan@sjtu.edu.cn
mohamad.rastan@yahoo.com

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

	// neighbor list and volumes (cheap to reuse each step)
	const labelListList& neigh = mesh.cellCells();
	const scalarField& vol = mesh.V();
	
    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        

        Info<< "Calculating the TKE terms \n" << endl;
// * * * * * * * * * * * * * * * * * * * * * * * * * * interFoamPlus additive * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

                		// * * * * * * * * * * * * compute Favre (mass-weighted filtered) velocity Utilda (1-ring box filter) * * * * * * * * * * * //
// small regularizer for safe division
static const scalar SMALL = 1.0e-30;

// Recompute instantaneous momentum for this time step
rhoU_inst.internalField() = (rho * U).internalField();  // update instantaneous rho*U

// 1-ring volume-weighted average (overwrite rhoU_bar / rho_bar)
forAll(mesh.cells(), cellI)
{
    // start with the cell itself
    vector sumRhoU = rhoU_inst.internalField()[cellI] * vol[cellI];
    scalar sumRho  = rho.internalField()[cellI] * vol[cellI];
    scalar volSum  = vol[cellI];

    const labelList& nb = neigh[cellI];
    forAll(nb, iNb)
    {
        label nbCell = nb[iNb];
        sumRhoU += rhoU_inst.internalField()[nbCell] * vol[nbCell];
        sumRho  += rho.internalField()[nbCell] * vol[nbCell];
        volSum  += vol[nbCell];
    }

    // store filtered (volume-average)
    rhoU_bar.internalField()[cellI] = sumRhoU / (volSum + SMALL);
    rho_bar.internalField()[cellI]  = sumRho  / (volSum + SMALL);
}

// Ensure filtered fields have consistent BCs before division
rhoU_bar.correctBoundaryConditions();
rho_bar.correctBoundaryConditions();

// compute Favre velocity: Utilda = rhoU_bar / rho_bar  (safe division)
Utilda.internalField() = rhoU_bar.internalField() / (rho_bar.internalField() + SMALL);
Utilda.correctBoundaryConditions();

                		// * * * * * * * * * * * * Parameter definition * * * * * * * * * * * //
           
        nul = mixture.nu(); 										
		
        UPrime = (Utilda - UtildaMean);                    			// Resolved velocity fluctuation vector
        pPrime = (p_rgh - p_rghMean);            					// Resolved pressure fluctuation
        
		strainTensor = symm(fvc::grad(UPrime));            	        // Strain rate tensor of resolved fluctuations
		SGSstress = (2.0*symm(fvc::grad(Utilda)))*(rho*nut);		// Tau_SGS			
        
        TMF = (rho*UPrime);                    						// Turbulent mass flux
        TKE = (0.5*rho*magSqr(UPrime));            					// Turbulent kinetic energy
        
        
                		// * * * * * * * * * * * * Resolved mass-weighted second moment * * * * * * * * * * * //
        
         ReyStr = symm(TMF * UPrime);         			        	// The time-mean of this term is used for turbulent production I calculation; 
	
        				// * * * * * * * * * * * * TKE correlation terms * * * * * * * * * * * //
     			
        turbDiff = -0.5*rho*(UPrime*magSqr(UPrime));          		// Turbulent diffusion +++ divergence operator is applied at the end of simulation through "interFoamPost"
        pressDiff = -pPrime*UPrime;         		     			// Pressure diffusion +++ divergence operator is applied at the end of simulation through "interFoamPost"
        pressDila = pPrime * fvc::div(UPrime);              		// Pressure dilatation
        pDiffDill = -UPrime & fvc::grad(pPrime);              		// Pressure diffusion + Pressure dilatation, i.e. pressure contribution/work
        viscDiff = 2.0*rho*nul*(strainTensor & UPrime);   		// Viscous diffusion +++ divergence operator is applied at the end of simulation through "interFoamPost"
        viscDiss = -2.0*rho*nul*(strainTensor && strainTensor);  // Viscous dissipation
        SGSDiff = -rho*(UPrime & SGSstress);                        // SGS diffusion +++ divergence operator is applied at the end of simulation through "interFoamPost"
        SGSDiss = rho*(strainTensor && SGSstress);                  // SGS dissipation
      
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //