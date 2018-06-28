/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    viscoelasticFluidDyMFoam

Description
    Transient solver for incompressible, laminar flow of viscoelastic fluids.

Author
    Alberto Serrano. All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "viscoelasticModel.H"
#include "extrapolatedCalculatedFvPatchFields.H"
//#include "extrapolatedCalculatedFvPatchField.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
// #   include "createMesh.H"
#   include "createDynamicFvMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

#       include "readControls.H"
#       include "checkTotalVolume.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;


        /// ADDED start

        bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());

#       include "volContinuity.H"

        if (correctPhi && (mesh.moving() || meshChanged))
        {
            // Fluxes will be corrected to absolute velocity
            // HJ, 6/Feb/2009
#           include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (mesh.moving() && checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }

// #       include "UEqn.H"

        /// end


        // Pressure-velocity SIMPLE corrector loop
        for (int corr = 0; corr < nCorr; corr++)
        {

#           include "UEqn.H"

            // Does the order of the next two lines matter?
            UEqn().relax(); //What's relax?

            p.boundaryField().updateCoeffs();

            rAU = 1.0/UEqn().A(); // from rUA -> rAU (why?)
            U = rAU*UEqn().H();
            UEqn.clear();
            phi = fvc::interpolate(U) & mesh.Sf();
            adjustPhi(phi, U, p);

            // Store pressure for under-relaxation
            p.storePrevIter();

            // Non-orthogonal pressure corrector loop
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                // pEqn.solve();

                // NEW start
                if( corr == nCorr - 1 && nonOrth == nNonOrthCorr ) {
                    pEqn.solve(mesh.solutionDict().solver(p.name() + "Final"));
                } else {
                    pEqn.solve(mesh.solutionDict().solver(p.name()));
                } // end

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Make the fluxes relative to the mesh motion (ADDED)
            fvc::makeRelative(phi, U);

            // Momentum corrector
            U -= rAU*fvc::grad(p);
            U.correctBoundaryConditions();

            visco.correct();
        }

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
