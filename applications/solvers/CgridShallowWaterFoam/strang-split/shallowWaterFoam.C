/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
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
    AdImExShallowWaterFoam

Description
    Transient solver for inviscid shallow-water equations with rotation with
    adaptive implicit-explicit advection.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.
    
    Adaptive implicit-explicit not implemented yet. 
    
    There is an option "opSplit" to use or not use operator splitting

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvMesh.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcLaplacian.H"
#include "fvcReconstruct.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "numericalParameters.H"
    #define dt runTime.deltaT()
    #define alpha num.alpha
    #include "readEarthProperties.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.name() << nl << endl;

        #include "CourantNo.H"
        
        // half Coriolis
        volVectorField U_prime = U - (dt/2 * F ^ U);
        U_prime = U - (dt/4 * F ^ (U + U_prime));
        U = U - (dt/4 * F^(U + U_prime));
        
        surfaceScalarField phi_C = hf*fvc::flux(U);
        volVectorField U_C = U;
        volVectorField hU_C = h*U;
        
        // half advection
        for(int outerCorr = 0; outerCorr < num.nOuterCorrs; outerCorr++)
        {
            fvScalarMatrix hEqnAdv
            (
                fvm::Sp(1,h) - h.oldTime()
              + 1/4 * dt * fvc::div(phi_C)
              + 1/4 * dt * fvm::div(fvc::flux(U),h)
            );
            hEqnAdv.solve();
            hf = fvc::interpolate(h);
            
            fvVectorMatrix UEqnAdv
            (
                fvm::Sp(1,hU) - hU_C
              + 1/4 * dt * fvc::div(fvc::flux(U_C),hU_C)
              + 1/4 * dt * fvm::div(fvc::flux(U),hU)
            );
            UEqnAdv.solve();
            U = hU/h;
        }
        volScalarField h_A = h;
        surfaceScalarField hf_A = fvc::interpolate(h_A);
        volVectorField U_A = U;
        volVectorField hU_A = hU;
            
        
        
        // full gravity
        for(int outerCorr = 0; outerCorr < num.nOuterCorrs; outerCorr++)
        {
            fvScalarMatrix hEqnGrav
            (
                fvm::Sp(1,h) - h_A
              + dt * fvc::div(hf_A*fvc::flux(U_A))
              - 1/4 * dt * fvc::laplacian(dt*magg*hf_A, h_A)
              - 1/4 * dt * fvm::laplacian(dt*magg*hf, h)
              - 1/4 * dt * fvc::laplacian(dt*magg*h0, h_A)
              - 1/4 * dt * fvm::laplacian(dt*magg*h0, h)
            );
            hEqnGrav.solve();
            hf = fvc::interpolate(h);
            
            fvVectorMatrix UEqnGrav
            (
                fvm::Sp(1,hU) - hU
              + 1/2 * dt * fvc::reconstruct(magg*hf*fvc::snGrad(h+h0)*mesh.magSf())
              + 1/2 * dt * fvc::reconstruct(magg*hf_A*fvc::snGrad(h_A+h0)*mesh.magSf())
            );
            UEqnGrav.solve();
            U = hU/h;
        }
        surfaceScalarField phi_G = hf*fvc::flux(U);
        volVectorField U_G = U;
        volVectorField hU_G = h*U;
        volScalarField h_G = h;
        
        // half advection
        for(int outerCorr = 0; outerCorr < num.nOuterCorrs; outerCorr++)
        {
            fvScalarMatrix hEqnAdv2
            (
                fvm::Sp(1,h) - h_G
              + 1/4 * dt * fvc::div(phi_G)
              + 1/4 * dt * fvm::div(fvc::flux(U),h)
            );
            hEqnAdv2.solve();
            hf = fvc::interpolate(h);
            
            fvVectorMatrix UEqnAdv2
            (
                fvm::Sp(1,hU) - hU_G
              + 1/4 * dt * fvc::div(fvc::flux(U_G),hU_G)
              + 1/4 * dt * fvm::div(fvc::flux(U),hU)
            );
            UEqnAdv2.solve();
            U = hU/h;
        }
        
        // half Coriolis
        U_prime = U - (dt/2 * F ^ U);
        U_prime = U - (dt/4 * F ^ (U + U_prime));
        U = U - (dt/4 * F^(U + U_prime));
        
        hU = h*U;

        E = 0.5*(h)*magSqr(U) + 0.5*magg*sqr(h);
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
