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
    CgridShallowWaterFoamSplit2809

Description
    Transient solver for inviscid shallow-water equations with rotation with
    an RK IMEX scheme.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.

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
    #include "readEarthProperties.H"
    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.name() << nl << endl;

        #include "CourantNo.H"
        double tolerance = 1e-07;
        
        // half Coriolis
        for(int it = 0; it < 2; it++)
        {
            U = U.oldTime() - 0.5*dt*(F ^ U);
        }
        
        // half gravity (expl)
        phi = fvc::flux(h*U);
        if (num.FB) h = h.oldTime() - 0.5*dt*fvc::div(phi);
        hf = fvc::interpolate(h);
        
        U = 
        (
            h.oldTime()*U
          - 0.5*dt*fvc::reconstruct(magg*hf*fvc::snGrad(h+h0)*mesh.magSf())
        )/h;
        
        if (!num.FB) h = h.oldTime() - 0.5*dt*fvc::div(phi);
        phi = fvc::flux(h*U);
        
        
        // Momentum advection with iterations for the non-linearity
        // First store the previous velocity (from after half a gravity)
        U.oldTimeRef() = U;
        
        // First implicit RK stage
        for(int its = 0; its < num.maxIts; its++)
        {
            fvVectorMatrix UEqn
            (
                1/dt * (fvm::Sp(h,U) - h*U.oldTime()) //fvm::ddt(h, U)
              + 0.25*fvm::div(phi, U)
            );
            UEqn.solve();
            phi = fvc::flux(h*U);
        };
        
        U.oldTimeRef() = U;
        volVectorField U_a = U;
        surfaceScalarField phi_a = phi;
        
        // Second implicit RK stage
        for(int its = 0; its < num.maxIts; its++)
        {
            fvVectorMatrix UEqn
            (
                1/dt * (fvm::Sp(h,U) - h*U.oldTime()) //fvm::ddt(h, U)
              + 2/3*fvc::div(phi_a, U_a, "div(phi,U)")
              + 1/3*fvm::div(phi, U)
            );
            UEqn.solve();
            phi = fvc::flux(h*U);
        };
        
        // Previous stages are no longer needed
        U.oldTimeRef() = U;
        h.oldTimeRef() = h;
        phi = fvc::flux(h*U);
        
        // Final gravity (implicit), with non-linear iterations
        
        for(int its = 0; its < num.maxIts; its++)
        {
            hf = fvc::interpolate(h);
            
            fvScalarMatrix hEqn
            (
                1/dt * (fvm::Sp(1,h) - h.oldTime()) //fvm::ddt(h)
              + 0.5*fvc::div(phi)
              - 0.25*fvm::laplacian(dt*magg*hf, h)
              - 0.25*fvc::laplacian(dt*magg*hf, h0)
            );
            hEqn.solve();
        };

        U = (h.oldTime()*U.oldTime() - 0.5*dt*h*magg*fvc::reconstruct(fvc::snGrad(h+h0)*mesh.magSf()))/h;
        phi = phi - 0.5*dt*magg*hf*fvc::snGrad(h+h0)*mesh.magSf();
        
        // Final half Coriolis
        U -= 0.5*dt*(F ^ U);
        

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
