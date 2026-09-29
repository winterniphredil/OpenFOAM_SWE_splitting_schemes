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
        
        
        volScalarField hprev = h;
        volVectorField Uprev = U;
        
        int its = 0;
        do
        {
            Uprev = U;
            
            // First implicit RK stage
            fvVectorMatrix UEqn
            (
                1/dt * (fvm::Sp(h,U) - h*U.oldTime()) //fvm::ddt(h, U)
              + 0.25*fvm::div(phi, U)
            );
            UEqn.solve();
            phi = fvc::flux(h*U);
            
            // Second implicit RK stage
            UEqn = fvVectorMatrix
            (
                1/dt * (fvm::Sp(h,U) - h*U.oldTime()) //fvm::ddt(h, U)
              + 2/3*fvc::div(phi, U)
              + 1/3*fvm::div(phi, U)
            );
            UEqn.solve();
            phi = fvc::flux(h*U);
            its++;
        }
        while (max(mag(U-Uprev)).value()>tolerance && its<num.maxIts);
        
        // Previous stages are no longer needed
        U.oldTimeRef() = U;
        h.oldTimeRef() = h;
        phi.oldTimeRef() = fvc::flux(h*U);
        
        // Final gravity (implicit), with non-linear iterations
        
        its = 0;
        do
        {
            hprev = h;
            Uprev = U;
            hf = fvc::interpolate(h);
            
            fvScalarMatrix hEqn
            (
                1/dt * (fvm::Sp(1,h) - h.oldTime())
              + 0.5*fvc::div(phi)
              - 0.25*fvm::laplacian(dt*magg*hf, h)
              - 0.25*fvc::laplacian(dt*magg*hf, h0)
            );
            hEqn.solve();
            U = U.oldTime() - 0.5*dt*magg*fvc::reconstruct(fvc::snGrad(h)*mesh.magSf());
            phi = phi.oldTime() - 0.5*dt*magg*hf*fvc::snGrad(h+h0)*mesh.magSf();
            its++;
        }
        while ((max(mag(h-hprev)).value()>tolerance || max(mag(U-Uprev)).value()>tolerance) && its<num.maxIts);
        
        // Final half Coriolis
        U -= 0.5*dt*(F ^ U);
        
        E = 0.5*h*magSqr(U) + 0.5*magg*sqr(h);
        Info<< "E registered? "
            << mesh.foundObject<volScalarField>("E")
            << nl;

        Info<< "E write option = "
            << E.writeOpt()
            << nl;

        Info<< "E name = "
            << E.name()
            << nl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
