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
    
    const int nIters = readLabel(mesh.solution().lookup("nIterations"));
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.name() << nl << endl;

        #include "CourantNo.H"
        
        surfaceScalarField Uf = fvc::interpolate(U)&mesh.Sf();
        
        // CORIOLIS OUTSIDE
        
        // half Coriolis
        
        volVectorField U_c1 = U - dt/2*(F ^ U);
        volVectorField U_c2 = U - dt/2*(F ^ U_c1);
        volVectorField U_c3 = U - dt/2*(F ^ U_c2);
        U = U_c3;
        volVectorField hU_C = h*U;
        volVectorField U_C = U;
        
        // half gravity (expl)
        
        volVectorField dhUdt = -1/2 * fvc::reconstruct(magg*hf*fvc::snGrad(h+h0)*mesh.magSf());
        h = h.oldTime() - dt/2 * h * fvc::div(U);
        hU += dhUdt * dt;
        
        volScalarField h_G = h;
        volVectorField hU_G = hU;
        U = hU/h;
        
        // first advection step
        
        volScalarField hprev = h;
        volVectorField Uprev = hU/h;
        volVectorField hUprev = hU;
        do 
        {
            hprev = h;
            Uprev = hU/h;
            hUprev = hU;
            
            fvScalarMatrix h_adv_1
            (
                fvm::Sp(1,h) - h_G
              + dt/2 * fvm::div(fvc::interpolate(Uprev)&mesh.Sf(),h,"div(U,h)")
              - dt/2 * fvm::Sp(fvc::div(Uprev,"div(U)"),h)
            );
            h_adv_1.solve();
            
            fvVectorMatrix hU_adv_1
            (
                fvm::Sp(1,hU) - hU_G
              + dt/4 * fvm::div(fvc::flux(hUprev/h),hU)
            );
            hU_adv_1.solve();
        }
        while (max(mag(h-hprev)).value()>1e-7 && max(mag(hU-hUprev)).value()>1e-7);
        
        volScalarField h_A = h;
        volVectorField hU_A = hU;
        volVectorField U_A = U;
        
        
        // second advection step
        
        do 
        {
            hprev = h;
            Uprev = hU/h;
            hUprev = hU;
            
            fvScalarMatrix h_adv_1
            (
                fvm::Sp(1,h) - h_G
              + dt/3 * fvm::div(fvc::interpolate(Uprev)&mesh.Sf(),h,"div(U,h)")
              - dt/3 * fvm::Sp(fvc::div(Uprev,"div(U)"),h)
              + 2*dt/3 * fvc::div(fvc::interpolate(U_A)&mesh.Sf(),h_A)
              - 2*dt/3 * h_A * fvc::div(U_A)
            );
            h_adv_1.solve();
            
            fvVectorMatrix hU_adv_1
            (
                fvm::Sp(1,hU) - hU_G
              + dt/3 * fvm::div(fvc::flux(hUprev/h),hU)
              + 2*dt/3 * fvc::div(fvc::flux(U_A),hU_A)
            );
            hU_adv_1.solve();
        }
        while (max(mag(h-hprev)).value()>1e-7 && max(mag(hU-hUprev)).value()>1e-7);
        
        h_A = h;
        U_A = U;
        hU_A = hU;
        
        
        // second gravity step 
        
        do 
        {
            hprev = h;
            Uprev = hU/h;
            hUprev = hU;
            
            fvScalarMatrix hEqnGrav
            (
                fvm::Sp(1,h) - h_A
              + dt/2 * fvc::div(U_A) * h_A
              - 1/4 * fvm::laplacian(dt*dt*magg*fvc::interpolate(hprev), h)
              - 1/4 * fvm::laplacian(dt*dt*magg*h0, h)
            );
            
            hEqnGrav.solve();
            
            dhUdt = fvc::reconstruct(hEqnGrav.flux()*2)/(dt*dt);
        }
        while (max(mag(h-hprev)).value()>1e-7);
        
        
        hU += dhUdt * dt;
        U = hU/h;
        
        // half Coriolis
        
        U -= dt/2 * (F ^ (U + U_c3 - U_c2));
        
        hU = U*h;
        
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
