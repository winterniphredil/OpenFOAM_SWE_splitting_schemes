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
        
        // CORIOLIS OUTSIDE
        
        // half Coriolis
        
        volVectorField U_c1 = U - dt/2*(F ^ U);
        U = U - dt/2*(F ^ U_c1);
        volVectorField hU_C = h*U;
        volVectorField U_C = U;
        
        // half advection (expl) - add iteration (also over whole thing?)
        
        volScalarField h_C = h;
        
        for (int iIt=0; iIt < nIters; iIt++)
        {
            volScalarField dhdt = U & fvc::grad(h);
            h = h_C - dt/2 * dhdt;
            hU = hU_C - dt/2 * fvc::div(fvc::flux(hU),U);
            U = hU/h;
        }
        
        volScalarField h_A = h;
        volVectorField U_A = U;
        
        hf = fvc::interpolate(h);
        surfaceScalarField hf_A = hf;
        
        // full gravity (impl)
        
        
        for (int iIt=0; iIt < nIters; iIt++)
        {
            fvScalarMatrix hEqnGrav
            (
                fvm::Sp(1,h) - h_A
              + dt * h_A *fvc::div(U)
              - 1/4 * fvc::laplacian(dt*dt*magg*hf_A, h_A)
              - 1/4 * fvm::laplacian(dt*dt*magg*hf, h)
              - 1/4 * fvc::laplacian(dt*dt*magg*h0, h_A)
              - 1/4 * fvm::laplacian(dt*dt*magg*h0, h)
            );
            
            hEqnGrav.solve();
            hf = fvc::interpolate(h);
            
            dhUdt = fvc::reconstruct(hEqnGrav.flux()*2)/(dt*dt);
        }
        
        hU += dhUdt * dt;
        
        U = hU/h;
        surfaceScalarField Uf = fvc::interpolate(U)&mesh.Sf();
        volScalarField h_G = h;
        volVectorField hU_G = hU;
        
        // full advection (impl)
        
        for (int iIt=0; iIt < nIters; iIt++)
        {
            fvScalarMatrix hEqnAdv
            (
                fvm::Sp(1,h) - h_G
              + dt/2 * fvm::div(Uf,h)
              - dt/2 * fvm::Sp(fvc::div(U),h)
            );
            hEqnAdv.solve();
            
            fvVectorMatrix UEqnAdv
            (
                fvm::Sp(h,U) - hU_G
              + dt/2 * fvm::div(fvc::flux(hU),U)
            );
            UEqnAdv.solve();
            
            hU = h*U;
            Uf = fvc::interpolate(U)&mesh.Sf();
        }
        
        // half Coriolis
        
        U -= dt/2 * (F ^ (U + U_A - U_C));
        
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
