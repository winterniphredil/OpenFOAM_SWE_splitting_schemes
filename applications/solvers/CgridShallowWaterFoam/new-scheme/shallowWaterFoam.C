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
    CgridShallowWaterFoam

Description
    Semi-implicit solver for inviscid non-linear shallow-water equations
    with rotation on a C-grid.

    If the geometry is 3D then it is assumed to be one layers of cells and
    the component of the velocity normal to gravity is removed.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvMesh.H"

#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcReconstruct.H"

#include "fvmDdt.H"
#include "fvmLaplacian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "readEarthProperties.H"
    #include "createFields.H"

    const int nIters = readLabel(mesh.solution().lookup("nIterations"));
    const scalar alpha =readScalar(mesh.solution().lookup("timeOffCentre"));
    const int e = readLabel(mesh.solution().lookup("e"));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.name() << endl;

        #include "CourantNo.H"

        // Initial explicit momentum update
        dhUdt = - hf*((F^Uf) & mesh.Sf())
                - fvc::flux(fvc::div(flux,U));
        flux = flux.oldTime() + dt*(dhUdt);
        h = h.oldTime() - dt*fvc::div(flux.oldTime()) + dt * h*fvc::div(U.oldTime(), "div(U)");
        hf = fvc::interpolate(h);
        U = fvc::reconstruct(flux/hf);
        Uf = fvc::interpolate(U);

        // Outer Iterations
        for (int iIt=0; iIt < nIters; iIt++)
        {
            // Solve momentum equation on faces without the pressure gradient
            dhUdt = - hf*((F^Uf) & mesh.Sf())
                   - hf*magg*fvc::snGrad(h0)*mesh.magSf() //+h
                   - fvc::flux(fvc::div(flux,U));
            flux = flux.oldTime()
                 + dt*((1-alpha)*dhUdt.oldTime() + alpha*dhUdt);
                 
            
            flux -= (1-alpha)*dt*magg*hf.oldTime()
                  *fvc::snGrad(h.oldTime())*mesh.magSf();

            // Create the pressure equation
            fvScalarMatrix hEqn
            (
                fvm::ddt(h)
              + fvc::div((1-alpha)*flux.oldTime() + alpha*flux)
              - fvm::laplacian(sqr(alpha)*dt*magg*hf, h)
            );
            hEqn.solve();
            hf = fvc::interpolate(h);
            
            // Back substitutions
            if (alpha > 0) 
            {
                flux += hEqn.flux()/alpha;
                dhUdt += hEqn.flux()/(sqr(alpha)*dt);
            }

            // Update the velocity field on cell centres and faces
            U = fvc::reconstruct(flux/hf);
            Uf = fvc::interpolate(U);
        }

        E = 0.5*h*magSqr(U) + 0.5*magg*sqr(h);
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
