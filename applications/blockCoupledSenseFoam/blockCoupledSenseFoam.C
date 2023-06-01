/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockCoupledSenseFoam

Description
    Pressure-based, steady-state solver for compressible flow of ideal
    gases. Solves Sambasivam's set of extended Navier-Stokes equations (SENSE)
    to include self-diffusion of mass. Solves a block coupled U-p-e system.

Authors
    Klas Jareteg, Chalmers University of Technology,
    Vuko Vukcevic, FMENA Zagreb.
    Johannes Schwarz, Fraunhofer EMFT

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrixFix.H"

#include "RASModel.H"
#include "basicPsiThermo.H"
#include "constrainmdT.H"


int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "initConvergenceCheck.H"
#   include "readNSEDict.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
#       include "readBlockSolverControls.H"
#       include "readFieldBounds.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // p.storePrevIter();

        // Initialize the Up block system (matrix, source and reference to Up)
        fvBlockMatrixFix<vector5> UpEEqn(UpE);

        const volScalarField K("K", 0.5*magSqr(U));
        const surfaceScalarField nmuByp = fvc::interpolate(-thermo.mu()/p);
        const surfaceScalarField muBy2e = fvc::interpolate(0.5*thermo.mu()/thermo.e());
        const surfaceScalarField phiPsi = fvc::interpolate(psi*U) & mesh.Sf();
        const surfaceScalarField rhof = fvc::interpolate(rho);

        // Assemble and insert momentum equations (including couplings)
#       include "UEqn.H"

        // Assemble and insert pressure equation (including couplings)
#       include "pEqn.H"

        // Assemble and insert energy equation (including couplings)
#       include "EEqn.H"

        // Solve the block matrix
        maxResidual = cmptMax(UpEEqn.solve().initialResidual());

        // Retrieve solution

        UpEEqn.retrieveSolution(0, U.internalField());
        UpEEqn.retrieveSolution(3, p.internalField());
        UpEEqn.retrieveSolution(4, thermo.e().internalField());

        thermo.e().correctBoundaryConditions();
        thermo.correct();
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

#       include "boundPU.H"
        //p.relax();
        rho = thermo.rho();

        phi = (fvc::interpolate(rho*U) & mesh.Sf());
        phimdp = epEqn.flux();
        phimdT = eInp.flux();
        phimd = phimdp + phimdT;

        mdT = constrainmdT(thermo.mu()*1/(2*thermo.T())*fvc::grad(thermo.T()));
        mdp = -thermo.mu()*(1/p*fvc::grad(p));
        md = mdp + mdT;
        Ud = md/rho;
        Ut = U + Ud;

#       include "contErrors.H"

        turbulence->correct();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

#       include "convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}
