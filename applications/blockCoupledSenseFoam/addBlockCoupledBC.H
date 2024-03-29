    // Hack block-coupled boundary conditions: due for rewrite
    const volScalarField nuEff = turbulence->muEff();

    forAll (U.boundaryField(), patchI)
    {
        if (U.boundaryField()[patchI].blockCoupled())
        {
            // Insert correcting fully implicit coupling coefficient

            const labelList fc = mesh.boundary()[patchI].faceCells();

            const fvPatchVectorField& Up = U.boundaryField()[patchI];

            // Warning: hacked for nuEff in viscosity
            const scalarField nutpMagSf =
                nuEff.boundaryField()[patchI]*
                mesh.magSf().boundaryField()[patchI];

            // Get boundary condition contribution to matrix diagonal
            tensorField patchDiag =
                -Up.blockGradientInternalCoeffs()().asSquare()*nutpMagSf;

            // Get matrix diagonal
            CoeffField<vector5>::squareTypeField& blockDiag =
                UpEEqn.diag().asSquare();

            forAll (fc, faceI)
            {
                blockDiag[fc[faceI]](0, 0) += patchDiag[faceI].xx();
                blockDiag[fc[faceI]](0, 1) += patchDiag[faceI].xy();
                blockDiag[fc[faceI]](0, 2) += patchDiag[faceI].xz();

                blockDiag[fc[faceI]](1, 0) += patchDiag[faceI].yx();
                blockDiag[fc[faceI]](1, 1) += patchDiag[faceI].yy();
                blockDiag[fc[faceI]](1, 2) += patchDiag[faceI].yz();

                blockDiag[fc[faceI]](2, 0) += patchDiag[faceI].zx();
                blockDiag[fc[faceI]](2, 1) += patchDiag[faceI].zy(); // (3,1) -> (2,1)
                blockDiag[fc[faceI]](2, 2) += patchDiag[faceI].zz(); // (3,2) -> (2,2)
            }

         // Get boundary condition contribution to matrix source
         vectorField patchSource =
             -Up.blockGradientBoundaryCoeffs()*nutpMagSf;

         // Get matrix source
         Field<vector5>& blockSource = UpEEqn.source();

            forAll (fc, faceI)
            {
                blockSource[fc[faceI]](0) -= patchSource[faceI](0);
                blockSource[fc[faceI]](1) -= patchSource[faceI](1);
                blockSource[fc[faceI]](2) -= patchSource[faceI](2);
            }
        }
    }
