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

\*---------------------------------------------------------------------------*/

#include "fvBlockMatrixFix.H"


template<class Type>
Foam::fvBlockMatrixFix<Type>::fvBlockMatrixFix
(
    GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    fvBlockMatrix<Type>(psi)
{}


template<class Type>
Foam::fvBlockMatrixFix<Type>::fvBlockMatrixFix
(
    const fvBlockMatrixFix<Type>& bxs
)
:
    fvBlockMatrix<Type>(bxs)
{}

template<class Type>
void Foam::fvBlockMatrixFix<Type>::insertCouplingDiagSource
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& matrix
)
{
    // Prepare the diagonal and source
    scalarField diag = matrix.diag();
    scalarField source = matrix.source();

    // Add boundary source contribution
    matrix.addBoundaryDiag(diag, 0);
    matrix.addBoundarySource(source, false);

    // Get reference to block diagonal of the block system
    typename CoeffField<Type>::squareTypeField& blockDiag =
        this->diag().asSquare();

    // Get reference to this source field of the block system
    Field<Type>& b = this->source();

    // Set off-diagonal coefficient
    forAll(diag, cellI)
    {
        blockDiag[cellI](dirI, dirJ) += diag[cellI];
        b[cellI](dirI) += source[cellI];
    }
}

template<class Type>
void Foam::fvBlockMatrixFix<Type>::insertCouplingUpperLower
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& matrix
)
{
    if (matrix.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    if (matrix.symmetric() && this->symmetric())
    {
        Info<< "Both fvScalarMatrix and block matrix are symmetric: " << nl
            << "inserting only upper triangle"
            << endl;
    }
    else
    {
        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = matrix.lower();

        typename CoeffField<Type>::squareTypeField& blockLower =
            this->lower().asSquare();

        forAll (lower, cellI)
        {
            blockLower[cellI](dirI, dirJ) += lower[cellI];
        }
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();

        typename CoeffField<Type>::squareTypeField& blockUpper =
            this->upper().asSquare();

        forAll (upper, cellI)
        {
            blockUpper[cellI](dirI, dirJ) += upper[cellI];
        }
    }
    else
    {
        FatalErrorIn
        (
            "void fvBlockMatrix<Type>::insertCouplingUpperLower\n"
            "(\n"
            "    const direction dirI,\n"
            "    const direction dirJ,\n"
            "    const fvScalarMatrix& m\n"
            ")"
        )   << "Error in matrix insertion: problem with block structure."
            << abort(FatalError);
    }

    // Insert block interface fields
    forAll (this->interfaces(), patchI)
    {
        if (this->interfaces().set(patchI))
        {
            // Couple upper and lower
            const scalarField& cUpper = matrix.boundaryCoeffs()[patchI];
            const scalarField& cLower = matrix.internalCoeffs()[patchI];

            typename CoeffField<Type>::squareTypeField& blockUpper =
                this->coupleUpper()[patchI].asSquare();

            typename CoeffField<Type>::squareTypeField& blockLower =
                this->coupleLower()[patchI].asSquare();

            forAll (cUpper, faceI)
            {
                blockUpper[faceI](dirI, dirJ) += cUpper[faceI];
                blockLower[faceI](dirI, dirJ) += cLower[faceI];
            }
        }
    }
}



template<class Type>
void Foam::fvBlockMatrixFix<Type>::insertEquationCoupling
(
    const direction dirI,
    const direction dirJ,
    fvScalarMatrix& matrix
)
{
    insertCouplingDiagSource(dirI, dirJ, matrix);
    insertCouplingUpperLower(dirI, dirJ, matrix);
}

/*
template<class Type>
void Foam::fvBlockMatrixFix<Type>::insertVectorEquationCoupling
(
    const direction dirI,
    const direction dirJ,
    fvVectorMatrix& matrix,
    const int32_t component
)
{
    if(component == -1)
    {
        const direction nCmpts = pTraits<vector>::nComponents;
        for (direction cmpt = 0; cmpt < nCmpts; cmpt++)
        {
            insertVectorEquationCoupling(dirI + cmpt, dirJ + cmpt, matrix, cmpt);
        }

    }
    else
    {
        fvScalarMatrix sMatrix(matrix.psi().component(component), matrix.dimensions());
        sMatrix.diag() = matrix.diag();
        sMatrix.upper() = matrix.upper();
        sMatrix.lower() = matrix.lower();

        forAll (matrix.boundaryCoeffs(), patchI)
        {
            sMatrix.boundaryCoeffs()[patchI] = matrix.boundaryCoeffs()[patchI].component(component);
            sMatrix.internalCoeffs()[patchI] = matrix.internalCoeffs()[patchI].component(component);
        }
        insertEquationCoupling(dirI, dirJ, sMatrix);
    }
}
*/


