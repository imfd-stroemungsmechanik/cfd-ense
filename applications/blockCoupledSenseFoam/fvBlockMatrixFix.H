/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
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

Class
    Foam::fvBlockMatrix

Description
    fvBlockMatrix is an extension of fvMatrix for block coupled types. It holds
    general insertion and retrieval tools for block coupling and specific
    functions for p-U coupling.

Author
    Vuko Vukcevic, FMENA Zagreb.
    Update by Hrvoje Jasak
    Minorest change by JS

SourceFiles
    fvBlockMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef fvBlockMatrixFix_H
#define fvBlockMatrixFix_H

#include "fvBlockMatrix.H"

namespace Foam
{


template<class Type>
class fvBlockMatrixFix;


template<class Type>
class fvBlockMatrixFix
:
    public fvBlockMatrix<Type>
{

private:
            //- Insert coupling matrix diag and source into this fvBlockMatrix
            void insertCouplingDiagSource
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );

            //- Insert coupling matrix lower and upper into this fvBlockMatrix
            void insertCouplingUpperLower
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );


public:

    // Constructors

        //- Construct given a field to solve for
        fvBlockMatrixFix(GeometricField<Type, fvPatchField, volMesh>&);

        //- Construct as copy
        fvBlockMatrixFix(const fvBlockMatrixFix<Type>&);




            void insertEquationCoupling
            (
                const direction dirI,
                const direction dirJ,
                fvScalarMatrix& matrix
            );

//            void insertVectorEquationCoupling
//            (
//                const direction dirI,
//                const direction dirJ,
//                fvVectorMatrix& matrix,
//                const int32_t component = -1
//            );


    //- Destructor
    virtual ~fvBlockMatrixFix()
    {};
};

}

//****

#ifdef NoRepository
    #include "fvBlockMatrixFix.C"
#endif

#endif