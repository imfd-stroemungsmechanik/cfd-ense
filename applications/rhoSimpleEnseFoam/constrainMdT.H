// Constrain face normal component of mdT on wall boundaries for corresponding matrix
template <class Type>
tmp<fvMatrix<Type>> constrainmdT(tmp<fvMatrix<Type>> tmatrix)
{
    fvMatrix<Type>& matrix = tmatrix.ref();

    const polyBoundaryMesh& patches = matrix.psi().mesh().boundaryMesh();

    forAll(patches, patchi)
    {
        if(patches[patchi].type() == "wall")
        {
            const label size = matrix.psi().boundaryField()[patchi].size();
            matrix.internalCoeffs()[patchi] = Field<Type>(size, pTraits<Type>::zero);
            matrix.boundaryCoeffs()[patchi] = Field<Type>(size, pTraits<Type>::zero);
        }
    }

    return tmatrix;
}

// Constrain face normal component of mdT on wall boundaries for mdT field
tmp<volVectorField> constrainmdT(tmp<volVectorField> tmdT)
{
    volVectorField& mdT = tmdT.ref();

    const polyBoundaryMesh& patches = mdT.mesh().boundaryMesh();

    forAll(patches, patchi)
    {
        if(patches[patchi].type() == "wall")
        {
            const fvPatchVectorField& patchField = mdT.boundaryField()[patchi];

            const vectorField pNf = patchField.patch().nf();
            vectorField& mdTwall = mdT.boundaryFieldRef()[patchi];

            mdTwall -= ((patchField & pNf) * pNf);
        }
    }

    return tmdT;
}
