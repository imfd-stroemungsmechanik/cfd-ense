// Constrain face normal component of mdT to zero on wall boundaries
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

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
}