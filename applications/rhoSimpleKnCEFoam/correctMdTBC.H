{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        if(patches[patchi].type() == "wall")
        {
            const vectorField pNf = mdT.boundaryField()[patchi].patch().nf(); 
            vectorField& mdTwall = mdT.boundaryFieldRef()[patchi];

            mdTwall -= ((mdT.boundaryField()[patchi] & pNf) * pNf);        
        }
    }    
}