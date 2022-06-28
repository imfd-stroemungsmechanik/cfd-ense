{
    const label wallId = mesh.boundaryMesh().findPatchID("WALL");
    const vectorField pNf = mdT.boundaryField()[wallId].patch().nf(); 
    vectorField& mdTwall = mdT.boundaryFieldRef()[wallId];   
    mdTwall -= ((mdT.boundaryField()[wallId] & pNf) * pNf);
}