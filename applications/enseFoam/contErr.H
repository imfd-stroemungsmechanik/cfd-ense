{

massRes = fvc::div(phi) + fvc::div(phiMd); 

scalar sumLocalmassRes = runTime.deltaTValue()*mag(massRes)().weightedAverage(mesh.V()).value();
scalar globalmassRes = runTime.deltaTValue()*massRes.weightedAverage(mesh.V()).value();

cumulativeContErr += globalmassRes;
Info<< "time step continuity errors : sum = " << gSum(massRes)
 << ", sum local = " << sumLocalmassRes
 << ", global = " << globalmassRes
 << ", cumulative = " <<  cumulativeContErr
 << endl;

}      