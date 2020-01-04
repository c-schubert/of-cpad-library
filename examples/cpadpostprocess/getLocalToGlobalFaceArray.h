Foam::labelList procFaceToGlobalFaceList; 

if (Foam::Pstream::parRun())
{
    Foam::labelIOList localFaceProcAddr
    (
        IOobject
        (
            "faceProcAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

        procFaceToGlobalFaceList = localFaceProcAddr;
}


