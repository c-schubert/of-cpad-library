/* 
cpad.C

License (GNU General Public License version 3 or later)
Copyright (c) 2019 Christian Schubert
*/


#include "cpad.H"
#include "IOstreamOption.H"

// Constructor Functions: 

FoamCpad::cpad::cpad
(
    Foam::immiscibleIncompressibleTwoPhaseMixture& mixture,
    Foam::Time& runTime
)
{
    IOdictionary cpadDict
    (
        IOobject
        (
            "cpad",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    alphaMin = readScalar(cpadDict.lookup("alphaMinThreshold"));
    reportFileStr = word(cpadDict.lookup("reportFileName"));
    targetAlpha =  readLabel(cpadDict.lookup("targetPhase"));
    Foam::word cellConnectionModeStr = word(cpadDict.lookup("detectOver"));


    Info << cellConnectionModeStr << endl;
    if (Foam::Pstream::parRun())
    {
        Info << "cpa detection for parallel run:" << endl;
        Info << cellConnectionModeStr << endl;
        // There are additional boundaries called processor where the exchange between mpi process is handeled ...
        reportFileStr = Foam::word(std::to_string(Foam::Pstream::myProcNo())) 
                        + "_" + reportFileStr;
    }

    if (cellConnectionModeStr == "faces")
    {
        cellConDectMode = 0;
    }
    else if (cellConnectionModeStr == "edges")
    {
        cellConDectMode = 1;
    }
    else if (cellConnectionModeStr == "points")
    {
        cellConDectMode = 2;
    }
    else
    {
        Info<< "Warning: No valid setting for detectOver setting given. "
            << "Falling back to face detection mode" << endl;
        cellConDectMode = 0;
    }

     Foam::volScalarField& alpha = (targetAlpha == 2) ? mixture.alpha2() : mixture.alpha1();

     const Foam::dimensionedScalar& rho = (targetAlpha == 2) ? mixture.rho2() : mixture.rho1();
    // default value constructor

    Info<<"\n----------- CPAD -----------" << endl;
    Info<< "Starting continous phase detection for phase " << targetAlpha 
        << " with alphaMin " << alphaMin << " using ";

    if (cellConDectMode == 0)
    {
        Info<<"cell face detection ";
    }
    else if (cellConDectMode == 1)
    {
        Info<<"cell edge detection ";
    }
    else if (cellConDectMode == 2)
    {
        Info<<"cell point detection ";
    }

    Info<<"at time " << runTime.value() << endl;

    if (!(targetAlpha == 1 || targetAlpha == 2))
    {
        Info<< "Warning: Value of targetPhase not corret (musst be 1 or 2). "
            << " setting targetPhase to 1 ...!" << endl;
    }

    Foam::boolList isCpadCellTemp(alpha.mesh().nCells(), false);
    isCpadCell = isCpadCellTemp;

    // Limit cells to check on range above metal bath
    forAll(alpha, c)
    {
        if (reduceCpadArea)
        {
            // TODO: Implentent choosing of cpa area (similar to setFields)...
        }
        else
        {
            noCpadCells++;
            isCpadCell[c] = true;
        }
    }

    FoamCpad::cpad::cpadSearch
    (
        alpha,
        alphaMin,
        rho,
        isCpadCell,
        cellConDectMode,
        cpaList
    );
}


void FoamCpad::cpad::cpadSearch
(
    Foam::volScalarField& alpha,
    Foam::scalar alphaMin,
    const Foam::dimensionedScalar& rho,
    Foam::boolList& isCpadCell,
    Foam::label cellConDectMode,
    std::vector<cpa>& cpaList
)
{
    Foam::boolList isCellChecked(alpha.size(), false);

    forAll(alpha, c)
    {
        if ((isCpadCell[c] == true) && !isCellChecked[c])
        {
            if(alpha[c] > alphaMin)
            {
                FoamCpad::cpa contPhaseArea
                (
                    c,
                    alpha,
                    alphaMin,
                    rho,
                    isCpadCell,
                    isCellChecked,
                    cellConDectMode
                );
                cpaList.push_back(contPhaseArea);
            }
            else
            {
                isCellChecked[c] = true;
            }
        }
    }

    Info<< "Found: " << cpaList.size() << " cont. phase areas!\n\n" << endl;
}


Foam::label FoamCpad::cpad::getNoCpas()
{
    return cpaList.size();
}


void FoamCpad::cpad::appendTimeStepReport
(
    Foam::scalar currentTime, 
    Foam::labelList procFaceToGlobalFaceList
)
{   
    bool append = true;

    Foam::IOstreamOption reportFsSettings;
    Foam::OFstream fs(reportFileStr, IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED, append);


    fs<<"{"<< endl;
    fs<<currentTime<<endl;
    fs<<"cell_count[com.(x)_in_m,com.(y),com.(z),mass_in_kg,volume_in_m³,"
    "cellAveragedAlpha,maxAlpha,[list_of_boundary_names],"
    "[list_of_processor_cells]]" << endl;

    for(auto &cpa : cpaList)
    {
        fs<<cpaReportString(cpa, procFaceToGlobalFaceList) << endl;
    }

    fs<<"}"<<endl;
}



Foam::word FoamCpad::cpad::cpaReportString
(
    FoamCpad::cpa mycpa, 
    Foam::labelList procFaceToGlobalFaceList
)
{
    Foam::word reportStr =  Foam::name(mycpa.get_noCells()) + "[" 
                            + Foam::name(mycpa.get_com().x()) + ","
                            + Foam::name(mycpa.get_com().y()) + ","
                            + Foam::name(mycpa.get_com().z()) + ","
                            + Foam::name(mycpa.get_mass()) + ","
                            + Foam::name(mycpa.get_volume()) + ","
                            + Foam::name(mycpa.get_cellAveragedAlpha()) + ","
                            + Foam::name(mycpa.get_maxAlpha()) + ",";

    reportStr += "[";

     Foam::wordList adjacentBoundariesNames 
        = mycpa.get_adjacentBoundariesNames();

    forAll(adjacentBoundariesNames, i)
    {   
        if (i < (adjacentBoundariesNames.size()-1))
        {
            reportStr += adjacentBoundariesNames[i] + ",";
        }
        else
        {
            reportStr += adjacentBoundariesNames[i];
        }
    }

    reportStr += "]";
    reportStr += ",";
    reportStr += Foam::name(mycpa.get_noProcessBoundaryCells()) + "[";

    Foam::labelList processorBoundaryFaceIDs 
        = mycpa.get_processorBoundaryFaceIDs();

    forAll(processorBoundaryFaceIDs, i)
    {
        Foam::label faceID;

        if (Foam::Pstream::parRun())
        {
            
            faceID = procFaceToGlobalFaceList[processorBoundaryFaceIDs[i]];
        }
        else
        {

            faceID = processorBoundaryFaceIDs[i];
        }

        if (i < (processorBoundaryFaceIDs.size()-1))
        {
            reportStr += Foam::name(faceID) + ",";
        }
        else
        {
            reportStr += Foam::name(faceID);
        }
    }
    reportStr += "]";
    reportStr += "]";


    return reportStr;
}


FoamCpad::cpad::~cpad()
{
    // cpaList.clear();
}

