/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
           dropletDetection | Copyright (C) 2018 IOB RWTH Aachen university
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your optio 
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    dropletDetection

Group
    

Description
    WIP:
    This code shall enable the postum detection (volume, centroid) of segregated
    droplets in two phase scalar fields.

    add: monitor mass cells detatched to electrode
    add: monitor mass cells detatched to outlet


    dropletDetection is developed by Christian, Schubert (2018)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


/* solver specific loads */
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "dynamicFvMesh.H"
#include "subCycle.H"
#include "fvOptions.H"
#include "cpad.H"

/*Ersetzen durch Detektion ob AMR benutzt wird und ob das Gitter sich geändert hat*/
#define AMR 0

void printVector(Foam::volScalarField &alpha1, label i)
{
       Info<<alpha1[i]<<endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    
    Info<< "-----------------------------------------------------------------"
        << endl;
    Info<< "Starting continous phase area detection for processed cases:" 
        << endl;
    Info<< "\tCase directory: " << runTime.path() << endl;
    Info<< "\tCase name: " << runTime.caseName() << endl;
    Info<< "-----------------------------------------------------------------\n"
        << endl;
    Foam::instantList availTime(runTime.times());

    {
        forAll(availTime, t)
        {
            runTime.setTime(availTime[t].value(), 1);

            #include "createDynamicFvMesh.H"
            #include "createFields.H"
            #include "getLocalToGlobalFaceArray.h"

            FoamCpad::cpad cpas(mixture, runTime);
            cpas.appendTimeStepReport
            (
                runTime.value(), 
                procFaceToGlobalFaceList
            );
        }
    }

    Info<< "Cpa detection finished\n" << endl;
    return 0;
}


// ************************************************************************* //
