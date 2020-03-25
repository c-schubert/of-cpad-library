/* 
cpad.H

License (GNU General Public License version 3 or later)
Copyright (c) 2019 Christian Schubert
*/
#ifndef CPAD_H
#define CPAD_H

#include "cpa.H"
#include <vector>
#include "twoPhaseMixture.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "twoPhaseMixtureThermo.H"

namespace FoamCpad
{
class cpad
{
private:

    // Private data
        Foam::fileName reportFileStr = "cpa";
        bool newReport = "true";

        Foam::label targetAlpha = 1;
        Foam::scalar alphaMin = 0.05;
        Foam::label cellConDectMode = 0;
        bool reduceCpadArea = false;

        std::vector<cpa> cpaList; // cpa muss ein Template sein? Damit das hier funktioniert?
        Foam::boolList isCpadCell;
        Foam::label noCpadCells=0;


    // Private member functions
        Foam::word cpaReportString
        (
            FoamCpad::cpa mycpa, 
            Foam::labelList procFaceToGlobalFaceList
        );
protected:

public:

    // Constructors
        cpad();

        cpad
        (
            Foam::immiscibleIncompressibleTwoPhaseMixture& mixture,
            Foam::Time& runTime
        );

        cpad
        (
            Foam::twoPhaseMixtureThermo& mixture,
            Foam::Time& runTime
        );

        // cpad
        // (
        //     Foam::immiscibleIncompressibleTwoPhaseMixture& mixture,
        //     Foam::IOdictionary
        // );

    // Member Functions
        // void setCpadArea();
        IOdictionary cpadGetDict(Foam::Time& runTime);

        void cpadInit
        (
            IOdictionary& cpadDict,
            Foam::volScalarField& alpha,
            Foam::Time& runTime
        );

        void cpadSearch
        (
            Foam::volScalarField& alpha,
            Foam::scalar alphaMin,
            const Foam::dimensionedScalar& rho,
            Foam::boolList& isCpadCell,
            Foam::label  cellConDectMode,
            std::vector<cpa>& cpaList
        );

        void cpadSearchThermo
        (
            Foam::volScalarField& alpha,
            Foam::scalar alphaMin,
            Foam::volScalarField& rho,
            Foam::boolList& isCpadCell,
            Foam::label cellConDectMode,
            std::vector<cpa>& cpaList
        );

        Foam::label getNoCpas();

        void appendTimeStepReport
        (
            Foam::scalar currentTime,
            Foam::labelList procFaceToGlobalFaceList
        );
        // void writeReport();
    
    // Decstructor
        ~cpad();
};
}

#endif
