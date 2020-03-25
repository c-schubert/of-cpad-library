/* 
cpa.H

License (GNU General Public License version 3 or later)
Copyright (c) 2019 Christian Schubert
*/
#ifndef CPA_H
#define CPA_H


#include "bool.H"
#include "error.H"
#include "fvCFD.H"

namespace FoamCpad
{

class cpa 
{
private:
    // Private data
        Foam::labelList cells;
        Foam::label noCells;
        Foam::vector com;
        Foam::scalar mass;
        Foam::scalar volume;
        Foam::scalar cellAveragedAlpha;
        Foam::scalar maxAlpha;
        bool hasCellsOnProcessorBoundary;
        Foam::label noProcessBoundaryCells = 0;
        Foam::labelList processorBoundaryFaceIDs;
        Foam::wordList adjacentBoundariesNames;

    // Private member functions
        void getProperties
        (
            Foam::volScalarField &alpha,
            const Foam::dimensionedScalar& rho
        );

        void getPropertiesThermo
        (
            Foam::volScalarField &alpha,
            Foam::volScalarField &rho
        );

        void getCpaBoundaryList
        (
            Foam::volScalarField& alpha
        );

        void cpaTreeSearch
        (
            Foam::label initialCell,
            Foam::volScalarField& alpha,
            Foam::scalar alphaMin,
            Foam::boolList& isCpadCell,
            Foam::boolList& isCellChecked,
            Foam::labelList& cpaCells,
            Foam::label cellConDectMode
        );


protected:


public:
    // Constructors
        cpa
        (
            Foam::label initialCell,
            Foam::volScalarField& alpha,
            Foam::scalar alphaMin,
            const Foam::dimensionedScalar& rho,
            Foam::boolList& isCpadCell,
            Foam::boolList& isCellChecked,
            Foam::label cellConDectMode
        );

        cpa
        (
            Foam::label initialCell,
            Foam::volScalarField& alpha,
            Foam::scalar alphaMin,
            Foam::volScalarField& rho,
            Foam::boolList& isCpadCell,
            Foam::boolList& isCellChecked,
            Foam::label cellConDectMode
        );

    // Member Functions


        void printProperties();

        Foam::labelList get_cells();
        Foam::label get_noCells();
        Foam::vector get_com();
        Foam::scalar get_mass();
        Foam::scalar get_volume();
        Foam::scalar get_cellAveragedAlpha();
        Foam::scalar get_maxAlpha();
        Foam::label get_noProcessBoundaryCells();
        Foam::labelList get_processorBoundaryFaceIDs();
        Foam::wordList get_adjacentBoundariesNames();

    // Destructor
        ~cpa();
};


Foam::labelList cellFaceNeighbors
(
    Foam::label cell,
    const Foam::fvMesh &mesh
);


Foam::labelList cellEdgeNeighbors
(
    Foam::label cell,
    const Foam::fvMesh &mesh
);

Foam::labelList cellPointNeighbors
(
    Foam::label cell,
    const Foam::fvMesh &mesh
);

}

#endif
