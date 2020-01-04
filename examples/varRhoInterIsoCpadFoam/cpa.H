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
            bool cellNodeNeighborDetectionActive
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
            bool cellNodeNeighborDetectionActive
        );

    // Member Functions


        void printProperties();

        Foam::labelList cellEdgeNeighbors
        (
            Foam::label cell,
            const Foam::fvMesh &mesh
        );

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
}
#endif
