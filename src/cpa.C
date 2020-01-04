/* 
cpa.C

License (GNU General Public License version 3 or later)
Copyright (c) 2019 Christian Schubert
*/
#include "cpa.H"
#include "error.H"


FoamCpad::cpa::cpa
(
    Foam::label initialCell,
    Foam::volScalarField& alpha,
    Foam::scalar alphaMin,
    const Foam::dimensionedScalar& rho,
    Foam::boolList& isCpadCell,
    Foam::boolList& isCellChecked,
    bool cellNodeNeighborDetectionActive
)
{

    com = {0,0,0};
    mass = 0;
    volume = 0;

    FoamCpad::cpa::cpaTreeSearch
    (
        initialCell,
        alpha,
        alphaMin,
        isCpadCell,
        isCellChecked,
        cells,
        cellNodeNeighborDetectionActive
    );


    FoamCpad::cpa::getProperties
    (
        alpha,
        rho
    );


    FoamCpad::cpa::getCpaBoundaryList
    (
        alpha
    );

    // Foam:word reportStr = FoamCpad::cpa::reportString();

    // Info << reportStr << endl;

    // FoamCpad::cpa::printProperties();
}

void FoamCpad::cpa::cpaTreeSearch
(
    Foam::label initialCell,
    Foam::volScalarField& alpha,
    Foam::scalar alphaMin,
    Foam::boolList& isCpadCell,
    Foam::boolList& isCellChecked,
    Foam::labelList& cpaCells,
    bool cellNodeNeighborDetectionActive
)
{
    const Foam::fvMesh &mesh = alpha.mesh();

    Foam::DynamicList<label> dynamicCpaCellList;
    dynamicCpaCellList.append(initialCell);
    isCellChecked[initialCell] = true;

    Info << "Cpa tree search from inital cell " << initialCell << endl;

    forAll(dynamicCpaCellList, i)
    {
        Foam::label c = dynamicCpaCellList[i]; 

        // Info << "Cecking Cell: " << c << endl;
        Foam::labelList cNeighboarCells;

        /* only face cell neighboars */
        if (!cellNodeNeighborDetectionActive)
        {
            /* All of cell c's - face neighboars */
            cNeighboarCells = mesh.cellCells()[c];
        }
        else
        {
            /* All of cell c's - cell neighboars (including over node point cell 
            neighboars)  */
            cNeighboarCells = FoamCpad::cpa::cellEdgeNeighbors(c, mesh);
        }

        Foam::labelList newCpaCells;

        forAll(cNeighboarCells, j)
        {
            Foam::label cn = cNeighboarCells[j]; // neighboar cell

            if
            (   
                !isCellChecked[cn]
            &&  isCpadCell[cn]
            &&  alpha[cn] > alphaMin
            )
            {
                newCpaCells.append(cn);
            }

            isCellChecked[cn] = true;
        }

        if (newCpaCells.size() > 0)
        {
            if ((dynamicCpaCellList.size() + newCpaCells.size()) > dynamicCpaCellList.capacity())
            {   
                dynamicCpaCellList.setCapacity(dynamicCpaCellList.capacity() + 4000);
            }

            dynamicCpaCellList.append(newCpaCells);
            newCpaCells.clear();
        }
    }

    if (dynamicCpaCellList.size() > 0)
    {  
        dynamicCpaCellList.shrink();

        cpaCells.setSize(dynamicCpaCellList.size());
        cpaCells.transfer(dynamicCpaCellList); // Entweder so oder andersherum...
    }
}

Foam::labelList FoamCpad::cpa::cellEdgeNeighbors
(
    Foam::label cell,
    const Foam::fvMesh &mesh
)
{
    Foam::labelList cEdgeNodes = mesh.cellPoints()[cell];
    Foam::label maxUniqueNeighbors = 40;
    Foam::DynamicList<label> uniqueCellEdgeNeighbors(maxUniqueNeighbors);

    // Info << "Edge nodes: " << cEdgeNodes.size() << endl;
   
     // musst be greater or equal than max possible length...

    forAll(cEdgeNodes, j)
    {
        Foam::labelList cNeighboarCells = mesh.pointCells()[cEdgeNodes[j]];

        // Info << "Edge node Cells: " << cNeighboarCells.size() << endl;
        
        forAll(cNeighboarCells, i)
        {
            if(uniqueCellEdgeNeighbors.size() < 1)
            {
                uniqueCellEdgeNeighbors.append(cNeighboarCells[i]);
            }
            else
            {
                bool CellIdIsUnique = true;

                forAll(uniqueCellEdgeNeighbors, k)
                {
                    if (uniqueCellEdgeNeighbors[k] == cNeighboarCells[i])
                    {
                        CellIdIsUnique = false;
                    }
                }

                if (CellIdIsUnique)
                {
                    if (uniqueCellEdgeNeighbors.size() >= maxUniqueNeighbors)
                    {   

                        Foam::error("Overflow of unique neighbor cell, "
                            "this should not happen!"
                            " Aborting ...");

                        Info<< "Error" << endl;
                    }
                    else
                    {
                        uniqueCellEdgeNeighbors.append(cNeighboarCells[i]);
                    }
                }
            }
        }
    }

    if (uniqueCellEdgeNeighbors.size() > 0)
    {
        uniqueCellEdgeNeighbors.shrink();
    }

    Foam::labelList cellEdgeNeighbors(uniqueCellEdgeNeighbors);

    return cellEdgeNeighbors;
}



void FoamCpad::cpa::getProperties
(
    Foam::volScalarField &alpha,
    const Foam::dimensionedScalar& rho
)
{
    Foam::scalar xG = 0;
    Foam::scalar yG = 0;
    Foam::scalar zG = 0;
    Foam::scalar alphaSum = 0;
    maxAlpha = 0;
    mass = 0;
    volume = 0;
    Foam::scalar denominator_sum = 0;
    Foam::scalar RHO = rho.value();
    const Foam::fvMesh &mesh = alpha.mesh();

    forAll(cells, i)
    {
        Foam::label c = cells[i];

        volume += alpha[c]*mesh.V()[c];
        mass += RHO*alpha[c]*mesh.V()[c];
        alphaSum += alpha[c];

        scalar denominator = RHO*alpha[c]*mesh.V()[c];
        denominator_sum += denominator;

        xG += mesh.cellCentres()[c].x() * denominator;
        yG += mesh.cellCentres()[c].y() * denominator;
        zG += mesh.cellCentres()[c].z() * denominator;

        maxAlpha = std::max(alpha[c], maxAlpha);
    }

    com.x() = xG/denominator_sum;
    com.y() = yG/denominator_sum;
    com.z() = zG/denominator_sum;

    noCells = cells.size();

    if (noCells >  0)
    {
        cellAveragedAlpha = alphaSum/Foam::scalar(noCells);
    }
}


void FoamCpad::cpa::getCpaBoundaryList
(
    Foam::volScalarField& alpha
)
{

    // -----------------------------------------------------------------------
    // TODO: In eigenen Datentyp auslagern einmal pro Case

    const Foam::fvMesh &mesh = alpha.mesh();
    const Foam::fvBoundaryMesh &bMesh = mesh.boundary();
    const Foam::polyBoundaryMesh &pbMesh = mesh.boundaryMesh();

    Foam::boolList isBoundaryCell(alpha.size(), false);
    Foam::boolList isProcessBoundaryCell(alpha.size(), false);
    Foam::labelListList procBoundaryFaceID;
    Foam::labelList boundaryID(alpha.size(), -1);
    Foam::label startProcFaceListId;
    Foam::label counter = 0;

    procBoundaryFaceID.setSize(alpha.size());

    forAll(pbMesh, i)
    {
        
        Foam::polyPatch boundaryPolyPatch = pbMesh[i];
        Foam::word BCtype(pbMesh.types()[i]);

        Foam::labelListList boundaryPatchFaceList 
            = boundaryPolyPatch.faceFaces();

        Foam::labelList boundaryPatchCellList = bMesh[i].faceCells();
        // boundaryPatchCellList.size() = nfaces in processor

        startProcFaceListId = boundaryPolyPatch.start();

        forAll(boundaryPatchCellList, j)
        {
            Foam::label c = boundaryPatchCellList[j];
            Foam::label pProcFaceIdx = startProcFaceListId + j;
            isBoundaryCell[c] = true;
            boundaryID[c] = i;

            if (BCtype == "processor")
            {
                isProcessBoundaryCell[c] = true;
                procBoundaryFaceID[c].append(pProcFaceIdx);
                counter++;
            }
        }
    }

    // Info<<"Cells near processor boundary " << counter << endl;

    /*TODO Move out of Object for better perfomance or use other function 
    (there realy should be a more elegant way ...*/
    // Foam::labelIOList localFaceProcAddr
    // (
    //     IOobject
    //     (
    //         "faceProcAddressing",
    //         mesh.facesInstance(),
    //         mesh.meshSubDir,
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE
    //     )
    // );

    noProcessBoundaryCells = 0;

    // ------------------------------------------------------------------------

    forAll(cells, i)
    {
        Foam::label c = cells[i];

        if (isBoundaryCell[c])
        {
            Foam::word bNameCell = bMesh[boundaryID[c]].name();
            bool isBoundaryUnique = true;

            if (isProcessBoundaryCell[c])
            {
                // is MPI boundary
                processorBoundaryFaceIDs.append
                (
                    // localFaceProcAddr[procBoundaryFaceID[c]]
                    procBoundaryFaceID[c]
                );

                noProcessBoundaryCells += procBoundaryFaceID[c].size();
            }

            if (adjacentBoundariesNames.size() > 0)
            {
                forAll(adjacentBoundariesNames,j)
                {
                    if (adjacentBoundariesNames[j] == bNameCell)
                    {
                        isBoundaryUnique = false;
                    }
                }
            }

            if (isBoundaryUnique)
            {
                adjacentBoundariesNames.append(bNameCell);
            }
        }
    }

}

void FoamCpad::cpa::printProperties()
{
    Info<< "\t Cpa Properties:" << endl;
    Info<<"\t\tMass: " << mass << endl;
    Info<<"\t\tVolume: " << volume << endl;
    Info<<"\t\tCOM: " << com << endl;
    Info<<"\t\tAlpha Mean (over cells): " << cellAveragedAlpha << endl;
    Info<<"\t\tAlpha Max: " << maxAlpha << endl;
    Info<<"\t\tAdjacent Boundaries: " << adjacentBoundariesNames << endl;
    Info<<"\t\tProcessor boundary cells: " << processorBoundaryFaceIDs << endl;
}


Foam::labelList FoamCpad::cpa::get_cells()
{
    return cells;
}
Foam::label FoamCpad::cpa::get_noCells()
{
    return noCells;
}
Foam::vector FoamCpad::cpa::get_com()
{
    return com;
}
Foam::scalar FoamCpad::cpa::get_mass()
{
    return mass;
}
Foam::scalar FoamCpad::cpa::get_volume()
{
    return volume;
}
Foam::scalar FoamCpad::cpa::get_cellAveragedAlpha()
{
    return cellAveragedAlpha;
}
Foam::scalar FoamCpad::cpa::get_maxAlpha()
{
    return maxAlpha;
}
Foam::label FoamCpad::cpa::get_noProcessBoundaryCells()
{
    return noProcessBoundaryCells;
}
Foam::labelList FoamCpad::cpa::get_processorBoundaryFaceIDs()
{
    return processorBoundaryFaceIDs;
}
Foam::wordList FoamCpad::cpa::get_adjacentBoundariesNames()
{
    return adjacentBoundariesNames;
}


FoamCpad::cpa::~cpa()
{
    // processorBoundaryFaceIDs.~List();
    // adjacentBoundariesNames.~List();
}
