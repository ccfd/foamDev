

#include "H12.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(H12, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        H12,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::H12::calc()
{
    if (foundObject<volVectorField>(fieldName_))
    {   
        // get U field and define interpolation on it
        const volVectorField& U = lookupObject<volVectorField>(fieldName_);
        interpolationCellPoint<vector> interpU = 
            interpolationCellPoint<vector>(U); 

        //define fields to be written
        tmp<volScalarField> deltaStar = scalar(0)*mag(U);
        tmp<volScalarField> theta = scalar(0)*mag(U);
        tmp<volScalarField> h12 = scalar(0)*mag(U);

        //extract boundary references
        label patchID = mesh_.boundaryMesh().findPatchID(patchName_);
        fvPatchScalarField& deltaStarBoundary = 
            deltaStar->boundaryFieldRef()[patchID];
        fvPatchScalarField& thetaBoundary = theta->boundaryFieldRef()[patchID];
        fvPatchScalarField& h12Boundary = h12->boundaryFieldRef()[patchID];
        
        
        forAll(normals_,  face)
        {   
            //get nuber of sampling points for this face
            label nSamplePoints = samplePoints_[face].size();

            // get local (for this face) sampling cell and ponts references
            const List<point>& sPoints_l = samplePoints_[face];
            const List<label>& sCells_l = sampleCells_[face];

            // get information about U on the sampling line
            vectorField U_values(nSamplePoints);
            
            vector U_local;
            forAll(sCells_l, samplePointI)
            {
                U_local = interpU.interpolate
                (
                    sPoints_l[samplePointI], sCells_l[samplePointI]
                );
                U_values[samplePointI] = 
                (
                    U_local - ((U_local & normals_[face]) * normals_[face])
                );
            }
            vector U_max = Foam::max(U_values);

            // compute densities (intagtated functions) of momentum and
            // displacement thickness 
            scalarField fraction = 
                  (U_values & U_max)
                / ( (U_max & U_max) + SMALL );

            scalarField deltaStar_d = scalar(1.0) - fraction;
            scalarField theta_d = (scalar(1.0) - fraction) * fraction;

            // integrate along sampling line
            for 
            (
                label samplePointI = 1; 
                samplePointI < (sCells_l).size(); 
                samplePointI++
            )
            {
                deltaStarBoundary[face] += trapz(
                    deltaStar_d[samplePointI-1],
                    deltaStar_d[samplePointI],
                    dy_
                );
                thetaBoundary[face] += trapz(
                    theta_d[samplePointI-1],
                    theta_d[samplePointI],
                    dy_
                );      
            }

            h12Boundary[face] = 
                  deltaStarBoundary[face] 
                / (thetaBoundary[face] + SMALL);

        } 


        if (write_deltaStar_) store(fieldName_deltaStar_, deltaStar);
        if (write_theta_)     store(fieldName_theta_, theta);

        return store(resultName_, h12);
    }
    else
    {
        return false;
    }
}

bool Foam::functionObjects::H12::findSamples
(
    const label face,
    meshSearch& searchEngine,
    label startingCell
)
{
    // wall distance
    scalar y = 0;
    vector n = normals_[face];
    vector p = faceCenters_[face];


    // find current cell
    label currentCell = searchEngine.findCell(p, startingCell);
    if (currentCell == -1 || currentCell != startingCell){
        FatalErrorInFunction
        << "Something went wrong! First cell was not found!"
        << abort(FatalError);
    }

    //init cell data with first value
    List<label> cellist = {currentCell};
    List<vector> pointList = {p};

    while(y <= limit_[face])
    {   
        // move in the wall normal direction
        p += n*dy_;
        y += dy_;

        //find cell
        currentCell = searchEngine.findCell(p, currentCell);
        if (currentCell == -1) break;

        cellist.append(currentCell);
        pointList.append(p);
    }

    samplePoints_[face] = pointList;
    sampleCells_[face] = cellist;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::H12::H12
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression
    (
        name,
        runTime,
        dict,
        dict.lookupOrDefault<word>("field","U")
    ),
    patchName_(dict.lookup("patch")),
    limitType_(dict.lookupOrDefault<word>("limitType","const")),
    dy_(dict.lookupType<scalar>("dy"))
{
    setResultName(typeName, fieldName_);

    Info << "H12 -> Vector field for shape factor calculation: " << fieldName_ << endl;
    
    Info << "H12 -> Selecting patch: " << patchName_ << endl;
    label patchID = mesh_.boundaryMesh().findPatchID(patchName_);

    if (patchID == -1)
    {
        FatalErrorInFunction
            << "Patch <" << patchName_ << "> not found."
            << abort(FatalError);
    }

    label nFaces = mesh_.boundary()[patchID].Cf().size();
    if (nFaces == 0){
        FatalErrorInFunction
        << "Something went wrong! Patch "<< patchName_ <<" has zero faces!"
        << abort(FatalError);
    }

    Info << "H12 -> Calculating normals... ";
    normals_ = scalar(-1) *
                mesh_.boundary()[patchID].Sf() / 
                mesh_.boundary()[patchID].magSf();
    Info << "Done." << endl;

    Info << "H12 -> Reading face centers... ";
    faceCenters_ = mesh_.boundary()[patchID].Cf();
    Info << "Done." << endl;

    if (limitType_=="const")
    {
        Info << "H12 -> Selecting constant limit" << endl;
        limit_ = 
              scalar( dict.subDict("constCoeff").lookupType<scalar>("limit"))
            * mag(normals_);
    }
    else if (limitType_=="linFuncX")
    {   
        Info << "H12 -> Selecting limit based on linear function of x " 
             << "coordinate of face centers" << endl;
        scalar A = dict.subDict("linFuncXCoeff").lookupType<scalar>("A");
        scalar B = dict.subDict("linFuncXCoeff").lookupType<scalar>("B");
        limit_ = A*faceCenters_.component(0)+ B;
            
    }
    else
    {
        FatalErrorInFunction
            << "Unknown limitType: "<< limitType_
            << abort(FatalError);
    }

    Info << "H12 -> Integration step dy = " << dy_ << endl;
    Info << "H12 -> Selecting auxiliary fields for writing:\n";
    M_CheckIfShouldWrite(deltaStar)
    M_CheckIfShouldWrite(theta)
    M_CheckIfShouldWrite(delta)
    
    meshSearch searchEngine(mesh_);
    samplePoints_ = List<List<point>>(nFaces);
    sampleCells_ = List<List<label>>(nFaces);

    Info << "H12 -> Generating sample points... ";
    forAll(mesh_.boundary()[patchID], faceI)
    {
        
        label faceGlobalIndex = faceI + mesh_.boundaryMesh()[patchID].start();
        label adjacentCell = mesh_.owner()[faceGlobalIndex];

        findSamples(faceI,searchEngine, adjacentCell);
    }
    Info << "Done." << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::H12::~H12()
{}


// ************************************************************************* //
