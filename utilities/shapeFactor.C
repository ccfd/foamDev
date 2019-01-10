
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Utility for determining the shape factor of turbomachinery blades.

Author: Robert Crow (r.tykocki.crow@gmail.com)


License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
//#include "fvPatchFields.H"
//#include "volFields.H"
//#include "meshTools.H"
//#include "fvOptions.H"
//#include "Vector.H"

using namespace Foam;

class shapeFactor
{
private:
    //word approximationType_;



public:
    int id;

    shapeFactor()
    {
        Info << "shapeFactor has object has been constructed" << endl;
        id = 1;

    }

    ~shapeFactor()
    {}
};


class boundaryLayer
{
private:
    word patchName;





};



int main(int argc, char *argv[])
{
    #include "setRootCase.H" // reads and checks OF folder structure
    #include "createTime.H" // creates a time format based on folder structure
    #include "createMesh.H" // creates mesh based on constant/polyMeshs

    Info << "Reading mesh data..." << endl;
    

    // Evaluating ID of the desired patch
    const word patch_name("plate");
    double dx = 0.001;
    label patch_ID(0);
    patch_ID = mesh.boundaryMesh().findPatchID(patch_name);

    // Loading mesh data belonging to the patch
    const polyPatch& patch_data = mesh.boundaryMesh()[patch_ID];
    int n_patch_faces = mesh.boundary()[patch_ID].Cf().size();
    
    //const vectorField& patch_faces = mesh.boundary()[patch_ID].Cf();
    const List<point>& patch_face_centres = mesh.boundaryMesh()[patch_ID].faceCentres();
    const List<point>& patch_face_nvectors = mesh.boundary()[patch_ID].Sf();
    
    Info << "finished" << endl;
    Info << "Patch " << patch_name << " has " << n_patch_faces << " faces." << endl;

    // Loading velocity and vorticity field
    Info << "Reading velocity field ..." ;
    volVectorField U // note that velocity is a vector field
	(
		IOobject
		(
		    "U",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

    Info << "...done." << endl;

    Info<< "Reading field vorticity...\n" << endl;
	volVectorField vorticity // note that velocity is a vector field
	(
		IOobject
		(
		    "vorticity",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

    Info << "...done." << endl;
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

		// Loop over all cells in the mesh and calculate the pressure value.
		for (label cellI=0; cellI<mesh.C().size(); cellI++)
		{
			// cellI describes a series of integers, each corresponding to an index of an individual cell in the grid.

			// Call the method and compute p.
			// Note how mesh.C() and p elements are accessed by using the [] operators, as in a regular C array.
			// .value() is also called to convert the time to a dim-less scalar
			p[cellI] = calculatePressure(runTime.time().value(), mesh.C()[cellI], originVector, rFarCell);

            // NOTE: it is also possbile to interact with boundary face values, but
            // this will be addressed in a separate tutorial.
		}
    }

    Info << "End\n" << endl;

}
