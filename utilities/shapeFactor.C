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
#include "IOobject.H"
#include "interpolationCellPoint.H"
#include "boundaryLayer.H"
//#include "fvPatchFields.H"
//#include "volFields.H"
//#include "meshTools.H"
//#include "fvOptions.H"
//#include "Vector.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    #include "setRootCase.H" // reads and checks OF folder structure
    #include "createTime.H" // creates a time format based on folder structure
    #include "createMesh.H" // creates mesh based on constant/polyMeshs
    #include "OFstream.H"
    Info << "Reading mesh data..." << endl;
    //const point sample = args.argRead<point>(1);
    // Evaluating ID of the desired patch
    const word patch_name("plate");
    label patch_ID = mesh.boundaryMesh().findPatchID(patch_name);

    // Loading mesh data belonging to the patch
    const polyPatch& patch_data = mesh.boundaryMesh()[patch_ID];
    int n_patch_faces = mesh.boundary()[patch_ID].Cf().size();
    
    //const vectorField& patch_faces = mesh.boundary()[patch_ID].Cf();
    const List<point>& patch_face_centres = mesh.boundaryMesh()[patch_ID].faceCentres();
    const List<point>& patch_face_nvectors = mesh.boundary()[patch_ID].Sf();
    
    Info << "finished" << endl;
    Info << "Patch " << patch_name << " has " << n_patch_faces << " faces." << endl;
//    fvMesh elo = mesh.clone() const;
//    boundaryLayer boundaryLayer(elo, patch_name);

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

    Info << "done." << endl;

    // Find maximum vorticity values

//    List<point> boundary_layer_edge = patch_face_centres;
//    Field<scalar> U_mean = mag(U.boundaryField()[patch_ID]);
//    List<scalar> momentum_thickness;
//    List<scalar> displacement_thickness;
//    scalar eps = 0.01; // scaling factor for rotation rate magnitude tensor
//    scalar dx = 0.001;
//    scalar l = 0.0;
//    vector position(0.0, 0.0, 0.0);
//    vector normal_vector(0.0, 0.0, 0.0);
//    label assigned_cell;
//    scalar vorticity_cell;
//    vector U_cell;
//    scalar U_mag = 0.0;

//    // Iterate
//    while (runTime.loop())
//    {
//        if(runTime.outputTime())
//        {
//            Info << "Time = " << runTime.timeName() << nl << endl;
//            fileName outputFile1("boundary.txt");
//            OFstream os(runTime.timeName()/outputFile1);
//            os << "This is the first line in the file.\n";

//            volTensorField gradU = fvc::grad(U);
//            volTensorField vorticity_tensor(skew(gradU));
//            volScalarField vorticity_mag = pow(2 * vorticity_tensor && vorticity_tensor, 0.5);
//            fvPatchField<double> vorticity_max = vorticity_mag.boundaryField()[patch_ID];
//            interpolationCellPoint<scalar> vorticity_interp(vorticity_mag);
//            interpolationCellPoint<vector> U_interp(U);

//            forAll(mesh.boundary()[patch_ID].Cf(), iter)
//            {
//                bool convergence = false;

//                // evaluating average velocity (assuming compressible flow)

//                while (!convergence)
//                {
//                    position = boundary_layer_edge[iter];
//                    normal_vector = patch_face_nvectors[iter];
//                    position = position - dx * normal_vector;
//                    assigned_cell = mesh.findCell(position);

//                    if (mesh.findCell(position))
//                    {
//                        boundary_layer_edge[iter] = patch_face_centres[iter];
//                        Info << "Application failed to evaluate average velocity" << endl;
//                        convergence = true;
//                        break;
//                    }
                    
                    
//                    U_cell = U_interp.interpolate(position, assigned_cell);
//                    U_mag = U_mag + mag(U_cell) * dx;
//                    l = l + dx;

//                    //vorticity_cell = vorticity_interp.interpolate(position, assigned_cell);
//                    //boundary_layer_edge[iter] = position;

//                    // if (vorticity_cell < eps * vorticity_max[iter])
//                    // {
                        
//                    //     os << boundary_layer_edge[iter];
//                    //     os << endl;
//                    //     convergence = true;
//                    // }


//                }

//                U_mean[iter] = U_mag / l;
//                U_mag = 0;
//                l = 0;
//                Info << U_mean[iter] << endl;

//            }

            // forAll(mesh.boundary()[patch_ID].Cf(), iter)
            // {
            //     bool convergence = false;

            //     while (!convergence)
            //     {
            //         position = boundary_layer_edge[iter];
            //         normal_vector = patch_face_nvectors[iter];
            //         position = position - dx * normal_vector;
            //         assigned_cell = mesh.findCell(position);

            //         if (assigned_cell == -1)
            //         {
            //             boundary_layer_edge[iter] = patch_face_centres[iter];
            //             Info << "Application failed to evaluate boundary layer edge" << endl;
            //             convergence = true;

            //         }

            //         U_cell = U_interp.interpolate(position, assigned_cell);
            //         boundary_layer_edge[iter] = position;

            //         if (vorticity_cell < eps * vorticity_max[iter])
            //         {
                        
            //             os << boundary_layer_edge[iter];
            //             os << endl;
            //             convergence = true;
            //         }


            //     }
                

            // }



            
//        }
//    }

    
    

    Info << "End\n" << endl;


}
