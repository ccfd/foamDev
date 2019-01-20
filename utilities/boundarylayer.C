#include "boundarylayer.h"

// ***** Constructors ***** //

boundaryLayer::boundaryLayer
(
    const word patch_name,
    const fvMesh& mesh,
    const volVectorField& U
)
:
    patch_name_ = patch_name;
    mesh_& = mesh;
    U_& = U;
    tolerance_ = 0.01;
    n_steps_ = 1000;
{}

boundaryLayer::boundaryLayer
(
    const word patch_name,
    const fvMesh& mesh,
    const volVectorField& U,
    const int n_iter
)
:
    patch_name_ = patch_name;
    mesh_& = mesh;
    U_& = U;
    tolerance_ = 0.01;
    n_steps_ = n_iter;
{}

void boundaryLayer::extractPatchData()
{

    Info << "Reading mesh data..." << endl;
    // Evaluating ID of the desired patch
    patch_ID_ = mesh.boundaryMesh().findPatchID(patch_name_);

    // Loading mesh data belonging to the patch
    n_patch_faces_ = mesh.boundary()[patch_ID].Cf().size();

    //const vectorField& patch_faces = mesh.boundary()[patch_ID].Cf();
    List<point>& patch_face_centres_ = mesh.boundaryMesh()[patch_ID].faceCentres();
    List<point>& patch_face_nvectors = mesh.boundary()[patch_ID].Sf();

    Info << "finished" << endl;
    Info << "Patch " << patch_name << " has " << n_patch_faces << " faces." << endl;

}
