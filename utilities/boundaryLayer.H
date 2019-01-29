#ifndef boundarylayer_H
#define boundarylayer_H

#include "fvCFD.H"
#include "fvMesh.H"

class boundaryLayer
{
private:
    fvMesh& mesh_;
    word patch_name_;
    volVectorField U_;
    scalar tolerance_;
    int n_steps_;
    label patch_ID_;
    int n_patch_faces_;
    List<point> patch_face_centres_;
    List<point> patch_face_nvectors_;

public:

    // Constructors
    boundaryLayer
    (
        fvMesh& mesh,
        word patch_name
    )

    : mesh_(mesh) {patch_name_ = patch_name}

    ~boundaryLayer();
//    boundaryLayer
//    (
//        const word;
//        const fvMesh&;
//        const volVectorField&;
//        const int;

//    );

    // Member functions

//    inline label get() const {return tolarance;};
//    inline void set(label patch_ID) { patch_ID_ = patch_ID;};
 //   inline void set(List<point>& patch_face_centres) { patch_face_centres_ = patch_face_centres;};

//inline void set(List<point>& patch_face_nvectors_) { patch_face_nvectors_ = patch_face_nvectors;};

virtual void extractPatchData()
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



};


};

#endif
