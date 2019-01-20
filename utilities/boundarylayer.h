#ifndef boundarylayer_H
#define boundarylayer_H

#include "fvCFD.H"

class boundaryLayer
{
private:

    word patch_name_;
    fvMesh mesh_;
    volVectorField U_;
    scalar tolerance_;
    int n_steps_;
    label patch_ID_;
    int n_patch_faces_;
    List<point>& patch_face_centres_;
    List<point>& patch_face_nvectors_;

public:

    // Constructors
    boundaryLayer
    (
        const word;
        const fvMesh&;
        const volVectorField&;
    );

    boundaryLayer
    (
        const word;
        const fvMesh&;
        const volVectorField&;
        const int;

    );

    // Member functions

    inline label get() const {return tolarance;};

    inline void set(int ud_n_steps_) {n_steps_ = ud_n_steps_;};

    virtual void extractPatchData();



};

#endif
