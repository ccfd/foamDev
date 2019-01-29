#ifndef BOUNDARYLAYERCONTAINER_H
#define BOUNDARYLAYERCONTAINER_H

#endif // BOUNDARYLAYERCONTAINER_H

class boundaryLayerContainer
{
private:
    word patch_name_;
    label patch_ID_;
    int n_patch_faces_;
    List<point> patch_face_centres_;
    List<point> patch_face_nvectors_;

public:

    boundaryLayerContainer
    (
        word patch_name
    )

    : mesh_(mesh) {patch_name_ = patch_name}

    ~boundaryLayer();










}
