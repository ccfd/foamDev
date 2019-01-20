#ifndef SHAPEFACTOR_H
#define SHAPEFACTOR_H
#include "fvCFD.H"

class shapeFactor
{
private:
    List<point> boundary_layer_edge_;
    Field<scalar> U_mag_;
    List<scalar> momentum_thickness_;
    List<scalar> displacement_thickness_;
    List<scalar> shape_factor_;



public:
    shapeFactor();
};

#endif // SHAPEFACTOR_H
