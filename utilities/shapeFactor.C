
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

    int meshSize_ = mesh.C().size;







    shapeFactor obj1;
    Info << "the id is:" << obj1.id << endl;


}
