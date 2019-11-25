/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
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

Description
    Macros for BaseAlgKOmega model.

\*---------------------------------------------------------------------------*/

#ifndef BASEALGKOMEGAM_H
#define BASEALGKOMEGAM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Macro helping in volScalarField initialisation in model initialisation list
// Variables: 
//      word name - name of the field (e.g. "psi")
//      dimensionSet diemsnion - physical dimensions of the field 
//      scalar value - initial value
#define M_scalarFieldInit(name, dimension, value)                         \
    IOobject                                                              \
    (                                                                     \
        IOobject::groupName(name, U.group()),                             \
        this->runTime_.timeName(),                                        \
        this->mesh_,                                                      \
        IOobject::NO_READ,                                                \
        IOobject::NO_WRITE                                                \
    ),                                                                    \
    this->mesh_,                                                          \
    dimensionedScalar                                                     \
    (                                                                     \
        IOobject::groupName(name, U.group()),                             \
        dimension,                                                        \
        value                                                             \
    ),                                                                    \
    "calculated"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif //BASEALGKOMEGAM_H

// ************************************************************************* //

