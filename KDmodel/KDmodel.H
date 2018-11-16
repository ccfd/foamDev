/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::RASModels::KDmodel

Group
    grpRASTurbulence

Description
    Standard high Reynolds-number k-omega turbulence model for
    incompressible and compressible flows.

    References:
    \verbatim
        Wilcox, D. C. (1998).
        Turbulence modeling for CFD
        (Vol. 2, pp. 103-217). La Canada, CA: DCW industries.
    \endverbatim

    The default model coefficients are
    \verbatim
        KDmodelCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            beta        0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    KDmodel.C

\*---------------------------------------------------------------------------*/

#ifndef KDmodel_H
#define KDmodel_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class KDmodel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class KDmodel
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

	    dimensionedScalar Agamma_;
        dimensionedScalar Cs_;
        dimensionedScalar Ca_;
        dimensionedScalar Ckh_;
        dimensionedScalar Ck_;
	    dimensionedScalar Csep_;
        dimensionedScalar Anu_;
        dimensionedScalar a1_;
        dimensionedScalar a2_;
        dimensionedScalar betaStar_;
        dimensionedScalar alfa_;
        dimensionedScalar beta0_;
        dimensionedScalar sigma_;
        dimensionedScalar sigmaStar_;
        dimensionedScalar sigmaD0_;
        dimensionedScalar Clim_;
	
	
        // Fields
	    volScalarField y_;
        volScalarField k_;
        volScalarField omega_;


    // Declarations of  Protected Member Functions

        virtual void correctNut();
        volScalarField S();
        volScalarField O();
        volScalarField P
        (
            volScalarField& S,
            volScalarField& O
        );

        volScalarField gamma();

        void valueCheck
        (
            volScalarField& fieldValue
        );

        volScalarField Psep
        (
            volScalarField& S
        );

        volScalarField beta();

        volScalarField crossDiff
        (
            volScalarField& beta
        );
        volScalarField gammaDif
        (
            volScalarField& Gamma
        );	



public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;	
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("KDmodel");


    // Constructors

        //- Construct from components
        KDmodel
        (
	    const alphaField& alpha,
	    const rhoField& rho,
            const volVectorField& U,
	    const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~KDmodel()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    sigmaStar_*this->nut_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    sigma_*this->nut_ + this->nu()
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "KDmodel.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
