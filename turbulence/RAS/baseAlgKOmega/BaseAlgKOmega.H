/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
  Updated implementation of algebraic intermittency model of 
  laminar-turbulent transition. Model is based on Wilcox 2008
  k-omega. Derivation of the model is described in the following
  papers:
  
  \verbatim
    Kubacki, Slawomir & Dick, Erik. (2016). 
    An algebraic intermittency model for bypass, separation-induced 
    and wake-induced transition. International Journal of Heat and Fluid 
    Flow. 62. 344-361. 10.1016/j.ijheatfluidflow.2016.09.013. 
  \endverbatim
  \verbatim
    Kubacki, Slawomir & Dick, Erik. (2016). 
    An algebraic model for prediction of bypass and separation-induced 
    transition in turbomachinery boundary layer flows. In D. Termini, 
    T. Tucciarelli, & E. Napoli (Eds.), 11th International ERCOFTAC 
    Symposium on Engineering Turbulence Modelling and Measurements (pp. 1–6). 
  \endverbatim

Authors:
    Robert Tykocki-Crow (r.tykocki.crow@gmail.com)
    Wojciech Sadowski   (wojciech1sadowski@gmail.com)

SourceFiles
    BaseAlgKOmega.C

SeeAlso
    BaseAlgKOmega.H
    BaseAlgKOmegaM.C
\*---------------------------------------------------------------------------*/


#ifndef BASEALGKOMEGA_H
#define BASEALGKOMEGA_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "BaseAlgKOmegaM.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class BaseAlgKOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class BaseAlgKOmega
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients
	    dimensionedScalar c_Agamma_;
        dimensionedScalar c_Cs_;
        dimensionedScalar c_CA_;
        dimensionedScalar c_Ckh_;
        dimensionedScalar c_Ck_;
	    dimensionedScalar c_Csep_;
        dimensionedScalar c_Av_;
        dimensionedScalar c_a1_;
        dimensionedScalar c_a2_;
        dimensionedScalar c_betaStar_;
        dimensionedScalar c_alfa_;
        dimensionedScalar c_beta0_;
        dimensionedScalar c_sigma_;
        dimensionedScalar c_sigmaStar_;
        dimensionedScalar c_sigmaD0_;
        dimensionedScalar c_Clim_;
	
        // wall distance field
        volScalarField y_;
	
        // Fields - modelled - always written
        volScalarField omega_;
        volScalarField k_;

	    // Fields - auxiliary - optionaly written
        volScalarField gamma_;
        volScalarField nu_s_;
        volScalarField nu_l_;
        volScalarField k_s_;
        volScalarField k_l_;

        volScalarField f_ss_;
        volScalarField C_ss_;
        volScalarField P_k_;
        volScalarField P_sep_;
        volScalarField F_sep_;
        volScalarField f_k_;
        volScalarField chi_;
        volScalarField Re_v_;
        volScalarField Re_y_;
        volScalarField f_d_;

    // Declarations of  Protected Member Functions

        void valueCheck(volScalarField& fieldValue);

        virtual void correctNut();

        volScalarField ShearRateMag();

        volScalarField VorticityMag();
        
        void f_d();

        volScalarField betaShield();
        
        void Re_y(volScalarField& S);

        void Re_v(volScalarField& S);

        void f_ss(volScalarField& S, volScalarField& O);

        volScalarField omegaTilda(volScalarField& S, dimensionedScalar& a);

        inline volScalarField omegaTilda_s(volScalarField& S){ 
            return omegaTilda(S, this->c_a1_);
        };

        inline volScalarField omegaTilda_l(volScalarField& S){ 
            return omegaTilda(S, this->c_a2_);
        };

        inline dimensionedScalar readCoeff(word name, scalar defaultValue)
        {
            return dimensioned<scalar>::lookupOrAddToDict
            (
                name,
                this->coeffDict_,
                defaultValue
            );
        };

        void gamma();

        void P_sep
        (
            volScalarField& S
        );

        volScalarField beta();

        volScalarField crossDiff();

        volScalarField nu_small(volScalarField& S, volScalarField& O);

        volScalarField nu_large(volScalarField& S, volScalarField& O);

        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField& S,
            const volScalarField& beta
        );


       

public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;	
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    // TypeName("BaseAlgKOmega");   //commented out since this is a base class
                                    //uncomment for testing of the base model

    // Constructors

        //- Construct from components
        BaseAlgKOmega
        (
	    const alphaField& alpha,
	    const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName,
        const word& type
        );


    //- Destructor
    virtual ~BaseAlgKOmega()
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
                    (c_sigmaStar_*k_ / omega_) + this->nu()
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
                    (c_sigma_ * k_ / omega_) + this->nu()
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
                    c_betaStar_*k_*omega_,
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
 #   include "BaseAlgKOmega.C"
 #endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // BASEALGKOMEGA_H

// ************************************************************************* //
