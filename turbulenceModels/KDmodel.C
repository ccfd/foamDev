
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

\*---------------------------------------------------------------------------*/

#include "KDmodel.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
#include "Vector.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //




template<class BasicTurbulenceModel>
void KDmodel<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
volScalarField KDmodel<BasicTurbulenceModel>::S() //strain rate magnitude
{

//	return sqrt((2*(magSqr(symm(fvc::grad(this->U_))))));
	return mag(2*(symm(fvc::grad(this->U_))));
	
}

template<class BasicTurbulenceModel>
volScalarField KDmodel<BasicTurbulenceModel>::O() //strain rate magnitude
{
 //	return (sqrt(2*(magSqr(skew(fvc::grad(this->U_))))));
 	return mag(2*(skew(fvc::grad(this->U_))));
}
template<class BasicTurbulenceModel>
volScalarField KDmodel<BasicTurbulenceModel>::gamma()
{
    volScalarField gamma
    (
        min
		(
			max
			(
				(sqrt(k_) * y_)/(Agamma_ * this->nu()) - 1.0,
				0.0
			),
			1.0
		)

    );
    return gamma;

}


template<class BasicTurbulenceModel>
volScalarField KDmodel<BasicTurbulenceModel>::P // production
(
    volScalarField& S,
    volScalarField& O
)
{
    tmp<volScalarField> f_k = 1.0 - tanh(k_/(Ck_ * this->nu() * omega_));
    tmp<volScalarField> tmp1
    (
        S-O
	    +dimensionedScalar("small",dimensionSet(0,0,-1,0,0), SMALL)
    );
    tmp<volScalarField> chi = tanh( (-tmp1 * O ) / (Ckh_ * pow( betaStar_* omega_,2.0)));
    tmp<volScalarField> tmp2
    (
        f_k * chi
		+dimensionedScalar("small", dimensionSet(0,0,0,0,0), SMALL)
    );
    tmp<volScalarField> f_ss = exp( -pow((Cs_* (1.0 + Ca_ * tmp2) * this->nu()/ (sqrt(k_) * y_ )),2));
    tmp<volScalarField> tmp3
    (
        max
		(
			omega_,
			Clim_ * S / a1_
		)
    );
    tmp<volScalarField> nu_s = f_ss * k_ / tmp3;
    
    return nu_s * pow(S,2.0);

}


template<class BasicTurbulenceModel> 
void KDmodel<BasicTurbulenceModel>::valueCheck
(
	volScalarField& fieldValue
)
{
	Foam::dimensioned<double> argmin = min(fieldValue);
	Info << "Minimum value is:" << endl;
	Info << argmin << endl;
	Foam::dimensioned<double> argmax = max(fieldValue);
	Info << "Maximum value is:" << endl;
	Info << argmax << endl;

}

template<class BasicTurbulenceModel>
volScalarField KDmodel<BasicTurbulenceModel>::Psep
(
    volScalarField& S
)
{
    tmp<volScalarField> rv = pow(y_,2) * S / this->nu();
    tmp<volScalarField> f_sep
    (   
        min
        (
            max
            (
                rv/(2.2 *Anu_) - 1.0,
                0.0
            ),
            1.0
        )
    );
    return  Csep_ * f_sep * this->nu() * S * S;

}

template<class BasicTurbulenceModel>
volScalarField KDmodel<BasicTurbulenceModel>::beta()	
{

	volSymmTensorField Sij = symm(fvc::grad(this->U_));
	tmp<volTensorField> Oij = skew(fvc::grad(this->U_));
    tmp<volTensorField> tmp1(Oij & Oij);
    tmp<volScalarField> tmp2(tmp1 && Sij);
    volScalarField chi_w = mag( tmp2 / pow(betaStar_ * omega_ , 3));
    tmp<volScalarField> f_beta = (1.0 + 85.0 * chi_w)/(1.0 + 100.0 * chi_w);
    tmp<volScalarField> corrected_beta = f_beta * beta0_;
    return corrected_beta;

}

template<class BasicTurbulenceModel>
volScalarField KDmodel<BasicTurbulenceModel>::crossDiff
(
	volScalarField& beta
)
{
	volScalarField sigmaD = beta;
	volVectorField gradK = fvc::grad(k_);
	volVectorField gradOmega = fvc::grad(omega_);
	volScalarField gradKO = gradK & gradOmega;

	forAll(gradKO, celli){
	if ( gradKO[celli] <=  0)
	{
		sigmaD[celli] = 0;
	}
	else
	{
		sigmaD[celli] = 0.125;
	}
	}
	return sigmaD * gradKO / omega_ ;
	
}


template<class BasicTurbulenceModel>
volScalarField KDmodel<BasicTurbulenceModel>::gammaDif
(
	volScalarField& Gamma
)
{
	return 1- Gamma;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
KDmodel<BasicTurbulenceModel>::KDmodel
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:

    eddyViscosity<RASModel<BasicTurbulenceModel>>
	(
		type,
		alpha,
		rho,
		U,
		alphaRhoPhi,
		phi,
		transport,
		propertiesName    

	),


    Agamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Agamma",
            this->coeffDict_,
            12.0
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            21.0
        )
    ),
    Ca_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ca",
            this->coeffDict_,
            1.0
        )
    ),
    Ckh_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ckh",
            this->coeffDict_,
            10.0
        )
    ),
    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            6.0
        )
    ),
	Csep_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Csep",
            this->coeffDict_,
            2.0
        )
    ),
	Anu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anu",
            this->coeffDict_,
            550.0
        )
    ),
	a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.3
        )
    ),
	a2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a2",
            this->coeffDict_,
            0.45
        )
    ),
	betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
	alfa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alfa",
            this->coeffDict_,
            0.52
        )
    ),
	beta0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaO",
            this->coeffDict_,
            0.0708
        )
    ),
	sigma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigma",
            this->coeffDict_,
            0.5
        )
    ),
	sigmaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaStar",
            this->coeffDict_,
            0.6
        )
    ),
	sigmaD0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaD0",
            this->coeffDict_,
            0.125
        )
    ),
	Clim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim",
            this->coeffDict_,
            0.875
        )
    ),
    
    

    y_(wallDist::New(this->mesh_).y()),



    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
    {
   	 bound(k_, this->kMin_);
    	bound(omega_, this->omegaMin_);

    	if (type == typeName)
    	{
        	this->printCoeffs(type);
    	}

    }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool KDmodel<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Agamma_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        Ca_.readIfPresent(this->coeffDict());
        Ckh_.readIfPresent(this->coeffDict());
        Ck_.readIfPresent(this->coeffDict());
	    Csep_.readIfPresent(this->coeffDict());
        Anu_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        a2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
	    alfa_.readIfPresent(this->coeffDict());
        beta0_.readIfPresent(this->coeffDict());
        sigma_.readIfPresent(this->coeffDict());
        sigmaStar_.readIfPresent(this->coeffDict());
        sigmaD0_.readIfPresent(this->coeffDict());
	    Clim_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void KDmodel<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references

    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    const rhoField& rho = this->rho_;
    const alphaField& alpha = this->alpha_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const surfaceScalarField&  phi = this->phi_;  
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));
    
	volScalarField S(this->S());
	volScalarField O(this->O());
    volScalarField gamma(this->gamma());
    volScalarField P(this->P(S,O));
	volScalarField Psep(this->Psep(S));
    volScalarField beta(this->beta());
    volScalarField crossDiff(this->crossDiff(beta));
	volScalarField gammaDif(this->gammaDif(gamma));

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G
    (
        this->GName(),
        nut*(tgradU() && dev(twoSymm(tgradU())))
    );
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        fvm::Sp(alfa_*alpha*rho*P/k_,omega_)
      - fvm::Sp(beta0_*alpha*rho*omega_, omega_)
      + crossDiff * alpha * rho
	
    );


    omegaEqn.ref().relax();
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    bound(omega_, this->omegaMin_);
    

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi,k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        P*alpha*rho*gamma
      + gammaDif* alpha*rho* Psep
      - fvm::Sp(betaStar_*alpha*rho*omega_, k_)
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, this->kMin_);

    correctNut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace RASModels
} // End namespace Foam
// ************************************************************************* //
