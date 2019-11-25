/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "BaseAlgKOmega.H"
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
void BaseAlgKOmega<BasicTurbulenceModel>::valueCheck
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

//- Update nut field
template<class BasicTurbulenceModel>    
void BaseAlgKOmega<BasicTurbulenceModel>::correctNut()
{   
    //correct nu_s and and nu_s
    this->nut_ = (nu_l_ + nu_s_);
    this->nut_.relax();
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    BasicTurbulenceModel::correctNut();
}

//- Strain rate tensor magnitude
template<class BasicTurbulenceModel>
volScalarField BaseAlgKOmega<BasicTurbulenceModel>::ShearRateMag() 
{
    volTensorField gradU = fvc::grad(this->U_);
    // S_ij = 1/2(u_i,j + u_j,i) -1/3 * u_k,k * delta_ij
    return sqrt
    (   
        2 * magSqr(symm( gradU ) - (1.0 / 3.0) * tr( gradU ) * tensor::I ) 
    );
}

//- Vorticity tensor magnitude
template<class BasicTurbulenceModel>
volScalarField BaseAlgKOmega<BasicTurbulenceModel>::VorticityMag() 
{
    //O_ij = 1/2*(u_i,j - u_j,i)
    return sqrt( 2 * ( magSqr( skew( fvc::grad(this->U_) ) ) ) );
}

//- Re_y = sqrt(k)*y/nu Reynolds number
template<class BasicTurbulenceModel>
void BaseAlgKOmega<BasicTurbulenceModel>::Re_y 
(
    volScalarField& S
)
{
    this->Re_y_ = (sqrt(k_) * y_ ) / this->nu();
}

//- Re_v = S*y^2/nu Reynolds number - estimates Re_theta
template<class BasicTurbulenceModel>
void BaseAlgKOmega<BasicTurbulenceModel>::Re_v 
(
    volScalarField& S
)
{
    this->Re_v_ = pow(y_,2) * S / this->nu();
} 

//- Calculate shear sheltering factor
template<class BasicTurbulenceModel>
void BaseAlgKOmega<BasicTurbulenceModel>::f_ss 
(
    volScalarField& S,
    volScalarField& O
)
{   
    this->f_k_ = 1.0 - tanh(k_/(c_Ck_ * this->nu() * omega_));

    this->chi_ = 
        tanh
        ( 
            ( -O *( S - O ) ) 
            /
            ((c_Ckh_ * pow( c_betaStar_* omega_,2.0)) )
        );
    
    this->C_ss_ = (c_Cs_* (1.0 + c_CA_ * f_k_ * chi_));
    this->f_ss_ = exp( -pow(C_ss_ / Re_y_,2));
}

//- Calculate shielding function f_d
template<class BasicTurbulenceModel>
void BaseAlgKOmega<BasicTurbulenceModel>::f_d()
{
    tmp<volScalarField> r_d = 
          ( this->nu() + this->nut_ ) 
        / ( sqrt( magSqr( fvc::grad(this->U_) ) ) * pow( y_ * 0.41 , 2.0) );

    this->f_d_ = 1 - tanh(pow( 8 * r_d ,3.0)) ;
}

//- Calculate beta function employing LES shielding function f_d
template<class BasicTurbulenceModel>
volScalarField BaseAlgKOmega<BasicTurbulenceModel>::betaShield()
{
    f_d();

    volScalarField f_dgamma = min( this->f_d_, gamma_ );
    
    tmp<volSymmTensorField> Sij = 
          symm
          ( 
              fvc::grad(this->U_) 
            - (1.0/3.0) * tr(fvc::grad(this->U_)) * tensor::I
          );
	
    volTensorField Oij = skew( fvc::grad(this->U_) );

    volScalarField chi_w = 
    mag(
          (Oij & Oij) && Sij
        / pow(c_betaStar_ * omega_ , 3)
    );

    tmp<volScalarField> f_beta = (1.0 + 85.0 * chi_w)/(1.0 + 100.0 * chi_w);
    tmp<volScalarField> f_betaS = f_dgamma * f_beta + (1 - f_dgamma);

    tmp<volScalarField> corrected_beta = f_betaS* c_beta0_;
    return corrected_beta;
}

//- Calculate standard beta function
template<class BasicTurbulenceModel>
volScalarField BaseAlgKOmega<BasicTurbulenceModel>::beta()	
{ 
    tmp<volSymmTensorField> Sij = symm( fvc::grad(this->U_) );
	volTensorField          Oij = skew( fvc::grad(this->U_) );
    
    volScalarField chi_w = mag(
          (Oij & Oij) && Sij
        / pow(c_betaStar_ * omega_ , 3)
    );
    
    tmp<volScalarField> f_beta = (1.0 + 85.0 * chi_w) / (1.0 + 100.0 * chi_w);
    tmp<volScalarField> corrected_beta = f_beta * c_beta0_;
    
    return corrected_beta;
}

//- omegaTilda helper function
template<class BasicTurbulenceModel>
volScalarField BaseAlgKOmega<BasicTurbulenceModel>::omegaTilda 
(
    volScalarField& S,
    dimensionedScalar& a
)
{
    tmp<volScalarField> omegaTilda
        (
            max
            (
                omega_,
                c_Clim_ * S / a
            )
        );
    return omegaTilda;
}

//- Intermitenncy factor update
template<class BasicTurbulenceModel>
void BaseAlgKOmega<BasicTurbulenceModel>::gamma()
{
    this->gamma_ =
        min
        (
            max
            (
                (sqrt(k_) * y_)/(c_Agamma_ * this->nu()) - 1.0,
                0.0
            ),
            1.0
        ); 
}

//- Separation-induced production
template<class BasicTurbulenceModel>
void BaseAlgKOmega<BasicTurbulenceModel>::P_sep
(
    volScalarField& S
)
{
    this->F_sep_ = 
        min
        (
            max
            (
                Re_v_ / (2.2 *c_Av_) - 1.0,
                0.0
            ),
            1.0
        );

    this->P_sep_ = c_Csep_ * F_sep_ * this->nu() * S * S ;
}


//- Cross-diffusion term
template<class BasicTurbulenceModel>
volScalarField BaseAlgKOmega<BasicTurbulenceModel>::crossDiff()
{
	// volScalarField sigmaD = beta;
    volScalarField sigmaD(
        IOobject 
        ( 
            IOobject::groupName("sigmaD", this->U_.group()), 
            this->runTime_.timeName(), 
            this->mesh_, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE 
        ), 
        this->mesh_, 
        dimensionedScalar 
        ( 
            IOobject::groupName("sigmaD", this->U_.group()), 
            dimensionSet(0,0,0,0,0,0,0), 
            0.0 
        ), 
        "calculated"
    );

	volVectorField gradK = fvc::grad(k_);
	volVectorField gradOmega = fvc::grad(omega_);
	volScalarField gradKO = gradK & gradOmega;

	forAll(gradKO, celli){
        if ( gradKO[celli] >  0)
        {
            sigmaD[celli] = c_sigmaD0_.value();
        }
	}
	return sigmaD * gradKO / omega_ ;
	
}

//- Calculate small scale turbulent viscosity 
template<class BasicTurbulenceModel>
volScalarField BaseAlgKOmega<BasicTurbulenceModel>::nu_small(
    volScalarField& S, 
    volScalarField& O
)
{
    this->k_s_ = k_ * this->f_ss_;

    return this->k_s_ / this->omegaTilda_s(S);
}

//- Calculate large scale turbulent viscosity
template<class BasicTurbulenceModel>
volScalarField BaseAlgKOmega<BasicTurbulenceModel>::nu_large
(
    volScalarField& S, 
    volScalarField& O
)
{
    this->k_l_ = k_ *(1 - this->f_ss_);
    return this->k_l_ / this->omegaTilda_l(S);
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> BaseAlgKOmega<BasicTurbulenceModel>::Qsas
(
    const volScalarField& S,
    const volScalarField& beta
)
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
BaseAlgKOmega<BasicTurbulenceModel>::BaseAlgKOmega
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
    eddyViscosity<RASModel<BasicTurbulenceModel> >
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
    // constants
    c_Agamma_(    readCoeff("Agamma",    12.0   ) ),
    c_Cs_(        readCoeff("Cs",        21.0   ) ),
    c_CA_(        readCoeff("Ca",        1.0    ) ),
    c_Ckh_(       readCoeff("Ckh",       10.0   ) ),
    c_Ck_(        readCoeff("Ck",        10.0   ) ),//6.0->10 
	c_Csep_(      readCoeff("Csep",      2.0    ) ),
	c_Av_(        readCoeff("Av",        550.0  ) ),
	c_a1_(        readCoeff("a1",        0.3    ) ),
	c_a2_(        readCoeff("a2",        0.6    ) ),//0.45->0.6
	c_betaStar_(  readCoeff("betaStar",  0.09   ) ),
	c_alfa_(      readCoeff("alfa",      0.52   ) ),
	c_beta0_(     readCoeff("betaO",     0.0708 ) ),
	c_sigma_(     readCoeff("sigma",     0.5    ) ),
	c_sigmaStar_( readCoeff("sigmaStar", 0.6    ) ),
	c_sigmaD0_(   readCoeff("sigmaD0",   0.125  ) ),
	c_Clim_(      readCoeff("Clim",      0.875  ) ),
    
    // wall distance field
    y_(wallDist::New(this->mesh_).y()),

    // Fields - modelled - always written
    omega_(IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_),
    k_(IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_),

    // Fields - auxiliary - optionaly written
    gamma_( M_scalarFieldInit("gamma", dimensionSet(0,0, 0,0,0,0,0), 0.0) ),
    nu_s_(  M_scalarFieldInit("nu_s",  dimensionSet(0,2,-1,0,0,0,0), 0.0) ),
    nu_l_(  M_scalarFieldInit("nu_l",  dimensionSet(0,2,-1,0,0,0,0), 0.0) ) ,
    k_s_(   M_scalarFieldInit("k_s",   dimensionSet(0,2,-2,0,0,0,0), 0.0) ) ,
    k_l_(   M_scalarFieldInit("k_l",   dimensionSet(0,2,-2,0,0,0,0), 0.0) ),
    f_ss_(  M_scalarFieldInit("f_ss",  dimensionSet(0,0, 0,0,0,0,0), 0.0) ),
    C_ss_(  M_scalarFieldInit("C_ss",  dimensionSet(0,0, 0,0,0,0,0), 0.0) ),
    P_k_(   M_scalarFieldInit("P_k",   dimensionSet(0,2,-3,0,0,0,0), 0.0) ),
    P_sep_( M_scalarFieldInit("P_sep", dimensionSet(0,2,-3,0,0,0,0), 0.0) ),
    F_sep_( M_scalarFieldInit("F_sep", dimensionSet(0,0, 0,0,0,0,0), 0.0) ),
    f_k_(   M_scalarFieldInit("f_k",   dimensionSet(0,0, 0,0,0,0,0), 0.0) ),
    chi_(   M_scalarFieldInit("chi",   dimensionSet(0,0, 0,0,0,0,0), 0.0) ),
    Re_v_(  M_scalarFieldInit("Re_v",  dimensionSet(0,0, 0,0,0,0,0), 0.0) ),
    Re_y_(  M_scalarFieldInit("Re_y",  dimensionSet(0,0, 0,0,0,0,0), 0.0) ),
    f_d_(   M_scalarFieldInit("f_d",   dimensionSet(0,0, 0,0,0,0,0), 0.0) )
        
    {
   	    bound(k_, this->kMin_);
    	bound(omega_, this->omegaMin_);

        // Uncomment below code for testing of the baseModel
        /* 
        if (type == typeName)
    	{
            this->printCoeffs(type);
    	}
        */
        
        this->nut_.storePrevIter();

    }
    

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool BaseAlgKOmega<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        c_Agamma_   .readIfPresent( this->coeffDict() );
        c_Cs_       .readIfPresent( this->coeffDict() );
        c_CA_       .readIfPresent( this->coeffDict() );
        c_Ckh_      .readIfPresent( this->coeffDict() );
        c_Ck_       .readIfPresent( this->coeffDict() );
	    c_Csep_     .readIfPresent( this->coeffDict() );
        c_Av_       .readIfPresent( this->coeffDict() );
        c_a1_       .readIfPresent( this->coeffDict() );
        c_a2_       .readIfPresent( this->coeffDict() );
        c_betaStar_ .readIfPresent( this->coeffDict() );
	    c_alfa_     .readIfPresent( this->coeffDict() );
        c_beta0_    .readIfPresent( this->coeffDict() );
        c_sigma_    .readIfPresent( this->coeffDict() );
        c_sigmaStar_.readIfPresent( this->coeffDict() );
        c_sigmaD0_  .readIfPresent( this->coeffDict() );
	    c_Clim_     .readIfPresent( this->coeffDict() );

        return true;
    }
    else
    {
        return false;
    }
}



template<class BasicTurbulenceModel>
void BaseAlgKOmega<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Local references
    const volVectorField& U = this->U_;
    const rhoField& rho = this->rho_;
    const alphaField& alpha = this->alpha_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    
    volScalarField& nut = this->nut_;

     
	volScalarField S( ShearRateMag() );
	volScalarField O( VorticityMag() );
    volScalarField nu_s( nu_small(S,O) );
    volScalarField beta( betaShield());


    // Compute algebraic relations
    Re_y(S);
    Re_v(S);
    f_ss(S,O);
    gamma();
    this->P_k_ = nu_s * S * S;
    P_sep(S);

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //Computing G is necessary for general OF implementation, this may be 
    //reworked in future version of the code
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G
    (
        this->GName(),
       
        nut*(2 * symm(tgradU()) && symm(tgradU()))
    );
    tgradU.clear();

    // Update omega at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian( (alpha*rho) * DomegaEff(), omega_)
     ==
        fvm::Sp( (alpha*rho) * c_alfa_*P_k_ /k_, omega_) 
        // below term was added to enchance stability 
      - fvm::SuSp((alpha*rho) * (2.0/3.0)*c_alfa_*divU, omega_)  
      - fvm::Sp( (alpha*rho) * beta*omega_, omega_)
      + (alpha*rho) * crossDiff()
      + Qsas(S, beta)
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
      - fvm::laplacian((alpha*rho) * DkEff(), k_)
     ==
        (alpha*rho) * gamma_ * P_k_
        // below term was added to enchance stability 
      - fvm::SuSp((alpha*rho) * (2.0/3.0)*divU, k_) 
      + (alpha*rho) * (1.0 - gamma_) * P_sep_
      - fvm::Sp( (alpha*rho) * c_betaStar_*omega_, k_)
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, this->kMin_);

    // save and update viscosity
    this->nu_s_ = this->nu_small(S, O);
    this->nu_l_ = this->nu_large(S, O);
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace RASModels
} // End namespace Foam
// ************************************************************************* //
