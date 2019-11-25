#include "AlgKOmegaSAS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> AlgKOmegaSAS<BasicTurbulenceModel>::Qsas
(
    const volScalarField& S,
    const volScalarField& beta
)
{
    this->L_ = sqrt(this->k_)/(pow025(this->c_betaStar_)*this->omega_);


    this->Lvk_ = 
    (
        max
        (
            kappa_*S
           /(
                mag(fvc::laplacian(this->U_))
              + dimensionedScalar
                (
                    "rootVSmall",
                    dimensionSet(0, -1, -1, 0, 0),
                    rootVSmall
                )
            ),
            Cs_*sqrt(kappa_*zeta2_/(beta/this->c_betaStar_ - this->c_alfa_))*delta()
        )
    );

    this->Qsas_ = 
    max
    (
          zeta2_*kappa_*sqr(S)*sqr(this->L_/this->Lvk_)
        - (2*C_/sigmaPhi_)*this->k_
        * max
        (
            magSqr(fvc::grad(this->omega_))/sqr(this->omega_),
            magSqr(fvc::grad(this->k_))/sqr(this->k_)
        ),    
        dimensionedScalar("0", dimensionSet(0, 0, -2, 0, 0), 0)
    );

    return fvm::Su
    (
        this->alpha_()*this->rho_()
        * min
        (
            this->Qsas_(),
            this->omega_()/(0.1*this->omega_.time().deltaT())  
        ),
        this->omega_
    );
}

template<class BasicTurbulenceModel>
AlgKOmegaSAS<BasicTurbulenceModel>::AlgKOmegaSAS
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
    BaseAlgKOmega<BasicTurbulenceModel>
        (
	    alpha,
	    rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
        ),

    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.11
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    zeta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "zeta2",
            this->coeffDict_,
            3.51
        )
    ),
    sigmaPhi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPhi",
            this->coeffDict_,
            2.0/3.0
        )
    ),
    C_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C",
            this->coeffDict_,
            2
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            *this,
            this->coeffDict_
        )
    ),
    L_(     M_scalarFieldInit("L",    dimensionSet(0,1, 0,0,0,0,0), 0.0)),
    Lvk_(   M_scalarFieldInit("Lvk",  dimensionSet(0,1, 0,0,0,0,0), 0.0)),
    Qsas_(  M_scalarFieldInit("Qsas", dimensionSet(0,0,-2,0,0,0,0), 0.0))
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

template<class BasicTurbulenceModel>
bool AlgKOmegaSAS<BasicTurbulenceModel>::read()
{
    if (BaseAlgKOmega<BasicTurbulenceModel>::read())
    {
        Cs_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        sigmaPhi_.readIfPresent(this->coeffDict());
        zeta2_.readIfPresent(this->coeffDict());
        C_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //