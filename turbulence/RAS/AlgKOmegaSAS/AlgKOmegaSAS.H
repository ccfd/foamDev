
#ifndef ALGKOMEGASAS_H
#define ALGKOMEGASAS_H

#include "BaseAlgKOmega.H"
#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class BaseAlgKOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class AlgKOmegaSAS
:
    public BaseAlgKOmega<BasicTurbulenceModel>
{
protected:
    // Protected data

        // Model constants

            dimensionedScalar Cs_;
            dimensionedScalar kappa_;
            dimensionedScalar zeta2_;
            dimensionedScalar sigmaPhi_;
            dimensionedScalar C_;


        // Fields

            //- Run-time selectable delta model
            autoPtr<Foam::LESdelta> delta_;
            volScalarField L_;
            volScalarField Lvk_;
            volScalarField Qsas_;

    // Protected Member Functions

        //- SAS omega source
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
    TypeName("AlgKOmegaSAS");

    // Constructors

        //- Construct from components
        AlgKOmegaSAS
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
    virtual ~AlgKOmegaSAS()
    {}


    // Member Functions
        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function to filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }
    


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 #ifdef NoRepository
 #   include "AlgKOmegaSAS.C"
 #endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ALGKOMEGASAS_H

// ************************************************************************* //
