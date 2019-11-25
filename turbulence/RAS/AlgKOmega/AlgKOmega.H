
#ifndef ALGKOMEGA_H
#define ALGKOMEGA_H

#include "BaseAlgKOmega.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class BaseAlgKOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class AlgKOmega
:
    public BaseAlgKOmega<BasicTurbulenceModel>
{

public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;	
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("AlgKOmega");

    // Constructors

        //- Construct from components
        AlgKOmega
        (
	    const alphaField& alpha,
	    const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
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
        )
        {
            if (type == typeName)
            {
                this->printCoeffs(type);
            }
        }


    //- Destructor
    virtual ~AlgKOmega()
    {}


    // Member Functions

       

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam


#endif // ALGKOMEGA_H

// ************************************************************************* //
