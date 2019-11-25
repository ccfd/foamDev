#ifndef functionObjects_H12_H
#define functionObjects_H12_H

#include "fieldExpression.H"
#include "fvCFD.H"
#include "meshTools.H"
#include "meshSearch.H"
#include "fvPatchFields.H"
#include "volFields.H"
#include "Vector.H"
#include "IOobject.H"
#include "interpolationCellPoint.H"

#define M_CheckIfShouldWrite(name)                                          \
    if (dict.lookupOrDefault<bool>(#name, false))                           \
    {                                                                       \
        write_##name##_ = true;                                             \
        fieldName_##name##_ = #name;                                        \
        Info << "    * " << #name << "\n";                                  \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class H12 Declaration
\*---------------------------------------------------------------------------*/

class H12
:
    public fieldExpression
{
    //Private members
        
        //- processed patch name
        const word patchName_;

        //- limit type
        const word limitType_;

        //- integration step
        const scalar dy_;

        //- normal vectors of selected boundary faces
        vectorField normals_;
        
        //- centers of selected boundary faces
        vectorField faceCenters_;

        //- scalar field of integration limit 
        //  each sampling line is integrated to limit_[face]
        scalarField limit_;

        //- List of sample point lists for each face
        List<List<point>> samplePoints_;
        
        //- List of sample cells lists for each face
        List<List<label>> sampleCells_;

        bool write_deltaStar_ = false;
        word fieldName_deltaStar_ = word();
        bool write_theta_ = false;
        word fieldName_theta_ = word();
        bool write_delta_ = false;
        word fieldName_delta_ = word();

        

    // Private Member Functions

        //- Calculate the H12 field and return true if successful
        virtual bool calc();

        //- Simple trapezoid rule integration
        inline scalar trapz(scalar a, scalar b, scalar h)
        {
            return (a + b)/2*h;
        }

        //- Find sampling points for boundary layer integral 
        //  parameters. Points are looked for in direction 
        //  to the wall. Returns true if is succesful
        bool findSamples
        (
            label face, 
            meshSearch& searchEngine,
            label startingCell
        );

public:

    //- Runtime type information
    TypeName("H12");


    // Constructors

        //- Construct from Time and dictionary
        H12
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~H12();

    inline bool writeObjectWithCheck(word objName, bool auxFields)
    {
        return (writeObject(objName) && auxFields);
    }

    //- 
    bool write()
    {
        bool auxFields = true;
        if (write_deltaStar_)
        {
            auxFields = writeObjectWithCheck(fieldName_deltaStar_, auxFields);
        }
        if (write_theta_)
        {
            auxFields = writeObjectWithCheck(fieldName_theta_, auxFields);
        }
        if (write_delta_)
        {
            auxFields = writeObjectWithCheck(fieldName_delta_, auxFields);
        }  

        return (writeObject(resultName_) && auxFields);

    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif //functionObjects_H12_H

// ************************************************************************* //
