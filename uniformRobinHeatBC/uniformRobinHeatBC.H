/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::uniformRobinHeatBC

Description
    Boundary condition specifies a wall, given heat flux value, the thermal conductivity of the fluid
   , heat transfer coefficient of an external fluid, and temperature of the external fluid.

    @verbatim
        inlet
        {
            type                    uniformRobinHeatBC;
            flux                    100.0;
            k                       1.1;
            ho                      0.1;
            To                      280.0;
        }
    @endverbatim

SourceFiles
    uniformRobinHeatBC.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef uniformRobinHeatBC_H
#define uniformRobinHeatBC_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformRobinFvPatch Declaration
\*---------------------------------------------------------------------------*/

class uniformRobinHeatBC
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Flux
        scalar flux_;

        //- thermal condoctivity
        scalar k_;

        //- convective heat transfer
        scalar ho_;

        //- Out temperature
        scalar To_;


public:

    //- Runtime type information
    TypeName("uniformRobinHeatBC");


    // Constructors

        //- Construct from patch and internal field
        uniformRobinHeatBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformRobinHeatBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformRobinHeatBC onto
        //  a new patch
        uniformRobinHeatBC
        (
            const uniformRobinHeatBC&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformRobinHeatBC
        (
            const uniformRobinHeatBC&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformRobinHeatBC(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uniformRobinHeatBC
        (
            const uniformRobinHeatBC&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformRobinHeatBC(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
