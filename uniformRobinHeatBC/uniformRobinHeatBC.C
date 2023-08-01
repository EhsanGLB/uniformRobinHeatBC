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

\*---------------------------------------------------------------------------*/

#include "uniformRobinHeatBC.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformRobinHeatBC::uniformRobinHeatBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    flux_(0.0),
    k_(0.0),
    ho_(0.0),
    To_(0.0)
{}


Foam::uniformRobinHeatBC::uniformRobinHeatBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    flux_(readScalar(dict.lookup("flux"))),
    k_(readScalar(dict.lookup("k"))),
    ho_(readScalar(dict.lookup("ho"))),
    To_(readScalar(dict.lookup("To")))
{

    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::uniformRobinHeatBC::uniformRobinHeatBC
(
    const uniformRobinHeatBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    flux_(ptf.flux_),
    k_(ptf.k_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


Foam::uniformRobinHeatBC::uniformRobinHeatBC
(
    const uniformRobinHeatBC& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    flux_(ptf.flux_),
    k_(ptf.k_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


Foam::uniformRobinHeatBC::uniformRobinHeatBC
(
    const uniformRobinHeatBC& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    flux_(ptf.flux_),
    k_(ptf.k_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformRobinHeatBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);

    gradient() = ( flux_ - ho_ * (Tp - To_) ) / k_;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::uniformRobinHeatBC::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("flux") << flux_ << token::END_STATEMENT << nl;
    os.writeKeyword("k") << k_ << token::END_STATEMENT << nl;
    os.writeKeyword("ho") << ho_ << token::END_STATEMENT << nl;
    os.writeKeyword("To") << To_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        uniformRobinHeatBC
    );
}

// ************************************************************************* //
