// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup OnePBoxModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_ROOTSYSTEM_PROPERTY_DEFAULTS_HH
#define DUMUX_ROOTSYSTEM_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/cellcentered/properties.hh>
#include <dumux/implicit/growth/gridgrowthproperties.hh>
#include <dumux/implicit/growth/gridgrowthindicatordefault.hh>
#include <dumux/implicit/growth/gridgrowthhelperdefault.hh>

#include <dumux/porousmediumflow/1p/implicit/localresidual.hh>

#include "model.hh"
#include "volumevariables.hh"
#include "indices.hh"
#include "spatialparams.hh"
#include "newtoncontroller.hh"
#include "fluxvariables.hh"
#include "properties.hh"

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1p.hh>
#include <dumux/implicit/1d/fvelementgeometry.hh>

namespace Dumux
{
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(Rootsystem, NumEq, 1); //!< set the number of equations to 1
SET_INT_PROP(Rootsystem, NumPhases, 1); //!< The number of phases in the rootsystem model is 1

//! The local residual function
SET_TYPE_PROP(Rootsystem, LocalResidual, OnePLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(Rootsystem, Model, RootsystemModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(Rootsystem, VolumeVariables, RootsystemVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(Rootsystem, FluxVariables, RootsystemFluxVariables<TypeTag>);

//! The class of the newton controller
SET_TYPE_PROP(Rootsystem, NewtonController, RootsystemNewtonController<TypeTag>);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(Rootsystem, Indices, RootsystemIndices);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParamsRootsystem by default.
SET_TYPE_PROP(Rootsystem, SpatialParams, RootsystemSpatialParams<TypeTag>);

//! the fv geometry for one-dimensional networks
SET_TYPE_PROP(Rootsystem, FVElementGeometry, OneDFVElementGeometry<TypeTag>);

//! The weight of the upwind control volume when calculating
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(Rootsystem, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(Rootsystem, ImplicitMobilityUpwindWeight, 1.0);

//! The fluid system to use by default
SET_TYPE_PROP(Rootsystem, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(Rootsystem, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

// enable gravity by default
SET_BOOL_PROP(Rootsystem, ProblemEnableGravity, false);

// no growing grid is the default
SET_BOOL_PROP(Rootsystem, GrowingGrid, false);

// standard setting
SET_INT_PROP(Rootsystem, GrowthInterval, 1);

//! Set the default indicator class models for growth
SET_TYPE_PROP(Rootsystem, GrowthIndicator, ImplicitGridGrowthIndicatorDefault<TypeTag>);

//! Set the default indicator class models for growth
SET_TYPE_PROP(Rootsystem, GrowthHelper, ImplicitGridGrowthHelperDefault<TypeTag>);

// \}
} // end namespace Properties

} // end namespace Dumux
#endif
