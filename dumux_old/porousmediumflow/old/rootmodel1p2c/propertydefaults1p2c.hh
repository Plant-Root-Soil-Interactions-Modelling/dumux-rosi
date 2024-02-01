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
#ifndef DUMUX_ROOTSYSTEM_PROPERTY_DEFAULTS_1P2C_HH
#define DUMUX_ROOTSYSTEM_PROPERTY_DEFAULTS_1P2C_HH

//#include <dumux/porousmediumflow/1p2c/implicit/properties.hh>
//#include <dumux/implicit/growth/gridgrowthproperties.hh>
#include <dumux/implicit/growth/gridgrowthindicatordefault.hh>
#include <dumux/implicit/growth/gridgrowthhelperdefault.hh>
//#include <dumux/porousmediumflow/1p/implicit/localresidual.hh>
//#include <dumux/porousmediumflow/1p2c/implicit/localresidual.hh>
#include "localresidual1p2c.hh"

//#include "model.hh"
#include "model1p2c.hh"

//#include "volumevariables.hh"
#include "dumux/porousmediumflow/1p2c/implicit/volumevariables.hh"

//#include "indices.hh"
#include "dumux/porousmediumflow/1p2c/implicit/indices.hh"

#include "dumux/porousmediumflow/1d/rootsystem/spatialparams.hh"
//#include "spatialparams.hh"
//#include "newtoncontroller.hh"
#include "dumux/porousmediumflow/1d/rootsystem/newtoncontroller.hh"

//#include "fluxvariables.hh"
//#include "dumux/porousmediumflow/1p2c/implicit/fluxvariables.hh"
#include "fluxvariables1p2c.hh"
#include "properties1p2c.hh"

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1p.hh>
#include <dumux/implicit/1d/fvelementgeometry.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttau.hh>

namespace Dumux
{
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(RootsystemOnePTwoC, NumEq, 2); //!< set the number of equations to 1
SET_INT_PROP(RootsystemOnePTwoC, NumPhases, 1); //!< The number of phases in the rootsystem model is 1
SET_INT_PROP(RootsystemOnePTwoC, NumComponents, 2); //!< The number of components in the 1p2c model is 2
SET_SCALAR_PROP(RootsystemOnePTwoC, Scaling, 1); //!< Scaling of the model is set to 1 by default ?!
SET_BOOL_PROP(RootsystemOnePTwoC, UseMoles, true); //!< Define that mole fractions are used in the balance equations

//! The local residual function
//SET_TYPE_PROP(Rootsystem, LocalResidual, OnePLocalResidual<TypeTag>); //OnePTwoCLocalResidual
SET_TYPE_PROP(RootsystemOnePTwoC, LocalResidual, RootSystemOnePTwoCLocalResidual<TypeTag>);
//SET_TYPE_PROP(RootsystemOnePTwoC, LocalResidual, OnePTwoCLocalResidual<TypeTag>);


//! the Model property
SET_TYPE_PROP(RootsystemOnePTwoC, Model, RootsystemOnePTwoCModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(RootsystemOnePTwoC, VolumeVariables, OnePTwoCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(RootsystemOnePTwoC, FluxVariables, RootsystemFluxOnePTwoCFluxVariables<TypeTag>);

//! The class of the newton controller
SET_TYPE_PROP(RootsystemOnePTwoC, NewtonController, RootsystemNewtonController<TypeTag>);

//! The indices required by the isothermal single-phase model

SET_TYPE_PROP(RootsystemOnePTwoC, Indices, OnePTwoCIndices<TypeTag>);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParamsRootsystem by default.
SET_TYPE_PROP(RootsystemOnePTwoC, SpatialParams, RootsystemSpatialParams<TypeTag>);

//! the fv geometry for one-dimensional networks
SET_TYPE_PROP(RootsystemOnePTwoC, FVElementGeometry, OneDFVElementGeometry<TypeTag>);

//! The weight of the upwind control volume when calculating
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(RootsystemOnePTwoC, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(RootsystemOnePTwoC, ImplicitMobilityUpwindWeight, 1.0);

//! The fluid system to use by default
//SET_TYPE_PROP(Rootsystem, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

//SET_PROP(Rootsystem, Fluid)
//{ private:
//    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//public:
//    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
//};

SET_PROP(RootsystemOnePTwoC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

////! The model after Millington (1961) is used for the effective diffusivity ?! porosity of root ?!?!?!?!
SET_PROP(RootsystemOnePTwoC, EffectiveDiffusivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef DiffusivityMillingtonQuirk<Scalar> type;
//    typedef DiffusivityConstantTau<TypeTag, Scalar> type;
};

// enable gravity by default
SET_BOOL_PROP(RootsystemOnePTwoC, ProblemEnableGravity, false);

// no growing grid is the default
SET_BOOL_PROP(RootsystemOnePTwoC, GrowingGrid, false);

// standard setting
SET_INT_PROP(RootsystemOnePTwoC, GrowthInterval, 1);

// set Phase idx
SET_INT_PROP(RootsystemOnePTwoC, PhaseIdx, 0);

//! Set the default indicator class models for growth
SET_TYPE_PROP(RootsystemOnePTwoC, GrowthIndicator, ImplicitGridGrowthIndicatorDefault<TypeTag>);

//! Set the default indicator class models for growth
SET_TYPE_PROP(RootsystemOnePTwoC, GrowthHelper, ImplicitGridGrowthHelperDefault<TypeTag>);

// \}
} // end namespace Properties

} // end namespace Dumux
#endif
