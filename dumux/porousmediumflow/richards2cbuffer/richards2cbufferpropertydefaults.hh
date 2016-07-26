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
 * \ingroup RichardsTwoCModel
 * \file
 *
 * \brief Contains the default definitions for the properties required
 *        by the RichardsTwoC fully implicit model.
 */
#ifndef DUMUX_RICHARDS_2C_BUFFER_PROPERTY_DEFAULTS_HH
#define DUMUX_RICHARDS_2C_BUFFER_PROPERTY_DEFAULTS_HH

#include "richards2cbufferfluxvariables.hh"
#include "richards2cbuffermodel.hh"
#include "richards2cbufferproblem.hh"
#include "richards2cbufferindices.hh"
#include "richards2cbuffervolumevariables.hh"
#include "richards2cbufferproperties.hh"
#include "richards2cbuffernewtoncontroller.hh"
#include "richards2cbufferlocalresidual.hh"

#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/spatialparams/implicit.hh>

namespace Dumux
{
// \{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////
//! Number of equations required by the model
SET_INT_PROP(RichardsTwoCBuffer, NumEq, GET_PROP_VALUE(TypeTag, IsothermalNumEq));
//! Number of fluid phases considered
SET_INT_PROP(RichardsTwoCBuffer, NumPhases, 2);
SET_INT_PROP(RichardsTwoCBuffer, NumComponents, 2); //!< The number of components in the 1p2c model is 2
SET_SCALAR_PROP(RichardsTwoCBuffer, Scaling, 1); //!< Scaling of the model is set to 1 by default
SET_BOOL_PROP(RichardsTwoCBuffer, UseMoles, false); //!< Define that mole fractions are used in the balance equations

//! Use pressure [Pa] by default
SET_BOOL_PROP(RichardsTwoCBuffer, UsePH, false);

//! The local residual operator
SET_TYPE_PROP(RichardsTwoCBuffer,
              LocalResidual,
              typename GET_PROP_TYPE(TypeTag, IsothermalLocalResidual));

//! The global model used
SET_TYPE_PROP(RichardsTwoCBuffer, Model, typename GET_PROP_TYPE(TypeTag, IsothermalModel));

//! The class for the volume averaged quantities
SET_TYPE_PROP(RichardsTwoCBuffer, VolumeVariables, typename GET_PROP_TYPE(TypeTag, IsothermalVolumeVariables));

//! The class for the quantities required for the flux calculation
SET_TYPE_PROP(RichardsTwoCBuffer, FluxVariables, typename GET_PROP_TYPE(TypeTag, IsothermalFluxVariables));

//! The class of the newton controller
SET_TYPE_PROP(RichardsTwoCBuffer, NewtonController, RichardsTwoCBufferNewtonController<TypeTag>);

//! The upwind weight for the mass conservation equations
SET_SCALAR_PROP(RichardsTwoCBuffer, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(RichardsTwoCBuffer, ImplicitMobilityUpwindWeight, 1.0);

//! The class with all index definitions for the model
SET_TYPE_PROP(RichardsTwoCBuffer, Indices, typename GET_PROP_TYPE(TypeTag, IsothermalIndices));

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(RichardsTwoCBuffer, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(RichardsTwoCBuffer, EffectiveDiffusivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef DiffusivityMillingtonQuirk<Scalar> type;
};

/*!
 * \brief Set type of the parameter objects for the material law
 *
 * By default this is just retrieved from the material law.
 */
SET_PROP(RichardsTwoCBuffer, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. Please be aware that you
 * should be careful to use the RichardsTwoC model in conjunction with
 * liquid non-wetting phases. This is only meaningful if the viscosity
 * of the liquid phase is _much_ lower than the viscosity of the
 * wetting phase.
 */
SET_PROP(RichardsTwoCBuffer, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. This doed not need to be
 * specified by the problem for the RichardsTwoC model to work because the
 * RichardsTwoC model does not conserve the non-wetting phase.
 */
SET_PROP(RichardsTwoCBuffer, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::GasPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the immiscible twophase fluid system. The
 * actual fluids used are specified using in the problem definition by
 * the WettingPhase and NonwettingPhase properties. Be aware that
 * using different fluid systems in conjunction with the RichardsTwoC
 * model only makes very limited sense.
 */
SET_PROP(RichardsTwoCBuffer, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Dumux::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

// disable velocity output by default
SET_BOOL_PROP(RichardsTwoCBuffer, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(RichardsTwoCBuffer, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(BoxModel, SpatialParamsForchCoeff, 0.55);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(NonIsothermal, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! temperature is already written by the isothermal model
SET_BOOL_PROP(RichardsTwoCBufferNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(RichardsTwoCBuffer, IsothermalModel, RichardsTwoCBufferModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(RichardsTwoCBuffer, IsothermalFluxVariables, RichardsTwoCBufferFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(RichardsTwoCBuffer, IsothermalVolumeVariables, RichardsTwoCBufferVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(RichardsTwoCBuffer, IsothermalLocalResidual, RichardsTwoCBufferLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(RichardsTwoCBuffer, IsothermalIndices, RichardsTwoCBufferIndices<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(RichardsTwoCBuffer, IsothermalNumEq, 2);

// \}
}

} // end namespace

#endif
