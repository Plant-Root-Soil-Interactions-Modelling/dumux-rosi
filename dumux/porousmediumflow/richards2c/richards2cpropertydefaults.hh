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
#ifndef DUMUX_RICHARDS_2C_PROPERTY_DEFAULTS_HH
#define DUMUX_RICHARDS_2C_PROPERTY_DEFAULTS_HH

#include "richards2cfluxvariables.hh"
#include "richards2cmodel.hh"
#include "richards2cproblem.hh"
#include "richards2cindices.hh"
#include "richards2cvolumevariables.hh"
#include "richards2cproperties.hh"
#include "richards2cnewtoncontroller.hh"
#include "richards2clocalresidual.hh"

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
SET_INT_PROP(RichardsTwoC, NumEq, GET_PROP_VALUE(TypeTag, IsothermalNumEq));
//! Number of fluid phases considered
SET_INT_PROP(RichardsTwoC, NumPhases, 2);
SET_INT_PROP(RichardsTwoC, NumComponents, 2); //!< The number of components in the 1p2c model is 2
SET_SCALAR_PROP(RichardsTwoC, Scaling, 1); //!< Scaling of the model is set to 1 by default
SET_BOOL_PROP(RichardsTwoC, UseMoles, false); //!< Define that mole fractions are used in the balance equations

//! Use pressure [Pa] by default
SET_BOOL_PROP(RichardsTwoC, UsePH, false);

//! The local residual operator
SET_TYPE_PROP(RichardsTwoC,
              LocalResidual,
              typename GET_PROP_TYPE(TypeTag, IsothermalLocalResidual));

//! The global model used
SET_TYPE_PROP(RichardsTwoC, Model, typename GET_PROP_TYPE(TypeTag, IsothermalModel));

//! The class for the volume averaged quantities
SET_TYPE_PROP(RichardsTwoC, VolumeVariables, typename GET_PROP_TYPE(TypeTag, IsothermalVolumeVariables));

//! The class for the quantities required for the flux calculation
SET_TYPE_PROP(RichardsTwoC, FluxVariables, typename GET_PROP_TYPE(TypeTag, IsothermalFluxVariables));

//! The class of the newton controller
SET_TYPE_PROP(RichardsTwoC, NewtonController, RichardsTwoCNewtonController<TypeTag>);

//! The upwind weight for the mass conservation equations
SET_SCALAR_PROP(RichardsTwoC, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(RichardsTwoC, ImplicitMobilityUpwindWeight, 1.0);

//! The class with all index definitions for the model
SET_TYPE_PROP(RichardsTwoC, Indices, typename GET_PROP_TYPE(TypeTag, IsothermalIndices));

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(RichardsTwoC, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(RichardsTwoC, EffectiveDiffusivityModel)
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
SET_PROP(RichardsTwoC, MaterialLawParams)
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
SET_PROP(RichardsTwoC, WettingPhase)
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
SET_PROP(RichardsTwoC, NonwettingPhase)
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
SET_PROP(RichardsTwoC, FluidSystem)
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
SET_BOOL_PROP(RichardsTwoC, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(RichardsTwoC, ProblemEnableGravity, true);

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
SET_BOOL_PROP(RichardsTwoCNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(RichardsTwoC, IsothermalModel, RichardsTwoCModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(RichardsTwoC, IsothermalFluxVariables, RichardsTwoCFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(RichardsTwoC, IsothermalVolumeVariables, RichardsTwoCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(RichardsTwoC, IsothermalLocalResidual, RichardsTwoCLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(RichardsTwoC, IsothermalIndices, RichardsTwoCIndices<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(RichardsTwoC, IsothermalNumEq, 2);

// \}
}

} // end namespace

#endif
