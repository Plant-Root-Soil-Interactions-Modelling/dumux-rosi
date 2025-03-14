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

//#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
//#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidstates/compositional.hh>
//#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/spatialparams/implicit1p.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>

namespace Dumux
{
// \{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////
//! Number of equations required by the model
SET_INT_PROP(RichardsTwoCBuffer, NumEq, 2);
//! Number of fluid phases considered
SET_INT_PROP(RichardsTwoCBuffer, NumPhases, 1);
SET_INT_PROP(RichardsTwoCBuffer, NumComponents, 2); //!< The number of components in the 1p2c model is 2
SET_SCALAR_PROP(RichardsTwoCBuffer, Scaling, 1); //!< Scaling of the model is set to 1 by default
SET_BOOL_PROP(RichardsTwoCBuffer, UseMoles, false); //!< Define that mole fractions are used in the balance equations

//! Use pressure [Pa] by default
SET_BOOL_PROP(RichardsTwoCBuffer, UsePH, false);

//! The local residual operator
SET_TYPE_PROP(RichardsTwoCBuffer,
              LocalResidual,
              RichardsTwoCBufferLocalResidual<TypeTag>);

//! The global model used
SET_TYPE_PROP(RichardsTwoCBuffer, Model, RichardsTwoCBufferModel<TypeTag>);

//! The class for the volume averaged quantities
SET_TYPE_PROP(RichardsTwoCBuffer, VolumeVariables, RichardsTwoCBufferVolumeVariables<TypeTag>);

//! The class for the quantities required for the flux calculation
SET_TYPE_PROP(RichardsTwoCBuffer, FluxVariables, RichardsTwoCBufferFluxVariables<TypeTag>);

//! The class of the newton controller
SET_TYPE_PROP(RichardsTwoCBuffer, NewtonController, RichardsTwoCBufferNewtonController<TypeTag>);

//! The upwind weight for the mass conservation equations
SET_SCALAR_PROP(RichardsTwoCBuffer, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(RichardsTwoCBuffer, ImplicitMobilityUpwindWeight, 1.0);

//! The class with all index definitions for the model
SET_TYPE_PROP(RichardsTwoCBuffer, Indices, RichardsTwoCBufferIndices<TypeTag>);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(RichardsTwoCBuffer, SpatialParams, ImplicitSpatialParamsOneP<TypeTag>);

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
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(RichardsTwoCBuffer, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Set the phaseIndex per default to zero (important for two-phase fluidsystems).
SET_INT_PROP(RichardsTwoCBuffer, PhaseIdx, 0);

// disable velocity output by default
SET_BOOL_PROP(RichardsTwoCBuffer, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(RichardsTwoCBuffer, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(RichardsTwoCBuffer, SpatialParamsForchCoeff, 0.55);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(RichardsTwoCBufferNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

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
