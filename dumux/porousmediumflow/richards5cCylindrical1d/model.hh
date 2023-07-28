// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \file
 * \ingroup Richards5CModel
 * \brief Base class for all models which use the Richards,
 *        n-component fully implicit model.
 *
 * In the unsaturated zone, Richards' equation
 *\f{eqnarray*}
 && \frac{\partial (\sum_w \varrho_w X_w^\kappa \phi S_w )}
 {\partial t}
 - \sum_w  \text{div} \left\{ \varrho_w X_w^\kappa
 \frac{k_{rw}}{\mu_w} \mbox{\bf K}
 (\text{grad}\, p_w - \varrho_{w}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_w \text{div} \left\{{\bf D_{w, pm}^\kappa} \varrho_{w} \text{grad}\, X^\kappa_{w} \right\}
 - \sum_w q_w^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 w \in \{w, g\}
 \f}
 * is frequently used to
 * approximate the water distribution above the groundwater level.
 *
 * In contrast to the full two-phase model, the Richards model assumes
 * gas as the non-wetting fluid and that it exhibits a much lower
 * viscosity than the (liquid) wetting phase. (For example at
 * atmospheric pressure and at room temperature, the viscosity of air
 * is only about \f$1\%\f$ of the viscosity of liquid water.) As a
 * consequence, the \f$\frac{k_{r\alpha}}{\mu_\alpha}\f$ term
 * typically is much larger for the gas phase than for the wetting
 * phase. For this reason, the Richards model assumes that
 * \f$\frac{k_{rn}}{\mu_n}\f$ is infinitely large. This implies that
 * the pressure of the gas phase is equivalent to the static pressure
 * distribution and that therefore, mass conservation only needs to be
 * considered for the wetting phase.
 *
 * The model thus chooses the absolute pressure of the wetting phase
 * \f$p_w\f$ as its only primary variable. The wetting phase
 * saturation is calculated using the inverse of the capillary
 * pressure, i.e.
 \f[
 S_w = p_c^{-1}(p_n - p_w)
 \f]
 * holds, where \f$p_n\f$ is a given reference pressure. Nota bene,
 * that the last step is assumes that the capillary
 * pressure-saturation curve can be uniquely inverted, so it is not
 * possible to set the capillary pressure to zero when using the
 * Richards model!
 */

#ifndef DUMUX_RICHARDS5C_MODEL_CYL_HH
#define DUMUX_RICHARDS5C_MODEL_CYL_HH

#include <dumux/common/properties.hh>


#include <dumux/material/spatialparams/fv1p.hh>

#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
//#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>

#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/soil.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/solidConstant.hh>
#include <dumux/material/fluidsystems/liquidphase3c.hh>
#include <dumux/material/fluidstates/compositionalnc.hh>

#include <dumux/material/solidstates/compositional3Csolidstate.hh>
#include <dumux/material/solidsystems/solidphase3c.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include <dumux/porousmediumflow/richardsCylindrical1d/model.hh>

#include <dumux/porousmediumflow/richardsnc/indices.hh>
#include <dumux/porousmediumflow/richardsnc/iofields.hh>
#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

/*!
 * \ingroup Richards5CModel
 * \brief Specifies a number properties of the Richards n-components model.
 *
 * \tparam nComp the number of components to be considered.
 * \tparam useMol whether to use mass or mole balances
 */
template<int nComp, bool useMol, int nCompS, int repCompEqIdx = nComp >
struct Richards5CModelTraits
{
    using Indices = RichardsNCIndices;

    static constexpr int numFluidPhases() { return 1; }
    static constexpr int numFluidComponents() { return nComp; }
    static constexpr int numSolidComponents() { return nCompS; }
    //static constexpr int soilIdx() { return nCompS; }
    static constexpr int numInertSolidComps() { return 1; }
    static constexpr int replaceCompEqIdx() { return repCompEqIdx; }
    static constexpr int numEq() { return nComp + nCompS - numInertSolidComps(); }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }

    static constexpr bool useMoles() { return useMol; }
};

namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
// Create new type tags
namespace TTag {
struct Richards5C { using InheritsFrom = std::tuple<PorousMediumFlow>; };
struct Richards5CNI { using InheritsFrom = std::tuple<Richards5C>; };
} // end namespace TTag
//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the model traits class
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::Richards5C>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
//,getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()
    using type = Richards5CModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>(),SolidSystem::numComponents>;
};
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Richards5C> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; };

//! Define that per default mole fractions are used in the balance equations
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Richards5C> { static constexpr bool value = true; };

//! Use the dedicated local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Richards5C> { using type = CompositionalLocalResidual<TypeTag>; };

//! We set the replaceCompIdx to 0, i.e. the first equation is substituted with
//! the total mass balance, i.e. the phase balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::Richards5C> { static constexpr int value = 0; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Richards5C>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");

    using Traits = RichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = Richards5CVolumeVariables<Traits>;
};

//! The default richardsnc model computes no diffusion in the air phase
//! Turning this on leads to the extended Richards equation (see e.g. Vanderborght et al. 2017)
template<class TypeTag>
struct EnableWaterDiffusionInAir<TypeTag, TTag::Richards5C> { static constexpr bool value = false; };

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the liquid phase fluid system with simple H2O.
 */
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Richards5C>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseThreeC<Scalar, Components::SimpleH2O<Scalar>, 
														 Components::Constant<1, Scalar>,
													     Components::Constant<2, Scalar>>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::Richards5C>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};



//! per default solid state is inert
template<class TypeTag>
struct SolidState<TypeTag, TTag::Richards5C>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = Compositional3CSolidState<Scalar, SolidSystem>;
};

// one main constant component
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Richards5C>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::Soil<Scalar>;
    //using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
	using type = SolidSystems::SolidPhaseThreeC<Scalar,
													 Components::solidConstant<0, Scalar>,
													 Components::solidConstant<1, Scalar>,
													 Components::solidConstant<2, Scalar>,
													 Components::solidConstant<3, Scalar>,
													 InertComponent>;
};


//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::Richards5C> { using type = RichardsNCIOFields; };

////! The model after Millington (1961) is used for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::Richards5C> { using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//template<class TypeTag>
//struct EffectiveDiffusivityModel<TypeTag, TTag::Richards5C> { using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };


//! average is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::Richards5CNI> { using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>; };

//////////////////////////////////////////////////////////////////
// Property values for non-isothermal Richards n-components model
//////////////////////////////////////////////////////////////////

//! set non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Richards5CNI>
{
private:
    using IsothermalTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};


} // end namespace Properties
} // end namespace Dumux

#endif
