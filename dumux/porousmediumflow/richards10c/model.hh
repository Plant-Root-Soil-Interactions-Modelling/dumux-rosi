// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Richards10CModel
 * \brief Base class for all models which use the Richards,
 *        n-component fully implicit model.
 *
 * This extension of Richards' equation, allows for
 * the wetting phase to consist of multiple components:
 *\f{eqnarray*}
 && \frac{\partial (\sum_w \varrho_w X_w^\kappa \phi S_w )}
 {\partial t}
 - \sum_w  \nabla \cdot \left\{ \varrho_w X_w^\kappa
 \frac{k_{rw}}{\mu_w} \mathbf{K}
 (\nabla  p_w - \varrho_{w}  \mathbf{g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_w \nabla \cdot \left\{{\bf D_{w, pm}^\kappa} \varrho_{w} \nabla  X^\kappa_{w} \right\}
 - \sum_w q_w^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 w \in \{w, g\},
 \f}
 * where:
 * * \f$ \phi \f$ is the porosity of the porous medium,
 * * \f$ S_w \f$ represents the saturation of the wetting phase,
 * * \f$ \varrho_w \f$ is the mass density of the wetting phase,
 * * \f$ k_{rw} \f$ is the relative permeability of the wetting phase,
 * * \f$ \mu_w \f$ is the dynamic viscosity of the wetting phase,
 * * \f$ \mathbf{K} \f$ is the intrinsic permeability tensor,
 * * \f$ p_w \f$ is the pressure of the wetting phase,
 * * \f$ \mathbf{g} \f$ is the gravitational acceleration vector,
 * * \f$ \bf D_{w,pm}^{k} \f$ is the effective diffusivity of component \f$ \kappa \f$ in the wetting phase,
 * * \f$ X_w^k \f$ is the mass fraction of component \f$ \kappa \f$ in the wetting phase,
 * * \f$ q_w \f$ is a source or sink term in the wetting phase.
 */

#ifndef DUMUX_RICHARDS10C_MODEL_HH
#define DUMUX_RICHARDS10C_MODEL_HH

#include <dumux/common/properties.hh>


#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
//#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>

#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/soil.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/solidConstant.hh>
#include <dumux/material/fluidsystems/liquidphase3c.hh>
#include <dumux/material/fluidstates/compositionalnc.hh>

//actially this works for any number of components.
#include <dumux/material/solidstates/compositionalNCsolidstate.hh> 
#include <dumux/material/solidsystems/solidphase7c.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/flux/referencesystemformulation.hh>
																
												   

#include <dumux/porousmediumflow/richardsCylindrical1d/model.hh>

#include <dumux/porousmediumflow/richardsnc/indices.hh>
#include <dumux/porousmediumflow/richardsnc/iofields.hh>
#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

/*!
 * \ingroup Richards10CModel
 * \brief Specifies a number properties of the Richards n-components model.
 *
 * \tparam nComp the number of components to be considered.
 * \tparam useMol whether to use mass or mole balances
 */
template<int nComp, bool useMol, int nCompS, int repCompEqIdx = nComp>
struct Richards10CModelTraits
{
    using Indices = RichardsNCIndices;

    static constexpr int numEq() { return nComp + nCompS - numInertSolidComps(); }
    static constexpr int numFluidPhases() { return 1; }
    static constexpr int numFluidComponents() { return nComp; }
    static constexpr int numSolidComponents() { return nCompS; }
    static constexpr int numInertSolidComps() { return 1; }
    static constexpr int replaceCompEqIdx() { return repCompEqIdx; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }
    static constexpr bool enableCompositionalDispersion() { return false; }
    static constexpr bool enableThermalDispersion() { return false; }

    static constexpr bool useMoles() { return true;}//useMol; }
    static constexpr ReferenceSystemFormulation getReferenceSystemFormulation() { return ReferenceSystemFormulation::molarAveraged;}
	
};

/*!
 * \ingroup Richards10CModel
 * \brief Traits class for the Richards n-components model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 * \tparam DT The diffusion type
 * \tparam EDM The effective diffusivity model
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT, class DT, class EDM>
struct Richards10CVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
    using DiffusionType = DT;
    using EffectiveDiffusivityModel = EDM;
};

namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
// Create new type tags
namespace TTag {
struct Richards10C { using InheritsFrom = std::tuple<PorousMediumFlow>; };
struct Richards10CNI { using InheritsFrom = std::tuple<Richards10C>; };
} // end namespace TTag
//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the model traits class
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::Richards10C>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = Richards10CModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>(), 
											SolidSystem::numComponents,getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
};
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Richards10C> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; };

//! Define that per default mole fractions are used in the balance equations
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Richards10C> { static constexpr bool value = true; };

//! Redefine Fick s law to use moles
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::Richards10C> 
{ 
	using type = FicksLaw<TypeTag, ReferenceSystemFormulation::molarAveraged >; 
};

//! Use the dedicated local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Richards10C> { using type = CompositionalLocalResidual<TypeTag>; };

//! We set the replaceCompIdx to 0, i.e. the first equation is substituted with
//! the total mass balance, i.e. the phase balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::Richards10C> { static constexpr int value = 0; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Richards10C>
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
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using NCTraits = Richards10CVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;

public:
    using type = Richards10CVolumeVariables<NCTraits>;
};

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the liquid phase fluid system with simple H2O.
 */
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Richards10C>
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
struct FluidState<TypeTag, TTag::Richards10C>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};


//! per default solid state is inert
template<class TypeTag>
struct SolidState<TypeTag, TTag::Richards10C>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = CompositionalNCSolidState<Scalar, SolidSystem>;
};

// one main constant component
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Richards10C>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::Soil<Scalar>;
    //using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
	using type = SolidSystems::SolidPhaseSevenC<Scalar,
													 Components::solidConstant<0, Scalar>,
													 Components::solidConstant<1, Scalar>,
													 Components::solidConstant<2, Scalar>,
													 Components::solidConstant<3, Scalar>,
													 Components::solidConstant<4, Scalar>,
													 Components::solidConstant<5, Scalar>,
													 InertComponent>;
};

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::Richards10C> { using type = RichardsNCIOFields; };

//! The model after Millington (1961) is used for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::Richards10C> { using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//! average is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::Richards10CNI> { using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>; };

//////////////////////////////////////////////////////////////////
// Property values for non-isothermal Richards n-components model
//////////////////////////////////////////////////////////////////

//! set non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Richards10CNI>
{
private:
    using IsothermalTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Richards10CNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using BaseTraits = Richards10CVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;

    using ETCM = GetPropType< TypeTag, Properties::ThermalConductivityModel>;
    template<class BaseTraits, class ETCM>
    struct NCNITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };
public:
    using type = Richards10CVolumeVariables<NCNITraits<BaseTraits, ETCM>>;
};


} // end namespace Properties
} // end namespace Dumux

#endif
