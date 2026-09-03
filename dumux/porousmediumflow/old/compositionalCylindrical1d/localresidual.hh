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
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 */

#ifndef DUMUX_COMPOSITIONAL_CYL_LOCAL_RESIDUAL_HH
#define DUMUX_COMPOSITIONAL_CYL_LOCAL_RESIDUAL_HH

#include <vector>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 */
template<class TypeTag>
class CompositionalLocalResidual: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
	using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
	using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = ModelTraits::useMoles();

    enum { conti0EqIdx = Indices::conti0EqIdx };

    //! The index of the component balance equation that gets replaced with the total mass balance
    static constexpr int replaceCompEqIdx = ModelTraits::replaceCompEqIdx();
    static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < numComponents;

public:

    using ParentType::ParentType;

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };
        // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
            	auto eqIdx = conti0EqIdx + compIdx;
            	Scalar b = problem.bufferPower(scv, volVars, compIdx);
                if (eqIdx != replaceCompEqIdx) {// all component except water
                    storage[eqIdx] += (volVars.porosity()*volVars.saturation(phaseIdx)+b)
                                      * massOrMoleDensity(volVars, phaseIdx)
                                      * massOrMoleFraction(volVars, phaseIdx, compIdx)*scv.center()[0];
                }
            }

            // for water
            if (useTotalMoleOrMassBalance) {
                storage[replaceCompEqIdx] += massOrMoleDensity(volVars, phaseIdx)
                                             *(volVars.porosity()*volVars.saturation(phaseIdx))*scv.center()[0];
            }

            //! The energy storage in the fluid phase with index phaseIdx
            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }

        //! The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        // get upwind weights into local scope
        NumEqVector flux(0.0);

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
			diffusiveFluxes *= scvf.center()[0];
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = conti0EqIdx + compIdx;

                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&massOrMoleDensity, &massOrMoleFraction, phaseIdx, compIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*massOrMoleFraction(volVars, phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                if (eqIdx != replaceCompEqIdx)
                    flux[eqIdx] += (fluxVars.advectiveFlux(phaseIdx, upwindTerm)*scvf.center()[0]); //

                // diffusive fluxes (only for the component balances)
                if(eqIdx != replaceCompEqIdx)
                    flux[eqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                            : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
            }

            // in case one balance is substituted by the total mole balance
            if (useTotalMoleOrMassBalance)
            {
                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&massOrMoleDensity, phaseIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*volVars.mobility(phaseIdx); };

                flux[replaceCompEqIdx] += (fluxVars.advectiveFlux(phaseIdx, upwindTerm)*scvf.center()[0]); // ;

                for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                    flux[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                                       : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
            }

            //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;
    }

    /*!
     * \brief Adds the storage derivative
     *
     * \param partialDerivatives The partial derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curVolVars The current volume variables
     * \param scv The sub control volume
     */
    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars,
                               const SubControlVolume& scv) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "richards/localresidual.hh: Analytic Jacobian only supports incompressible fluids!");
        // static_assert(useMoles,
                      // "richards/localresidual.hh: need to use moles!");

        const auto poreVolume = Extrusion::volume(fvGeometry, scv)*curVolVars.porosity()*curVolVars.extrusionFactor();
        static const auto rho = useMoles ? curVolVars.molarDensity(0) : curVolVars.density(0);

        // partial derivative of storage term w.r.t. p_w
        // d(Sw*rho*phi*V/dt)/dpw = rho*phi*V/dt*dsw/dpw = rho*phi*V/dt*dsw/dpc*dpc/dpw = -rho*phi*V/dt*dsw/dpc
        const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, InvalidElemSol{});
        partialDerivatives[conti0EqIdx][0] += -rho*poreVolume/this->timeLoop().timeStepSize()*fluidMatrixInteraction.dsw_dpc(curVolVars.capillaryPressure());
    }

    /*!
     * \brief Adds source derivatives for wetting and nonwetting phase.
     *
     * \param partialDerivatives The partial derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curVolVars The current volume variables
     * \param scv The sub control volume
     *
     * \todo Maybe forward to problem for the user to implement the source derivatives?
     */
    template<class PartialDerivativeMatrix>
    void addSourceDerivatives(PartialDerivativeMatrix& partialDerivatives,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& curVolVars,
                              const SubControlVolume& scv) const
    { /* TODO maybe forward to problem for the user to implement the source derivatives?*/ }

    /*!
     * \brief Adds flux derivatives for wetting and nonwetting phase for cell-centered FVM using TPFA
     *
     * Compute derivatives for the wetting and the nonwetting phase flux with respect to \f$p_w\f$
     * and \f$S_n\f$.
     *
     * \param derivativeMatrices The partial derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethods::cctpfa, void>
    addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "richards/localresidual.hh: Analytic Jacobian only supports incompressible fluids!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "richards/localresidual.hh: Analytic Jacobian only supports fluids with constant viscosity!");
        // static_assert(useMoles,
                      // "richards/localresidual.hh: need to use moles!");

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto outsideScvIdx = scvf.outsideScvIdx();
        const auto outsideElement = fvGeometry.gridGeometry().element(outsideScvIdx);
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[outsideScvIdx];

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho = useMoles ? insideVolVars.molarDensity(0) : insideVolVars.density(0); //insideVolVars.density(0);
        static const auto mu = insideVolVars.viscosity(0);
        static const auto rho_mu = rho/mu;

        // upwind term
        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");
        const auto flux = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 0, elemFluxVarsCache);
        const auto insideWeight = std::signbit(flux) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight = 1.0 - insideWeight;
        const auto upwindTerm = rho*insideVolVars.mobility(0)*insideWeight + rho*outsideVolVars.mobility(0)*outsideWeight;

        const auto insideFluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, insideScv, InvalidElemSol{});
        const auto outsideFluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(outsideElement, outsideScv, InvalidElemSol{});

        // material law derivatives
        const auto insideSw = insideVolVars.saturation(0);
        const auto outsideSw = outsideVolVars.saturation(0);
        const auto insidePc = insideVolVars.capillaryPressure();
        const auto outsidePc = outsideVolVars.capillaryPressure();
        const auto dkrw_dsw_inside = insideFluidMatrixInteraction.dkrw_dsw(insideSw);
        const auto dkrw_dsw_outside = outsideFluidMatrixInteraction.dkrw_dsw(outsideSw);
        const auto dsw_dpw_inside = -insideFluidMatrixInteraction.dsw_dpc(insidePc);
        const auto dsw_dpw_outside = -outsideFluidMatrixInteraction.dsw_dpc(outsidePc);

        // the transmissibility
        const auto tij = elemFluxVarsCache[scvf].advectionTij();

        // get references to the two participating derivative matrices
        auto& dI_dI = derivativeMatrices[insideScvIdx];
        auto& dI_dJ = derivativeMatrices[outsideScvIdx];

        // partial derivative of the wetting phase flux w.r.t. pw
        dI_dI[conti0EqIdx][0] += tij*upwindTerm + rho_mu*flux*insideWeight*dkrw_dsw_inside*dsw_dpw_inside;
        dI_dJ[conti0EqIdx][0] += -tij*upwindTerm + rho_mu*flux*outsideWeight*dkrw_dsw_outside*dsw_dpw_outside;
    }

    /*!
     * \brief Adds flux derivatives for box method
     *
     * \param A The Jacobian Matrix
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
    template<class JacobianMatrix, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethods::box, void>
    addFluxDerivatives(JacobianMatrix& A,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "richards/localresidual.hh: Analytic Jacobian only supports incompressible fluids!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "richards/localresidual.hh: Analytic Jacobian only supports fluids with constant viscosity!");
        // static_assert(useMoles,
                      // "richards/localresidual.hh: need to use moles!");

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto outsideScvIdx = scvf.outsideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[outsideScvIdx];

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho =  useMoles ? insideVolVars.molarDensity(0) : insideVolVars.density(0); //insideVolVars.density(0);
        static const auto mu = insideVolVars.viscosity(0);
        static const auto rho_mu = rho/mu;

        // upwind term
        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");
        const auto flux = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 0, elemFluxVarsCache);
        const auto insideWeight = std::signbit(flux) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight = 1.0 - insideWeight;
        const auto upwindTerm = rho*insideVolVars.mobility(0)*insideWeight + rho*outsideVolVars.mobility(0)*outsideWeight;

        const auto insideFluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, insideScv, InvalidElemSol{});
        const auto outsideFluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, outsideScv, InvalidElemSol{});

        // material law derivatives
        const auto insideSw = insideVolVars.saturation(0);
        const auto outsideSw = outsideVolVars.saturation(0);
        const auto insidePc = insideVolVars.capillaryPressure();
        const auto outsidePc = outsideVolVars.capillaryPressure();
        const auto dkrw_dsw_inside = insideFluidMatrixInteraction.dkrw_dsw(insideSw);
        const auto dkrw_dsw_outside = outsideFluidMatrixInteraction.dkrw_dsw(outsideSw);
        const auto dsw_dpw_inside = -insideFluidMatrixInteraction.dsw_dpc(insidePc);
        const auto dsw_dpw_outside = -outsideFluidMatrixInteraction.dsw_dpc(outsidePc);

        // so far it was the same as for tpfa
        // the transmissibilities (flux derivatives with respect to all pw-dofs on the element)
        const auto ti = AdvectionType::calculateTransmissibilities(
            problem, element, fvGeometry, curElemVolVars, scvf, elemFluxVarsCache[scvf]
        );

        // get the rows of the jacobian matrix for the inside/outside scv
        auto& dI_dJ_inside = A[insideScv.dofIndex()];
        auto& dI_dJ_outside = A[outsideScv.dofIndex()];

        // add the partial derivatives w.r.t all scvs in the element
        for (const auto& scvJ : scvs(fvGeometry))
        {
            // the transmissibily associated with the scvJ
            const auto& tj = ti[scvJ.indexInElement()];
            const auto globalJ = scvJ.dofIndex();

            // partial derivative of the wetting phase flux w.r.t. p_w
            dI_dJ_inside[globalJ][conti0EqIdx][0] += tj*upwindTerm;
            dI_dJ_outside[globalJ][conti0EqIdx][0] += -tj*upwindTerm;
        }

        const auto upwindContributionInside = rho_mu*flux*insideWeight*dkrw_dsw_inside*dsw_dpw_inside;
        const auto upwindContributionOutside = rho_mu*flux*outsideWeight*dkrw_dsw_outside*dsw_dpw_outside;

        // additional contribution of the upwind term only for inside and outside dof
        A[insideScv.dofIndex()][insideScv.dofIndex()][conti0EqIdx][0] += upwindContributionInside;
        A[insideScv.dofIndex()][outsideScv.dofIndex()][conti0EqIdx][0] += upwindContributionOutside;

        A[outsideScv.dofIndex()][insideScv.dofIndex()][conti0EqIdx][0] -= upwindContributionInside;
        A[outsideScv.dofIndex()][outsideScv.dofIndex()][conti0EqIdx][0] -= upwindContributionOutside;
    }

    /*!
     * \brief Adds cell-centered Dirichlet flux derivatives for wetting and nonwetting phase
     *
     * Compute derivatives for the wetting and the nonwetting phase flux with respect to \f$p_w\f$
     * and \f$S_n\f$.
     *
     * \param derivativeMatrices The matrices containing the derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
    template<class PartialDerivativeMatrices>
    void addCCDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                       const Problem& problem,
                                       const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       const ElementVolumeVariables& curElemVolVars,
                                       const ElementFluxVariablesCache& elemFluxVarsCache,
                                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "richards/localresidual.hh: Analytic Jacobian only supports incompressible fluids!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "richards/localresidual.hh: Analytic Jacobian only supports fluids with constant viscosity!");
		// static_assert(useMoles, "should use useMoles");


        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[scvf.outsideScvIdx()];
        const auto insideFluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, insideScv, InvalidElemSol{});

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho =  useMoles ? insideVolVars.molarDensity(0) : insideVolVars.density(0); //insideVolVars.density(0);
        static const auto mu = insideVolVars.viscosity(0);
        static const auto rho_mu = rho/mu;

        // upwind term
        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");
        const auto flux = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 0, elemFluxVarsCache);
        const auto insideWeight = std::signbit(flux) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight = 1.0 - insideWeight;
        const auto upwindTerm = rho*insideVolVars.mobility(0)*insideWeight + rho*outsideVolVars.mobility(0)*outsideWeight;

        // material law derivatives
        const auto insideSw = insideVolVars.saturation(0);
        const auto insidePc = insideVolVars.capillaryPressure();
        const auto dkrw_dsw_inside = insideFluidMatrixInteraction.dkrw_dsw(insideSw);
        const auto dsw_dpw_inside = -insideFluidMatrixInteraction.dsw_dpc(insidePc);

        // the transmissibility
        const auto tij = elemFluxVarsCache[scvf].advectionTij();

        // partial derivative of the wetting phase flux w.r.t. pw
        derivativeMatrices[insideScvIdx][conti0EqIdx][0] += tij*upwindTerm + rho_mu*flux*insideWeight*dkrw_dsw_inside*dsw_dpw_inside;
    }

    /*!
     * \brief Adds Robin flux derivatives for wetting and nonwetting phase
     *
     * \param derivativeMatrices The matrices containing the derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
    template<class PartialDerivativeMatrices>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    {}

protected:

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }
};

} // end namespace Dumux

#endif
