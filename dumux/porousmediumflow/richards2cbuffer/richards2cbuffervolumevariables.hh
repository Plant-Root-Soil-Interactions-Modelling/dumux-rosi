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
 * \file
 *
 * \brief Volume averaged quantities required by the RichardsTwoC model.
 */
#ifndef DUMUX_RICHARDS_2C_BUFFER_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDS_2C_BUFFER_VOLUME_VARIABLES_HH

#include "richards2cbufferproperties.hh"

#include <dumux/implicit/volumevariables.hh>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsTwoCModel
 * \ingroup ImplicitVolumeVariables
 * \brief Volume averaged quantities required by the RichardsTwoC model.
 *
 * This contains the quantities which are are constant within a finite
 * volume in the RichardsTwoCBuffer model
 */
template <class TypeTag>
class RichardsTwoCBufferVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffectiveDiffusivityModel;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    static const bool usePH = GET_PROP_VALUE(TypeTag, UsePH);
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    //! The type returned by the fluidState() method
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    //indices of primary variables
    enum{
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx
    };

    enum{
         phaseIdx = Indices::phaseIdx,
         phaseCompIdx = Indices::phaseCompIdx,
         transportCompIdx = Indices::transportCompIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> DimVector;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:


    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx,
                const bool isOldSol)
    {
        ParentType::update(priVars, problem, element, fvGeometry, scvIdx, isOldSol);

        //calculate all secondary variables from the primary variables and store results in fluidstate
        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, fluidState_);

        const MaterialLawParams &matParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);
        relativePermeabilityWetting_ =
            MaterialLaw::krw(matParams,
                             fluidState_.saturation(phaseIdx));

        porosity_ = problem.spatialParams().porosity(element, fvGeometry, scvIdx);

        dispersivity_ = problem.spatialParams().dispersivity(element, fvGeometry, scvIdx);
        buffer_ = problem.spatialParams().buffer(element, fvGeometry, scvIdx);


        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, phaseIdx);

        diffCoeff_ = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                             paramCache,
                                                             phaseIdx,
                                                             phaseCompIdx,
                                                             transportCompIdx);
        effDiffCoeff_=EffectiveDiffusivityModel::effectiveDiffusivity(porosity_,
                                                                     fluidState_.saturation(phaseIdx),
                                                                     diffCoeff_)
                        /(fluidState_.saturation(phaseIdx)* porosity_+buffer_);
        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(dispersivity_);
        Valgrind::CheckDefined(diffCoeff_);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const int scvIdx,
                                   FluidState& fluidState)
    {
        // temperature
        Scalar t = Implementation::temperature_(priVars, problem, element,
                                                fvGeometry, scvIdx);
        fluidState.setTemperature(t);

        const MaterialLawParams &matParams =
                problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

        fluidState.setPressure(phaseIdx, priVars[pressureIdx]);

        //Scalar minPc = MaterialLaw::pc(matParams, 1.0);

        //fluidState.setPressure(nPhaseIdx, std::max(pnRef, priVars[pwIdx] + minPc));

        // saturations
        Scalar pnRef = problem.referencePressure(element, fvGeometry, scvIdx);
        Scalar sw = MaterialLaw::sw(matParams, pnRef - fluidState.pressure(phaseIdx));
        fluidState.setSaturation(phaseIdx, sw);

        if(useMoles)
        {
            fluidState.setMoleFraction(phaseIdx, phaseCompIdx, 1 - priVars[massOrMoleFracIdx]);
            fluidState.setMoleFraction(phaseIdx, transportCompIdx, priVars[massOrMoleFracIdx]);
        }
        else
        {
            // setMassFraction() has only to be called 1-numComponents times
            fluidState.setMassFraction(phaseIdx, phaseCompIdx, 1 - priVars[massOrMoleFracIdx]);
            //fluidState.setMassFraction(phaseIdx, transportCompIdx, priVars[massOrMoleFracIdx]);
        }

        // density and viscosity
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        Scalar value;
        value = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, value);
        value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
        fluidState.setViscosity(phaseIdx, value);
        // compute and set the enthalpy
        Scalar h = Implementation::enthalpy_(fluidState, paramCache, phaseIdx);
        fluidState.setEnthalpy(phaseIdx, h);

    //    Scalar x1 = priVars[massOrMoleFracIdx]; //mole or mass fraction of component 1
    //    if(!useMoles) //mass-fraction formulation
    //    {
    //        // convert mass to mole fractions
    //        Scalar M0 = FluidSystem::molarMass(phaseCompIdx);
    //        Scalar M1 = FluidSystem::molarMass(transportCompIdx);
    //        //meanMolarMass if x1_ is a massfraction
    //        Scalar meanMolarMass = M0*M1/(M1 + x1*(M0 - M1));
    //        x1 *= meanMolarMass/M1;
    //    }
    //    fluidState.setMoleFraction(phaseIdx, phaseCompIdx, 1 - x1);
    //    fluidState.setMoleFraction(phaseIdx, transportCompIdx, x1);
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     */
    Scalar molarDensity() const
    { return fluidState_.molarDensity(phaseIdx);}

    /*!
     * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar moleFraction(int compIdx) const
    { return fluidState_.moleFraction(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar massFraction(int compIdx) const
    { return fluidState_.massFraction(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar molarity(int compIdx) const
    { return fluidState_.molarity(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */
    Scalar diffCoeff() const
    { return diffCoeff_; }

    Scalar effDiffCoeff() const
    {
        return effDiffCoeff_;
    }

    /*!
     * \brief Returns the dispersivity of the fluid's streamlines.
     */
    const GlobalPosition &dispersivity() const
    { return dispersivity_; }

    /*!
     * \brief Returns the average porosity [] within the control volume.
     *
     * The porosity is defined as the ratio of the pore space to the
     * total volume, i.e. \f[ \Phi := \frac{V_{pore}}{V_{pore} + V_{rock}} \f]
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the average absolute saturation [] of a given
     *        fluid phase within the finite volume.
     *
     * The saturation of a fluid phase is defined as the fraction of
     * the pore volume filled by it, i.e.
     * \f[ S_\alpha := \frac{V_\alpha}{V_{pore}} = \phi \frac{V_\alpha}{V} \f]
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the average mass density \f$\mathrm{[kg/m^3]}\f$ of a given
     *        fluid phase within the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar density() const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * For the non-wetting phase (i.e. the gas phase), we assume
     * infinite mobility, which implies that the non-wetting phase
     * pressure is equal to the finite volume's reference pressure
     * defined by the problem.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar pressure() const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns average temperature \f$\mathrm{[K]}\f$ inside the control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Returns the effective mobility \f$\mathrm{[1/(Pa*s)]}\f$ of a given phase within
     *        the control volume.
     *
     * The mobility of a fluid phase is defined as the relative
     * permeability of the phase (given by the chosen material law)
     * divided by the dynamic viscosity of the fluid, i.e.
     * \f[ \lambda_\alpha := \frac{k_{r\alpha}}{\mu_\alpha} \f]
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar mobility(const int phaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of a given phase
     *        within the control volume.
     */
    Scalar viscosity() const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns relative permeability [-] of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar relativePermeability(const int phaseIdx) const
    {
        return relativePermeabilityWetting_;
    }

    /*!
     * \brief Returns the water content
     *        fluid phase within the finite volume.
     *
     * The water content is defined as the fraction of
     * the saturation devided by the porosity

     * \param phaseIdx The index of the fluid phase
     */
    Scalar waterContent (const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx)* porosity_; }

    /*!
     * \brief Returns the bufferpower of component within the control volume.
     *
     */
    Scalar buffer() const
    { return buffer_; }

protected:
    static Scalar temperature_(const PrimaryVariables &primaryVariables,
                            const Problem& problem,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const int scvIdx)
    {
        return problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            int phaseIdx)
    {
        return 0;
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &priVars,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       const bool isOldSol)
    { }

    FluidState fluidState_;
    Scalar relativePermeabilityWetting_;
    Scalar porosity_;
    GlobalPosition dispersivity_;
    Scalar diffCoeff_;
    Scalar buffer_;
    Scalar effDiffCoeff_;


private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
