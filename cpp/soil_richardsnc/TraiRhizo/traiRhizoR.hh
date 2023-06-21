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
 * \ingroup Chemistry
 * \brief The source and sink terms due to reactions are calculated in this class.
 */
#ifndef DUMUX_BIOMIN_REACTIONS_HH
#define DUMUX_BIOMIN_REACTIONS_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup Chemistry
 * \brief The source and sink terms due to reactions are calculated in this class.
 */
template< class NumEqVector, class VolumeVariables >
class SimpleBiominReactions
{
    using SolidSystem = typename VolumeVariables::SolidSystem;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using Scalar = typename FluidSystem::Scalar;

public:

    SimpleBiominReactions()
    {   //ureolysis kinetic parameters
        kub_ = getParam<Scalar>("UreolysisCoefficients.Kub");
        kurease_ = getParam<Scalar>("UreolysisCoefficients.Kurease");
        ku_ = getParam<Scalar>("UreolysisCoefficients.Ku");
    }

    static constexpr int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
    static constexpr int numComponents = FluidSystem::numComponents;

    static constexpr int H2OIdx = FluidSystem::H2OIdx;
    static constexpr int CO2Idx = FluidSystem::CO2Idx;
    static constexpr int CaIdx = FluidSystem::CaIdx;
    static constexpr int UreaIdx = FluidSystem::UreaIdx;

    // phase indices when used for the SolidSystem context (e.g. density)
    static constexpr int BiofilmPhaseIdx = SolidSystem::BiofilmIdx;
    // overall indices when used in the problem context (source term)
    static constexpr int BiofilmIdx = SolidSystem::BiofilmIdx+numComponents;
    static constexpr int CalciteIdx = SolidSystem::CalciteIdx+numComponents;

    /*!
     * \brief Returns the molality of a component x (mol x / kg solvent) for a given
     * mole fraction (mol x / mol solution)
     * The salinity and the mole Fraction of CO2 are considered
     *
     */
    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        Scalar molalityX = moleFracX / (1 - moleFracSalinity - moleFracCTot) / FluidSystem::molarMass(H2OIdx);
        return molalityX;
    }

    /*!
     * \brief Calculate the source/sink term due to reactions.
     *
     * \param Source The source
     * \param volVars The volume variables
     */
    void reactionSource(NumEqVector &q, const VolumeVariables &volVars)
    {
        //   define and compute some parameters for convenience:
        Scalar xwCa = volVars.moleFraction(liquidPhaseIdx,CaIdx);
        Scalar densityBiofilm = volVars.solidComponentDensity(BiofilmPhaseIdx);

        Scalar volFracBiofilm = volVars.solidVolumeFraction(BiofilmPhaseIdx);

        if (volFracBiofilm < 0)
            volFracBiofilm = 0;

        // TODO: dumux-course-task
        // implement mass of biofilm
        Scalar massBiofilm = densityBiofilm * volFracBiofilm;

        Scalar molalityUrea = moleFracToMolality(volVars.moleFraction(liquidPhaseIdx,UreaIdx),
                                                 xwCa,
                                                 volVars.moleFraction(liquidPhaseIdx,CO2Idx));  // [mol_urea/kg_H2O]

        // TODO: dumux-course-task
        // compute rate of ureolysis by implementing Z_urease,biofilm and r_urea
        Scalar Zub = kub_ *  massBiofilm; // [kg urease/m³]
        Scalar rurea = kurease_ * Zub * molalityUrea / (ku_ + molalityUrea); // [mol/m³s]

        // compute/set dissolution and precipitation rate of calcite
        Scalar rprec = 0.0;
        rprec = rurea;

        // set source terms
        // TODO: dumux-course-task
        // update terms according to stochiometry
        q[H2OIdx]     += 0.0;
        q[CO2Idx]     += rurea - rprec ;
        q[CaIdx]      += - rprec;
        q[UreaIdx]    += - rurea;
        q[BiofilmIdx] += 0.0;
        q[CalciteIdx] += rprec;
    }

private:
    // urease parameters
    Scalar kub_;
    Scalar kurease_;
    Scalar ku_;
};

} //end namespace Dumux

#endif
