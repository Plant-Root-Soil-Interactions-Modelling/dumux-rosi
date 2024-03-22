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
 * \ingroup SolidSystems
 * \brief A solid phase consisting of multiple inert solid components.
 */
#ifndef DUMUX_SOLIDSYSTEMS_COMPOSITIONAL_SOLID_PHASE7C_HH
#define DUMUX_SOLIDSYSTEMS_COMPOSITIONAL_SOLID_PHASE7C_HH

#include <string>
#include <dune/common/exceptions.hh>

namespace Dumux {
namespace SolidSystems {

/*!
 * \ingroup SolidSystems
 * \brief A solid phase consisting of multiple inert solid components + main component.
 * \note a solid is considered inert if it cannot dissolve in a liquid and
 *       and cannot increase its mass by precipitation from a fluid phase.
 * \note inert components have to come after all non-inert components
 */
template <class Scalar, class Component1, class Component2, class Component3,
			class Component4, class Component5, class Component6, class Component7>
class SolidPhaseSevenC
{
public:
    using ComponentOne = Component1;
    using ComponentTwo = Component2;
    using ComponentThree = Component3;
    using ComponentFour = Component4;
    using ComponentFive = Component5;
    using ComponentSix = Component6;
    using ComponentSeven = Component7;
    using MainComponent = Component7;


    /****************************************
     * Solid phase related static parameters
     ****************************************/
    static constexpr int numComponents = 7;
    static constexpr int numInertComponents = 1;
    static constexpr int mainCompIdx = 6;
    static constexpr int comp0Idx = 0;
    static constexpr int comp1Idx = 1;
    static constexpr int comp2Idx = 2;
    static constexpr int comp3Idx = 3;
    static constexpr int comp4Idx = 4;
    static constexpr int comp5Idx = 5;
    static constexpr int comp6Idx = 6;


    /*!
     * \brief Return the human readable name of a solid phase
     *
     * \param compIdx The index of the solid phase to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return ComponentOne::name();
            case comp1Idx: return ComponentTwo::name();
            case comp2Idx: return ComponentThree::name();
            case comp3Idx: return ComponentFour::name();
            case comp4Idx: return ComponentFive::name();
            case comp5Idx: return ComponentSix::name();
            case comp6Idx: return ComponentSeven::name();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief A human readable name for the solid system.
     */
    static std::string name()
    { return "s"; }

    /*!
     * \brief Returns whether the phase is incompressible
     */
    static constexpr bool isCompressible(int compIdx)
    { return false; }

    /*!
     * \brief Returns whether the component is inert (doesn't react)
     */
    static constexpr bool isInert()
    { return (numComponents == numInertComponents); }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
            case mainCompIdx: return MainComponent::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "molarMass, Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState& solidState)
    {
		return MainComponent::solidDensity();
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState& solidState, const int compIdx)
    {
        switch (compIdx)
        {
            case mainCompIdx: return MainComponent::solidDensity(solidState.temperature());
            default: DUNE_THROW(Dune::InvalidStateException, "density, Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief The molar density of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar molarDensity(const SolidState& solidState, const int compIdx)
    {
        switch (compIdx)
        {
            case mainCompIdx: return MainComponent::solidDensity(solidState.temperature())/MainComponent::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "molarDensity, Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief Thermal conductivity of the solid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class SolidState>
    static Scalar thermalConductivity(const SolidState &solidState)
    {
        return MainComponent::solidThermalConductivity(solidState.temperature());
    }

    /*!
     * \brief Specific isobaric heat capacity of the pure solids \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class SolidState>
    static Scalar heatCapacity(const SolidState &solidState)
    {
        return MainComponent::solidHeatCapacity(solidState.temperature());
    }

};

} // end namespace SolidSystems
} // end namespace Dumux

#endif
