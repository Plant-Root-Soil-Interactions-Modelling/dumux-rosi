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
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::FluidSystems::H2OSolute
 */
#ifndef DUMUX_H2O_SOLUTE_FLUID_SYSTEM_HH
#define DUMUX_H2O_SOLUTE_FLUID_SYSTEM_HH

//#include <dumux/material/idealgas.hh>
#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/richardsNCsolute.hh>

#include <dumux/material/binarycoefficients/richardsNCh2oSolute.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/material/fluidsystems/base.hh>

#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A compositional fluid system with water and solute
 *        components in both the liquid and the gas phase.
 */
template <class Scalar,
         // class H2OType = Dumux::Components::TabulatedComponent<Dumux::Components::H2O<Scalar> > >
            class H2OType = Dumux::Components::SimpleH2O<Scalar> >
class H2OSOLUTE
    : public Base<Scalar, H2OSOLUTE<Scalar, H2OType> >
{
    using ThisType = H2OSOLUTE<Scalar, H2OType>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

public:
    using Solute = Dumux::Components::SOLUTE<Scalar>;
    using H2O = H2OType;


    static const int numPhases = 1;
    static const int numComponents = 2;

    static const int wPhaseIdx = 0; // index of the water phase
    //static const int nPhaseIdx = 1; // index of the NAPL phase
    //static const int gPhaseIdx = 2; // index of the gas phase

    static const int H2OIdx = 0;
    static const int SoluteIdx = 1;

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static const int wCompIdx = H2OIdx;
    static const int nCompIdx = SoluteIdx;

    /*!
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
       // if (H2O::isTabulated)
       // {
       //     H2O::init(tempMin, tempMax, nTemp,
       //               pressMin, pressMax, nPress);
       // }
    }

    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        switch (phaseIdx)
        {
            case wPhaseIdx: return IOName::aqueousPhase();
           // case nPhaseIdx: return IOName::naplPhase();
           // case gPhaseIdx: return IOName::gaseousPhase();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the human readable name of a component (used in indices)
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx) {
            case H2OIdx: return H2O::name();
            case SoluteIdx: return Solute::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
            case H2OIdx: return H2O::molarMass();
            case SoluteIdx: return Solute::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState, int phaseIdx)
    { return 1000;   }

    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    { return 1e-3;  }

  /*  using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState, int phaseIdx)
    { 
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        return BinaryCoeff::H2O_SOLUTE::liquidDiffCoeff(T, p);   }*/
};
} // end namespace FluidSystems
} // end namespace Dumux

#endif
