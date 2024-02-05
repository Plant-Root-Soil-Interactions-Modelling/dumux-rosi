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
/**
 * \file
 * \ingroup TracerTests
 * \brief Definition of a problem for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */

#ifndef DUMUX_TRACER_TEST_PROBLEM_HH
#define DUMUX_TRACER_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "tracerparams.hh"

namespace Dumux {
/**
 * \ingroup TracerTests
 * \brief Definition of a problem for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
template <class TypeTag>
class TracerTest;

namespace Properties {

// Create new type tags
namespace TTag {
struct TracerTest { using InheritsFrom = std::tuple<Tracer>; };
struct TracerTestTpfa { using InheritsFrom = std::tuple<TracerTest, CCTpfaModel>; };
struct TracerTestMpfa { using InheritsFrom = std::tuple<TracerTest, CCMpfaModel>; };
struct TracerTestBox { using InheritsFrom = std::tuple<TracerTest, BoxModel>; };
} // end namespace TTag

// enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTest> { using type = TracerTest<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTest>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerTestSpatialParams<FVGridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTest> { static constexpr bool value = false; };

//! A simple fluid system with one tracer component
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>,
                                                               TracerFluidSystem<TypeTag>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    static constexpr bool isTracerFluidSystem() { return true; }

    //! None of the components are the main component of the phase
    static constexpr int getMainComponent(int phaseIdx) { return -1; }

    //! The number of components
    static constexpr int numComponents = 1;
    static constexpr int numPhases = 1;

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx) { return "tracer_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0) { return "Water"; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx) { return 0.300; }

    //! Binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv) { return problem.diffusionCoefficient; }

    static constexpr bool isCompressible(int phaseIdx) { return false; } ///< \copydoc Dumux::FluidSystems::Base::isCompressible

    static constexpr bool viscosityIsConstant(int phaseIdx) { return true; } ///< \copydoc  Dumux::FluidSystems::Base::viscosityIsConstant

};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTest> { using type = TracerFluidSystem<TypeTag>; };

} // end namespace Properties


/*!
 * \ingroup TracerTests
 *
 * \brief Definition of a problem for the tracer problem:
 * A lens of contaminant tracer is diluted by diffusion and a base groundwater flow
 *
 * This problem uses the \ref TracerModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxtracer -ParameterFile ./test_boxtracer.input</tt> or
 * <tt>./test_cctracer -ParameterFile ./test_cctracer.input</tt>
 */
template <class TypeTag>
class TracerTest : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename FVGridGeometry::GridView;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
//    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>; // wrong todo
//    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>; // wrong todo
//    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;

    enum {
        dimWorld = GridView::dimensionworld, // world dimension
        isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::Box // discretization method
    };

    enum BCTypes {
        constantConcentration = 1,
        constantFlux = 2,
		linear = 3,
		michaelisMenten = 4,
		outflow = 5
    };

public:

    TracerTest(std::shared_ptr<const FVGridGeometry> fvGridGeom)
	: PorousMediumFlowProblem<TypeTag>(fvGridGeom) {

    	// BC
        bcTopType_ = getParam<int>("Soil.BC.Top.Type"); // todo type as a string might be nicer
        bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
        bcTopValue_ = getParam<Scalar>("Soil.BC.Top.Value",0.);
        bcBotValue_ = getParam<Scalar>("Soil.BC.Bot.Value",0.);

        // IC
        initialSoil_ = InputFileFunction("Soil.IC", "P", "Z", 0.); // [cm]([m]) pressure head, conversions hard coded

        std::cout << "TracerProblem constructed: bcTopType " << bcTopType_ << ", " << bcTopValue_ << "; bcBotType "
            <<  bcBotType_ << ", " << bcBotValue_ << "\n" << std::flush;
    }

    /*!
     * \copydoc FVProblem::initial
     *
     * called by FVProblem::applyInitialSolution(...)
     */
    template<class Entity>
    PrimaryVariables initial(const Entity& entity) const {
        auto eIdx = this->gridGeometry().elementMapper().index(entity);
        Scalar z = entity.geometry().center()[dimWorld - 1];
        // std::cout << "tracer initial " << z << ", " << initialSoil_.f(z,eIdx) << " \n";
        return PrimaryVariables(toPa_(initialSoil_.f(z, eIdx)));
    }


    /*!
     * @copydoc FVProblem::boundaryTypesAtPos
     *
     * discretization dependent, e.g. called by BoxElementBoundaryTypes::boundaryTypes(...)
     * when?
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const {
        BoundaryTypes bcTypes;
        if (onUpperBoundary_(globalPos)) { // top or outer bc
            switch (bcTopType_) {
            case constantConcentration:
                bcTypes.setAllDirichlet();
                break;
            case constantFlux:
                bcTypes.setAllNeumann();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Top or outer boundary type not implemented");
            }
        } else if (onLowerBoundary_(globalPos)) { // bot bc
            switch (bcBotType_) {
            case constantConcentration:
                bcTypes.setAllDirichlet();
                break;
            case constantFlux:
                bcTypes.setAllNeumann();
                break;
            case michaelisMenten:
                bcTypes.setAllNeumann();
                break;
            case linear:
                bcTypes.setAllNeumann();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Bottom or inner boundary type not implemented");
            }
        } else {
            bcTypes.setAllNeumann(); // no top not bottom is no flux
        }
        return bcTypes;
    }

    /*!
     * \copydoc FVProblem::dirichletAtPos
     *
     * dirchlet(...) is called by the local assembler, e.g. BoxLocalAssembler::evalDirichletBoundaries
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const {
        PrimaryVariables values;
        if (onUpperBoundary_(globalPos)) { // top bc
            switch (bcTopType_) {
            case constantConcentration:
                values[0] = toPa_(bcTopValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Top boundary or outer type Dirichlet: unknown boundary type");
            }
        } else if (onLowerBoundary_(globalPos)) { // bot bc
            switch (bcBotType_) {
            case constantConcentration:
                values[0] = toPa_(bcBotValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Bottom or inner boundary type Dirichlet: unknown boundary type");
            }
        }
        return values;
    }

    /*!
     * \copydoc FVProblem::neumann // [kg/(m²*s)]
     *
     * called by BoxLocalResidual::evalFlux,
     *  mass flux in \f$ [ kg / (m^2 \cdot s)] \f$
     */
    NumEqVector neumann(const Element& element,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf) const {

        NumEqVector flux;
        GlobalPosition pos = scvf.center();
        if (onUpperBoundary_(pos)) { // top or outer bc
            switch (bcTopType_) {
            case constantFlux: { // with switch for maximum in- or outflow
                Scalar constflux = -bcTopValue_*rho_/(24.*60.*60.) / 100.; // cm/day -> kg/(m²*s)
                flux[0] = constflux;
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann: unknown error");
            }
        } else if (onLowerBoundary_(pos)) { // bot or inner bc
            switch (bcBotType_) {
            case constantFlux: { // with switch for maximum in- or outflow
                Scalar constflux = -bcBotValue_*rho_/(24.*60.*60.) / 100.; // cm/day -> kg/(m²*s)
                flux[0] = constflux;
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann: unknown error");
            }
        } else {
            flux[0] = 0.;
        }
        return flux;
    }

    /*!
     * Sets the current simulation time (within the simulation loop) for atmospheric look up [s]
     *
     * eventually, called in the main file (example specific, richards.cc)
     */
    void setTime(Scalar t, Scalar dt) {
        time_ = t;
        dt_ = dt; // currently unused
    }

    Scalar diffusionCoefficient = 1.; // the diffusion coefficient of the tracer

//    void setVelocity(std::shared_ptr<RichardsProblem> problem) { // TODO
//			spatialParams().setVelocityFun(problem.velocity) ??? / todoo
//    }

//    /**
//     * Callback function for spatial parameters (bad design)
//     *
//     * spatialparams has no type tag (that is good)
//     * BUT for PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>,
//     * I need GridVariables, and FluxVariables again.
//     *
//     * so that must be good enough... TODO move to Richards, otherwise types are wrong e.g. SolutionVector
//     */
//    GlobalPosition velocity(const Element &element) const {
//        GlobalPosition vel; // return value
//        auto fvGeometry = localView(this->gridGeometry());
//        fvGeometry.bindElement(element);
//        auto elemVolVars = localView(gridVars.curGridVolVars());
//        elemVolVars.bindElement(element, fvGeometry, currSol);
//    	pmVelocity->calculateVelocity(vel, elemVolVars, fvGeometry, element, 0); // o: there is only one phase, i assume pase = 0
//    	return vel;
//    }
    //    SolutionVector& currSol;
    //    GridVariables& gridVars;
    //    PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>* pmVelocity = nullptr;


private:

    //! cm pressure head -> Pascal
    Scalar toPa_(Scalar ph) const {
        return pRef_ + ph / 100. * rho_ * g_;
    }

    //! Pascal -> cm pressure head
    Scalar toHead_(Scalar p) const {
        return (p - pRef_) * 100. / rho_ / g_;
    }

    //! true if on the point lies on the upper boundary
    bool onUpperBoundary_(const GlobalPosition &globalPos) const {
        return globalPos[dimWorld - 1] > this->gridGeometry().bBoxMax()[dimWorld - 1] - eps_;
    }

    //! true if on the point lies on the upper boundary
    bool onLowerBoundary_(const GlobalPosition &globalPos) const {
        return globalPos[dimWorld - 1] < this->gridGeometry().bBoxMin()[dimWorld - 1] + eps_;
    }

    // Initial
    InputFileFunction initialSoil_; // initial concentration kg/m3

    // BC
    int bcTopType_;
    int bcBotType_;
    Scalar bcTopValue_;
    Scalar bcBotValue_;

    // additional BC params (inner)
    // Michaelis Menten maximal flux == bcBotValue_
    Scalar mmK  =1.; // Michaelis Menten half concentration

    Scalar time_ = 0.;
    Scalar dt_ = 0.;

    static constexpr Scalar eps_ = 1.e-7;
    static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
    static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)
    static constexpr Scalar pRef_ = 1.e5; // Pa

};

} // end namespace Dumux

#endif
