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
 * \ingroup Richards Equation Solver
 * \brief Uses Dumux as an easy to use Richards equation solver, where most parameters can be set dynamically.
 */
#ifndef DUMUX_RICHARDS_PROBLEM_HH
#define DUMUX_RICHARDS_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/problem.hh> // base class

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <RootSystem.h>

#include "richardsparams.hh"

namespace Dumux {

template <class TypeTag>
class RichardsProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct RichardsTT { using InheritsFrom = std::tuple<Richards>; };
struct RichardsBox { using InheritsFrom = std::tuple<RichardsTT, BoxModel>; };
struct RichardsCC { using InheritsFrom = std::tuple<RichardsTT, CCTpfaModel>; };
} // end namespace TTag

// Set grid type
#ifndef GRIDTYPE
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsTT> { using type = Dune::YaspGrid<3>; }; // Use 3d YaspGrid per default
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsTT> { using type = GRIDTYPE; };  // Use GRIDTYPE from CMakeLists.txt
#endif

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsTT> { using type = RichardsProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsTT> {
    using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>;
};
} // end namespace properties

/*!
 *
 * \ingroup RichardsModel
 *
 */
template <class TypeTag>
class RichardsProblem : public PorousMediumFlowProblem<TypeTag>
{
public:

    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using MaterialLaw = typename GetPropType<TypeTag, Properties::SpatialParams>::MaterialLaw;
    using MaterialLawParams = typename MaterialLaw::Params;

    enum {
        // copy some indices for convenience
        pressureIdx = Indices::pressureIdx,
        conti0EqIdx = Indices::conti0EqIdx,
        bothPhases = Indices::bothPhases,
        // world dimension
        dimWorld = GridView::dimensionworld
    };

    enum BCTypes {
        constantPressure = 1,
        constantFlux = 2,
        atmospheric = 4,
        freeDrainage = 5
    };

    enum GridParameterIndex {
        materialLayerNumber = 0
    };

    /*!
     * \brief Constructor: constructed in the main file
     */
    RichardsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {
        // BC
        bcTopType_ = getParam<int>("Soil.BC.Top.Type"); // todo type as a string might be nicer
        bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
        bcTopValue_ = getParam<Scalar>("Soil.BC.Top.Value",0.);
        bcBotValue_ = getParam<Scalar>("Soil.BC.Bot.Value",0.);
        // precipitation
        if (bcTopType_==atmospheric) {
            precipitation_ = InputFileFunction("Climate", "Precipitation", "Time"); // cm/day (day)
            precipitation_.setVariableScale(1./(24.*60.*60.)); // s -> day
            precipitation_.setFunctionScale(1.e3/(24.*60.*60.)/100); // cm/day -> kg/(m²*s)
            std::string filestr = this->name() + ".csv"; // output file
            myfile_.open(filestr.c_str());
        }
        // IC
        initialSoil_ = InputFileFunction("Soil.IC", "P", "Z", this->spatialParams().layerIFF()); // [cm]([m]) pressure head, conversions hard coded
    }

    /**
     * \brief Eventually, closes output file
     */
    ~RichardsProblem() {
        if (bcTopType_==atmospheric) {
            std::cout << "closing file \n";
            myfile_.close();
        }
    }

    /*!
     * \brief Temperature [K] within a finite volume. This problem assumes a temperature of 10 degrees Celsius.
     *
     * called EnergyVolumeVariablesImplementation::updateTemperature(...) in porousmediumflow/nonisothermal/volumevariables.hh,
     * included by porousmediumflow/volumevariables.hh,
     *
     * todo this makes very little sense for isothermal !
     *
     * overwrites PorousMediumFlowProblem::temperature (compiles without, throws exception of base class)
     */
    Scalar temperature() const {
        // std::cout << "\n\n dont dont dont \n\n";
        return 273.15 + 10; // -> 10°C
    }

    /*!
     * \brief Reference pressure [Pa] of the non-wetting. This problem assumes a constant reference pressure of 1 bar.
     *
     * called by porousmediumflow/richards/volumevariables.hh
     */
    Scalar nonWettingReferencePressure() const {
        return pRef_;
    }

    /*!
     * \copydoc FVProblem::boundaryTypesAtPos
     *
     * discretization dependent, e.g. called by BoxElementBoundaryTypes::boundaryTypes(...)
     * when?
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (onUpperBoundary_(globalPos)) { // top bc
            switch (bcTopType_) {
            case constantPressure:
                bcTypes.setAllDirichlet();
                break;
            case constantFlux:
                bcTypes.setAllNeumann();
                break;
            case atmospheric:
                bcTypes.setAllNeumann();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Top boundary type not implemented");
            }
        } else if (onLowerBoundary_(globalPos)) { // bot bc
            switch (bcBotType_) {
            case constantPressure:
                bcTypes.setAllDirichlet();
                break;
            case constantFlux:
                bcTypes.setAllNeumann();
                break;
            case freeDrainage:
                bcTypes.setAllNeumann();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type not implemented");
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
            case constantPressure:
                values[Indices::pressureIdx] = toPa_(bcTopValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,
                    "Top boundary type Dirichlet: unknown boundary type");
            }
        } else if (onLowerBoundary_(globalPos)) { // bot bc
            switch (bcBotType_) {
            case constantPressure:
                values[Indices::pressureIdx] = toPa_(bcBotValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,
                    "Bottom boundary type Dirichlet: unknown boundary type");
            }
        }
        values.setState(Indices::bothPhases);
        return values;
    }

    /*!
     * \copydoc FVProblem::neumann // [kg/(m²*s)]
     *
     * called by BoxLocalResidual::evalFlux
     */
    NumEqVector neumann(const Element& element,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf) const {

        NumEqVector values;
        GlobalPosition pos = scvf.center();
        if (onUpperBoundary_(pos)) { // top bc
            switch (bcTopType_) {
            case constantFlux: {
                values[conti0EqIdx] = -bcTopValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                break;
            }
            case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
                Scalar s = elemVolVars[scvf.insideScvIdx()].saturation(0);
                Scalar Kc = this->spatialParams().hydraulicConductivity(element); //  [m/s]
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar p = MaterialLaw::pc(params, s) + pRef_;
                Scalar h = -toHead_(p); // todo why minus -pc?
                GlobalPosition ePos = element.geometry().center();
                Scalar dz = 100 * 2 * std::abs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // cm
                Scalar prec = -precipitation_.f(time_);
                if (prec < 0) { // precipitation
                    Scalar imax = rho_ * Kc * ((h - 0.) / dz - 1.); // maximal infiltration
                    Scalar v = std::max(prec, imax);
                    values[conti0EqIdx] = v;
                } else { // evaporation
                    Scalar krw = MaterialLaw::krw(params, s);
                    Scalar emax = rho_ * krw * Kc * ((h - criticalPressure_) / dz - 1.); // maximal evaporation
                    Scalar v = std::min(prec, emax);
                    // std::cout << prec << ", " << emax << ", " << h << "\n";
                    values[conti0EqIdx] = v;
                }
                // hack for benchmark 4 TODO some better concept for output
                if (time_ > last_time_) { // once per time step
                    myfile_ << time_ << ", "; //
                    myfile_ << values[conti0EqIdx] << "\n";
                    last_time_ = time_;
                }
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException,
                    "Top boundary type Neumann: unknown error");
            }
        } else if (onLowerBoundary_(pos)) { // bot bc
            switch (bcBotType_) {
            case constantFlux: {
                values[conti0EqIdx] = -bcBotValue_*rho_/(24.*60.*60.)/100; // cm/day -> kg/(m²*s)
                break;
            }
            case freeDrainage: {
                Scalar Kc = this->spatialParams().hydraulicConductivity(element);
                Scalar s = elemVolVars[scvf.insideScvIdx()].saturation(0);
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar krw = MaterialLaw::krw(params, s);
                values[conti0EqIdx] = krw * Kc * rho_; // * 1 [m]
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException,
                    "Bottom boundary type Neumann: unknown error");
            }
        } else {
            values[conti0EqIdx] = 0.;
        }
        return values;
    }

    /*!
     * \copydoc FVProblem::source
     *
     * called by FVLocalResidual:computeSource(...)
     */
    NumEqVector source(const Element &element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
        const SubControlVolume &scv) const {
        const auto eIdx = this->spatialParams().fvGridGeometry().elementMapper().index(element);
        if (source_ != nullptr) {
            return source_->at(eIdx);
        } else {
            return 0.;
        }
    }

    /*!
     * \copydoc FVProblem::initial
     *
     * called by FVProblem::applyInitialSolution(...)
     */
    template<class Entity>
    PrimaryVariables initial(const Entity& entity) const {
        auto eIdx = this->fvGridGeometry().elementMapper().index(entity);
        Scalar z = entity.geometry().center()[dimWorld - 1];
        PrimaryVariables v(0.0);
        v[pressureIdx] = toPa_(initialSoil_.f(z,eIdx));
        v.setState(bothPhases);
        return v;
    }

    /*!
     * Sets the current simulation time (within the simulation loop) for atmospheric look up [s]
     *
     * eventually, called in the main file (example specific, richards.cc)
     */
    void setTime(Scalar t) {
        time_ = t;
    }

    /*!
     * Source per element index \f$ [ kg / (m^3 \cdot s)] \f$
     *
     * eventually, called in the main file (example specific, richards.cc)
     */
    void setSource(std::vector<double>* s) {
        source_ = s;
    }

    /*!
     * sets the critical pressure for evaporation [cm] (default = -10000 cm)
     *
     *  eventually, called in the main file (example specific, richards.cc)
     */
    void criticalPressure(Scalar p) {
        criticalPressure_ = p;
    }

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
        return globalPos[dimWorld - 1]
                         > this->fvGridGeometry().bBoxMax()[dimWorld - 1] - eps_;
    }

    //! true if on the point lies on the upper boundary
    bool onLowerBoundary_(const GlobalPosition &globalPos) const {
        return globalPos[dimWorld - 1]
                         < this->fvGridGeometry().bBoxMin()[dimWorld - 1] + eps_;
    }

    // Initial
    InputFileFunction initialSoil_;

    // BC
    int bcTopType_;
    int bcBotType_;
    Scalar bcTopValue_;
    Scalar bcBotValue_;

    // Source
    std::vector<double>* source_ = nullptr;

    InputFileFunction precipitation_;
    Scalar criticalPressure_ = -1.e4; // cm
    Scalar time_ = 0.;

    mutable std::ofstream myfile_;
    mutable Scalar last_time_ = -1.;

    static constexpr Scalar eps_ = 1.e-7;
    static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
    static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)
    static constexpr Scalar pRef_ = 1.e5; // Pa

};

} //end namespace Dumux

#endif
