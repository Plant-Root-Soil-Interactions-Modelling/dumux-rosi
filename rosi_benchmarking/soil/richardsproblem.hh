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
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <RootSystem.h>

#include "richardsparams.hh"

namespace Dumux {

/*!
 * \ingroup RichardsTests
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards box model.
 */
template <class TypeTag>
class RichardsProblem;

// Specify the properties for the lens problem
namespace Properties {
NEW_TYPE_TAG(RichardsTypeTag, INHERITS_FROM(Richards));
NEW_TYPE_TAG(RichardsBoxTypeTag, INHERITS_FROM(BoxModel, RichardsTypeTag));
NEW_TYPE_TAG(RichardsCCTypeTag, INHERITS_FROM(CCTpfaModel, RichardsTypeTag));

#ifndef GRIDTYPE
// Use 3d YaspGrid per default
SET_TYPE_PROP(RichardsTypeTag, Grid, Dune::YaspGrid<3>);
#else
// Use GRIDTYPE from CMakeLists.txt
SET_TYPE_PROP(RichardsTypeTag, Grid, GRIDTYPE);
#endif

// Set the physical problem to be solved
SET_TYPE_PROP(RichardsTypeTag, Problem, RichardsProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(RichardsTypeTag, SpatialParams, RichardsParams<typename GET_PROP_TYPE(TypeTag, Grid), typename GET_PROP_TYPE(TypeTag, FVGridGeometry), typename GET_PROP_TYPE(TypeTag, Scalar)>);
} // end namespace Dumux


/*!
 * \ingroup RichardsModel
 *
 */
template <class TypeTag>
class RichardsProblem : public PorousMediumFlowProblem<TypeTag>
{
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

public:

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
     * \brief Constructor
     */
    RichardsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, GridManager<Grid>* gridManager)
    : ParentType(fvGridGeometry) {
        name_ = getParam<std::string>("Problem.Name");
        std::cout << "NAAMMEEE" << name_ << ", "<< this->name() << "\n";
        // BC
        bcTopType_ = getParam<int>("Soil.BC.Top.Type"); // todo type as a string might be nicer
        bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
        bcTopValue_ = getParam<Scalar>("Soil.BC.Top.Value",0.);
        bcBotValue_ = getParam<Scalar>("Soil.BC.Bot.Value",0.);
        // precipitation
        if (bcTopType_==atmospheric) {
            std::string filestr = name_ + ".csv";
            myfile_.open(filestr.c_str());
            precipitation_ = InputFileFunction("Climate.Precipitation", "Climate.Time"); // in [kg/(s m²)] , todo better [cm/day]?
        }
        // IC
        initialSoil_ = InputFileFunction("Soil.IC.P","Soil.IC.Z", this->spatialParams().layerIFF());
    }

    /**
     * eventually close file
     */
    ~RichardsProblem() {
        if (bcTopType_==atmospheric) {
            std::cout << "closing file \n";
            myfile_.close();
        }
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name2() const {
        return name_;
    }

    //! change problem name
    void setName(std::string str) {
        name_= str;
    }

    /*!
     * \brief Temperature [K] within a finite volume. This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const {
        return 273.15 + 10; // -> 10°C
    };

    /*!
     * \brief Reference pressure [Pa] of the non-wetting. This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonWettingReferencePressure() const {
        return 1.0e5;
    };

    /*!
     * \copydoc FVProblem::boundaryTypesAtPos
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
                values[conti0EqIdx] = -10 * bcTopValue_ / (24. * 60. * 60.); // [kg/(m²*s)] = 1/10 [cm/s] * rho
                break;
            }
            case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
                Scalar Kc = this->spatialParams().hydraulicConductivity(element);
                Scalar mS = 0;
                auto numScv = fvGeometry.numScv();
                for (auto i = 0; i < numScv; i++) {
                    mS += (elemVolVars[i].saturation() / numScv);
                }
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar krw = MaterialLaw::krw(params, mS);
                Scalar p = MaterialLaw::pc(params, mS) + nonWettingReferencePressure();
                Scalar h = -toHead_(p) / 100.; // from Pa -> m pressure head
                GlobalPosition ePos = element.geometry().center();
                Scalar dz = 2 * std::abs(ePos[dimWorld - 1] - pos[dimWorld - 1]); // 0.01; // m // fvGeometry.geometry().volume()?; TODO
                Scalar prec = precipitation_.f(time_);

                if (prec < 0) { // precipitation
                    Scalar imax = rho_ * Kc * ((h - 0.) / dz - 1.); // maximal infiltration
                    Scalar v = std::max(prec, imax);
                    values[conti0EqIdx] = v;
                } else { // evaporation
                    Scalar emax = rho_ * krw * Kc * ((h - (-100)) / dz - 1.); // maximal evaporation (-100 m = -10.000 cm) // TODO make a parameter
                    Scalar v = std::min(prec, emax);
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
                values[conti0EqIdx] = -10 * bcBotValue_ / (24. * 60. * 60.); // [kg/(m²*s)] = 1/10 [cm/s] *rho
                break;
            }
            case freeDrainage: { // TODO needs improvement
                Scalar Kc = this->spatialParams().hydraulicConductivity(
                    element);
                Scalar mS = 0; // mean saturation
                auto numScv = fvGeometry.numScv();
                for (auto i = 0; i < numScv; i++) {
                    mS += (elemVolVars[i].saturation() / numScv);
                }
                MaterialLawParams params =
                    this->spatialParams().materialLawParams(element);
                Scalar krw = MaterialLaw::krw(params, mS);
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
     * \copydoc FVProblem::initial
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

    /**
     * Sets the current simulation time (within the simulation loop) for atmospheric look up
     */
    void setTime(Scalar t) {
        time_ = t;
    }

private:

    //! cm pressure head -> Pascal
    Scalar toPa_(Scalar ph) const {
        return nonWettingReferencePressure() + ph / 100. * rho_ * g_;
    }

    //! Pascal -> cm pressure head
    Scalar toHead_(Scalar p) const {
        return (p - nonWettingReferencePressure()) * 100. / rho_ / g_;
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

    InputFileFunction initialSoil_;
    InputFileFunction precipitation_;

    std::string name_;
    Scalar time_ = 0.;

    // BC
    int bcTopType_;
    int bcBotType_;
    Scalar bcTopValue_;
    Scalar bcBotValue_;

    mutable std::ofstream myfile_;
    mutable Scalar last_time_ = -1.;

    static constexpr Scalar eps_ = 1.e-7;
    static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
    static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)

};

} //end namespace Dumux

#endif
