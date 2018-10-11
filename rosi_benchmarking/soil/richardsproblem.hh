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

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

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
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, SpatialParams)::MaterialLaw;
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
    : ParentType(fvGridGeometry)
    {
        gridManager_ = gridManager;
        this->spatialParams().setGridManager(gridManager); // SpatialParams constructor is hidden in parent class, and can not be modified

        name_ = getParam<std::string>("Problem.Name");
        // BC
        bcTopType_ = getParam<int>("Soil.BC.Top.Type"); // todo type as a string might be nicer
        bcBotType_ = getParam<int>("Soil.BC.Bot.Type");
        bcTopValue_ = getParam<Scalar>("Soil.BC.Top.Value",0.);
        bcBotValue_ = getParam<Scalar>("Soil.BC.Bot.Value",0.);
        if (bcTopType_==atmospheric) {
            std::string filestr = name_ + ".csv";
            myfile_.open(filestr.c_str());
            precData_ = getParam<std::vector<Scalar>>("Climate.Precipitation"); // in [kg/(s m²)] , todo better [cm/day]?
            precTime_ = getParam<std::vector<Scalar>>("Climate.Time");
        } else {
            precData_ = std::vector<Scalar>(0); // not used
            precTime_ = std::vector<Scalar>(0);
        }
        // IC
        initialPressure_ = getParam<std::vector<Scalar>>("Soil.IC.Pressure", std::vector<Scalar>{-100.});
        if (initialPressure_.size() == 1) { // constant pressure
            constInitial_ = true;
            gridInitial_ = false;
            initialZ_ = std::vector<Scalar>(0); // not used
        } else {
            constInitial_ = false;
            try {
                initialZ_ = getParam<std::vector<Scalar>>("Soil.IC.Z"); // table look up
            } catch (std::exception& e) { // grid parameter look up
                initialZ_ = std::vector<Scalar>(0); // not used
                gridInitial_ = true;
            }
        }
        //        std::cout << "Grid boundarys " << this->fvGridGeometry().bBoxMax()
        //            << "\n";
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
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief Temperature [K] within a finite volume. This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 273.15 + 10; // -> 10°C
    };

    /*!
     * \brief Reference pressure [Pa] of the non-wetting. This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonWettingReferencePressure() const
    {
        return 1.0e5;
    };

    /*!readGridData
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
        } else { // bot bc
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
        }
        return bcTypes;
    }

    /*!
     * \copydoc FVProblem::dirichletAtPos
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        if (onUpperBoundary_(globalPos)) { // top bc
            switch (bcTopType_) {
            case constantPressure:
                values[Indices::pressureIdx] =toPa_(bcTopValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Top boundary type Dirichlet: unknown boundary type");
            }
        } else { // bot bc
            switch (bcBotType_) {
            case constantPressure:
                values[Indices::pressureIdx] = toPa_(bcBotValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Dirichlet: unknown boundary type");
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
        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values;
        GlobalPosition pos = scvf.center();

        if (onUpperBoundary_(pos)) { // top bc
            switch (bcTopType_) {
            case constantFlux: {
                values[conti0EqIdx] = -10*bcTopValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] * rho
                break;
            }
            case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
                Scalar Kc = this->spatialParams().hydraulicConductivity(element);
                Scalar mS = 0;
                auto numScv = fvGeometry.numScv();
                for (auto i = 0; i<numScv; i++) {
                    mS += (elemVolVars[i].saturation()/numScv);
                }
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar krw = MaterialLaw::krw(params, mS);
                Scalar p = MaterialLaw::pc(params, mS)+nonWettingReferencePressure();
                Scalar h = -toHead_(p)/100.; // from Pa -> m pressure head

                GlobalPosition ePos = element.geometry().center();
                Scalar dz = 2 * std::abs(ePos[dimWorld-1]- pos[dimWorld-1]); // 0.01; // m // fvGeometry.geometry().volume()?; TODO
                Scalar prec = getPrec_(time_); // precipitation or evaporation TODO

                if (prec<0) { // precipitation
                    Scalar imax = rho_*Kc*((h-0.)/dz -1.); // maximal infiltration
                    Scalar v = std::max(prec,imax);
                    values[conti0EqIdx] = v;
                } else { // evaporation
                    Scalar emax = rho_*krw*Kc*((h-(-100))/dz -1.); // maximal evaporation (-100 m = -10.000 cm) // TODO make a parameter
                    Scalar v  = std::min(prec,emax);
                    values[conti0EqIdx] = v;
                }
                // hack for benchmark 4 TODO some better concept for output
                myfile_ << time_ << ", "; //
                myfile_ << values[conti0EqIdx] << "\n";
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException,"Top boundary type Neumann: unknown error");
            }
        } else  { // bot bc
            switch (bcBotType_) {
            case constantFlux: {
                values[conti0EqIdx] = -10*bcBotValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] *rho
                break;
            }
            case freeDrainage: { // TODO needs improvement
                Scalar Kc = this->spatialParams().hydraulicConductivity(element);
                Scalar mS = 0; // mean saturation
                auto numScv = fvGeometry.numScv();
                for (auto i = 0; i<numScv; i++) {
                    mS += (elemVolVars[i].saturation()/numScv);
                }
                MaterialLawParams params = this->spatialParams().materialLawParams(element);
                Scalar krw = MaterialLaw::krw(params, mS);
                values[conti0EqIdx] = krw*Kc*rho_; // * 1 [m]
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Neumann: unknown error");
            }
        }
        return values;
    }

    /*!
     * \copydoc FVProblem::initial
     */
    template<class Entity>
    PrimaryVariables initial(const Entity& entity) const
    {
        PrimaryVariables values;
        if (constInitial_) { // constant initial data
            values[pressureIdx] = toPa_(initialPressure_[0]);
        } else {
            if (gridInitial_) { // obtain layer number from grid data
                size_t i = size_t(
                    gridManager_->getGridData()->parameters(entity).at(
                        materialLayerNumber));
                values[pressureIdx] = toPa_(initialPressure_.at(i));
            } else { // obtain pressure by table look up and linear interpolation
                Scalar z = entity.geometry().center()[dimWorld-1];
                values[pressureIdx] = toPa_(this->spatialParams().interp1(z,initialPressure_, initialZ_));
            }
        }
        values.setState(bothPhases);
        return values;
    }

    /**
     * Sets the current simulation time (within the simulation loop) for atmospheric look up
     */
    void setTime(Scalar t)
    {
        time_ = t;
    }

private:

    //! cm pressure head -> Pascal
    Scalar toPa_(Scalar ph) const
    {
        return nonWettingReferencePressure() +  ph/100.*rho_*g_;
    }

    //! Pascal -> cm pressure head
    Scalar toHead_(Scalar p) const
    {
        return (p-nonWettingReferencePressure()) * 100./rho_/g_;
    }

    //! true if on the point lies on the upper boundary
    bool onUpperBoundary_(const GlobalPosition &globalPos) const {
        return globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_;
    }

    //! returns dynamic precipitation from table from input time
    Scalar getPrec_(Scalar t) const {
        return this->spatialParams().interp1(t, precData_, precTime_);
    }

    std::string name_;
    mutable std::ofstream myfile_;

    GridManager<Grid>* gridManager_;

    Scalar time_ = 0.;

    // IC
    bool constInitial_;
    bool gridInitial_;
    std::vector<Scalar> initialPressure_;
    std::vector<Scalar> initialZ_;

    // BC
    int bcTopType_;
    int bcBotType_;
    Scalar bcTopValue_;
    Scalar bcBotValue_;
    std::vector<Scalar> precTime_; // climatic data for atmospheric bc
    std::vector<Scalar> precData_; // climatic data for atmospheric bc

    static constexpr Scalar eps_ = 1.e-7;
    static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)
    static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)

};

} //end namespace Dumux

#endif
