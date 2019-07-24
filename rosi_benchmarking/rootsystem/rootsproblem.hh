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
 * \ingroup OnePTests
 * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
#ifndef ROOTS_PROBLEM_HH
#define ROOTS_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dumux/common/reorderingdofmapper.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/growth/soillookup.hh>

#include <math.h>
#include <map>

#include "rootspatialparams_dgf.hh"
#include "rootspatialparams_rb.hh"

namespace Dumux {

template<class TypeTag>
class RootsProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct Roots {
    using InheritsFrom = std::tuple<OneP>;
};
struct RootsCCTpfa {
    using InheritsFrom = std::tuple<Roots, CCTpfaModel>;
};
struct RootsBox {
    using InheritsFrom = std::tuple<Roots, BoxModel>;
};
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_FOAMGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::Roots> {using type = Dune::FoamGrid<1, 3>;};
#endif

// if we have pt scotch use the reordering dof mapper to optimally sort the dofs (cc)
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::RootsCCTpfa> {
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>; // ReorderingDofMapper
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, CCTpfaDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

// if we have pt scotch use the reordering dof mapper to optimally sort the dofs (box)
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::RootsBox> {
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>; //ReorderingDofMapper
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, BoxDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Roots> {
    using type = RootsProblem<TypeTag>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Roots> {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
};

} // end namespace Properties


/*!
 * \ingroup RootsProblem
 * \brief A test problem for roots
 */
template<class TypeTag>
class RootsProblem: public PorousMediumFlowProblem<TypeTag> {

    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    enum {
        conti0EqIdx = Indices::conti0EqIdx, // indices of the primary variables
        pressureIdx = Indices::pressureIdx
    };
    enum {
        isBox = GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box
    };
    enum {
        bcDirichlet = 0,
        bcNeumann = 1
    };

    static const int dimWorld = GridView::dimensionworld;

public:

    RootsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry) :
        ParentType(fvGridGeometry) {
        auto sf = InputFileFunction("Soil.IC.P", "Soil.IC.Z");
        soil_ = new GrowthModule::SoilLookUpTable(sf);
        try {
            collar_ = InputFileFunction("RootSystem.Collar", "Transpiration", "TranspirationT");
            bcType_ = bcNeumann;
        } catch (...) {
            collar_ = InputFileFunction("RootSystem.Collar", "P", "RootSystem.Collar.PT");
            bcType_ = bcDirichlet;
        }
        file_at_.open(this->name() + "_actual_transpiration.txt");
    }

    virtual ~RootsProblem() {
        std::cout << "closing file \n" << std::flush;
        file_at_.close();
    }

    //! evaluates user defined data for vtk fields
    void userData(std::string name, const SolutionVector& sol) {
        const auto& gridView = this->fvGridGeometry().gridView();
        userData_[name] = std::vector<Scalar>(gridView.size(0));
        auto eMapper = this->fvGridGeometry().elementMapper();
        auto vMapper = this->fvGridGeometry().vertexMapper();
        for (const auto& e : elements(gridView)) {
            auto eIdx = eMapper.index(e);
            double d = 0;
            if (name=="kr") {
                d = this->spatialParams().kr(eIdx);
            }
            if (name=="kx") {
                d = this->spatialParams().kx(eIdx);
            }
            if (name=="age") {
                d = this->spatialParams().age(eIdx);
            }
            if (name=="order") {
                d = this->spatialParams().order(eIdx);
            }
            if (name=="radius") {
                d = this->spatialParams().radius(eIdx);
            }
            if (name=="initialPressure") {
                d = initialAtPos(e.geometry().center());
            }
            if (name=="radialFlux") {
                auto geo = e.geometry();
                auto length = geo.volume();
                auto kr = this->spatialParams().kr(eIdx);
                auto a = this->spatialParams().radius(eIdx);
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                auto p = geo.center();
                d =  2 * a * M_PI * length* kr * (soil(p) - (sol[i1] + sol[i0]) / 2); // m^3 / s
            }
            if (name=="axialFlux") {
                auto geo = e.geometry();
                auto length = geo.volume();
                auto kx = this->spatialParams().kx(eIdx);
                auto i0 = vMapper.subIndex(e, 0, 1);
                auto i1 = vMapper.subIndex(e, 1, 1);
                d = kx * ((sol[i1] - sol[i0]) / length - rho_ * g_); // m^3 / s
            }
            userData_[name][eIdx] = d;
        }
    }

    //! vtk fields call back functions
    std::vector<Scalar>& radialFlux() { return userData_["radialFlux"]; } // [m3/s]
    std::vector<Scalar>& axialFlux() { return userData_["axialFlux"]; } // [m3/s]
    std::vector<Scalar>& kr() { return userData_["kr"]; } // [m/Pa/s]
    std::vector<Scalar>& kx() { return userData_["kx"]; } // [m4/Pa/s]
    std::vector<Scalar>& age() { return userData_["age"]; } // [s]
    std::vector<Scalar>& order() { return userData_["order"]; } // [1]
    std::vector<Scalar>& radius() { return userData_["radius"]; } // [m]
    std::vector<Scalar>& initialPressure() { return userData_["initialPressure"]; } // [Pa]

    //! calculates transpiraton, as the sum of radial fluxes (slow but accurate)
    Scalar transpiration(const SolutionVector& sol) {
        userData("radialFlux", sol);
        return std::accumulate(userData_["radialFlux"].begin(), userData_["radialFlux"].end(), 0.);
    }

    /*
     * \brief Return the temperature within the domain in [K]. (actually needed? why?)
     */
    Scalar temperature() const {
        return 273.15 + 10; // 10C
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &pos) const {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann(); // default
        if (onUpperBoundary_(pos)) { // root collar
            if (bcType_ == bcDirichlet) {
                bcTypes.setAllDirichlet();
            } else {
                if (!critical_) {
                    bcTypes.setAllNeumann();
                } else {
                    bcTypes.setAllDirichlet();
                }
            }
        } else { // for all other (i.e. root tips)
            bcTypes.setAllNeumann();
        }
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &pos) const {
        if (critical_) {
            return criticalCollarPressure_;
        } else {
            return PrimaryVariables(collar());
        }
    }

    /*
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
        const SubControlVolumeFace& scvf) const {
        const auto globalPos = scvf.center();
        if (onUpperBoundary_(globalPos)) {
            auto& volVars = elemVolVars[scvf.insideScvIdx()];
            Scalar p = volVars.pressure();
            auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            Scalar kx = this->spatialParams().kx(eIdx);
            auto dist = (globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
            Scalar maxTrans = volVars.density(0) * kx * (p - criticalCollarPressure_) / dist;  // todo!!!!
            Scalar trans = collar(); // kg/s
            Scalar v = std::min(trans, maxTrans);
            lastActualTrans_ = v; // the one we return
            lastTrans_ = trans;  // potential transpiration
            lastMaxTrans_ = maxTrans; // maximal transpiration at this saturation
            lastP_ = p;
            v /= volVars.extrusionFactor(); // convert from kg/s to kg/(s*m^2)
            return NumEqVector(v);
        } else {
            return NumEqVector(0.); // no flux at root tips
        }
    }

    /*!
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    NumEqVector source(const Element &element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
        const SubControlVolume &scv) const {

        NumEqVector values;
        auto params = this->spatialParams();
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        Scalar a = params.radius(eIdx); // root radius (m)
        Scalar kr = params.kr(eIdx); //  radial conductivity (m^2 s / kg)
        Scalar phx = elemVolVars[scv.localDofIndex()].pressure(); // kg/m/s^2
        Scalar phs = soil(scv.center()); // kg/m/s^2
        values[conti0EqIdx] = kr * 2 * a * M_PI * (phs - phx); // m^3/s
        values[conti0EqIdx] /= (a * a * M_PI); // 1/s
        values[conti0EqIdx] *= rho_; // (kg/s/m^3)
        return values;

    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& p) const {
        return PrimaryVariables(soil(p)); // soil(p)
    }

    void setSoil(CRootBox::SoilLookUp* s) {
        std::cout << "setSoil(...): manually changed soil to " << s->toString() << "\n";
        soil_ = s;
    }

    //! soil pressure (called by initial, and source term)
    Scalar soil(const GlobalPosition& p) const {
        auto p2 = CRootBox::Vector3d(p[0] * 100, p[1] * 100, p[2] * 100);
        double d = soil_->getValue(p2);
        // std::cout << "rootsproblem::soil() " << p2.toString() << ", " << d << "\n";
        return toPa_(d);
    }

    //! sets the current simulation time [s] (within the simulation loop) for collar boundary look up
    void setTime(Scalar t) {
        this->spatialParams().setTime(t);
        time_ = t;
    }

    /*!
     * writes the actual transpiration into a text file:
     * 0 time, 1 actual transpiration, 2 potential transpiration, 3 maximal transpiration, 4 collar pressure, 5 calculated actual transpiration,
     * 1 - 4 work only for neuman bc
     */
    void writeTranspirationRate(const SolutionVector& sol) {
        Scalar trans = this->transpiration(sol);
        Scalar p = lastP_;
        file_at_ << time_ << ", " << lastActualTrans_ << ", " << lastTrans_ << ", " << lastMaxTrans_ << ", " << p << ", " << trans << "\n"; // << std::setprecision(17)
    }

    //! pressure or transpiration rate at the root collar (called by dirichletor neumann, respectively)
    Scalar collar() const {
        Scalar t = time_/24./3600; // s -> day
        if (bcType_ == bcDirichlet) {
            return toPa_(collar_.f(t)); // Pa
        } else {
            return collar_.f(t)/24./3600.; // kg/day -> kg/s
        }
    }

    //! if true, sets bc to dirichlet at criticalCollarPressure (false per default)
    void setCritical(bool b) {
        critical_ = b;
    }

    //! sets the criticalCollarPressure [Pa]
    void criticalCollarPressure(Scalar p) {
        criticalCollarPressure_ = p;
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here makes extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     *
     * called by volumevariables (why there?), no compilation error if you remove it, just wrong results
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
        const SubControlVolume &scv,
        const ElementSolution& elemSol) const
    {
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        const auto radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
    }


private:

    CRootBox::SoilLookUp* soil_;
    InputFileFunction collar_;
    size_t bcType_;
    Scalar time_ = 0.;
    Scalar criticalCollarPressure_ = -1.4e6;

    bool critical_ = false;

    Scalar toPa_(Scalar ph) const {     // cm -> Pa
        return pRef_ + ph / 100. * rho_ * g_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const {  // on root collar
        return globalPos[dimWorld - 1] > this->fvGridGeometry().bBoxMax()[dimWorld - 1] - eps_;
    }

    static constexpr Scalar g_ = 9.81; // cm / s^2
    static constexpr Scalar rho_ = 1.e3; // kg / m^3
    static constexpr Scalar pRef_ = 1.e5; // Pa
    static constexpr Scalar eps_ = 1e-8;

    std::ofstream file_at_; // file for actual transpiration
    mutable Scalar lastActualTrans_ = 0;
    mutable Scalar lastTrans_ = 0.;
    mutable Scalar lastMaxTrans_ = 0.;
    mutable Scalar lastP_ = 0.;

    std::map<std::string, std::vector<Scalar>> userData_;

};

} //end namespace Dumux

#endif
