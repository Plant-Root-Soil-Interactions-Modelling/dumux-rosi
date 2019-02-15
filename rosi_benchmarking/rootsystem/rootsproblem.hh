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

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        conti0EqIdx = Indices::conti0EqIdx, // indices of the primary variables
        pressureIdx = Indices::pressureIdx
    };

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
        isBox = GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box
    };

    enum {
        bcDirichlet = 0, bcNeumann = 1
    };

public:

    RootsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry) :
        ParentType(fvGridGeometry) {
        auto sf = InputFileFunction("Soil.IC.P", "Soil.IC.Z");
        soil_ = new GrowthModule::SoilLookUpTable<FVGridGeometry>(sf, fvGridGeometry);
        try {
            collar_ = InputFileFunction("RootSystem.Collar.Transpiration", "RootSystem.Collar.TranspirationT", true);
            bcType_ = bcNeumann;
        } catch (...) {
            collar_ = InputFileFunction("RootSystem.Collar.P", "RootSystem.Collar.PT", true);
            bcType_ = bcDirichlet;
        }
        file_at_.open(this->name() + "_actual_transpiration.txt");
    }

    ~RootsProblem() {
        std::cout << "closing file \n";
        file_at_.close();
    }

    //! calculates axial fluxes from a given solution (for vtk output) [m^3 / s]
    void axialFlux(const SolutionVector& sol) {
        const auto& gridView = this->fvGridGeometry().gridView();
        axialFlux_ = std::vector<Scalar>(gridView.size(0));
        auto eMapper = this->fvGridGeometry().elementMapper();
        auto vMapper = this->fvGridGeometry().vertexMapper();
        for (const auto& e : elements(gridView)) {
            const auto eIdx = eMapper.index(e);
            auto geo = e.geometry();
            auto length = geo.volume();
            auto kx = this->spatialParams().kx(eIdx);
            auto i0 = vMapper.subIndex(e, 0, 1);
            auto i1 = vMapper.subIndex(e, 1, 1);
            axialFlux_[eIdx] = rho_ * kx * ((sol[i1] - sol[i0]) / length - rho_ * g_); // m^3 / s
        }
    }

    //! calculates the radial fluxes from a given solution (for vtk output) [m^3 / s]
    void radialFlux(const SolutionVector& sol) {
        const auto& gridView = this->fvGridGeometry().gridView();
        radialFlux_ = std::vector<Scalar>(gridView.size(0));
        auto eMapper = this->fvGridGeometry().elementMapper();
        auto vMapper = this->fvGridGeometry().vertexMapper();
        for (const auto& e : elements(gridView)) {
            auto eIdx = eMapper.index(e);
            auto kr = this->spatialParams().kr(eIdx);
            auto i0 = vMapper.subIndex(e, 0, 1);
            auto i1 = vMapper.subIndex(e, 1, 1);
            auto p = e.geometry().center();
            radialFlux_[eIdx] = kr * (soil(p) - (sol[i1] + sol[i0]) / 2); // m^3 / s
        }
    }

    std::vector<Scalar>& radialFlux() {
        return radialFlux_;
    }
    std::vector<Scalar>& axialFlux() {
        return axialFlux_;
    }

    //! calculates transpiraton, as the netflux of first element (m^3 /s), assuming first element is collar
    Scalar transpiration(const SolutionVector& sol) {
//        auto vMapper = this->fvGridGeometry().vertexMapper();
//        auto eIdx = 0; // collar index
//        const auto& e = this->fvGridGeometry().element(eIdx);
//        auto i0 = vMapper.subIndex(e, 0, 1);
//        auto i1 = vMapper.subIndex(e, 1, 1);
//        // std::cout << "vertex index" << i0 << ", " << i1 << "\n";
//        auto length = e.geometry().volume();
//        auto kx = this->spatialParams().kx(eIdx);
//        return rho_ * kx * ((sol[i1] - sol[i0]) / length); // kg / s ; + rho_ * g_
        radialFlux(sol);
        return std::accumulate(radialFlux_.begin(), radialFlux_.end(), 0.);
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
//            auto p = volVars.pressure(0);
//            auto eIdx = this->fvGridGeometry().elementMapper().index(element);
//            Scalar kx = this->spatialParams().kx(eIdx);
//            auto dist = (globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
//            Scalar maxTrans = volVars.density(0) * kx * (p - criticalCollarPressure_) / (2*dist); // / volVars.viscosity(0)
//            Scalar trans = collar();
////            std::cout << trans << " kg/s, " << maxTrans << " kg/s, " << p << " Pa " << ", diff " << (p - criticalCollarPressure_) << " scale "
////                << volVars.density(0) * kx / (2 * dist) << " crit " << criticalCollarPressure_ << ", p= " << (p - pRef_) * 100 / rho_ / g_ << "\n";
//            Scalar v = std::min(trans, maxTrans);
//            lastActualTrans_ = v; // the one we return
//            lastTrans_ = trans;  // potential transpiration
//            lastMaxTrans_ = maxTrans; // maximal transpiration at this saturation
//            lastP_ = p;
//            v /= volVars.extrusionFactor(); // convert from kg/s to kg/(s*m^2)
//            // std::cout << volVars.extrusionFactor() << " cm2, " << scvf.area() << " cm2 "; // scvf.area() == 1
//            return NumEqVector(v);
            Scalar trans = collar();
            lastTrans_ = trans;
            return NumEqVector(trans/volVars.extrusionFactor());
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
        const auto eIdx = params.fvGridGeometry().elementMapper().index(element);
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
     * \brief Return how much the domain is extruded at a given sub-control volume.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element, const SubControlVolume& scv, const ElementSolution& elemSol) const {
        const auto eIdx = this->spatialParams_->fvGridGeometry().elementMapper().index(element);
        Scalar r = this->spatialParams_->radius(eIdx); // root radius (m)
        return M_PI * r * r;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& p) const {
        // std::cout << "initial pos " << p[2] << ": " << soil(p) << "\n";
        return PrimaryVariables(soil(p)); // soil(p)
    }

    void setSoil(CRootBox::SoilLookUp* s) {
        std::cout << "manually changed soil to " << s->toString() << "\n";
        soil_ = s;
    }

    //! soil pressure (called by initial, and source term)
    Scalar soil(const GlobalPosition& p) const {
        return toPa_(soil_->getValue(CRootBox::Vector3d(p[0] * 100, p[1] * 100, p[2] * 100)));
    }

    //! sets the current simulation time (within the simulation loop) for collar boundary look up
    void setTime(Scalar t) {
        this->spatialParams().setTime(t);
        time_ = t;
    }

    //! writes the transpiration file
    void writeTranspirationRate(const SolutionVector& sol) {
        Scalar trans = this->transpiration(sol);
        Scalar p = lastP_;
        Scalar dp = lastP_ - criticalCollarPressure_;
        file_at_ << std::setprecision(17) << time_ << ", " << lastActualTrans_ << ", " << lastTrans_ << ", " << lastMaxTrans_ << ", " << p << ", " << dp << ", "
            << std::setprecision(17) << sol[0] << ", " << std::setprecision(17) << sol[1] << ", " << trans << "\n";
    }

    //! pressure or transpiration rate at the root collar (called by dirichletor neumann, respectively)
    Scalar collar() const {
        if (bcType_ == bcDirichlet) {
            return toPa_(collar_.f(time_)); // Pa
        } else {
            return collar_.f(time_); // kg/s (?)
        }
    }

    void setCritical(bool b) {
        critical_ = b;
    }

private:

    CRootBox::SoilLookUp* soil_;
    InputFileFunction collar_;
    int bcType_;
    Scalar time_ = 0.;
    Scalar criticalCollarPressure_ = toPa_(-1e4);

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

    // vtk fields
    std::vector<Scalar> axialFlux_;
    std::vector<Scalar> radialFlux_;

    //    void initialPressure_() {
    //        const auto& gridView = this->fvGridGeometry().gridView();
    //        initialP_ = std::vector<Scalar>(gridView.size(0));
    //        auto eMapper = this->fvGridGeometry().elementMapper();
    //        for (const auto& element : elements(gridView)) {
    //            auto eIdx = eMapper.index(element);
    //            initialP_[eIdx] = initialAtPos(element.geometry().center());
    //        }
    //    }
    // std::vector<Scalar> initialP_;

};

} //end namespace Dumux

#endif
