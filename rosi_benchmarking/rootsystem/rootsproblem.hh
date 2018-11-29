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

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/common/reorderingdofmapper.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <math.h>

#include "rootspatialparams_dgf.hh"

namespace Dumux {

template <class TypeTag>
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

    using ElementMapper = ReorderingDofMapper<GridView>;
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
    using VertexMapper = ReorderingDofMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, BoxDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Roots> {
    using type = RootsProblem<TypeTag>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Roots> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsDGF<FVGridGeometry, Scalar>;
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
template <class TypeTag>
class RootsProblem : public PorousMediumFlowProblem<TypeTag>
{
	using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

	static const int dim = GridView::dimension;
	static const int dimWorld = GridView::dimensionworld;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        // indices of the primary variables
		conti0EqIdx = Indices::conti0EqIdx,
		pressureIdx = Indices::pressureIdx
	};

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    enum {
        isBox = GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box
    };

public:

	RootsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry) {

        soilP_ = toPa_(getParam<std::vector<Scalar>>("RootSystem.Soil.P"));
		soilTable_ = soilP_.size()>1;
		if (soilTable_) {
            soilZ_ = getParam<std::vector<Scalar>>("RootSystem.Soil.Z");
		}

        p0_ = toPa_(getParam<std::vector<Scalar>>("RootSystem.Collar.P"));
		collarTable_ = p0_.size()>2;
		collarSine_ = p0_.size() == 2;
		if (collarTable_) {
            pT_ = getParam<std::vector<Scalar>>("RootSystem.Collar.PT");
		}

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
	 *
	 * \param globalPos The position of the center of the finite volume
	 */
	BoundaryTypes boundaryTypesAtPos(const GlobalPosition &pos) const
	{
	    BoundaryTypes bcTypes;
	    bcTypes.setAllNeumann(); // default
	    if (onUpperBoundary_(pos)) { // root collar
	        bcTypes.setAllDirichlet();
	    } else { // for all other (i.e. root tips)

	    	bcTypes.setAllNeumann();
	    }
	    return bcTypes;
	}

	/*!
	 * \brief Evaluate the boundary conditions for a dirichlet
	 *        control volume.
	 *
	 * \param globalPos The center of the finite volume which ought to be set.
	 *
	 * For this method, the \a values parameter stores primary variables.
	 */
	PrimaryVariables dirichletAtPos(const GlobalPosition &pos) const
	{
	    Scalar p = getCollarP_(0.); // TODO time is an issue
	    return PrimaryVariables(p);
	}

	/*
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
	{
        NumEqVector values;
	    values[conti0EqIdx] = 0.;
	    return values[conti0EqIdx];
	}

	/*!
	 * For this method, the return parameter stores the conserved quantity rate
	 * generated or annihilate per volume unit. Positive values mean
	 * that the conserved quantity is created, negative ones mean that it vanishes.
	 * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
	 */
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
	{
        NumEqVector values;
        Scalar l = element.geometry().volume(); // length of element (m)
        auto params = this->spatialParams();
        const auto& elementMapper = params.fvGridGeometry().elementMapper();
        const auto eIdx = elementMapper.index(element);
        Scalar r = params.radius(eIdx); // root radius (m)
        Scalar kr = params.Kr(eIdx); //  radial conductivity (m^2 s / kg)
        Scalar phx = elemVolVars[0].pressure(); // kg/m/s^2
		Scalar z = 0.; // TODO how to I get the (mid) z coordinate from the element?
		Scalar phs = getSoilP_(z); // kg/m/s^2
		values[conti0EqIdx] = kr * 2*r*M_PI*l * (phs - phx); // m^3/s
		values[conti0EqIdx] /= (r*r*M_PI)*l; // 1/s
		values[conti0EqIdx] *= rho_; // (kg/s/m^3)
		return values;
	}

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        const auto& elementMapper = this->spatialParams_->fvGridGeometry().elementMapper();
        const auto eIdx = elementMapper.index(element);
        Scalar r = this->spatialParams_->radius(eIdx); // root radius (m)
        return M_PI * r * r;
    }

	/*!
	 * \brief Evaluate the initial value for a control volume.
	 */
	PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const {
		return PrimaryVariables(getSoilP_(globalPos[2]));
	}

    /**
     * Sets the current simulation time (within the simulation loop) for atmospheric look up
     */
    void setTime(Scalar t) {
        time_ = t;
    }

private:

	// cm -> Pa
	Scalar toPa_(Scalar ph) const {
		return pRef_ +  ph/100.*rho_*g_;
	}

	std::vector<Scalar> toPa_(const std::vector<Scalar>& ph) const {
		std::vector<Scalar> np = ph;
		for (auto& p : np) {
			p = toPa_(p);
		}
		return np;
	}

	// on root collar
	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_;
	}

	// table look up
	size_t map_(double x, const std::vector<double>& x_) const {
		unsigned int jr,jm,jl;
		jl = 0;
		jr = x_.size();
		while (jr-jl > 1) {
			jm=(jr+jl) >> 1; // thats a divided by two
			if (x >= x_[jm])
				jl=jm;
			else
				jr=jm;
		}
		return jl; // left index
	}

	// soil table look up
	Scalar getSoilP_(const Scalar z) const {
		if (soilTable_) {
			size_t i = map_(z, soilZ_);
			return soilP_.at(i);
		} else {
			return soilP_[0];
		}
	}

	Scalar getCollarP_(const Scalar t) const {
		if (collarTable_) { // table
			size_t i = map_(t, pT_);
			return p0_.at(i);
		} else  if (collarSine_) { // daily sine
			auto d = t/24/3600; // s -> days
			d = d - floor(d);
			d *= 3.1415/2.; // day -> pi/2
			return p0_[0] + cos(d)*(p0_[1]-p0_[0]);
		} else { // const
			return p0_[0];
		}
	}
    static constexpr Scalar g_ = 9.81; // cm / s^2
    static constexpr Scalar rho_ = 1.e3; // kg / m^3
    static constexpr Scalar pRef_ = 1.e5; // Pa

	static constexpr Scalar eps_ = 1e-8;

	std::vector<Scalar> soilP_;
	std::vector<Scalar> soilZ_ = std::vector<Scalar> (0);
	bool soilTable_;

	std::vector<Scalar> p0_;
	std::vector<Scalar> pT_= std::vector<Scalar> (0);;
	bool collarTable_;
	bool collarSine_;

    Scalar time_ = 0.;

};

} //end namespace Dumux

#endif
