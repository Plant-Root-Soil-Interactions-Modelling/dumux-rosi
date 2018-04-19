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

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "rootsparams.hh"

namespace Dumux {

template <class TypeTag>
class RootsProblem;

namespace Properties {

NEW_TYPE_TAG(RootsTypeTag, INHERITS_FROM(OneP));
NEW_TYPE_TAG(RootsCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, RootsTypeTag));
NEW_TYPE_TAG(RootsBoxTypeTag, INHERITS_FROM(BoxModel, RootsTypeTag));

// Set the grid type
SET_TYPE_PROP(RootsTypeTag, Grid, Dune::FoamGrid<1, 3>);

// Set the problem property
SET_TYPE_PROP(RootsTypeTag, Problem, RootsProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(RootsTypeTag, SpatialParams, RootsParams<TypeTag>);

// the fluid system
SET_PROP(RootsTypeTag, FluidSystem)
{
	using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
	using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
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
	using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
	using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;

	static const int dim = GridView::dimension;
	static const int dimWorld = GridView::dimensionworld;

    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
	enum { // indices of the primary variables
		conti0EqIdx = Indices::conti0EqIdx,
		pressureIdx = Indices::pressureIdx
	};

	using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
	using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
	using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Element = typename GridView::template Codim<0>::Entity;
	using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
	using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
	using SubControlVolume = typename FVElementGeometry::SubControlVolume;
	using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
	using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
	using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:

	RootsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry) {

		scenario_ = getParam<Scalar>("Parameter.Scenario");
		soilP_ = toPa_(getParam<Scalar>("Parameter.SoilP"));
		if (scenario_==1) { // dirichlet top and at tips
			p0_ = toPa_(getParam<Scalar>("Parameter.P0"));
			pL_ = toPa_(getParam<Scalar>("Parameter.PL"));
		}
		if (scenario_==2) { // dirichlet top, neumann tips
			p0_ = toPa_(getParam<Scalar>("Parameter.P0"));
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
	        if (scenario_==1) {
	            bcTypes.setAllDirichlet();
	        } else {
	            bcTypes.setAllNeumann();
	        }
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
	    Scalar p = 0.;
	    if (onUpperBoundary_(pos)) { // root collar
	        p = p0_;
	    } else { // for all other (i.e. root tips)
	        p = pL_;
	    }
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
	    const RootsParams<TypeTag>& params = this->spatialParams();
		Scalar r = params.radius(SubControlVolume()); // root radius (m)
	 	Scalar kz = params.axialConductivity(SubControlVolume()); // (m^5 s / kg) == ( m^4 / (s Pa) )
	    ResidualVector values;
	    values[conti0EqIdx] = rho_*g_*kz; // m^3 / s
	    values[conti0EqIdx] /= (r*r*M_PI);
	    values[conti0EqIdx] *= rho_;
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
		ResidualVector values;
		const RootsParams<TypeTag>& params = this->spatialParams();
		Scalar l = element.geometry().volume(); // length of element (m)
		Scalar r = params.radius(scv); // root radius (m)
		Scalar kr = params.radialConductivity(scv); //  radial conductivity (m^2 s / kg)
		Scalar phx = elemVolVars[0].pressure(); // kg/m/s^2
		Scalar phs = soilP_; // kg/m/s^2
		values[conti0EqIdx] = kr * 2*r*M_PI*l * (phs - phx); // m^3/s
		values[conti0EqIdx] /= (r*r*M_PI)*l; // 1/s <-- volume?
		values[conti0EqIdx] *= rho_; // (kg/s/m^3)
		return values;
	}

	/*!
	 * \brief Evaluate the initial value for a control volume.
	 */
	PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const {
		return PrimaryVariables(soilP_);
	}

private:

	// cm -> Pa
	Scalar toPa_(Scalar ph) const {
		return pRef_ +  ph/100.*rho_*g_;
	}

	// on root collar
	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_;
	}

    static constexpr Scalar g_ = 9.81; // cm / s^2
    static constexpr Scalar rho_ = 1.e3; // kg / m^3
    static constexpr Scalar pRef_ = 1.e5; // Pa

	static constexpr Scalar eps_ = 1e-8;

	Scalar scenario_;
	Scalar soilP_;
	Scalar p0_;
	Scalar pL_;

};

} //end namespace Dumux

#endif
