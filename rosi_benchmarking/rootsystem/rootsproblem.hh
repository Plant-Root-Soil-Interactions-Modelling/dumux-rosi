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
#ifndef DUMUX_ONEP_TUBES_TEST_PROBLEM_HH
#define DUMUX_ONEP_TUBES_TEST_PROBLEM_HH

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

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
	using type = FluidSystems::LiquidPhase<Scalar, Components::Constant<1, Scalar> >;
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
	using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

	// Grid and world dimension
	static const int dim = GridView::dimension;
	static const int dimWorld = GridView::dimensionworld;

	using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
	enum {
		// indices of the primary variables
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

	enum { isBox = GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod == DiscretizationMethod::box };

public:

	RootsProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry) {

		name_ = getParam<std::string>("Problem.Name");
		soilP_ = 1.e5 - toPa_(getParam<Scalar>("Parameter.SoilP"));

		//get hMax_ of the grid
		hMax_ = 0.0;
		for (const auto& element : elements(fvGridGeometry->gridView())) {
			hMax_ = std::max(element.geometry().volume(), hMax_);
		}

		std::cout << "hMax : " << hMax_ << "\n";

	}

	/*!
	 * \name Problem parameters
	 */
	// \{

	/*!
	 * \brief The problem name.
	 *
	 * This is used as a prefix for files generated by the simulation.
	 */
	const std::string& name() const {
		return name_;
	}

	/*!
	 * \brief Return the temperature within the domain in [K].
	 *
	 */
	Scalar temperature() const {
		return 273.15 + 10.0;  // temperature
	}

	/*!
	 * \brief Return how much the domain is extruded at a given sub-control volume.
	 *
	 * This means the factor by which a lower-dimensional (1D or 2D)
	 * entity needs to be expanded to get a full dimensional cell. The
	 * default is 1.0 which means that 1D problems are actually
	 * thought as pipes with a cross section of 1 m^2 and 2D problems
	 * are assumed to extend 1 m to the back.
	 */
	template<class ElementSolution>
	Scalar extrusionFactor(const Element &element,
			const SubControlVolume &scv,
			const ElementSolution& elemSol) const
	{
		// std::cout << " Someone actually called me \n";
		auto radius = this->spatialParams().radius(scv);
		return M_PI*radius*radius;
	}

	// \}
	/*!
	 * \name Boundary conditions
	 */
	// \{

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
		if (onUpperBoundary_(pos)) { // top bc
			std::cout << "top bc " << pos << " dirichlet"<< "\n";
			bcTypes.setAllDirichlet();
		} else if(onLowerBoundary_(pos)) {
			std::cout << "bot bc " << pos << " neumann"<< "\n";
			bcTypes.setAllDirichlet();
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
		Scalar p0 = 0.;
		std::cout << "Dirichlet " << p0 << " [Pa] \n";
		if (onUpperBoundary_(pos)) { // top bc
			p0 = 1.e5 - toPa_(-1000);
		} else if(onLowerBoundary_(pos)) {
			p0 = 1.e5 - toPa_(-500);
		}
		return PrimaryVariables(p0);
	}

	// \}

	/*!
	 * \name Volume terms
	 */
	// \{
	/*!
	 * \brief Evaluate the source term for all phases within a given
	 *        sub-control-volume.
	 *
	 * This is the method for the case where the source term is
	 * potentially solution dependent and requires some quantities that
	 * are specific to the fully-implicit method.
	 *
	 * \param element The finite element
	 * \param fvGeometry The finite-volume geometry
	 * \param elemVolVars All volume variables for the element
	 * \param scv The sub control volume
	 *
	 * For this method, the return parameter stores the conserved quantity rate
	 * generated or annihilate per volume unit. Positive values mean
	 * that the conserved quantity is created, negative ones mean that it vanishes.
	 * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
	 */
	ResidualVector source(const Element &element,
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
		values[conti0EqIdx] /= 1000; // rho (kg/s/m^3)

		std::cout << "l " << l << ", r " << r << " eq " << ( kr * 2*r*M_PI*l * (phs - phx) ) << " vol " << (r*r*M_PI*l) << ", final " << values[conti0EqIdx] << "\n";

		return values;
	}

	/*!
	 * \brief Evaluate the initial value for a control volume.
	 *
	 * For this method, the \a priVars parameter stores primary
	 * variables.
	 */
	PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
	{
		// std::cout << "Initial: " << soilP_ << "\n";
		return PrimaryVariables(soilP_);
	}

private:

	/**
	 *  pressure head to pascal
	 *
	 *  @param ph           pressure head [cm]
	 */
	Scalar toPa_(Scalar ph) const {
		const Scalar g = 9.81; // TODO stupid gravity bug
		return -ph*10.*g;
	}

	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_;
	}

	bool onLowerBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld-1] < this->fvGridGeometry().bBoxMin()[dimWorld-1] + eps_;
	}

	static constexpr Scalar eps_ = 1e-8;
	std::string name_;
	Scalar hMax_;
	Scalar soilP_;

};

} //end namespace Dumux

#endif
