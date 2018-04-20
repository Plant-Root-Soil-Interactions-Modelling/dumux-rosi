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
#ifndef ROOTS_PARAMS_HH
#define ROOTS_PARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

#include <dumux/material/components/simpleh2o.hh>

namespace Dumux {

/*!
 * Managing the root parameters
 */
template<class TypeTag>
class RootsParams : public FVSpatialParamsOneP<TypeTag>
{
	using ParentType = FVSpatialParamsOneP<TypeTag>;
	using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
	using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
	using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
	using ElementMapper = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::ElementMapper;
	using SubControlVolume = typename FVElementGeometry::SubControlVolume;
	using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
	using Element = typename GridView::template Codim<0>::Entity;
	using Water = Components::SimpleH2O<Scalar>;
	using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);

public:

	using PermeabilityType = Scalar;

	RootsParams(const Problem& problem) : ParentType(problem) {
		try { // parameters are constant for all elements
			radius_ = getParam<Scalar>("Parameter.Radius");
			kr_ = getParam<Scalar>("Parameter.Kr");
			kz_ = getParam<Scalar>("Parameter.Kz");
		} catch  (const std::exception& e) {
			params_=true; // parameters are set by element
	        const auto& gridView = this->problem().fvGridGeometry().gridView();
	        radii_.resize(gridView.size(0));
	        krs_.resize(gridView.size(0));
	        kzs_.resize(gridView.size(0));
	        for (const auto& element : elements(gridView)) {
	            const auto eIdx = elementMapper_.index(element);
	            auto p = GridCreator::parameters(element);
	            radii_[eIdx] = p.at(0);
	            krs_[eIdx] = p.at(1);
	            kzs_[eIdx] = p.at(2);
	        }
		}
	}

	/**
	 * Root radius (m)
	 */
	Scalar radius(const Element &element) const {
		if (params_) {
			auto eIdx = elementMapper_.index(element);
			return radii_[eIdx];
		} else {
			return radius_;
		}
	}

	/**
	 * Root radial conductivity (m^2 s / kg)
	 */
	Scalar radialConductivity(const Element &element) const {
		if (params_) {
			auto eIdx = elementMapper_.index(element);
			return krs_[eIdx];
		} else {
			return kr_;
		}
	}

	/**
	 * Root axial conductivity (m^5 s / kg)
	 */
	Scalar axialConductivity(const Element &element) const {
		if (params_) {
			auto eIdx = elementMapper_.index(element);
			return kzs_[eIdx];
		} else {
			return kz_;
		}
	}

	/*!
	 * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
	 */
	template<class ElementSolution>
	Scalar permeability(const Element& element,
			const SubControlVolume& scv,
			const ElementSolution& elemSol) const {
		Scalar mu = Water::liquidViscosity(0.,0.); // temperature, pressure
		Scalar a = this->radius(element);
		Scalar kz = this->axialConductivity(element);
		return kz*mu/(a*a*M_PI); 		// a^2 * k / mu = kz  --> k = kz/a^2*mu
	}

	/*!
	 * \brief Returns the porosity \f$[-]\f$
	 */
	template<class ElementSolution>
	Scalar porosity(const Element& element,
			const SubControlVolume& scv,
			const ElementSolution& elemSol) const {
		return 1.0;
	}

private:

	Scalar radius_;
	Scalar kr_;
	Scalar kz_;

	bool params_ = false;
    const ElementMapper& elementMapper_ = this->problem().fvGridGeometry().elementMapper();
    std::vector<Scalar> radii_ = std::vector<Scalar>(0);
	std::vector<Scalar> krs_ = std::vector<Scalar>(0);
	std::vector<Scalar> kzs_ = std::vector<Scalar>(0);

};

} // end namespace Dumux

#endif
