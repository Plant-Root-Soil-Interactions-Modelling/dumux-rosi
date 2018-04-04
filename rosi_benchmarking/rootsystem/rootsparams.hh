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
#ifndef DUMUX_ONEP_TUBES_TEST_SPATIALPARAMS_HH
#define DUMUX_ONEP_TUBES_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

#include <dumux/material/components/simpleh2o.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
template<class TypeTag>
class RootsParams : public FVSpatialParamsOneP<TypeTag>
{
    using ParentType = FVSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Water = SimpleH2O<Scalar>;

public:

    using PermeabilityType = Scalar;

    RootsParams(const Problem& problem) : ParentType(problem) {
        radius_ = getParam<Scalar>("Parameter.Radius");
        kr_ = getParam<Scalar>("Parameter.Kr");
        kz_ = getParam<Scalar>("Parameter.Kz");
    }

    /*!
     * \brief Return the radius of the circular pipe for the current sub-control volume in [m].
     *
     * \param scv The sub control volume
     */
    Scalar radius(const SubControlVolume &scv) const {
    	return radius_;
    }

    /**
     * The root radial conductivity (m^2 s / kg)
     */
    Scalar radialConductivity(const SubControlVolume &scv) const {
    	return kr_;
    }


    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    template<class ElementSolution>
    Scalar permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const {

    	// a^2 * k / mu = kz  --> k = kz/a^2*mu
    	Scalar mu = Water::liquidViscosity(0.,0.); // temperature, pressure
        return kz_*mu/(radius_*radius_*M_PI);
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the porosity
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

};

} // end namespace Dumux

#endif
