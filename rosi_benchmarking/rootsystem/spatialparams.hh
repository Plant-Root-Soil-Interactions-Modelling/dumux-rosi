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

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv1p.hh>

#include <dumux/material/components/simpleh2o.hh>

namespace Dumux {

/*!
 * minimal class for root systems, constant radius, constant conductivities from input file
 */
template<class FVGridGeometry, class Scalar>
class TubesTestSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                             TubesTestSpatialParams<FVGridGeometry, Scalar>>
{
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                           TubesTestSpatialParams<FVGridGeometry, Scalar>>;
    using Water = Components::SimpleH2O<Scalar>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    TubesTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        radius_ = getParam<Scalar>("RootSystem.Radius") / 100.; // cm-> m
        kr_ = getParam<Scalar>("RootSystem.Conductivity.Kr");
        kx_ = getParam<Scalar>("RootSystem.Conductivity.Kx");
    }

    /*!
     * \brief Return the radius of the circular pipe for the current sub-control volume in [m].
     *
     * \param scv The sub control volume
     */
    Scalar radius(const Element& element) const
    {
        return radius_;
    }

    Scalar axialConductivity(const Element& element) const {
        return kx_;
    }

    Scalar radialConductivity(const Element&element) const {
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
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {

        Scalar mu = Water::liquidViscosity(285.15, 1e5); // temperature, pressure
        Scalar a = this->radius(element);
        Scalar kx = this->axialConductivity(element);
        // std::cout << "params " << kx * 1e13 << ", " << a << ", " << mu << "\n";
        return kx * mu / (a * a * M_PI);        // a^2 * k / mu = kz  --> k = kz/a^2*mu
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 1;
    }

private:

    Scalar radius_;
    Scalar kr_;
    Scalar kx_;

};

} // end namespace Dumux

#endif
