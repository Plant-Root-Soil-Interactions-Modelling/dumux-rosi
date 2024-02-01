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
 *
 * \brief A croot box soillookup implementation for dumux
 */
#ifndef DUMUX_SOIL_LOOKUP_BBOXTREE_HH
#define DUMUX_SOIL_LOOKUP_BBOXTREE_HH

#include <cmath>
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <soil.h>

#include <dumux/io/inputfilefunction.hh>

namespace Dumux {

namespace GrowthModule {

/*!
 * A CPlantBox soillookup implementation for Dumux
 * using a BoundingBoxTree to obtain the soil element,
 * and linear finite elements for interpolation of the point
 *
 * todo: do i have to update the mappers??? better use fvGridGeometry?
 */
template<class FVGridGeometry>
class SoilLookUpBBoxTree : public CPlantBox::SoilLookUp
{
    using Grid = typename FVGridGeometry::Grid;
    using FeCache = Dune::PQkLocalFiniteElementCache<typename Grid::ctype, double, 3, 1>;
    using ShapeValue = typename Dune::FieldVector<double, 1>;
    using BBoxTree = typename FVGridGeometry::BoundingBoxTree;

public:

    SoilLookUpBBoxTree(const FVGridGeometry& fvGridGeometry, const std::vector<double>& sat, bool periodic = true) :
        fvGridGeometry_(fvGridGeometry),
        sat_(sat),
        bBoxTree_(fvGridGeometry.boundingBoxTree())    {

        if (periodic) {
            const auto size = fvGridGeometry_.bBoxMax() - fvGridGeometry_.bBoxMin();
            // setPeriodicDomain(100.*-size[0]/2.,100.*size[0]/2., 100.*-size[1]/2., 100.*size[1]/2. );
            setPeriodicDomain(0., 100.*size[0], 0., 100.*size[1]);
        }

    }

    /**
     *  Returns the interpolated saturation, pos [cm]
     */
    double getValue(const CPlantBox::Vector3d& pos,
    		const std::shared_ptr<CPlantBox::Organ> organ = nullptr) const final {

        auto p = periodic(pos.plus(shiftRB)); // periodic mapping
        const auto globalPos = Dune::FieldVector<double, 3>( { p.x * 0.01, p.y * 0.01, p.z * 0.01 });

        // std::cout << "getValue() "<< pos.toString() << "->" << p.toString() << " [cm] \n" << std::flush;

        const auto entities = intersectingEntities(globalPos, bBoxTree_); // function from <dumux/common/geometry/intersectingentities.hh>

        if (entities.empty()) {
            return 0.0;
        }

        const auto element = bBoxTree_.entitySet().entity(entities[0]);
        const auto geo = element.geometry();
        double sat = 0.0;

        const auto& localBasis = feCache_.get(geo.type()).localBasis();
        const auto ipLocal = geo.local(globalPos);
        localBasis.evaluateFunction(ipLocal, shapeValues_);


        auto vMapper = fvGridGeometry_.vertexMapper();
        for (int i = 0; i < geo.corners(); ++i) {
            sat += shapeValues_[i][0]*sat_[vMapper.subIndex(element, i, Grid::dimension)];
        }

        return sat;
    }

    /*
     * Returns the element index, -1 if no element was found in tree
     */
    int pick(const Dune::FieldVector<double, 3>& pos) {

        auto p = pos + shift; // shift
        auto pp = periodic(CPlantBox::Vector3d(100*p[0], 100*p[1], 100*p[2])); // periodic mapping
        p[0] = pp.x/100.; p[1] = pp.y/100.; p[2] = pp.z/100.;

        // std::cout << "pick() "<< pos << "->" << p << " [m] \n" << std::flush;

        const auto entities = intersectingEntities(p, bBoxTree_);
        if (entities.empty()) {
            return -1;
        }
        const auto element = bBoxTree_.entitySet().entity(entities[0]);
        int eIdx = fvGridGeometry_.elementMapper().index(element);
        return eIdx;
    }

    /**
     * shifts incoming picks
     */
    void setShift(Dune::FieldVector<double, 3> p) {
        shift = p;
        shiftRB = CPlantBox::Vector3d(100.*shift[0], 100.*shift[1], 100.*shift[2]);
    }

    std::string toString() const final {
        return "SoilLookUpBBoxTree, linear interpolation using a bounding box tree";
    }

private:

    const FeCache feCache_;
    mutable std::vector<ShapeValue> shapeValues_;

    const FVGridGeometry& fvGridGeometry_;
    const std::vector<double>& sat_;
    const BBoxTree& bBoxTree_;

    Dune::FieldVector<double, 3> shift = { 0., 0., 0. };
    CPlantBox::Vector3d shiftRB = CPlantBox::Vector3d();
};







/*!
 * A SoilLookUp implementation for Dumux using a static soil, passed via the input file [or grid file todo]
 * see, InputFileFunction
 *
 * todo periodicity
 */
class SoilLookUpTable: public CPlantBox::SoilLookUp {

public:

    SoilLookUpTable(const InputFileFunction iff):
        iff_(iff) {

    }

    //! Returns the saturation, pos [cm]
    double getValue(const CPlantBox::Vector3d& pos,
    		const std::shared_ptr<CPlantBox::Organ> organ = nullptr) const final {
        auto p = Dune::FieldVector<double, 3>( { pos.x * 0.01, pos.y * 0.01, pos.z * 0.01 });
        double v = iff_.f(p[2]);
        return v;
    }

    std::string toString() const final {
        return "SoilLookUpTable";
    }

private:

    const InputFileFunction iff_;

};

} // namespace Growth module

} // end namespace Dumux

#endif
