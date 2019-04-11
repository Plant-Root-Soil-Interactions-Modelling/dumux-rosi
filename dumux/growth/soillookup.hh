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
#include <RootSystem.h>

#include <dumux/io/inputfilefunction.hh>

namespace Dumux {

namespace GrowthModule {

/*!
 * A CRootBox soillookup implementation for Dumux
 * using a BoundingBoxTree to obtain the soil element,
 * and linear finite elements for interpolation of the point
 */
template<class Grid>
class SoilLookUpBBoxTree : public CRootBox::SoilLookUp
{
    using FeCache = Dune::PQkLocalFiniteElementCache<typename Grid::ctype, double, 3, 1>;
    using ShapeValue = typename Dune::FieldVector<double, 1>;
    using GridView = typename Grid::LeafGridView;
    using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using BBoxTree = BoundingBoxTree< GridViewGeometricEntitySet<GridView, 0> >;

public:

    SoilLookUpBBoxTree(const GridView& gridView, const BBoxTree& tree, const std::vector<double>& sat) :
        mapper_(gridView, Dune::mcmgVertexLayout()), gridView_(gridView), bBoxTree_(tree), sat_(sat)
    {
    }

    //! Returns the interpolated saturation, pos [cm]
    double getValue(const CRootBox::Vector3d& pos, const CRootBox::Root* root = nullptr) const final {
        const auto globalPos = Dune::FieldVector<double, 3>( { pos.x * 0.01, pos.y * 0.01, pos.z * 0.01 });
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

        for (int i = 0; i < geo.corners(); ++i) {
            sat += shapeValues_[i][0]*sat_[mapper_.subIndex(element, i, Grid::dimension)];
        }

        return sat;
    }

    //! Returns the element index, -1 if no element was found in tree
    int pick(const CRootBox::Vector3d& pos) {
        return pick(Dune::FieldVector<double, 3>( { pos.x * 0.01, pos.y * 0.01, pos.z * 0.01 }));
    }

    //! Returns the element index, -1 if no element was found in tree
    int pick(const Dune::FieldVector<double, 3>& pos) {
        const auto entities = intersectingEntities(pos, bBoxTree_);
        if (entities.empty()) {
            return -1;
        }
        const auto element = bBoxTree_.entitySet().entity(entities[0]);
        return element.geometry(); //mapper_.index(element);
    }

    std::string toString() const final {
        return "SoilLookUpBBoxTree, linear interpolation using a bounding box tree";
    }

private:
    const FeCache feCache_;
    mutable std::vector<ShapeValue> shapeValues_;
    const Mapper mapper_;
    const GridView gridView_;
    const BBoxTree& bBoxTree_;
    const std::vector<double>& sat_;

};




/*!
 * (TODO) A soillookup implementation for Dumux using a static soil, passed via the input file or grid file
 * see, InputFileFunction
 */
template<class FVGridGeometry>
class SoilLookUpTable: public CRootBox::SoilLookUp {

    //using BBoxTree = BoundingBoxTree< GridViewGeometricEntitySet<GridView, 0> >;

public:

    SoilLookUpTable(const InputFileFunction iff, std::shared_ptr<const FVGridGeometry> fvGridGeometry) : // , const BBoxTree* tree = nullptr
        iff_(iff), fvGridGeometry_(fvGridGeometry) { // , bBoxTree_(tree)
//        if (((iff_.type() == iff_.data) || (iff_.type() == iff_.perType)) && (tree == nullptr)) {
//            std::cout << "SoilLookTable: grid data is only available if box tree is set ";
//        } // TODO should be an assertion
    }

    //! Returns the saturation, pos [cm]
    double getValue(const CRootBox::Vector3d& pos, const CRootBox::Root* root = nullptr) const final {
        const auto p = Dune::FieldVector<double, 3>( { pos.x * 0.01, pos.y * 0.01, pos.z * 0.01 });
        size_t eIdx = 0;
//        if (bBoxTree_ != nullptr) {
//            const auto entities = intersectingEntities(p, bBoxTree_); // function from <dumux/common/geometry/intersectingentities.hh>
//            if (!entities.empty()) {
//                eIdx = fvGridGeometry_().elementMapper().index(entities);
//            } else {
//                std::cout << "SoilLookTable: no element found at " << p;
//            }
//        }
        double v = iff_.f(p[2], eIdx);
        // std::cout << "soil " << p[2] << ", " << v << "\n";
        return v;
    }

    std::string toString() const final {
        return "SoilLookUpTable";
    }

private:

    // const BBoxTree* bBoxTree_;
    const InputFileFunction iff_;
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

};

} // namespace Growth module

} // end namespace Dumux

#endif
