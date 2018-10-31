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

namespace Dumux {
namespace GrowthModule {

/*!
 * \brief A croot box soillookup implementation for dumux using vertex data
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
    SoilLookUpBBoxTree(const GridView& gridView, const BBoxTree& tree, const std::vector<double>& sat)
    : mapper_(gridView, Dune::mcmgVertexLayout())
    , gridView_(gridView)
    , bBoxTree_(tree)
    , sat_(sat)
    , weight_(getParam<double>("SoilLookup.Weight", 1.0))
    , power_(getParam<double>("SoilLookup.Power", 1.0))
    {}

    //! Returns a scalar property of the soil scaled from 0..1
    double getValue(const CRootBox::Vector3d& pos, const CRootBox::Root* root = nullptr) const final
    {
        const auto globalPos = GrowthModule::CRootBoxInterface::convert(pos);
        const auto entities = intersectingEntities(globalPos, bBoxTree_);
        if (entities.empty())
            return 0.0;

        const auto element = bBoxTree_.entitySet().entity(entities[0]);
        const auto geo = element.geometry();
        double sat = 0.0;

        const auto& localBasis = feCache_.get(geo.type()).localBasis();
        const auto ipLocal = geo.local(globalPos);
        localBasis.evaluateFunction(ipLocal, shapeValues_);

        for (int i = 0; i < geo.corners(); ++i)
            sat += shapeValues_[i][0]*sat_[mapper_.subIndex(element, i, Grid::dimension)];

        return std::pow(sat, power_)*weight_;
    }

    std::string toString() const final
    { return "SoilLookUp with bounding box tree"; }

private:
    const FeCache feCache_;
    mutable std::vector<ShapeValue> shapeValues_;
    const Mapper mapper_;
    const GridView gridView_;
    const BBoxTree& bBoxTree_;
    const std::vector<double>& sat_;
    double weight_;
    double power_;
};

} // end namespace GrowthModule
} // end namespace Dumux

#endif
