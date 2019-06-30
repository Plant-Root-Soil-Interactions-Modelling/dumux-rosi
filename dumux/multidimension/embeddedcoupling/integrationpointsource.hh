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
 * \ingroup EmbeddedCoupling
 * \brief An integration point source class,
 *        i.e. sources located at a single point in space
 *        associated with a quadrature point
 */

#ifndef DUMUX_INTEGRATION_POINTSOURCE_HH
#define DUMUX_INTEGRATION_POINTSOURCE_HH

#include <dumux/common/pointsource.hh>

namespace Dumux
{

/*!
 * \ingroup EmbeddedCoupling
 * \brief An integration point source class with an identifier to attach data
 *        and a quadrature weight and integration element
 */
template<class TypeTag, typename IdType>
class IntegrationPointSource : public Dumux::IdPointSource<TypeTag, IdType>
{
    typedef typename Dumux::IdPointSource<TypeTag, IdType> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    static const int dimworld = GridView::dimensionworld;
    typedef typename Dune::FieldVector<Scalar, dimworld> GlobalPosition;

public:
    //! Constructor for integration point sources
    IntegrationPointSource(GlobalPosition pos, PrimaryVariables values, IdType id,
                           Scalar qpweight, Scalar integrationElement,
                           const std::vector<unsigned int>& elementIndices)
      : ParentType(pos, values, id),
        qpweight_(qpweight), integrationElement_(integrationElement),
        elementIndices_(elementIndices) {}

    //! Constructor for integration point sources, when there is no
    // value known at the time of initialization
    IntegrationPointSource(GlobalPosition pos, IdType id,
                           Scalar qpweight, Scalar integrationElement,
                           const std::vector<unsigned int>& elementIndices)
      : ParentType(pos, id),
        qpweight_(qpweight), integrationElement_(integrationElement),
        elementIndices_(elementIndices) {}


    Scalar quadratureWeight() const
    {
        return qpweight_;
    }

    Scalar integrationElement() const
    {
        return integrationElement_;
    }

    const std::vector<unsigned int>& elementIndices() const
    {
        return elementIndices_;
    }

    //! Convenience = operator overload modifying only the values
    IntegrationPointSource& operator= (const PrimaryVariables& values)
    {
        ParentType::operator=(values);
        return *this;
    }

    //! Convenience = operator overload modifying only the values
    IntegrationPointSource& operator= (Scalar s)
    {
        ParentType::operator=(s);
        return *this;
    }

private:
    Scalar qpweight_;
    Scalar integrationElement_;
    std::vector<unsigned int> elementIndices_;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief A helper class calculating a DOF-index to point source map
 */
template<class TypeTag>
class IntegrationPointSourceHelper
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;

    static const int dim = GridView::dimension;
    static const int dimworld = GridView::dimensionworld;

    typedef Dumux::BoundingBoxTree<GridView> BoundingBoxTree;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    //! calculate a DOF index to point source map from given vector of point sources
    static void computePointSourceMap(const Problem& problem,
                                      const BoundingBoxTree& boundingBoxTree,
                                      std::vector<PointSource>& sources,
                                      std::map<std::pair<unsigned int, unsigned int>, std::vector<PointSource> >& pointSourceMap)
    {
        for (auto&& source : sources)
        {
            // compute in which elements the point source falls
            std::vector<unsigned int> entities = source.elementIndices();
            // split the source values equally among all concerned entities
            source.setEmbeddings(source.embeddings()*entities.size());
            // loop over all concernes elements
            for (unsigned int eIdx : entities)
            {
                if(isBox)
                {
                    // check in which subcontrolvolume(s) we are
                    const auto element = boundingBoxTree.entity(eIdx);
                    FVElementGeometry fvGeometry;
                    fvGeometry.update(problem.gridView(), element);
                    const auto globalPos = source.position();
                    // loop over all sub control volumes and check if the point source is inside
                    std::vector<unsigned int> scvs;
                    for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                    {
                        auto geometry = fvGeometry.subContVolGeometries[scvIdx];
                        if (BoundingBoxTreeHelper<dimworld>::pointInGeometry(geometry, globalPos))
                            scvs.push_back(scvIdx);
                    }
                    // for all scvs that where tested positiv add the point sources
                    // to the element/scv to point source map
                    for (unsigned int scvIdx : scvs)
                    {
                        const auto key = std::make_pair(eIdx, scvIdx);
                        if (pointSourceMap.count(key))
                            pointSourceMap.at(key).push_back(source);
                        else
                            pointSourceMap.insert({key, {source}});
                        // split equally on the number of matched scvs
                        auto& s = pointSourceMap.at(key).back();
                        s.setEmbeddings(s.embeddings()*scvs.size());
                    }
                }
                else
                {
                    // add the pointsource to the DOF map
                    const auto key = std::make_pair(eIdx, /*scvIdx=*/ 0);
                    if (pointSourceMap.count(key))
                        pointSourceMap.at(key).push_back(source);
                    else
                        pointSourceMap.insert({key, {source}});
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
