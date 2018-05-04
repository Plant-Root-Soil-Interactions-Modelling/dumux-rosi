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
 * \brief Represents the finite volume geometry of a single element in
 *        the cell-centered fv scheme.
 */
#ifndef DUMUX_MPFA_FV_ELEMENTGEOMETRY_HH
#define DUMUX_MPFA_FV_ELEMENTGEOMETRY_HH

#include <dune/common/version.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/intersectioniterator.hh>

#include <dumux/common/propertysystem.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
}

/*!
 * \ingroup CCModel
 * \brief Represents the finite volume geometry of a single element in
 *        the cell-centered fv-mpfa scheme.
 */
template<class TypeTag>
class MpfaFVElementGeometry
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolume) InteractionVolume;
    typedef typename GET_PROP_TYPE(TypeTag, MpfaInteractionVolumeManager) InteractionVolumeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum{dim = GridView::dimension};
    enum{dimWorld = GridView::dimensionworld};

    enum{maxNE = InteractionVolume::Properties::StandardStencilSize}; //! Maximum number of elements in stencil
    enum{maxFACES = dim < 3 ? 4 : 6}; //! maximum number of faces of element

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:

    struct SubControlVolume //! FV intersected with element
    {
        LocalPosition local; //!< local position
        GlobalPosition global; //!< global position
        Scalar volume; //!< volume of scv
        bool inner;
    };

    struct SubControlVolumeFace //! interior face of a sub control volume
    {
        int i,j; //!< scvf seperates corner i and j of elem
    };

    typedef SubControlVolumeFace BoundaryFace; //!< compatibility typedef

    LocalPosition elementLocal; //!< local coordinate of element center
    GlobalPosition elementGlobal; //!< global coordinate of element center
    Scalar elementVolume; //!< element volume
    SubControlVolume subContVol[1]; //!< data of the sub control volumes
    SubControlVolumeFace subContVolFace[maxFACES]; //!< data of the sub control volume faces

    int numScv; //!< number of subcontrol volumes
    int numScvf; //!< number of inner-domain subcontrolvolume faces
    int numNeighbors; //!< number of neighboring elements including the element itself
    std::vector<Element> neighbors; //!< stores the neighboring elements
    std::map<int, int> globalIdxToStencilIdx;

    int findElementInStencil(int globalElIdx) const
    {
        std::map<int, int>::const_iterator it = globalIdxToStencilIdx.find(globalElIdx);
        if (it != globalIdxToStencilIdx.end())
            return it->second;
        else
            DUNE_THROW(Dune::NotImplemented,
                      "Element not found in stencil!! See mpfafvelementgeometry.hh, l. 116");
    }

    void updateInner(const Element& element)
    {
        const Geometry geometry = element.geometry();

        elementVolume = geometry.volume();
        elementGlobal = geometry.center();
        elementLocal = geometry.local(elementGlobal);

        numScv = 1;
        numScvf = 0;

        subContVol[0].local = elementLocal;
        subContVol[0].global = elementGlobal;
        subContVol[0].inner = true;
        subContVol[0].volume = elementVolume;

        // initialize neighbors list with self:
        numNeighbors = 1;
        neighbors.clear();
        neighbors.reserve(maxNE);
        neighbors.push_back(element);
        globalIdxToStencilIdx.clear();
    }

    void update(const GridView& gridView, const Element& element, const Problem &problem)
    {
        updateInner(element);

        // get informations on the elements of the stencil
        numNeighbors = InteractionVolumeManager::getElementsOfStencil(neighbors, problem);

        // now fill the mapping of the global element idx to idx in stencil
        for(int i = 0; i < numNeighbors; i++)
        {
            // map element to get global index
            int globalIdx = problem.elementMapper().index(neighbors[i]);
            std::pair<int, int> indexPair(globalIdx,i);
            globalIdxToStencilIdx.insert(indexPair);
        }

        int scvfIdx = 0;

        // fill control volume face data:
        IntersectionIterator isEndIt = gridView.iend(element);
        for (IntersectionIterator isIt = gridView.ibegin(element); isIt != isEndIt; ++isIt)
        {
            scvfIdx++;
            const auto isGeometry = isIt->geometry();
            SubControlVolumeFace& scvFace = subContVolFace[isIt->indexInInside()];

            // set the indices of the neighbouring volume variables
            scvFace.i = 0;
            if(isIt->neighbor())
            {
                int globalIdx = problem.elementMapper().index(isIt->outside());
                scvFace.j = findElementInStencil(globalIdx);
            }
            else
                scvFace.j = numNeighbors + isIt->indexInInside();
        }
        // set the number of inner-domain subcontrolvolume faces
        numScvf = scvfIdx;
    }
};

}

#endif

