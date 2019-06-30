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
 *        the for network grids.
 */
#ifndef DUMUX_ONE_D_FV_ELEMENTGEOMETRY_HH
#define DUMUX_ONE_D_FV_ELEMENTGEOMETRY_HH

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dumux/common/propertysystem.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(ImplicitIsBox);
}

/*!
 * \ingroup OneDModel
 * \brief Represents the finite volume geometry of a single element in
 *        for network grids.
 */
template<class TypeTag>
class OneDFVElementGeometry
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
    enum{dimWorld = GridView::dimensionworld};

    enum{maxBF = 2*dim}; //! maximum number of faces

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef typename Dune::AffineGeometry<CoordScalar, dim, dimWorld> AffineGeometry;
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { maxNumScv = isBox ? 2 : 1};

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
        SubControlVolumeFace()
        {
            //default is two-point flux
            grad.reserve(2);
            shapeValue.reserve(2);
            fapIndices.reserve(2);
            fapDistances.reserve(2);
            numFap = 2;
        }

        int i,j; //!< scvf seperates element i (this scv) and j (neighbor)
        LocalPosition ipLocal; //!< integration point in local coords
        GlobalPosition ipGlobal; //!< integration point in global coords
        GlobalPosition normal; //!< normal on face pointing to CV j or outward of the domain with length equal to |scvf|
        Scalar area; //!< area of face
        std::vector<GlobalPosition> grad; //!< derivatives of shape functions at ip
        std::vector<Scalar> shapeValue; //!< value of shape functions at ip
        std::vector<int> fapIndices; //!< indices w.r.t.neighbors of the flux approximation points
        std::vector<Scalar> fapDistances; //!< distance of the flux approximation points to the integration point
        unsigned int numFap; //!< number of flux approximation points
        unsigned int fIdx; //!< the index (w.r.t. the element) of the face (codim 1 entity) that the scvf is part of
    };

    typedef SubControlVolumeFace BoundaryFace; //!< compatibility typedef

    LocalPosition elementLocal; //!< local coordinate of element center
    GlobalPosition elementGlobal; //!< global coordinate of element center
    Scalar elementVolume; //!< element volume
    SubControlVolume subContVol[maxNumScv]; //!< data of the sub control volumes
    std::vector<AffineGeometry> subContVolGeometries; //!< geometries of the subcontrol volumes
    std::vector<SubControlVolumeFace> subContVolFace; //!< data of the sub control volume faces
    BoundaryFace boundaryFace[maxBF]; //!< data of the boundary faces
    int numScv; //!< number of subcontrol volumes
    int numScvf; //!< number of inner-domain subcontrolvolume faces
    int numNeighbors; //!< number of neighboring elements including the element itself
    std::vector<Element> neighbors; //!< stores pointers for the neighboring elements

    void updateInner(const Element& element)
    {
        const Geometry geometry = element.geometry();

        elementVolume = geometry.volume();
        elementGlobal = geometry.center();
        elementLocal = geometry.local(elementGlobal);

        numScv = maxNumScv;
        numScvf = isBox ? 1 : element.subEntities(1); //initialization
        if(!isBox)
        {
            subContVol[0].local = elementLocal;
            subContVol[0].global = elementGlobal;
            subContVol[0].inner = true;
            subContVol[0].volume = elementVolume;
        }
        else
        {
            subContVol[0].local = 0.25;
            subContVol[0].global = geometry.global(0.25);
            subContVol[0].inner = true;
            subContVol[0].volume = elementVolume/2.0;

            subContVol[1].local = 0.75;
            subContVol[1].global = geometry.global(0.75);
            subContVol[1].inner = true;
            subContVol[1].volume = elementVolume/2.0;

            subContVolGeometries.clear();
            std::array<GlobalPosition, 2> corners = {{geometry.corner(0), geometry.global(0.5)}};
            subContVolGeometries.push_back(AffineGeometry(geometry.type(), corners));
            corners = {{geometry.global(0.5), geometry.corner(1)}};
            subContVolGeometries.push_back(AffineGeometry(geometry.type(), corners));
        }

        // initialize neighbors list with self:
        numNeighbors = 1;
        neighbors.clear();
        neighbors.reserve(numScvf+1);
        neighbors.push_back(element);
    }

    void update(const GridView& gridView, const Element& element, const Problem& problem)
    {
        update(gridView, element);
    }

    void update(const GridView& gridView, const Element& element)
    {
        updateInner(element);

        const Geometry geometry = element.geometry();

        subContVolFace.clear();
        subContVolFace.reserve(numScvf); //reserve memory

        if(isBox)
        {
            SubControlVolumeFace scvFace;
            scvFace.i = 0;
            scvFace.j = 1;

            scvFace.ipGlobal = geometry.global(0.5);
            scvFace.ipLocal =  geometry.local(scvFace.ipGlobal);

            scvFace.numFap = 2;
            scvFace.fapIndices.push_back(scvFace.i);
            scvFace.fapIndices.push_back(scvFace.j);
            scvFace.fapDistances.push_back(0.5*geometry.volume());
            scvFace.fapDistances.push_back(0.5*geometry.volume());

            subContVolFace.push_back(scvFace);

            // treat boundaries
            for (auto&& intersection : intersections(gridView, element))
            {
                if(intersection.boundary())
                {
                    int bfIdx = intersection.indexInInside();
                    SubControlVolumeFace& bFace = boundaryFace[bfIdx];

                    bFace.grad.clear();
                    bFace.shapeValue.clear();
                    bFace.fapIndices.clear();
                    bFace.fapDistances.clear();

                    bFace.ipGlobal = intersection.geometry().center();
                    bFace.ipLocal =  geometry.local(bFace.ipGlobal);
                    bFace.area = 1.0;

                    bFace.i = bfIdx;
                    bFace.j = bfIdx;

                    scvFace.numFap = 2;
                    bFace.fapIndices.push_back(bFace.i);
                    bFace.fapIndices.push_back(bFace.j);
                    bFace.fapDistances.push_back(0.5*geometry.volume());
                    bFace.fapDistances.push_back(0.5*geometry.volume());
                }
            }
        }
        else
        {
            bool onBoundary = false;

            // fill neighbor information and control volume face data:
            for (const auto& intersection : intersections(gridView, element))
            {
                const auto isGeometry = intersection.geometry();

                // neighbor information and inner control volume face data:
                if (intersection.neighbor())
                {
                    numNeighbors++;
                    neighbors.push_back(intersection.outside());

                    int scvfIdx = numNeighbors - 2;
                    SubControlVolumeFace scvFace;

                    scvFace.i = 0;
                    scvFace.j = scvfIdx + 1;

                    scvFace.ipGlobal = isGeometry.center();
                    scvFace.ipLocal =  geometry.local(scvFace.ipGlobal);
                    scvFace.normal = intersection.centerUnitOuterNormal();
                    Scalar volume = isGeometry.volume();
                    scvFace.normal *= volume;
                    scvFace.area = volume;

                    GlobalPosition distVec = elementGlobal
                                           - neighbors[scvfIdx+1].geometry().center();
                    distVec /= distVec.two_norm2();

                    // gradients using a two-point flux approximation
                    for (unsigned int fapIdx = 0; fapIdx < scvFace.numFap; fapIdx++)
                    {
                        scvFace.grad.push_back(distVec);
                        scvFace.shapeValue.push_back(0.5);
                    }
                    scvFace.grad[1] *= -1.0;

                    scvFace.fapIndices.push_back(scvFace.i);
                    scvFace.fapIndices.push_back(scvFace.j);

                    //calculate the distance from element midpoint to integration point
                    GlobalPosition fapDistance = scvFace.ipGlobal;
                    fapDistance -= elementGlobal;
                    scvFace.fapDistances.push_back(fapDistance.two_norm());

                    fapDistance = scvFace.ipGlobal;
                    fapDistance -= neighbors[scvfIdx+1].geometry().center();
                    scvFace.fapDistances.push_back(fapDistance.two_norm());

                    scvFace.fIdx = intersection.indexInInside();

                    //pushback the sub control volume face
                    subContVolFace.push_back(scvFace);
                }

                // boundary control volume face data
                if (intersection.boundary())
                {
                    onBoundary = true;
                    int bfIdx = intersection.indexInInside();
                    SubControlVolumeFace& bFace = boundaryFace[bfIdx];

                    bFace.grad.clear();
                    bFace.shapeValue.clear();
                    bFace.fapIndices.clear();
                    bFace.fapDistances.clear();

                    bFace.ipGlobal = intersection.geometry().center();
                    bFace.ipLocal =  geometry.local(bFace.ipGlobal);
                    bFace.normal = intersection.centerUnitOuterNormal();
                    Scalar volume = isGeometry.volume();
                    bFace.normal *= volume;
                    bFace.area = volume;
                    bFace.i = 0;
                    bFace.j = 0; //initialize: will be overwritten later

                    GlobalPosition distVec = elementGlobal - bFace.ipGlobal;
                    distVec /= distVec.two_norm2();

                    bFace.numFap = 2;
                    // gradients using a two-point flux approximation
                    for (unsigned int fapIdx = 0; fapIdx < bFace.numFap; fapIdx++)
                    {
                        bFace.grad.push_back(distVec);
                        bFace.shapeValue.push_back(0.5);
                    }
                    bFace.grad[1] *= -1.0;

                    bFace.fapIndices.push_back(bFace.i);
                    bFace.fapIndices.push_back(bFace.j);

                    GlobalPosition fapDistance = bFace.ipGlobal;
                    fapDistance -= elementGlobal;
                    bFace.fapDistances.push_back(fapDistance.two_norm());
                    bFace.fapDistances.push_back(fapDistance.two_norm());
                }
            }

            // treat elements on the boundary
            if (onBoundary)
            {
                for (int bfIdx = 0; bfIdx < element.subEntities(1); bfIdx++)
                {
                    SubControlVolumeFace& bFace = boundaryFace[bfIdx];
                    bFace.j = numNeighbors + bfIdx;
                    bFace.fapIndices[1] = bFace.j;
                    neighbors.push_back(element);
                }
            }

            // set the number of inner-domain subcontrolvolume faces
            // this update is only needed if the number of neighbors exceeds the number of faces
            numScvf = numNeighbors - 1;

            //Special case for multiple neighbors on the same global integration point which
            //occurs in multi-dimensional grids. Then, multipoint flux approximation at branching point is necessary.
            //TODO how to avoid this type conversion from enum to int here
            if(int(GridView::dimension) != int(GridView::dimensionworld)) //only works for quad elements
            {
                if ((numScvf > 2 && !onBoundary) || (numScvf > 1 && onBoundary))
                {
                    for(unsigned int scvfIdx = 0; scvfIdx < numScvf; scvfIdx++)
                    {
                        subContVolFace[scvfIdx].fapIndices.reserve(numNeighbors); //allocate more memory
                        subContVolFace[scvfIdx].fapDistances.reserve(numNeighbors); //allocate for memory

                        for(unsigned int scvfIdxOther = 0; scvfIdxOther < numScvf; scvfIdxOther++)
                        {
                            if(scvfIdx == scvfIdxOther)
                                continue;

                            //if two intersections on the same facet (but different neighbor)
                            if(subContVolFace[scvfIdx].fIdx == subContVolFace[scvfIdxOther].fIdx)
                            {
                                subContVolFace[scvfIdx].numFap++; //multipoint flux approximation necessary
                                subContVolFace[scvfIdx].fapIndices.push_back(scvfIdxOther+1);
                            }
                        }

                        //do this only for sub control volume faces that have multiple neighbors (numFap > 2)
                        for (unsigned int fapIdx = 2; fapIdx < subContVolFace[scvfIdx].numFap; fapIdx++)
                        {
                            //calculate the distance from element midpoint to integration point for the addtional neighbors
                            GlobalPosition fapDistance = subContVolFace[scvfIdx].ipGlobal;
                            fapDistance -= neighbors[subContVolFace[scvfIdx].fapIndices[fapIdx]].geometry().center();
                            subContVolFace[scvfIdx].fapDistances.push_back(fapDistance.two_norm());
                        }
                    }
                }
            }
        }
    }
    int boundaryFaceIndex(const int fIdx, const int vIdxInFace) const
    {
        return fIdx;
    }
};

}

#endif
