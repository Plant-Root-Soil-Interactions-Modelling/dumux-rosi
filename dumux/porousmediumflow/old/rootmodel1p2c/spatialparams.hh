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
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of rootsystem problems
 * using a fully implicit discretization method.
 */
#ifndef DUMUX_ROOTSYSTEM_SPATIAL_PARAMS_HH
#define DUMUX_ROOTSYSTEM_SPATIAL_PARAMS_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/math.hh>

//#include <dumux/implicit/properties.hh>

#include <dune/common/fmatrix.hh>

namespace Dumux {
// forward declaration of property tags
namespace Properties {
NEW_PROP_TAG(SpatialParams);
}

/*!
 * \ingroup SpatialParameters
 */


/**
 * \brief The base class for spatial parameters of one-phase problems
 * using a fully implicit discretization method.
 */
template<class TypeTag>
class RootsystemSpatialParams
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) Implementation;

    enum { dimWorld = GridView::dimensionworld };
    enum { dim = GridView::dimension};

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView:: template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

public:
    RootsystemSpatialParams(const GridView &gridView)
        : gridView_(gridView)
    { }

    ~RootsystemSpatialParams()
    {}

    /*!
     * \brief Averages the intrinsic permeability (Scalar).
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first node
     * \param K2 intrinsic permeability of the second node
     */
    void meanK(DimMatrix &result,
               Scalar K1,
               Scalar K2) const
    {
        const Scalar K = Dumux::harmonicMean(K1, K2);
        for (int i = 0; i < dimWorld; ++i) {
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = 0;
            result[i][i] = K;
        }
    }

    /*!
     * \brief Averages the intrinsic permeability (Tensor).
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first node
     * \param K2 intrinsic permeability of the second node
     */
    void meanK(DimMatrix &result,
               const DimMatrix &K1,
               const DimMatrix &K2) const
    {
        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = harmonicMean(K1[i][j], K2[i][j]);
    }
    /*!
     * \brief Function for defining the root radius.
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return rootRadius
     */
    //Scalar rootRadius(const Element &element,
    //        const FVElementGeometry &fvGeometry,
    //        int scvIdx) const
    //{
    //    return GET_RUNTIME_PARAM(TypeTag,
    //                              Scalar,
    //                              SpatialParams.rootRadius);
    //}

    Scalar rootSurface(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        return  GET_RUNTIME_PARAM(TypeTag,
                                  Scalar,
                                  SpatialParams.rootSurface);
    }

    /*!
     * \brief Function for defining the axial conductance (K_x).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return axial conductance (K_x)
     */
    Scalar Kx(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        return GET_RUNTIME_PARAM(TypeTag,
                                 Scalar,
                                 SpatialParams.Kx);
    }

    /*!
     * \brief Function for defining the radial conductanctivity (K_r^*).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return radial conductanctivity (K_r^*)
     */
    Scalar Kr(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        return  GET_RUNTIME_PARAM(TypeTag,
                                  Scalar,
                                  SpatialParams.Kr) ;
    }

    /*
      Doc me!
    */
    Scalar rootOrder(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        return 1;
    }

    /*
      Doc me!
    */
    Scalar rootBranch(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        return 1;
    }

    /*
      Doc me!
    */
    Scalar rootMass(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        return 0.00;
    }

    std::vector<double> returnRootParams(const Element &element)
    {
        std::vector<double> rootParameter;
        rootParameter.resize(7);
        rootParameter[0] = GET_RUNTIME_PARAM(TypeTag,
                                              Scalar,
                                              SpatialParams.rootSurface);; // surface
        rootParameter[1] = element.geometry().volume(); // length
        rootParameter[2] = 1; // order
        rootParameter[3] = 1; // branch
        rootParameter[4] = 0.00; // mass
        rootParameter[5] = GET_RUNTIME_PARAM(TypeTag,
                                              Scalar,
                                              SpatialParams.Kx); // Kx
        rootParameter[6] = GET_RUNTIME_PARAM(TypeTag,
                                              Scalar,
                                              SpatialParams.Kr); // Kr

        return rootParameter;
    }

    void insertRootParams(const Element &element, std::vector<double> rootParams)
    {
        int elemIdx = gridView_.indexSet().index(element);
        if (elemIdx > rootParameter_.size())
            rootParameter_.resize(elemIdx);

        std::vector< std::vector<double> >::iterator it = rootParameter_.begin();
        rootParameter_.insert(it+elemIdx, rootParams);
    }

private:

    const GridView gridView_;
    std::vector< std::vector<double> > rootParameter_;

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
