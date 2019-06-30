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
#ifndef DUMUX_ROOTSYSTEM_SPATIAL_PARAMS_DGF_HH
#define DUMUX_ROOTSYSTEM_SPATIAL_PARAMS_DGF_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/math.hh>

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
class RootsystemSpatialParamsDGF
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

    enum { dimWorld = GridView::dimensionworld };
    enum { dim = GridView::dimension};

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::GridPtr<Grid> GridPointer;
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

    enum Root {
        surfaceIdx = 0,
        radiusIdx = 1,
        orderIdx = 2,
        branchIdx = 3,
        massIdx = 4,
        axialPermIdx = 5,
        radialPermIdx = 6
    };

public:
    RootsystemSpatialParamsDGF(const GridView &gridView)
        : gridView_(gridView)
    {
        Kx_ = tree.template get<Scalar>("RootSystem.SpatialParams.Kx");
        Kr_ = tree.template get<Scalar>("RootSystem.SpatialParams.Kr");

    }

    ~RootsystemSpatialParamsDGF()
    {}

    void setParams()
    {
        int numElems = gridView_.size(0);
        rootParameter_.resize(numElems);

        for (const auto& element : elements(gridView_))
        {

            int elemIdx = gridView_.indexSet().index(element);
            rootParameter_[elemIdx].resize(7);

            Element level0eIt(element);
            for(int levelIdx = element.level(); levelIdx != 0; levelIdx--)
                level0eIt = level0eIt.father();

            Scalar rootLength = element.geometry().volume();
            Scalar rootSurface = GridCreator::parameters(level0eIt)[2]/(1 << element.level());
            // root radius
            rootParameter_[elemIdx][Root::radiusIdx] = rootSurface / rootLength / 2.0 / M_PI;

            // root order
            rootParameter_[elemIdx][Root::orderIdx] =GridCreator::parameters(level0eIt)[0];
            // root branch  -> count from 0!!
            rootParameter_[elemIdx][Root::branchIdx] = GridCreator::parameters(level0eIt)[1] -1;
            // root surface
            rootParameter_[elemIdx][Root::surfaceIdx] = rootSurface;
            // root mass
            rootParameter_[elemIdx][Root::massIdx] = GridCreator::parameters(level0eIt)[3];

            if ((int)rootParameter_[elemIdx][Root::orderIdx] == 1){
                rootParameter_[elemIdx][Root::axialPermIdx] = Kx_; //Kx
                rootParameter_[elemIdx][Root::radialPermIdx] = Kr_; //Kr
            }
            else if  ((int)rootParameter_[elemIdx][Root::orderIdx] == 2){
                rootParameter_[elemIdx][Root::axialPermIdx] = Kx_;  //Kx
                rootParameter_[elemIdx][Root::radialPermIdx] = Kr_;  //Kr
            }
            else //order >= 3
                rootParameter_[elemIdx][Root::axialPermIdx] = Kx_; //Kx
                rootParameter_[elemIdx][Root::radialPermIdx] = Kr_; //Kr

        }
    }

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
    Scalar rootRadius(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][Root::radiusIdx];
    }

    Scalar radius(int idx) const
    {
        return rootParameter_[idx][Root::radiusIdx];
    }

    Scalar rootSurface(int eIdx) const
    {
        return rootParameter_[eIdx][Root::surfaceIdx];
    }

    Scalar rootSurface(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][Root::surfaceIdx];
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
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][Root::axialPermIdx];
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
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][Root::radialPermIdx];
    }

    Scalar Kr(int eIdx) const
    {
        return rootParameter_[eIdx][Root::radialPermIdx];
    }

    /*
      Doc me!
    */
    Scalar rootOrder(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][Root::orderIdx];
    }

    /*
      Doc me!
    */
    Scalar rootBranch(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][Root::branchIdx];
    }

    /*
      Doc me!
    */
    Scalar rootMass(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][Root::massIdx];
    }

private:

    std::vector< std::vector<double> > rootParameter_;
    const GridView gridView_;
    Scalar Kx_;
    Scalar Kr_;

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif

