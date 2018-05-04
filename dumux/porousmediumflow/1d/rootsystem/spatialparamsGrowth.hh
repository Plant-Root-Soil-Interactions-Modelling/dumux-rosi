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
#ifndef DUMUX_ROOTSYSTEM_SPATIAL_PARAMS_GROWTH_HH
#define DUMUX_ROOTSYSTEM_SPATIAL_PARAMS_GROWTH_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/math.hh>

//#include <dumux/implicit/common/implicitproperties.hh>

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
class RootsystemSpatialParamsGrowth
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    enum { dimWorld = GridView::dimensionworld };
    enum { dim = GridView::dimension};

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename Grid::template Codim<0>::Entity::EntitySeed ElementSeed;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef typename std::multimap<int, ElementSeed> BranchMap;

public:
    RootsystemSpatialParamsGrowth(const GridView &gridView)
        : gridView_(gridView)
    { }

    ~RootsystemSpatialParamsGrowth()
    {}

    void setParams()
    {
        int numElems = gridView_.size(0);
        rootParameter_.resize(numElems);

        ElementIterator eIt = gridView_.template begin<0>();
        const ElementIterator eItEnd = gridView_.template end<0>();
        for (; eIt != eItEnd; ++eIt){

            int elemIdx = gridView_.indexSet().index(*eIt);
            rootParameter_[elemIdx].resize(8);

            auto level0eIt = *eIt;
            for(int levelIdx = eIt->level(); levelIdx != 0; levelIdx--)
                level0eIt = level0eIt.father();

            double rootLength0eIt_ = level0eIt.geometry().volume();
            double rootSurface0eIt_ = GridCreator::parameters(level0eIt)[2];
            // root radius
            rootParameter_[elemIdx][1] = rootSurface0eIt_ / rootLength0eIt_ / 2.0 / 3.1415;

            // root order
            rootParameter_[elemIdx][2] = GridCreator::parameters(level0eIt)[0];
            // root branch  -> count from 0!!
            rootParameter_[elemIdx][3] = GridCreator::parameters(level0eIt)[1] -1;
            // root surface
            rootParameter_[elemIdx][0] = GridCreator::parameters(level0eIt)[2];
            // root mass
            rootParameter_[elemIdx][4] = GridCreator::parameters(level0eIt)[3];

            if ((int)rootParameter_[elemIdx][2] == 1){
                rootParameter_[elemIdx][5] =  5.0968e-9; //Kx
                rootParameter_[elemIdx][6] =  2.04e-7;   //Kr
            }
            else if  ((int)rootParameter_[elemIdx][2] == 2){
                rootParameter_[elemIdx][5] =  5.0968e-9; //Kx
                rootParameter_[elemIdx][6] =  2.04e-7;   //Kr
            }
            else if  ((int)rootParameter_[elemIdx][2] == 3){
                rootParameter_[elemIdx][5] = 1e6 * 0.3058e-13;  //Kx
                rootParameter_[elemIdx][6] = 1e6 * 0.102e-11;   //Kr
            }
            else {

            }
            //age
            rootParameter_[elemIdx][7] = 0.0;

        }
    }

    void createBranches()
    {
        branches_.clear();

        ElementIterator eIt = gridView_.template begin<0>();
        const ElementIterator eItEnd = gridView_.template end<0>();
        for (; eIt != eItEnd; ++eIt){

            int elemIdx = gridView_.indexSet().index(*eIt);
            int bIdx = (int)rootParameter_[elemIdx][3];

            ElementSeed eSeed = eIt->seed();
            branches_.insert(std::pair<int,ElementSeed>(bIdx,eSeed));
        }
    }

    void printBranches()
    {
        int numElems = gridView_.size(0);
        int numBranches = (int)rootParameter_[numElems-1][3];

        typedef typename BranchMap::iterator Iterator;
        for (int i = 0; i< numBranches; ++i){
            std::cout <<  "Branch no: " << i<< '\n' ;
            for (Iterator it=branches_.equal_range(i).first; it!=branches_.equal_range(i).second; ++it){
                auto element = GridCreator::grid().entity((*it).second);
                int elemIdx = gridView_.indexSet().index(element);
                std::cout << " => " << elemIdx;
            }
            std::cout << '\n';
         }
         std::cout<<std::endl;

    }

    void returnBranches() const
    {
         return branches_;
    }

    void insertBranch(int branchIdx, const Element &newElement)
    {

        ElementSeed eSeed = newElement.seed();
        branches_.insert(std::pair<int,ElementSeed>(branchIdx, eSeed));


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
        return rootParameter_[elemIdx][1];
    }

    Scalar rootSurface(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][0];
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
        return rootParameter_[elemIdx][5];
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
        return rootParameter_[elemIdx][6];
    }

    /*
      Doc me!
    */
    Scalar rootOrder(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][2];
    }

    /*
      Doc me!
    */
    Scalar rootBranch(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][3];
    }

    /*
      Doc me!
    */
    Scalar rootMass(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        int elemIdx = gridView_.indexSet().index(element);
        return rootParameter_[elemIdx][4];
    }

    void getRootParams(const Element &element, std::vector<double> &rootParams)
    {
        int elemIdx = gridView_.indexSet().index(element);
        rootParams = rootParameter_[elemIdx];
    }

    void insertRootParams(const Element &element, const std::vector<double> &rootParams)
    {
        int elemIdx = gridView_.indexSet().index(element);
        if (elemIdx > rootParameter_.size())
            rootParameter_.resize(elemIdx);

        std::vector< std::vector<double> >::iterator it = rootParameter_.begin();
        rootParameter_.insert(it+elemIdx, rootParams);
    }
private:

    std::vector< std::vector<double> > rootParameter_;
    BranchMap branches_;

    const GridView gridView_;

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
