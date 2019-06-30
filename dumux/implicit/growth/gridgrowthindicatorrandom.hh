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
#ifndef DUMUX_GRIDGROWTH_INDICATOR_RANDOM_HH
#define DUMUX_GRIDGROWTH_INDICATOR_RANDOM_HH

#include <ctime>
#include <math.h>
#include <stdlib.h> // srand, rand
/**
 * @file
 * @brief  Class defining a boundary indicator for grid growth
 */
namespace Dumux
{
/*!\ingroup GridGrowthIndicatorRandom
 * @brief  Class defining a boundary indicator for grid growth calculating semi-random root growth
 *
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridGrowthIndicatorRandom
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum {dim = GridView::dimension };
    enum {dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    /*! \brief Constructs a GridAdaptionIndicator for initialization of an adaptive grid
     *
     * Default implementation
     *
     * \param problem The problem object
     * \param adaptionIndicator Indicator whether a be adapted
     */
    GridGrowthIndicatorRandom(Problem& problem) : problem_(problem)
    {  }

     /*! \brief Initializes the growth indicator class*/
    void init() {}


    void calculateIndicator()
    {
        markedForGrowthElments_.clear();
        growthPoints_.clear();
        facetIndex_.clear();

        srand (time(NULL));
       // vector with all element indices and mark if they shall grow
       markedForGrowthElments_.resize(problem_.gridView().size(0));

       for (const auto& element : elements(problem_.gridView()))
       {
           bool grow = false;
           int eIdx = problem_.elementMapper().index(element);

           std::vector<double > spatialParams;
           spatialParams.resize(8);
           problem_.spatialParams().getRootParams(element, spatialParams);
           const Scalar time = problem_.timeManager().time() +  problem_.timeManager().timeStepSize();

           IntersectionIterator isEndIt = problem_.gridView().iend(element);
           for (const auto& intersection : intersections(problem_.gridView(), element))
           {
               const Scalar age = spatialParams[7];
               double pos = intersection.geometry().center()[2];
               if(!( pos + 1e-6 > 0 ) && ((int)(time - age/2.0)% 3 == 0))
                   grow = true;
            }
           markedForGrowthElments_[eIdx] = grow;

           // if elelemnt is marked for growth -> calculate new growth point
           if (grow){

               Scalar length = element.geometry().volume();
               Scalar t = problem_.timeManager().time() +  problem_.timeManager().timeStepSize();
               GlobalPosition growthPoint = element.geometry().center();

               int tIdex = problem_.timeManager().timeStepIndex(); // time index

               GlobalPosition lengthCoords(0.0);

               for(auto&& intersection : intersections(problem_.gridView(), element))
                   lengthCoords -= intersection.geometry().center();

               for (int i = 0; i<lengthCoords.size(); ++i)
                   lengthCoords[i] *= -1.0;

               for(auto&& intersection : intersections(problem_.gridView(), element)) {
                   int intIndex = intersection.indexInInside();

                   if ( intIndex == 1) { // facet index == 1
                       if(!intersection.neighbor() ) // at inner elements
                       {
                           growthPoint = intersection.geometry().center();

                           int ran0 = rand() % 20;
                           double ran = (double)ran0 / 10.0;
                           Scalar z = t/2.0/3.1415*(-1.0);
                           z = ran*3.1415;

                           if ((growthPoint[0] == 0.0) && (growthPoint[1] == 0.0))
                               growthPoint[2] -= 1.01*length;
                           else
                          {
                              growthPoint[0] += lengthCoords[0]*0.1;
                              growthPoint[1] += lengthCoords[1]*0.1;
                              growthPoint[2] += lengthCoords[2]*0.01;
                          }
                       }
                       else if (tIdex%2 == 0 ) // at boundary elements
                       {
                           growthPoint = intersection.geometry().center();
                           int ran0 = rand() % 20;
                           double ran = (double)ran0 / 10.0;
                           Scalar z = t/2.0/3.1415*(-1.0);
                           z = ran*3.1415;
                           growthPoint[0] -= cos(z)*.00005;
                           growthPoint[1] -= sin(z)*.00005;
                           growthPoint[2] -= z*0.000005;
                       }
                   }
               }
               //store  growthPoint;
               growthPoints_.insert({eIdx, growthPoint});
               // store facet index at old element for the new growth point (here always 1)
               facetIndex_.insert({eIdx, 1});
           }
        }
    }

    /*! \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool markedForGrowth(const Element& element) const
    {
        int eIdx = problem_.elementMapper().index(element);
        return markedForGrowthElments_[eIdx];
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool markedForRemoval(const Element& element) const
    {
        return false;
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool markedForMerging(const Element& element) const
    {
        return false;
    }

    /*! \brief returns growth facet index where new element should be added*/
    int getGrowthFacetIndex(const Element& element) const
    {
        int eIdx = problem_.elementMapper().index(element);
        return facetIndex_.at(eIdx);
    }

    /*! \brief returns new point to be connected to the grid*/
    GlobalPosition getNewPoint(const Element& element) const
    {
        int eIdx = problem_.elementMapper().index(element);
        return growthPoints_.at(eIdx);
    }

protected:
    Problem& problem_;
    std::vector<bool> markedForGrowthElments_;
    std::map<int, GlobalPosition> growthPoints_;
    std::map<int, int> facetIndex_;

};
}

#endif
