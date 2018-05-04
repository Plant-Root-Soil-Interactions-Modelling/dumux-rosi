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
#ifndef DUMUX_IMPLICIT_GRIDGROWTHINDICATORDEFAULT_HH
#define DUMUX_IMPLICIT_GRIDGROWTHINDICATORDEFAULT_HH

/**
 * @file
 * @brief  Class defining a default indicator for grid growth
 */
namespace Dumux
{
/*!\ingroup ImplicitGridGrowthIndicator
 * @brief  Class defining a default indicator for grid growth
 *
 * Default implementation
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class ImplicitGridGrowthIndicatorDefault
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum {dimWorld = GridView::dimensionworld };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     */
    void calculateIndicator()
    {}

    /*! \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool markedForGrowth(const Element& element)
    {
        return false;
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool markedForRemoval(const Element& element)
    {
        return false;
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool markedForMerging(const Element& element)
    {
        return false;
    }

    /*! \brief Initializes the adaption indicator class*/
    void init()
    {};

    /*! \brief returns growth facet index where new element should be added*/
    int getGrowthFacetIndex(const Element& element)
    {
        return -1;
    }

    /*! \brief returns new point to be connected to the grid*/
    GlobalPosition getNewPoint(const Element& element)
    {
        return GlobalPosition();
    }

    /*! \brief Constructs a GridAdaptionIndicator for initialization of an adaptive grid
     *
     * Default implementation
     *
     * \param problem The problem object
     * \param adaptionIndicator Indicator whether a be adapted
     */
    ImplicitGridGrowthIndicatorDefault(Problem& problem)
    {}
};
}

#endif
