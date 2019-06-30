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
 * \brief Base class for grid growth models.
 */
#ifndef DUMUX_IMPLICIT_GRIDGROWTH_HH
#define DUMUX_IMPLICIT_GRIDGROWTH_HH

#include "gridgrowthproperties.hh"
#include "gridgrowthhelper.hh"
#include <unordered_map>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

namespace Dumux
{

/*!\ingroup ImplicitGridGrowth
 * @brief Standard Module for grid growth models.
 *
 * This class is created by the problem class with the template
 * parameters <TypeTag, true> and provides basic functionality
 * for grid growth methods:
 *
 * A standard implementation growGrid() will prepare everything
 * to calculate the next primary variables vector on the new grid.
 */
template<class TypeTag, bool growingGrid>
class ImplicitGridGrowth
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem)  Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, GrowthIndicator) GrowthIndicator;
    typedef Dumux::GrowthHelper<TypeTag> GrowthHelper;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {dim = GridView::dimension};
    enum {dimworld = GridView::dimensionworld};
    typedef Dune::FieldVector<Scalar, dimworld> GlobalPosition;

public:
    /*!
     * Constructor
     * @param problem The problem
     */
    ImplicitGridGrowth (Problem& problem)
        : growthHelper_(problem.gridView()),
          problem_(problem),
          growthIndicator_(problem),
          grew_(0),
          removed_(0)
    {}

    /*!
     * @brief Initalization of grid growth
     *
     * Prepares the grid for simulation after the initialization of the
     * problem. The applied indicator is selectable via the property
     * GrowthInitializationIndicator
     */
    void init()
    {
        growthIndicator_.init();
        problem_.model().init(problem_);
    }

    /*!
     * @brief Standard method for a growth step
     *
     * uses a standard procedure for adaptivity:
     * 1) Determine the refinement indicator
     * 2) Store primary variables in a map
     * 3) Adapt the grid, adapt variables sizes, update mappers
     * 4) Reconstruct primary variables, regain secondary variables
     */
    void growGrid()
    {
        Dune::Timer watch;
        growGrid(growthIndicator_);
        std::cout << "Added " << grew_ << " new elements in " << watch.elapsed() << std::endl;
    }

    /*!
     * @brief Method to adapt the grid with individual indicator vector
     *
     * @param indicator The refinement indicator that is applied
     *
     */
    template<class Indicator>
    void growGrid(Indicator& indicator)
    {
        // reset internal counter for marked elements
        grew_ = removed_ = 0;


        /**** 1) insert elements according to element-based growth indicator ***/
         indicator.calculateIndicator();
         insertElements(indicator);

        // abort if nothing in grid is marked
        if (grew_ == 0 && removed_ == 0)
            return;

        /****  2) Put primary variables in a map ***/
        growthHelper_.storePrimVars(problem_);

        /****  3) Grow Grid and size of variable vectors ***/
        problem_.grid().preGrow();
        bool insertedElements = problem_.grid().grow();
        if (!insertedElements) std::cout << "None of the elements could be inserted!" << std::endl;

        // update mapper to new leaf index set
        problem_.elementMapper().update();
        problem_.vertexMapper().update();

        // adapt size of vectors
        problem_.model().adaptVariableSize();

        /****  4) (Re-)construct primary variables to new grid ***/
        growthHelper_.reconstructPrimVars(problem_);

        // delete markers in grid
        problem_.grid().postGrow();
    }

    /*!
     * Mark Elements for grid refinement according to applied Indicator
     * @param indicator The refinement indicator that is applied
     * TODO This has to be more general with the possibilities
     * of the new grid growth interface in FoamGrid
     */
    template<class Indicator>
    void insertElements(const Indicator& indicator)
    {
        for (const auto& element : elements(problem_.gridView()))
            if (indicator.markedForGrowth(element))
            {
                unsigned int vIdx0 = problem_.grid().insertVertex(indicator.getNewPoint(element));
                unsigned int vIdx1 = problem_.vertexMapper().index(element.template subEntity<dim>(indicator.getGrowthFacetIndex(element)));
                std::vector<unsigned int> vertices = {vIdx1, vIdx0};
                problem_.grid().insertElement(element.type(), vertices);
                ++grew_;
            }
    }

    /*!
     * @brief Returns true if grid cells have been marked for growth
     */
    bool wasGrown()
    {
        return (grew_ != 0 || removed_ != 0);
    }

    const bool wasGrown() const
    {
        return (grew_ != 0 || removed_ != 0);
    }

    //! The growth indicator
    GrowthIndicator& growthIndicator()
    {
        return growthIndicator_;
    }

    //! The constant growth indicator
    GrowthIndicator& growthIndicator() const
    {
        return growthIndicator_;
    }

private:
    GrowthHelper growthHelper_;
    Problem& problem_;
    GrowthIndicator growthIndicator_;
    int grew_;
    int removed_;
    int growthInterval_;
};

/*!
 * @brief Class for simulations without grid growth
 *
 * This class provides empty methods for simulations without grid growth
 * for compilation reasons. If grid growth is desired, create the
 * class with template arguments <TypeTag, true> instead.
 */
template<class TypeTag>
class ImplicitGridGrowth<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:
    ImplicitGridGrowth (Problem& problem) {}
    void init() {}
    void growGrid() {}
    bool wasGrown() { return false; }
};

} // end namespace Dumux
#endif
