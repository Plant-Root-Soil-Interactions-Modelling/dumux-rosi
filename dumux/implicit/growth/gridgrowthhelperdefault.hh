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
#ifndef DUMUX_GROWTHHELPER_DEFAULT_HH
#define DUMUX_GROWTHHELPER_DEFAULT_HH

#include <dune/grid/utility/persistentcontainer.hh>

/**
 * @file
 * @brief Default implementation of user data transfer helper
 */

namespace Dumux
{

template<class TypeTag>
class ImplicitGridGrowthHelperDefault
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    struct GrowthValues
    {
        PrimaryVariables priVars;
        std::size_t count;
        GrowthValues() : count(0.0) {}
    };

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::PersistentContainer<Grid, GrowthValues> PersistentContainer;

    const GridView gridView_;
    const Grid& grid_;
    PersistentContainer growthMap_;

public:
    //! Constructs an adaptive helper object
    /**
     * In addition to providing a storage object for cell-centered Methods, this class provides
     * mapping functionality to adapt the grid.
     */
    ImplicitGridGrowthHelperDefault(const GridView& gridView)
    : gridView_(gridView),
      grid_(gridView.grid()),
      growthMap_(grid_, dofCodim) {}

    /*!
     * Store primary variables
     * @param problem The current problem
     */
    void storePrimVars(Problem& problem)
    {
        // size the growthMap to the current size
        growthMap_.resize();

        if(isBox)
            DUNE_THROW(Dune::NotImplemented, "Growth for box discretization!");
        else
        {
            for (const auto& element : elements(problem.gridView()))
            {
                // get reference to map entry
                GrowthValues &growthVars = growthMap_[element];

                // put your value in the map
                const unsigned int eIdx = problem.elementMapper().index(element);
                growthVars.priVars = problem.model().curSol()[eIdx];
            }
        }
    }

    /*!
     * Reconstruct missing primary variables (where elements are created/deleted)
     * @param problem The current problem
     */
    void reconstructPrimVars(Problem& problem)
    {
        growthMap_.resize();

        if(isBox)
            DUNE_THROW(Dune::NotImplemented, "Growth for box discretization!");
        else
        {
            for (const auto& element : elements(problem.gridView()))
            {
                // old elements get their old variables assigned
                if(!element.isNew())
                {
                    GrowthValues &growthVars = growthMap_[element];
                    const unsigned int newEIdx = problem.elementMapper().index(element);
                    problem.model().curSol()[newEIdx] = growthVars.priVars;
                }
                // initialize new elements with some value
                else
                {
                    const unsigned int newEIdx = problem.elementMapper().index(element);
                    problem.model().curSol()[newEIdx] = PrimaryVariables(0.0);
                }

            }
        }

        // reset entries in restrictionmap
        growthMap_.resize(typename PersistentContainer::Value());
        growthMap_.shrinkToFit();
        growthMap_.fill(typename PersistentContainer::Value());
    }
};
}
#endif
