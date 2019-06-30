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
#ifndef DUMUX_GROWTHHELPER_HH
#define DUMUX_GROWTHHELPER_HH

#include <dune/grid/utility/persistentcontainer.hh>

/**
 * @file
 * @brief  Base class holding the variables for implicit models.
 */

namespace Dumux
{

template<class TypeTag>
class GrowthHelper
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    struct GrowthValues
    {
        PrimaryVariables u;
        std::vector<double> spatialVars;
        int count;
        GrowthValues()
        {
            count = 0;
            spatialVars.resize(8);
        }
    };

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GridView::Grid Grid;
    typedef typename Grid::LevelGridView LevelGridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::PersistentContainer<Grid, GrowthValues> PersistentContainer;

    const GridView gridView_;
    const Grid& grid_;
    PersistentContainer growthMap_;
    PrimaryVariables priVars_;

public:
    //! Constructs an adaptive helper object
    /**
     * In addition to providing a storage object for cell-centered Methods, this class provides
     * mapping functionality to adapt the grid.
     *
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */
    GrowthHelper(const GridView& gridView) :
        gridView_(gridView), grid_(gridView.grid()), growthMap_(grid_, dofCodim)
    {}


    /*!
     * Store primary variables
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     * From upper level on downwards, the old solution is stored into an container
     * object, before the grid is adapted. Father elements hold averaged information
     * from the son cells for the case of the sons being coarsened.
     *
     * @param problem The current problem
     */
    void storePrimVars(Problem& problem)
    {
        growthMap_.resize();

        // loop over all levels of the grid
        for (int level = grid_.maxLevel(); level >= 0; level--)
        {
            //get grid view on level grid
            const auto& levelView = grid_.levelGridView(level);

            if(!isBox)
            {
                for (const auto& element : elements(levelView))
                {
                    //get your map entry
                    GrowthValues &growthVars = growthMap_[element];

                    // put your value in the map
                    if (element.isLeaf())
                    {
                        // get index
                        int eIdx = problem.elementMapper().index(element);
                        std::vector<double> spatialParams(8,0);
                        problem.spatialParams().getRootParams(element, spatialParams);

                        storeGrowthValues(growthVars, problem.model().curSol()[eIdx], spatialParams);
                    }
                }
            }
            else
            {
                DUNE_THROW(Dune::NotImplemented, "Growth for box discretization!");
            }
        }
    }

    /*!
     * Reconstruct missing primary variables (where elements are created/deleted)
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     * Starting from the lowest level, the old solution is mapped on the new grid:
     * Where coarsened, new cells get information from old father element.
     * Where refined, a new solution is reconstructed from the old father cell,
     * and then a new son is created. That is then stored into the general data
     * structure (CellData).
     *
     * @param problem The current problem
     */
    void reconstructPrimVars(Problem& problem)
    {
        growthMap_.resize();

        // map the old indices to the new ones for old entites
        for (int level = 0; level <= grid_.maxLevel(); level++)
        {
            const auto& levelView = grid_.levelGridView(level);
            for (const auto& element : elements(levelView))
            {
                if (!element.isNew() && element.isLeaf())
                {
                    if(!isBox)
                    {
                        GrowthValues &growthVars = growthMap_[element];
                        int newEIdx = problem.elementMapper().index(element);
                        std::vector<double> spatialParams(8,0);
                        setGrowthValues(growthVars, problem.model().curSol()[newEIdx], spatialParams);
                        problem.spatialParams().insertRootParams(element, spatialParams);
                    }
                    else
                    {
                        DUNE_THROW(Dune::NotImplemented, "Growth for box discretization!");
                    }
                }
                else if(element.isNew() && element.isLeaf())
                {
                    if(!isBox)
                    {
                        // new element are always on the leaf?
                        // get the index of the new element
                        int newEIdx = problem.elementMapper().index(element);

                        PrimaryVariables priVars = PrimaryVariables(0.0);
                        std::vector<double> spatialParams(8, 0.0);

                        // the new privars are averages over the neighbouring elements' privars
                        int count = 0;
                        for(const auto& intersection : intersections(gridView_, element))
                        {
                            if (intersection.neighbor())
                            {
                                GrowthValues &growthVars = growthMap_[intersection.outside()];
                                priVars += growthVars.u;
                                for (unsigned int i = 0; i < spatialParams.size(); i++)
                                {
                                    if (i == 2) //order
                                        spatialParams[i] += 2;
                                    else
                                        spatialParams[i] += growthVars.spatialVars[i];
                                }
                                count++;
                            }
                        }
                        priVars /= count;

                        for (unsigned int i = 0; i < spatialParams.size(); i++)
                            spatialParams[i] /= count;

                        spatialParams[7] = problem.timeManager().time() + problem.timeManager().timeStepSize();
                        // set current solution for new element in solution vector
                        problem.model().curSol()[newEIdx] = priVars;

                        // set new spatial params for new root element
                        problem.spatialParams().insertRootParams(element, spatialParams);

                        int branchIdx = (int)spatialParams[3];
                        problem.spatialParams().insertBranch(branchIdx, element);
                    }
                    else
                    {
                        DUNE_THROW(Dune::NotImplemented, "Growth for box discretization!");
                    }
                }
            }
        }
        // reset entries in restrictionmap
        growthMap_.resize( typename PersistentContainer::Value() );
        growthMap_.shrinkToFit();
        growthMap_.fill( typename PersistentContainer::Value() );
    }

    //! Stores values to be adapted in an adaptedValues container
    /**
     * Stores values to be adapted from the current primary variables objects into
     * the persistent container in order to be mapped on a new grid.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element to be stored
     */
    static void storeGrowthValues(GrowthValues& growthVars, const PrimaryVariables& u, const std::vector<double > & spatialVars)
    {
        growthVars.u = u;
        growthVars.spatialVars = spatialVars;
    }

    //! Set adapted values in CellData
    /**
     * This methods stores reconstructed values into the cellData object, by
     * this setting a newly mapped solution to the storage container of the
     * decoupled models.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element where things are stored.
     */
    static void setGrowthValues(GrowthValues& growthVars, PrimaryVariables& u, std::vector<double > & spatialVars)
    {
        u = growthVars.u;
        spatialVars = growthVars.spatialVars;
    }
};
}
#endif
