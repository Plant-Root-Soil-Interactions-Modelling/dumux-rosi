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
 * \brief A test problem for the rootsystem model:
 */
#ifndef DUMUX_ROOTSYSTEMTEST_PROBLEM_HH
#define DUMUX_ROOTSYSTEMTEST_PROBLEM_HH

#include <dumux/porousmediumflow/1d/rootsystem/model.hh>
#include <dumux/porousmediumflow/1d/rootsystem/problem.hh>
#include <dumux/linear/seqsolverbackend.hh>

// fluidsystem
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
// growth
#include <dumux/implicit/growth/gridgrowthindicatorrandom.hh>
// spatial params
#include "rootsystemtestspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class RootsystemTestProblem;

namespace Properties
{
NEW_TYPE_TAG(RootsystemTestProblem, INHERITS_FROM(Rootsystem));
NEW_TYPE_TAG(RootsystemTestBoxProblem, INHERITS_FROM(BoxModel, RootsystemTestProblem));
NEW_TYPE_TAG(RootsystemTestCCProblem, INHERITS_FROM(CCModel, RootsystemTestProblem));

SET_PROP(RootsystemTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(RootsystemTestProblem, Grid,Dune::FoamGrid<1, 3>);

// Set the problem property
SET_TYPE_PROP(RootsystemTestProblem, Problem, Dumux::RootsystemTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(RootsystemTestProblem, SpatialParams, Dumux::RootsystemTestSpatialParams<TypeTag>);

// Disable velocity output gravity
SET_BOOL_PROP(RootsystemTestProblem, VtkAddVelocity, false);

// Enable gravity
SET_BOOL_PROP(RootsystemTestProblem, ProblemEnableGravity, false);

// growth properties
SET_BOOL_PROP(RootsystemTestProblem, GrowingGrid, true);
SET_TYPE_PROP(RootsystemTestProblem, GrowthIndicator, GridGrowthIndicatorRandom<TypeTag>);
}

/*!
 * \ingroup RootsystemBoxModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * In the middle of the domain, a lens with low permeability (\f$K=10e-12\f$)
 * compared to the surrounding material (\f$ K=10e-10\f$) is defined.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxrootsystem -parameterFile test_boxrootsystem.input</tt> or
 * <tt>./test_ccrootsystem -parameterFile test_ccrootsystem.input</tt>
 *
 * The same parameter file can be also used for 3d simulation but you need to change line
 * <tt>typedef Dune::SGrid<2,2> type;</tt> to
 * <tt>typedef Dune::SGrid<3,3> type;</tt> in the problem file
 * and use <tt>1p_3d.dgf</tt> in the parameter file.
 */
template <class TypeTag>
class RootsystemTestProblem : public RootsystemProblem<TypeTag>
{
    typedef RootsystemProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pIdx = Indices::pIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    RootsystemTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        this->spatialParams().setParams();
        this->spatialParams().createBranches();
    }

    /*!
     * \name Problem parameters
     */
    //
    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return name_.c_str();
    }

    /*!
     * \brief Return the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    // \}
    /*!
     * \name Boundary conditions
     */
    void boundaryTypes(BoundaryTypes &values,
                        const Intersection &intersection) const
    {
        double pos = intersection.geometry().center()[2];
        if ( pos +eps_ > 0 )
            values.setAllDirichlet();
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values[pIdx] = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.CriticalCollarPressure);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.*/
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        values[conti0EqIdx] = 0.0;
    }


    void solDependentSource(PrimaryVariables &values,
                     const Element &element,
                     const FVElementGeometry &fvGeometry,
                     const int scvIdx,
                     const ElementVolumeVariables &elemVolVars) const
    {
        const SpatialParams &spatialParams = this->spatialParams();
        Scalar Kr = spatialParams.Kr(element, fvGeometry, scvIdx);
        Scalar rootSurface = spatialParams.rootSurface(element, fvGeometry, scvIdx);

        Scalar phx = elemVolVars[0].pressure();
        Scalar phs = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.SoilPressure);
        values =  Kr * rootSurface * ( phs - phx) / element.geometry().volume();
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        priVars[pIdx] = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.CriticalCollarPressure);
    }

private:
    std::string name_;
    const double eps_ = 1e-8;

};
} //end namespace

#endif
