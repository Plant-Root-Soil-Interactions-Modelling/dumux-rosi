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
 * \brief A test problem the rootsystem coupled with richards in the bulk
 */
#ifndef DUMUX_ROSI_1P2C_BUFFER_TEST_PROBLEM_HH
#define DUMUX_ROSI_1P2C_BUFFER_TEST_PROBLEM_HH

#include "rootsystem1p2ctestproblem.hh"
#include "soilRichards2cbuffertestproblem.hh"

#include <dumux/multidimension/problem.hh>
//#include <dumux/multidimension/embeddedcoupling/cellcentered/couplingmanager.hh>
#include <dumux/multidimension/embeddedcoupling/cellcentered/bboxtreecouplingmanager.hh>
#include <dumux/multidimension/embeddedcoupling/integrationpointsource.hh>
#include <dumux/multidimension/linear/directsolverbackend.hh>

namespace Dumux
{
template <class TypeTag>
class RosiRichardsTwoCBufferTestProblem;

namespace Properties
{
NEW_TYPE_TAG(RosiRichardsTwoCBufferTestProblem, INHERITS_FROM(MultiDimension));

// Set the problem property
SET_TYPE_PROP(RosiRichardsTwoCBufferTestProblem, Problem, Dumux::RosiRichardsTwoCBufferTestProblem<TypeTag>);

// Set the coupling manager
SET_TYPE_PROP(RosiRichardsTwoCBufferTestProblem, CouplingManager, Dumux::CCBBoxTreeEmbeddedCouplingManager<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(RosiRichardsTwoCBufferTestProblem, LowDimProblemTypeTag, TTAG(RootsystemOnePTwoCTestCCProblem));
SET_TYPE_PROP(RosiRichardsTwoCBufferTestProblem, BulkProblemTypeTag, TTAG(SoilRichardsTwoCBufferTestCCProblem));

// publish this problem in the sub problems
SET_TYPE_PROP(RootsystemOnePTwoCTestCCProblem, GlobalProblemTypeTag, TTAG(RosiRichardsTwoCBufferTestProblem));
SET_TYPE_PROP(SoilRichardsTwoCBufferTestCCProblem, GlobalProblemTypeTag, TTAG(RosiRichardsTwoCBufferTestProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(RootsystemOnePTwoCTestCCProblem, ParameterTree) : GET_PROP(TTAG(RosiRichardsTwoCBufferTestProblem), ParameterTree) {};
SET_PROP(SoilRichardsTwoCBufferTestCCProblem, ParameterTree) : GET_PROP(TTAG(RosiRichardsTwoCBufferTestProblem), ParameterTree) {};

// Set the point source type of the subproblems to an integration point source
SET_TYPE_PROP(RootsystemOnePTwoCTestCCProblem, PointSource, Dumux::IntegrationPointSource<TTAG(RootsystemOnePTwoCTestCCProblem), unsigned int>);
SET_TYPE_PROP(RootsystemOnePTwoCTestCCProblem, PointSourceHelper, Dumux::IntegrationPointSourceHelper<TTAG(RootsystemOnePTwoCTestCCProblem)>);
SET_TYPE_PROP(SoilRichardsTwoCBufferTestCCProblem, PointSource, Dumux::IntegrationPointSource<TTAG(SoilRichardsTwoCBufferTestCCProblem), unsigned int>);
SET_TYPE_PROP(SoilRichardsTwoCBufferTestCCProblem, PointSourceHelper, Dumux::IntegrationPointSourceHelper<TTAG(SoilRichardsTwoCBufferTestCCProblem)>);

//#if HAVE_UMFPACK
//SET_TYPE_PROP(RosiRichardsTwoCBufferTestProblem, LinearSolver, Dumux::MultiDimensionUMFPackBackend<TypeTag>);
//SET_TYPE_PROP(RosiRichardsTwoCBufferTestProblem, LinearSolver, ILU0BiCGSTABBackend<TypeTag>);
SET_TYPE_PROP(RosiRichardsTwoCBufferTestProblem, LinearSolver, UMFPackBackend<TypeTag>);
//#endif

}//end namespace properties

template <class TypeTag>
class RosiRichardsTwoCBufferTestProblem : public MultiDimensionProblem<TypeTag>
{
    typedef MultiDimensionProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // obtain the type tags of the sub problems
    typedef typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag) BulkProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag) LowDimProblemTypeTag;

    // obtain types from the sub problem type tags
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, Problem) BulkProblem;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem) LowDimProblem;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, GridView) BulkGridView;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView) LowDimGridView;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingManager) CouplingManager;

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

public:
    RosiRichardsTwoCBufferTestProblem(TimeManager &timeManager, const BulkGridView &bulkGridView, const LowDimGridView &lowDimgridView)
    : ParentType(timeManager, bulkGridView, lowDimgridView)
    {}
    void init()
    {
        ParentType::init();
        ParentType::couplingManager().initializeCouplingSources();

    }
    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behavior is to write out the solution for
     * every time step. This function is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteOutput() const
    {
        if (this->timeManager().time()<0) return true;
        if (tree.template get<bool>("TimeManager.OutputEveryTimeStep", false))
            return true;
        else
            return (this->timeManager().episodeWillBeFinished() || this->timeManager().willBeFinished());
    }
};

} //end namespace

#endif
