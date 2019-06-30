// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup MultiDimension
 * \brief Reference implementation of a controller class for the Newton solver.
 *
 * Usually this controller should be sufficient.
 */
#ifndef DUMUX_MULTIDIMENSION_NEWTON_CONTROLLER_HH
#define DUMUX_MULTIDIMENSION_NEWTON_CONTROLLER_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include "linearsolveracceptsmultitypematrix.hh"

namespace Dumux
{
template <class TypeTag>
class MultiDimensionNewtonController;

namespace Properties
{
//! Specifies the implementation of the Newton controller
NEW_PROP_TAG(NewtonController);

//! Specifies the type of the actual Newton method
NEW_PROP_TAG(NewtonMethod);

//! Specifies the type of a solution
NEW_PROP_TAG(SolutionVector);

//! Specifies the type of a global Jacobian matrix
NEW_PROP_TAG(JacobianMatrix);

//! specifies the type of the time manager
NEW_PROP_TAG(TimeManager);

/*!
 * \brief Specifies whether the update should be done using the line search
 *        method instead of the plain Newton method.
 *
 * Whether this property has any effect depends on whether the line
 * search method is implemented for the actual model's Newton
 * controller's update() method. By default line search is not used.
 */
NEW_PROP_TAG(NewtonUseLineSearch);

//! indicate whether the shift criterion should be used
NEW_PROP_TAG(NewtonEnableShiftCriterion);

//! the value for the maximum relative shift below which convergence is declared
NEW_PROP_TAG(NewtonMaxRelativeShift);

//! indicate whether the residual criterion should be used
NEW_PROP_TAG(NewtonEnableResidualCriterion);

//! the value for the residual reduction below which convergence is declared
NEW_PROP_TAG(NewtonResidualReduction);

//! indicate whether both of the criteria should be satisfied to declare convergence
NEW_PROP_TAG(NewtonSatisfyResidualAndShiftCriterion);

/*!
 * \brief The number of iterations at which the Newton method
 *        should aim at.
 *
 * This is used to control the time-step size. The heuristic used
 * is to scale the last time-step size by the deviation of the
 * number of iterations used from the target steps.
 */
NEW_PROP_TAG(NewtonTargetSteps);

//! Number of maximum iterations for the Newton method.
NEW_PROP_TAG(NewtonMaxSteps);

} // end namespace Properties

/*!
 * \ingroup MultiDimension
 * \brief A reference implementation of a Newton controller specific
 *        for the box scheme.
 *
 * If you want to specialize only some methods but are happy with the
 * defaults of the reference controller, derive your controller from
 * this class and simply overload the required methods.
 */
template <class TypeTag>
class MultiDimensionNewtonController
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) LinearSolver;
    typedef typename GET_PROP(TypeTag, SubProblemBlockIndices) SubProblemBlockIndices;

    typename SubProblemBlockIndices::BulkIdx bulkIdx;
    typename SubProblemBlockIndices::LowDimIdx lowDimIdx;

    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, Problem) BulkProblem;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem) LowDimProblem;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, SolutionVector) BulkSolutionVector;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, SolutionVector) LowDimSolutionVector;

    enum {
        numEqBulk = GET_PROP_VALUE(BulkProblemTypeTag, NumEq),
        numEqLowDim = GET_PROP_VALUE(LowDimProblemTypeTag, NumEq)
    };
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

public:
    /*!
     * \brief Constructor
     */
    MultiDimensionNewtonController(const Problem &problem)
    : endIterMsgStream_(std::ostringstream::out), linearSolver_(problem)
    {
        useLineSearch_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, UseLineSearch);
        enableShiftCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableShiftCriterion);
        enableResidualCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableResidualCriterion);
        satisfyResidualAndShiftCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, SatisfyResidualAndShiftCriterion);
        if (!enableShiftCriterion_ && !enableResidualCriterion_)
        {
            DUNE_THROW(Dune::NotImplemented,
                       "at least one of NewtonEnableShiftCriterion or "
                       << "NewtonEnableResidualCriterion has to be set to true");
        }

        setMaxRelativeShift(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, MaxRelativeShift));
        setResidualReduction(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, ResidualReduction));
        setTargetSteps(GET_PARAM_FROM_GROUP(TypeTag, int, Newton, TargetSteps));
        setMaxSteps(GET_PARAM_FROM_GROUP(TypeTag, int, Newton, MaxSteps));
        avoidNegativeValueMacroScale_ = tree.template get<bool>("Newton.AvoidNegativeValueMacroScale", true);
        fixedNewtonIterationSteps_ = tree.template get<Scalar>("Newton.FixedNewtonIterationSteps", 0);
        decouplingEquations_ = tree.template get<bool>("MultiDimension.DecouplingEquations", true);
        decomposeMatrixSolver_ = tree.template get<bool>("MultiDimension.DecomposeMatrixSolver", true);

        ratioSuggestedTimeStepSize_ = tree.template get<Scalar>("TimeManager.RatioSuggestedTimeStepSize", 0);
        verbose_ = true;
        numSteps_ = 0;
    }

    /*!
     * \brief Set the maximum acceptable difference of any primary variable
     * between two iterations for declaring convergence.
     *
     * \param tolerance The maximum relative shift between two Newton
     *                  iterations at which the scheme is considered finished
     */
    void setMaxRelativeShift(Scalar tolerance)
    { shiftTolerance_ = tolerance; }

    /*!
     * \brief Set the maximum acceptable residual norm reduction.
     *
     * \param tolerance The maximum reduction of the residual norm
     *                  at which the scheme is considered finished
     */
    void setResidualReduction(Scalar tolerance)
    { reductionTolerance_ = tolerance; }

    /*!
     * \brief Set the number of iterations at which the Newton method
     *        should aim at.
     *
     * This is used to control the time-step size. The heuristic used
     * is to scale the last time-step size by the deviation of the
     * number of iterations used from the target steps.
     *
     * \param targetSteps Number of iterations which are considered "optimal"
     */
    void setTargetSteps(int targetSteps)
    { targetSteps_ = targetSteps; }

    /*!
     * \brief Set the number of iterations after which the Newton
     *        method gives up.
     *
     * \param maxSteps Number of iterations after we give up
     */
    void setMaxSteps(int maxSteps)
    { maxSteps_ = maxSteps; }

    /*!
     * \brief Returns true if another iteration should be done.
     *
     * \param uCurrentIter The solution of the current Newton iteration
     */
    bool newtonProceed(const SolutionVector &uCurrentIter)
    {
    	if (avoidNegativeValueMacroScale_)
    	{
    		haveNegValue_ = false;
    		for (const auto sol : uCurrentIter[bulkIdx])
    		    if (sol[1] < 0)
    		        {
    		            haveNegValue_ = true;
                        std::cout<<"Negative value in soil: "<<sol[1]<<"\n";
    		            if (numSteps_ >= maxSteps_)
    		                return false;
    		            else return true;
    		        }
    		//for (const auto sol : uCurrentIter[lowDimIdx])
    		//    if (sol[1] < -1e-25)
    		//        {
            //            std::cout<<"Negative value in root: "<<sol[1]<<"\n";
    		//            haveNegValue_ = true;
    		//            if (numSteps_ >= maxSteps_)
    		//                return false;
    		//            else return true;
    		//        }
    	}
        if (numSteps_ < 2)
            return true; // we always do at least two iterations
        else if (asImp_().newtonConverged()) {
            return false; // we are below the desired tolerance
        }
        else if (numSteps_ >= maxSteps_) {
            // We have exceeded the allowed number of steps. If the
            // maximum relative shift was reduced by a factor of at least 4,
            // we proceed even if we are above the maximum number of steps.
            if (enableShiftCriterion_)
                return shift_*4.0 < lastShift_;
            else
                return reduction_*4.0 < lastReduction_;
        }

        return true;
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const
    {
        if (avoidNegativeValueMacroScale_ and haveNegValue_)
        {   std::cout << "\nNegative concentrations in solution, set Newton return to be not converged !!! \n \n";
            return false;
        }
        if (fixedNewtonIterationSteps_ > 0)
        {
            return numSteps_ > fixedNewtonIterationSteps_;
        }
        if (enableShiftCriterion_ && !enableResidualCriterion_)
        {
            return shift_ <= shiftTolerance_;
        }
        else if (!enableShiftCriterion_ && enableResidualCriterion_)
        {
            return reduction_ <= reductionTolerance_;
        }
        else if (satisfyResidualAndShiftCriterion_)
        {
            return shift_ <= shiftTolerance_
                    && reduction_ <= reductionTolerance_;
        }
        else
        {
            return shift_ <= shiftTolerance_
                    || reduction_ <= reductionTolerance_;
        }

        return false;
    }

    /*!
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param method The object where the NewtonMethod is executed
     * \param u The initial solution
     */
    void newtonBegin(NewtonMethod &method, const SolutionVector &u)
    {
        method_ = &method;
        numSteps_ = 0;
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void newtonBeginStep()
    {
        lastShift_ = shift_;
        if (numSteps_ == 0)
        {
            lastReduction_ = 1.0;
        }
        else
        {
            lastReduction_ = reduction_;
        }
    }

    /*!
     * \brief Returns the number of steps done since newtonBegin() was
     *        called.
     */
    int newtonNumSteps()
    { return numSteps_; }

    /*!
     * \brief Update the maximum relative shift of the solution compared to
     *        the previous iteration.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    void newtonUpdateShift(const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        shift_ = 0;

        for (std::size_t i = 0; i < uLastIter[bulkIdx].size(); ++i)
        {
            auto uNewI = uLastIter[bulkIdx][i];
            uNewI -= deltaU[bulkIdx][i];

            const Scalar shiftAtDof = model_().relativeShiftAtDof(uLastIter[bulkIdx][i],
                                                                  uNewI);
            shift_ = std::max(shift_, shiftAtDof);
        }

        for (std::size_t i = 0; i < uLastIter[lowDimIdx].size(); ++i)
        {
            auto uNewI = uLastIter[lowDimIdx][i];
            uNewI -= deltaU[lowDimIdx][i];

            const Scalar shiftAtDof = model_().relativeShiftAtDof(uLastIter[lowDimIdx][i],
                                                                  uNewI);
            shift_ = std::max(shift_, shiftAtDof);
        }
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * If the linear solver doesn't accept multitype matrices we copy the matrix
     * into a 1x1 block BCRS matrix for solving.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<typename T = TypeTag>
    typename std::enable_if<!LinearSolverAcceptsMultiTypeMatrix<T>::value, void>::type
    newtonSolveLinear(JacobianMatrix &A,
                      SolutionVector &x,
                      SolutionVector &b)
    {
        try
        {
            if (numSteps_ == 0)
                initialResidual_ = b.two_norm();

            // copy the matrix and the vector to types the IterativeSolverBackend can handle
            using MatrixBlock = typename Dune::FieldMatrix<Scalar, 1, 1>;
            using SparseMatrix = typename Dune::BCRSMatrix<MatrixBlock>;

            if ((!decomposeMatrixSolver_) and (!decouplingEquations_))
            {
                // get the new matrix sizes
                std::size_t numRows = numEqBulk*A[bulkIdx][bulkIdx].N() + numEqLowDim*A[lowDimIdx][bulkIdx].N();
                std::size_t numCols = numEqBulk*A[bulkIdx][bulkIdx].M() + numEqLowDim*A[bulkIdx][lowDimIdx].M();

                // check matrix sizes
                assert(A[bulkIdx][bulkIdx].N() == A[bulkIdx][lowDimIdx].N());
                assert(A[lowDimIdx][bulkIdx].N() == A[lowDimIdx][lowDimIdx].N());
                assert(numRows == numCols);

                // create the bcrs matrix the IterativeSolver backend can handle
                auto M = SparseMatrix(numRows, numCols, SparseMatrix::random);

                // set the rowsizes
                // A11 and A12
                for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                    for (std::size_t i = 0; i < numEqBulk; ++i)
                        M.setrowsize(numEqBulk*row.index() + i, row->size()*numEqBulk);
                for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                    for (std::size_t i = 0; i < numEqBulk; ++i)
                        M.setrowsize(numEqBulk*row.index() + i, M.getrowsize(numEqBulk*row.index() + i) + row->size()*numEqLowDim);
                // A21 and A22
                for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                    for (std::size_t i = 0; i < numEqLowDim; ++i)
                        M.setrowsize(numEqLowDim*row.index() + i + A[bulkIdx][bulkIdx].N()*numEqBulk, row->size()*numEqBulk);
                for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                    for (std::size_t i = 0; i < numEqLowDim; ++i)
                        M.setrowsize(numEqLowDim*row.index() + i + A[bulkIdx][bulkIdx].N()*numEqBulk, M.getrowsize(numEqLowDim*row.index() + i + A[bulkIdx][bulkIdx].N()*numEqBulk) + row->size()*numEqLowDim);
                M.endrowsizes();

                // set the indices
                for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                    for (auto col = row->begin(); col != row->end(); ++col)
                        for (std::size_t i = 0; i < numEqBulk; ++i)
                            for (std::size_t j = 0; j < numEqBulk; ++j)
                                M.addindex(row.index()*numEqBulk + i, col.index()*numEqBulk + j);

                for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                    for (auto col = row->begin(); col != row->end(); ++col)
                        for (std::size_t i = 0; i < numEqBulk; ++i)
                            for (std::size_t j = 0; j < numEqLowDim; ++j)
                                M.addindex(row.index()*numEqBulk + i, col.index()*numEqLowDim + j + A[bulkIdx][bulkIdx].M()*numEqBulk);

                for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                    for (auto col = row->begin(); col != row->end(); ++col)
                        for (std::size_t i = 0; i < numEqLowDim; ++i)
                            for (std::size_t j = 0; j < numEqBulk; ++j)
                                M.addindex(row.index()*numEqLowDim + i + A[bulkIdx][bulkIdx].N()*numEqBulk, col.index()*numEqBulk + j);

                for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                    for (auto col = row->begin(); col != row->end(); ++col)
                        for (std::size_t i = 0; i < numEqLowDim; ++i)
                            for (std::size_t j = 0; j < numEqLowDim; ++j)
                                M.addindex(row.index()*numEqLowDim + i + A[bulkIdx][bulkIdx].N()*numEqBulk, col.index()*numEqLowDim + j + A[bulkIdx][bulkIdx].M()*numEqBulk);
                M.endindices();

                // copy values
                for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                    for (auto col = row->begin(); col != row->end(); ++col)
                        for (std::size_t i = 0; i < numEqBulk; ++i)
                            for (std::size_t j = 0; j < numEqBulk; ++j)
                                M[row.index()*numEqBulk + i][col.index()*numEqBulk + j] = A[bulkIdx][bulkIdx][row.index()][col.index()][i][j];

                for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                    for (auto col = row->begin(); col != row->end(); ++col)
                        for (std::size_t i = 0; i < numEqBulk; ++i)
                            for (std::size_t j = 0; j < numEqLowDim; ++j)
                                M[row.index()*numEqBulk + i][col.index()*numEqLowDim + j + A[bulkIdx][bulkIdx].M()*numEqBulk] = A[bulkIdx][lowDimIdx][row.index()][col.index()][i][j];

                for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                    for (auto col = row->begin(); col != row->end(); ++col)
                        for (std::size_t i = 0; i < numEqLowDim; ++i)
                            for (std::size_t j = 0; j < numEqBulk; ++j)
                                M[row.index()*numEqLowDim + i + A[bulkIdx][bulkIdx].N()*numEqBulk][col.index()*numEqBulk + j] = A[lowDimIdx][bulkIdx][row.index()][col.index()][i][j];

                for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                    for (auto col = row->begin(); col != row->end(); ++col)
                        for (std::size_t i = 0; i < numEqLowDim; ++i)
                            for (std::size_t j = 0; j < numEqLowDim; ++j)
                                M[row.index()*numEqLowDim + i + A[bulkIdx][bulkIdx].N()*numEqBulk][col.index()*numEqLowDim + j + A[bulkIdx][bulkIdx].M()*numEqBulk] = A[lowDimIdx][lowDimIdx][row.index()][col.index()][i][j];


                    //std::cout<<"A[bulkIdx][bulkIdx].N(): "<<A[bulkIdx][bulkIdx].N()<<"\n";
                    //std::cout<<"A[bulkIdx][bulkIdx].M(): "<<A[bulkIdx][bulkIdx].M()<<"\n";
                    //Dune::printmatrix(std::cout,A[bulkIdx][bulkIdx],"A11","m",10,2);
                    //std::cout<<"A[bulkIdx][lowDimIdx].N(): "<<A[bulkIdx][lowDimIdx].N()<<"\n";
                    //std::cout<<"A[bulkIdx][lowDimIdx].M(): "<<A[bulkIdx][lowDimIdx].M()<<"\n";
                    //Dune::printmatrix(std::cout,A[bulkIdx][lowDimIdx],"A12","m",10,2);
                    //std::cout<<"A[lowDimIdx][bulkIdx].N(): "<<A[lowDimIdx][bulkIdx].N()<<"\n";
                    //std::cout<<"A[lowDimIdx][bulkIdx].M(): "<<A[lowDimIdx][bulkIdx].M()<<"\n";
                    //Dune::printmatrix(std::cout,A[lowDimIdx][bulkIdx],"A21","m",10,2);
                    //std::cout<<"A[lowDimIdx][lowDimIdx].N(): "<<A[lowDimIdx][lowDimIdx].N()<<"\n";
                    //std::cout<<"A[lowDimIdx][lowDimIdx].M(): "<<A[lowDimIdx][lowDimIdx].M()<<"\n";
                    //Dune::printmatrix(std::cout,A[lowDimIdx][lowDimIdx],"A22","m",10,2);
                    //std::cout<<" -- Matrix M: \n";
                    //Dune::printmatrix(std::cout,M,"M","m",10,2);

                // create the vector the IterativeSolver backend can handle
                using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
                using BlockVector = typename Dune::BlockVector<VectorBlock>;

                BlockVector y, bTmp;
                y.resize(numRows);
                bTmp.resize(numCols);
                for (std::size_t i = 0; i < b[bulkIdx].N(); ++i)
                    for (std::size_t j = 0; j < numEqBulk; ++j)
                        bTmp[i*numEqBulk + j] = b[bulkIdx][i][j];
                for (std::size_t i = 0; i < b[lowDimIdx].N(); ++i)
                    for (std::size_t j = 0; j < numEqLowDim; ++j)
                        bTmp[i*numEqLowDim + j + b[bulkIdx].N()*numEqBulk] = b[lowDimIdx][i][j];

                    //std::cout<<" ---- Vector b[bulkIdx]: \n";
                    //Dune::printvector(std::cout,b[bulkIdx],"b[bulkIdx]","m",10,2);
                    //std::cout<<" ---- Vector b[lowDimIdx]: \n";
                    //Dune::printvector(std::cout,b[lowDimIdx],"b[lowDimIdx]","m",10,2);

                //std::cout<<" --- Vector bTmp: \n";
                //Dune::printvector(std::cout,bTmp,"bTmp","m",10,2);
                // solve
                bool converged = linearSolver_.solve(M, y, bTmp);
                    //std::cout<<" --- Vector y: \n";
                    //Dune::printvector(std::cout,y,"y","m",10,2);
                // copy back the result y into x
                for (std::size_t i = 0; i < x[bulkIdx].N(); ++i)
                    for (std::size_t j = 0; j < numEqBulk; ++j)
                        x[bulkIdx][i][j] = y[i*numEqBulk + j];
                for (std::size_t i = 0; i < x[lowDimIdx].N(); ++i)
                    for (std::size_t j = 0; j < numEqLowDim; ++j)
                        x[lowDimIdx][i][j] = y[i*numEqLowDim + j + x[bulkIdx].N()*numEqBulk];

                    //std::cout<<"\n ++++ Vector x[bulkIdx]: \n";
                    //Dune::printvector(std::cout,x[bulkIdx],"x[bulkIdx]","m",10,2);
                    //std::cout<<" ++++ Vector x[lowDimIdx]: \n";
                    //Dune::printvector(std::cout,x[lowDimIdx],"x[lowDimIdx]","m",10,2);
                    //std::cin.ignore();
                if (!converged)
                    DUNE_THROW(NumericalProblem, "Linear solver did not converge");
            }
            else
            {
                int eqNumberMax = std::max(numEqBulk, numEqLowDim);
                for ( int eqIdx = 0; eqIdx < eqNumberMax; ++eqIdx)
                {
                    // get the new matrix sizes
                    bool isSmallerThanNumEqBulk = (eqIdx < numEqBulk);
                    bool isSmallerThanNumEqLowDim = (eqIdx < numEqLowDim);
                    std::size_t numRows = isSmallerThanNumEqBulk*A[bulkIdx][bulkIdx].N() + isSmallerThanNumEqLowDim*A[lowDimIdx][bulkIdx].N();
                    std::size_t numCols = isSmallerThanNumEqBulk*A[bulkIdx][bulkIdx].M() + isSmallerThanNumEqLowDim*A[bulkIdx][lowDimIdx].M();

                    // check matrix sizes
                    assert(A[bulkIdx][bulkIdx].N() == A[bulkIdx][lowDimIdx].N());
                    assert(A[lowDimIdx][bulkIdx].N() == A[lowDimIdx][lowDimIdx].N());
                    assert(numRows == numCols);

                    // create the bcrs matrix the IterativeSolver backend can handle
                    auto M = SparseMatrix(numRows, numCols, SparseMatrix::random);

                    // set the rowsizes
                    // A11 and A12
                    for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                        for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                            M.setrowsize(isSmallerThanNumEqBulk*row.index() + i, row->size()*isSmallerThanNumEqBulk);
                    for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                        for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                            M.setrowsize(isSmallerThanNumEqBulk*row.index() + i, M.getrowsize(isSmallerThanNumEqBulk*row.index() + i) + row->size()*isSmallerThanNumEqLowDim);
                    // A21 and A22
                    for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                        for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                            M.setrowsize(isSmallerThanNumEqLowDim*row.index() + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk, row->size()*isSmallerThanNumEqBulk);
                    for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                        for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                            M.setrowsize(isSmallerThanNumEqLowDim*row.index() + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk, M.getrowsize(isSmallerThanNumEqLowDim*row.index() + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk) + row->size()*isSmallerThanNumEqLowDim);
                    M.endrowsizes();

                    // set the indices
                    for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                        for (auto col = row->begin(); col != row->end(); ++col)
                            for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                                    M.addindex(row.index()*isSmallerThanNumEqBulk + i, col.index()*isSmallerThanNumEqBulk + j);

                    for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                        for (auto col = row->begin(); col != row->end(); ++col)
                            for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                for (std::size_t j = 0; j < isSmallerThanNumEqLowDim; ++j)
                                    M.addindex(row.index()*isSmallerThanNumEqBulk + i, col.index()*isSmallerThanNumEqLowDim + j + A[bulkIdx][bulkIdx].M()*isSmallerThanNumEqBulk);

                    for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                        for (auto col = row->begin(); col != row->end(); ++col)
                            for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                                    M.addindex(row.index()*isSmallerThanNumEqLowDim + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk, col.index()*isSmallerThanNumEqBulk + j);

                    for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                        for (auto col = row->begin(); col != row->end(); ++col)
                            for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                for (std::size_t j = 0; j < isSmallerThanNumEqLowDim; ++j)
                                    M.addindex(row.index()*isSmallerThanNumEqLowDim + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk, col.index()*isSmallerThanNumEqLowDim + j + A[bulkIdx][bulkIdx].M()*isSmallerThanNumEqBulk);
                    M.endindices();

                    // copy values
                    for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                        for (auto col = row->begin(); col != row->end(); ++col)
                            for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                                    M[row.index()*isSmallerThanNumEqBulk + i][col.index()*isSmallerThanNumEqBulk + j] = A[bulkIdx][bulkIdx][row.index()][col.index()][i+eqIdx][j+eqIdx];

                    for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                        for (auto col = row->begin(); col != row->end(); ++col)
                            for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                for (std::size_t j = 0; j < isSmallerThanNumEqLowDim; ++j)
                                    M[row.index()*isSmallerThanNumEqBulk + i][col.index()*isSmallerThanNumEqLowDim + j + A[bulkIdx][bulkIdx].M()*isSmallerThanNumEqBulk] = A[bulkIdx][lowDimIdx][row.index()][col.index()][i+eqIdx][j+eqIdx];

                    for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                        for (auto col = row->begin(); col != row->end(); ++col)
                            for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                                    M[row.index()*isSmallerThanNumEqLowDim + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk][col.index()*isSmallerThanNumEqBulk + j] = A[lowDimIdx][bulkIdx][row.index()][col.index()][i+eqIdx][j+eqIdx];

                    for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                        for (auto col = row->begin(); col != row->end(); ++col)
                            for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                for (std::size_t j = 0; j < isSmallerThanNumEqLowDim; ++j)
                                    M[row.index()*isSmallerThanNumEqLowDim + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk][col.index()*isSmallerThanNumEqLowDim + j + A[bulkIdx][bulkIdx].M()*isSmallerThanNumEqBulk] = A[lowDimIdx][lowDimIdx][row.index()][col.index()][i+eqIdx][j+eqIdx];

                    //std::cout<<"A[bulkIdx][bulkIdx].N(): "<<A[bulkIdx][bulkIdx].N()<<"\n";
                    //std::cout<<"A[bulkIdx][bulkIdx].M(): "<<A[bulkIdx][bulkIdx].M()<<"\n";
                    //Dune::printmatrix(std::cout,A[bulkIdx][bulkIdx],"A11","m",10,2);
                    //std::cout<<"A[bulkIdx][lowDimIdx].N(): "<<A[bulkIdx][lowDimIdx].N()<<"\n";
                    //std::cout<<"A[bulkIdx][lowDimIdx].M(): "<<A[bulkIdx][lowDimIdx].M()<<"\n";
                    //Dune::printmatrix(std::cout,A[bulkIdx][lowDimIdx],"A12","m",10,2);
                    //std::cout<<"A[lowDimIdx][bulkIdx].N(): "<<A[lowDimIdx][bulkIdx].N()<<"\n";
                    //std::cout<<"A[lowDimIdx][bulkIdx].M(): "<<A[lowDimIdx][bulkIdx].M()<<"\n";
                    //Dune::printmatrix(std::cout,A[lowDimIdx][bulkIdx],"A21","m",10,2);
                    //std::cout<<"A[lowDimIdx][lowDimIdx].N(): "<<A[lowDimIdx][lowDimIdx].N()<<"\n";
                    //std::cout<<"A[lowDimIdx][lowDimIdx].M(): "<<A[lowDimIdx][lowDimIdx].M()<<"\n";
                    //Dune::printmatrix(std::cout,A[lowDimIdx][lowDimIdx],"A22","m",10,2);
                    //std::cout<<" -- Matrix M (eq "<<eqIdx<<"): \n";
                    //Dune::printmatrix(std::cout,M,"M","m",10,2);

                    // create the vector the IterativeSolver backend can handle
                    using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
                    using BlockVector = typename Dune::BlockVector<VectorBlock>;

                    BlockVector y, bTmp;
                    y.resize(numRows);
                    bTmp.resize(numCols);
                    for (std::size_t i = 0; i < b[bulkIdx].N(); ++i)
                        for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                            bTmp[i*isSmallerThanNumEqBulk + j] = b[bulkIdx][i][j+eqIdx];
                    for (std::size_t i = 0; i < b[lowDimIdx].N(); ++i)
                        for (std::size_t j = 0; j < isSmallerThanNumEqLowDim; ++j)
                            bTmp[i*isSmallerThanNumEqLowDim + j + b[bulkIdx].N()*isSmallerThanNumEqBulk] = b[lowDimIdx][i][j+eqIdx];

                    //std::cout<<" ---- Vector b[bulkIdx]: \n";
                    //Dune::printvector(std::cout,b[bulkIdx],"b[bulkIdx]","m",10,2);
                    //std::cout<<" ---- Vector b[lowDimIdx]: \n";
                    //Dune::printvector(std::cout,b[lowDimIdx],"b[lowDimIdx]","m",10,2);
                    //std::cout<<" ---- Vector bTmp: \n";
                    //Dune::printvector(std::cout,bTmp,"bTmp","m",10,2);
                    if ((!decouplingEquations_) and (eqIdx > 0))
                    {
                        std::size_t numColsK = isSmallerThanNumEqBulk*A[bulkIdx][bulkIdx].M() + 1*A[bulkIdx][lowDimIdx].M();
                        // create the bcrs matrix the IterativeSolver backend can handle
                        auto K = SparseMatrix(numRows, numColsK, SparseMatrix::random);

                        // set the rowsizes
                        // A11 and A12
                        for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                            for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                K.setrowsize(isSmallerThanNumEqBulk*row.index() + i, row->size()*isSmallerThanNumEqBulk);
                        for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                            for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                K.setrowsize(isSmallerThanNumEqBulk*row.index() + i, K.getrowsize(isSmallerThanNumEqBulk*row.index() + i) + row->size()*1);
                        // A21 and A22
                        for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                            for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                K.setrowsize(isSmallerThanNumEqLowDim*row.index() + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk, row->size()*isSmallerThanNumEqBulk);
                        for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                            for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                K.setrowsize(isSmallerThanNumEqLowDim*row.index() + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk, K.getrowsize(1*row.index() + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk) + row->size()*1);
                        K.endrowsizes();

                        // set the indices
                        for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                            for (auto col = row->begin(); col != row->end(); ++col)
                                for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                    for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                                        K.addindex(row.index()*isSmallerThanNumEqBulk + i, col.index()*isSmallerThanNumEqBulk + j);

                        for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                            for (auto col = row->begin(); col != row->end(); ++col)
                                for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                    for (std::size_t j = 0; j < 1; ++j)
                                        K.addindex(row.index()*isSmallerThanNumEqBulk + i, col.index()*1 + j + A[bulkIdx][bulkIdx].M()*isSmallerThanNumEqBulk);

                        for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                            for (auto col = row->begin(); col != row->end(); ++col)
                                for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                    for (std::size_t j = 0; j < 1; ++j)
                                        K.addindex(row.index()*isSmallerThanNumEqLowDim + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk, col.index()*isSmallerThanNumEqBulk + j);

                        for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                            for (auto col = row->begin(); col != row->end(); ++col)
                                for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                    for (std::size_t j = 0; j < 1; ++j)
                                        K.addindex(row.index()*isSmallerThanNumEqLowDim + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk, col.index()*1 + j + A[bulkIdx][bulkIdx].M()*isSmallerThanNumEqBulk);
                        K.endindices();

                        // copy values
                        for (auto row = A[bulkIdx][bulkIdx].begin(); row != A[bulkIdx][bulkIdx].end(); ++row)
                            for (auto col = row->begin(); col != row->end(); ++col)
                                for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                    for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                                        K[row.index()*isSmallerThanNumEqBulk + i][col.index()*isSmallerThanNumEqBulk + j]
                                                                    = A[bulkIdx][bulkIdx][row.index()][col.index()][i+eqIdx][j];

                        for (auto row = A[bulkIdx][lowDimIdx].begin(); row != A[bulkIdx][lowDimIdx].end(); ++row)
                            for (auto col = row->begin(); col != row->end(); ++col)
                                for (std::size_t i = 0; i < isSmallerThanNumEqBulk; ++i)
                                    for (std::size_t j = 0; j < 1; ++j)
                                        K[row.index()*isSmallerThanNumEqBulk + i][col.index()*1 + j + A[bulkIdx][bulkIdx].M()*isSmallerThanNumEqBulk]
                                                                    = A[bulkIdx][lowDimIdx][row.index()][col.index()][i+eqIdx][j];

                        for (auto row = A[lowDimIdx][bulkIdx].begin(); row != A[lowDimIdx][bulkIdx].end(); ++row)
                            for (auto col = row->begin(); col != row->end(); ++col)
                                for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                    for (std::size_t j = 0; j < 1; ++j)
                                        K[row.index()*isSmallerThanNumEqLowDim + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk][col.index()*isSmallerThanNumEqBulk + j]
                                                                    = A[lowDimIdx][bulkIdx][row.index()][col.index()][i+eqIdx][j];

                        for (auto row = A[lowDimIdx][lowDimIdx].begin(); row != A[lowDimIdx][lowDimIdx].end(); ++row)
                            for (auto col = row->begin(); col != row->end(); ++col)
                                for (std::size_t i = 0; i < isSmallerThanNumEqLowDim; ++i)
                                    for (std::size_t j = 0; j < 1; ++j)
                                        K[row.index()*isSmallerThanNumEqLowDim + i + A[bulkIdx][bulkIdx].N()*isSmallerThanNumEqBulk][col.index()*1 + j + A[bulkIdx][bulkIdx].M()*isSmallerThanNumEqBulk]
                                                                    = A[lowDimIdx][lowDimIdx][row.index()][col.index()][i+eqIdx][j];

                        //std::cout<<" -- Matrix K (eq "<<eqIdx<<"): \n";
                        //Dune::printmatrix(std::cout,K,"K","m",10,2);

                        using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
                        using BlockVector = typename Dune::BlockVector<VectorBlock>;
                        BlockVector xEqIdx0;
                        xEqIdx0.resize(numColsK);
                        for (std::size_t i = 0; i < b[bulkIdx].N(); ++i)
                            for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                                xEqIdx0[i*isSmallerThanNumEqBulk + j] = x[bulkIdx][i][j];
                                //xEqIdx0[i*isSmallerThanNumEqBulk + j] = problem_().bulkProblem().model().curSol()[i][j];
                        for (std::size_t i = 0; i < b[lowDimIdx].N(); ++i)
                            for (std::size_t j = 0; j < 1; ++j)
                                xEqIdx0[i*1 + j + b[bulkIdx].N()*isSmallerThanNumEqBulk] = x[lowDimIdx][i][j];
                                //xEqIdx0[i*1 + j + b[bulkIdx].N()*isSmallerThanNumEqBulk] = problem_().lowDimProblem().model().curSol()[i][j];

                        //std::cout<<" ---- Vector xEqIdx0: \n";
                        //Dune::printvector(std::cout,xEqIdx0,"xEqIdx0","m",10,2);

                        BlockVector b_;
                        b_.resize(numCols);
                        K.mv(xEqIdx0, b_);
                        //std::cout<<" ---- Vector b_: \n";
                        //Dune::printvector(std::cout,b_,"b_","m",10,2);

                        bTmp -= b_;
                        //std::cout<<" ---- Vector bTmp: \n";
                        //Dune::printvector(std::cout,bTmp,"bTmp","m",10,2);
                    }
                    // solve
                    bool converged = linearSolver_.solve(M, y, bTmp);
                    //std::cout<<" ---- Vector y: \n";
                    //Dune::printvector(std::cout,y,"y","m",10,2);
                    // copy back the result y into x
                    for (std::size_t i = 0; i < x[bulkIdx].N(); ++i)
                        for (std::size_t j = 0; j < isSmallerThanNumEqBulk; ++j)
                            x[bulkIdx][i][j+eqIdx] = y[i*isSmallerThanNumEqBulk + j];
                    for (std::size_t i = 0; i < x[lowDimIdx].N(); ++i)
                        for (std::size_t j = 0; j < isSmallerThanNumEqLowDim; ++j)
                            x[lowDimIdx][i][j+eqIdx] = y[i*isSmallerThanNumEqLowDim + j + x[bulkIdx].N()*isSmallerThanNumEqBulk];

                    //std::cout<<"\n ++++ Vector x[bulkIdx]: \n";
                    //Dune::printvector(std::cout,x[bulkIdx],"x[bulkIdx]","m",10,2);
                    //std::cout<<" ++++ Vector x[lowDimIdx]: \n";
                    //Dune::printvector(std::cout,x[lowDimIdx],"x[lowDimIdx]","m",10,2);
                    //std::cin.ignore();
                    if (!converged)
                        DUNE_THROW(NumericalProblem, "Linear solver did not converge");
                }
            }
        }
        catch (const Dune::Exception &e)
        {
            Dumux::NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<typename T = TypeTag>
    typename std::enable_if<LinearSolverAcceptsMultiTypeMatrix<T>::value, void>::type
    newtonSolveLinear(JacobianMatrix &A,
                      SolutionVector &x,
                      SolutionVector &b)
    {
        try
        {
            if (numSteps_ == 0)
                initialResidual_ = b.two_norm();

            bool converged = linearSolver_.solve(A, x, b);

            if (!converged)
                DUNE_THROW(NumericalProblem, "Linear solver did not converge");
        }
        catch (const Dune::Exception &e)
        {
            Dumux::NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the newtonConverged() and
     * newtonProceed() methods should be updated inside this method.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param uCurrentIter The solution vector after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        if (enableShiftCriterion_)
            newtonUpdateShift(uLastIter, deltaU);

        if (useLineSearch_)
            lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);

        else
        {
            for (std::size_t i = 0; i < uLastIter[bulkIdx].size(); ++i)
            {
                uCurrentIter[bulkIdx][i] = uLastIter[bulkIdx][i];
                uCurrentIter[bulkIdx][i] -= deltaU[bulkIdx][i];
            }
            for (std::size_t i = 0; i < uLastIter[lowDimIdx].size(); ++i)
            {
                uCurrentIter[lowDimIdx][i] = uLastIter[lowDimIdx][i];
                uCurrentIter[lowDimIdx][i] -= deltaU[lowDimIdx][i];
            }

            if (enableResidualCriterion_)
            {
                SolutionVector tmp(uLastIter);
                reduction_ = model_().globalResidual(tmp, uCurrentIter);
                reduction_ /= initialResidual_;
            }
        }

        // copy the global solution to the sub problems
        model_().copySolutionToSubProblems();
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void newtonEndStep(const SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        ++numSteps_;

        if (verbose())
        {
            std::cout << "\rNewton iteration " << numSteps_ << " done";
            if (enableShiftCriterion_)
                std::cout << ", maximum relative shift = " << shift_;
            if (enableResidualCriterion_)
                std::cout << ", residual reduction = " << reduction_;
            std::cout << endIterMsg().str() << "\n";
        }
        endIterMsgStream_.str("");

        // When the newton iteration is done: ask the model to check if it makes sense.
        model_().checkPlausibility();
    }

    /*!
     * \brief Indicates that we're done solving the non-linear system
     *        of equations.
     */
    void newtonEnd()
    {}

    /*!
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonFail()
    {
        numSteps_ = targetSteps_*2;
    }

    /*!
     * \brief Called when the Newton method was successful.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonSucceed()
    {}

    /*!
     * \brief Suggest a new time-step size based on the old time-step
     *        size.
     *
     * The default behavior is to suggest the old time-step size
     * scaled by the ratio between the target iterations and the
     * iterations required to actually solve the last time-step.
     */
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        // be aggressive reducing the time-step size but
        // conservative when increasing it. the rationale is
        // that we want to avoid failing in the next Newton
        // iteration which would require another linearization
        // of the problem.
        //Scalar suggestedTimeStepCourantNumber = problem_().bulkProblem().suggestedTimeStepCourantNumber();
        //Scalar suggestedTimeStepCourantNumber = problem_().bulkProblem().suggestedTimeStepModifiedCourantNumber();

        if (ratioSuggestedTimeStepSize_ >= 1)
            //if (oldTimeStep*ratioSuggestedTimeStepSize_ < suggestedTimeStepCourantNumber)
                return oldTimeStep*ratioSuggestedTimeStepSize_;
            //else
            //    return suggestedTimeStepCourantNumber;

        if (numSteps_ > targetSteps_) {
            Scalar percent = Scalar(numSteps_ - targetSteps_)/targetSteps_;
            return oldTimeStep/(1.0 + percent);
        }

        Scalar percent = Scalar(targetSteps_ - numSteps_)/targetSteps_;
        return oldTimeStep*(1.0 + percent/1.2);
    }

    /*!
     * \brief Returns a reference to the current Newton method
     *        which is controlled by this controller.
     */
    NewtonMethod &method()
    { return *method_; }

    /*!
     * \brief Returns a reference to the current Newton method
     *        which is controlled by this controller.
     */
    const NewtonMethod &method() const
    { return *method_; }

    std::ostringstream &endIterMsg()
    { return endIterMsgStream_; }

    /*!
     * \brief Specifies if the Newton method ought to be chatty.
     */
    void setVerbose(bool val)
    { verbose_ = val; }

    /*!
     * \brief Returns true if the Newton method ought to be chatty.
     */
    bool verbose() const
    { return verbose_; }

protected:

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem_()
    { return method_->problem(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    { return method_->problem(); }

    /*!
     * \brief Returns a reference to the time manager.
     */
    TimeManager &timeManager_()
    { return problem_().timeManager(); }

    /*!
     * \brief Returns a reference to the time manager.
     */
    const TimeManager &timeManager_() const
    { return problem_().timeManager(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    Model &model_()
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Model &model_() const
    { return problem_().model(); }

    // returns the actual implementation for the controller we do
    // it this way in order to allow "poor man's virtual methods",
    // i.e. methods of subclasses which can be called by the base
    // class.
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    void lineSearchUpdate_(SolutionVector &uCurrentIter,
                           const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        Scalar lambda = 1.0;
        SolutionVector tmp(uLastIter);

        while (true)
        {
           uCurrentIter = deltaU;
           uCurrentIter *= -lambda;
           uCurrentIter += uLastIter;

           // calculate the residual of the current solution
           reduction_ = this->method().model().globalResidual(tmp, uCurrentIter);
           reduction_ /= initialResidual_;

           if (reduction_ < lastReduction_ || lambda <= 0.125)
           {
               this->endIterMsg() << ", residual reduction " << lastReduction_ << "->"  << reduction_ << "@lambda=" << lambda;
               return;
           }

           // try with a smaller update
           lambda /= 2.0;
        }
    }

    std::ostringstream endIterMsgStream_;

    NewtonMethod *method_;
    bool verbose_;

    // shift criterion variables
    Scalar shift_;
    Scalar lastShift_;
    Scalar shiftTolerance_;

    // residual criterion variables
    Scalar reduction_;
    Scalar lastReduction_;
    Scalar initialResidual_;
    Scalar reductionTolerance_;

    // optimal number of iterations we want to achieve
    int targetSteps_;
    // maximum number of iterations we do before giving up
    int maxSteps_;
    // actual number of steps done so far
    int numSteps_;

    // the linear solver
    LinearSolver linearSolver_;

    bool useLineSearch_;
    bool enableShiftCriterion_;
    bool enableResidualCriterion_;
    bool satisfyResidualAndShiftCriterion_;
    bool haveNegValue_, avoidNegativeValueMacroScale_, fixedNewtonIterationSteps_, decouplingEquations_, decomposeMatrixSolver_;
    Scalar ratioSuggestedTimeStepSize_;
};

} // namespace Dumux

#endif
