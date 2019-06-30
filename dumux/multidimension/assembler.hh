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
 * \ingroup MultiDimension
 * \brief An assembler for the global Jacobian matrix for multidimension models.
 *        We assume a bulk domain and a lower dimensional domain on separate grids.
 */
#ifndef DUMUX_MULTIDIMENSION_ASSEMBLER_HH
#define DUMUX_MULTIDIMENSION_ASSEMBLER_HH

#include <dumux/common/exceptions.hh>
#include <dumux/multidimension/properties.hh>

namespace Dumux {

/*!
 * \ingroup MultiDimension
 * \brief An assembler for the global Jacobian matrix for multidimension models.
 *        We assume a bulk domain and a lower dimensional domain on separate grids.
 */
template<class TypeTag>
class MultiDimensionAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // obtain the type tags of the sub problems
    typedef typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag) BulkProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag) LowDimProblemTypeTag;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, GridView) BulkGridView;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView) LowDimGridView;

    typedef typename BulkGridView::template Codim<0>::Entity BulkElement;
    typedef typename LowDimGridView::template Codim<0>::Entity LowDimElement;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP(TypeTag, SubProblemBlockIndices) SubProblemBlockIndices;

    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockBulk BulkMatrixBlock;
    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockBulkCoupling BulkCouplingMatrixBlock;
    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockLowDim LowDimMatrixBlock;
    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockLowDimCoupling LowDimCouplingMatrixBlock;

    typename SubProblemBlockIndices::BulkIdx bulkIdx;
    typename SubProblemBlockIndices::LowDimIdx lowDimIdx;

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    enum {
        bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox),
        lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox)
    };

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

    // copying the jacobian assembler is not a good idea
    MultiDimensionAssembler(const MultiDimensionAssembler &);

public:

    MultiDimensionAssembler() : problemPtr_(nullptr) {}

    /*!
     * \brief Initialize the jacobian assembler.
     *
     * At this point we can assume that all objects in the problem and
     * the model have been allocated :. We can not assume that they are
     * fully initialized, though.
     *
     * \param problem The problem object
     */
    void init(Problem& problem)
    {
        // save problem pointer
        problemPtr_ = &problem;

        // initialize the multitype matrix
        asImp_().createMatrix_();

        // initialize the jacobian matrix with zeros
        *matrix_ = 0.0;

        // allocate the residual vector (right-hand-side)
        residual_[bulkIdx].resize(problem.model().bulkNumDofs());
        residual_[lowDimIdx].resize(problem.model().lowDimNumDofs());

        printMatrix_ = tree.template get<bool>("MultiDimension.PrintMatrix", false);
    }

    /*!
     * \brief Assemble the global Jacobian of the residual and the residual for the current solution.
     *
     * The current state of affairs (esp. the previous and the current
     * solutions) is represented by the model object.
     */
    void assemble()
    {
        bool succeeded;
        try
        {
            asImp_().assemble_();
            succeeded = true;
        }
        catch (Dumux::NumericalProblem &e)
        {
            std::cout << " caught an exception while assembling:" << e.what()
                      << std::endl;
            succeeded = false;
        }

        if (!succeeded) {
            DUNE_THROW(NumericalProblem,
                       "A process did not succeed in linearizing the system");
        }

        if (printMatrix_)
        {
            std::cout<<"matrix_[bulkIdx][bulkIdx].N(): "<<(*matrix_)[bulkIdx][bulkIdx].N()<<"\n";
            std::cout<<"matrix_[bulkIdx][bulkIdx].M(): "<<(*matrix_)[bulkIdx][bulkIdx].M()<<"\n";
            Dune::printmatrix(std::cout,(*matrix_)[bulkIdx][bulkIdx],"A11","m",10,2);
            std::cout<<"matrix_[bulkIdx][lowDimIdx].N(): "<<(*matrix_)[bulkIdx][lowDimIdx].N()<<"\n";
            std::cout<<"matrix_[bulkIdx][lowDimIdx].M(): "<<(*matrix_)[bulkIdx][lowDimIdx].M()<<"\n";
            Dune::printmatrix(std::cout,(*matrix_)[bulkIdx][lowDimIdx],"A12","m",10,2);
            std::cout<<"matrix_[lowDimIdx][bulkIdx].N(): "<<(*matrix_)[lowDimIdx][bulkIdx].N()<<"\n";
            std::cout<<"matrix_[lowDimIdx][bulkIdx].M(): "<<(*matrix_)[lowDimIdx][bulkIdx].M()<<"\n";
            Dune::printmatrix(std::cout,(*matrix_)[lowDimIdx][bulkIdx],"A21","m",10,2);
            std::cout<<"matrix_[lowDimIdx][lowDimIdx].N(): "<<(*matrix_)[lowDimIdx][lowDimIdx].N()<<"\n";
            std::cout<<"matrix_[lowDimIdx][lowDimIdx].M(): "<<(*matrix_)[lowDimIdx][lowDimIdx].M()<<"\n";
            Dune::printmatrix(std::cout,(*matrix_)[lowDimIdx][lowDimIdx],"A22","m",10,2);
            //Dune::printmatrix(std::cout,residual_,"R","row", 200,1,3);
            //Dune::printmatrix(std::cout,(residual_)[bulkIdx][0],"R1","m",10,2);
            //Dune::printmatrix(std::cout,(residual_)[lowDimIdx][0],"R2","m",10,2);
        }
    }

    /*!
     * \brief Return constant reference to global Jacobian matrix.
     */
    const JacobianMatrix& matrix() const
    { return *matrix_; }
    JacobianMatrix& matrix()
    { return *matrix_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const SolutionVector& residual() const
    { return residual_; }
    SolutionVector& residual()
    { return residual_; }

protected:
    // reset the global linear system of equations.
    void resetSystem_()
    {
        // reset the right hand side.
        residual_ = 0.0;

        // reset the matrix
        *matrix_ = 0.0;
    }

    // linearize the whole system
    void assemble_()
    {
        resetSystem_();

        // assemble the elements
        for (const auto& element : elements(problem_().lowDimGridView()))
        {
            asImp_().assembleElement_(element);
        }

        for (const auto& element : elements(problem_().bulkGridView()))
        {
            asImp_().assembleElement_(element);
        }
    }

    void buildBulkMatrixBlocksCC_(BulkMatrixBlock& m, BulkCouplingMatrixBlock& cm)
    {
        for (const auto& element : elements(problem_().bulkGridView()))
        {
            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);
            const unsigned int eIdx = problem_().model().bulkElementMapper().index(element);
            m.setrowsize(eIdx, stencil.size());
            cm.setrowsize(eIdx, couplingStencil.size());
        }
        m.endrowsizes();
        cm.endrowsizes();

        for (const auto& element : elements(problem_().bulkGridView()))
        {
            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);
            const unsigned int globalI = problem_().model().bulkElementMapper().index(element);

            for (auto&& globalJ : stencil)
                m.addindex(globalI, globalJ);

            for (auto&& globalJ : couplingStencil)
                cm.addindex(globalI, globalJ);
        }
        m.endindices();
        cm.endindices();
    }

    void buildLowDimMatrixBlocksCC_(LowDimMatrixBlock& m, LowDimCouplingMatrixBlock& cm)
    {
        for (const auto& element : elements(problem_().lowDimGridView()))
        {
            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);
            const unsigned int eIdx = problem_().model().lowDimElementMapper().index(element);
            m.setrowsize(eIdx, stencil.size());
            cm.setrowsize(eIdx, couplingStencil.size());
        }
        m.endrowsizes();
        cm.endrowsizes();

        for (const auto& element : elements(problem_().lowDimGridView()))
        {
            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);
            const unsigned int globalI = problem_().model().lowDimElementMapper().index(element);

            for (auto&& globalJ : stencil)
                m.addindex(globalI, globalJ);

            for (auto&& globalJ : couplingStencil)
                cm.addindex(globalI, globalJ);
        }
        m.endindices();
        cm.endindices();
    }

    void buildBulkMatrixBlocksBox_(BulkMatrixBlock& m, BulkCouplingMatrixBlock& cm)
    {
        std::vector<std::set<unsigned int>> dofIndexSet;
        std::vector<std::set<unsigned int>> dofIndexSetCoupling;
        dofIndexSet.resize(problem_().bulkGridView().size(bulkDim));
        dofIndexSetCoupling.resize(problem_().bulkGridView().size(bulkDim));
        for (const auto& element : elements(problem_().bulkGridView()))
        {
            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

            for (unsigned int scvIdx = 0; scvIdx < element.subEntities(bulkDim); ++scvIdx)
            {
                auto globalI = problem_().model().bulkVertexMapper().subIndex(element, scvIdx, bulkDim);
                dofIndexSet[globalI].insert(std::begin(stencil), std::end(stencil));
                dofIndexSetCoupling[globalI].insert(std::begin(couplingStencil), std::end(couplingStencil));
            }
        }

        for (unsigned int globalI = 0; globalI < dofIndexSet.size(); ++globalI)
        {
            m.setrowsize(globalI, dofIndexSet[globalI].size());
            cm.setrowsize(globalI, dofIndexSetCoupling[globalI].size());
        }
        m.endrowsizes();
        cm.endrowsizes();

        for (unsigned int globalI = 0; globalI < dofIndexSet.size(); ++globalI)
        {
            for (auto&& globalJ : dofIndexSet[globalI])
                m.addindex(globalI, globalJ);

            for (auto&& globalJ : dofIndexSetCoupling[globalI])
                cm.addindex(globalI, globalJ);
        }
        m.endindices();
        cm.endindices();
    }

    void buildLowDimMatrixBlocksBox_(LowDimMatrixBlock& m, LowDimCouplingMatrixBlock& cm)
    {
        std::vector<std::set<unsigned int>> dofIndexSet;
        std::vector<std::set<unsigned int>> dofIndexSetCoupling;
        dofIndexSet.resize(problem_().lowDimGridView().size(lowDimDim));
        dofIndexSetCoupling.resize(problem_().lowDimGridView().size(lowDimDim));
        for (const auto& element : elements(problem_().lowDimGridView()))
        {
            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

            for (unsigned int scvIdx = 0; scvIdx < element.subEntities(lowDimDim); ++scvIdx)
            {
                auto globalI = problem_().model().lowDimVertexMapper().subIndex(element, scvIdx, lowDimDim);
                dofIndexSet[globalI].insert(std::begin(stencil), std::end(stencil));
                dofIndexSetCoupling[globalI].insert(std::begin(couplingStencil), std::end(couplingStencil));
            }
        }

        for (unsigned int globalI = 0; globalI < dofIndexSet.size(); ++globalI)
        {
            m.setrowsize(globalI, dofIndexSet[globalI].size());
            cm.setrowsize(globalI, dofIndexSetCoupling[globalI].size());
        }
        m.endrowsizes();
        cm.endrowsizes();

        for (unsigned int globalI = 0; globalI < dofIndexSet.size(); ++globalI)
        {
            for (auto&& globalJ : dofIndexSet[globalI])
                m.addindex(globalI, globalJ);

            for (auto&& globalJ : dofIndexSetCoupling[globalI])
                cm.addindex(globalI, globalJ);
        }
        m.endindices();
        cm.endindices();
    }

    // Construct the multitype matrix for the global jacobian
    void createMatrix_()
    {
        // create the multitype matrix
        matrix_ = std::make_shared<JacobianMatrix>();

        // get sub matrix sizes
        const auto bulkSize = problem_().model().bulkNumDofs();
        const auto lowDimSize = problem_().model().lowDimNumDofs();

        // allocate the sub matrices
        auto A11 = BulkMatrixBlock(bulkSize, bulkSize, BulkMatrixBlock::random);
        auto A12 = BulkCouplingMatrixBlock(bulkSize, lowDimSize, BulkCouplingMatrixBlock::random);
        auto A22 = LowDimMatrixBlock(lowDimSize, lowDimSize, LowDimMatrixBlock::random);
        auto A21 = LowDimCouplingMatrixBlock(lowDimSize, bulkSize, LowDimCouplingMatrixBlock::random);

        // cell-centered
        if (!lowDimIsBox)
            buildLowDimMatrixBlocksCC_(A22, A21);
        else
            buildLowDimMatrixBlocksBox_(A22, A21);

        if (!bulkIsBox)
            buildBulkMatrixBlocksCC_(A11, A12);
        else
            buildBulkMatrixBlocksBox_(A11, A12);

        (*matrix_)[bulkIdx][bulkIdx] = A11;
        (*matrix_)[bulkIdx][lowDimIdx] = A12;
        (*matrix_)[lowDimIdx][bulkIdx] = A21;
        (*matrix_)[lowDimIdx][lowDimIdx] = A22;
    }

    // assemble a bulk element
    void assembleElement_(const BulkElement &element)
    {
        problem_().model().bulkLocalJacobian().assemble(element);

        if (!bulkIsBox)
        {
            int globalI = problem_().model().bulkElementMapper().index(element);

            // update the right hand side
            residual_[bulkIdx][globalI] =  problem_().model().bulkLocalJacobian().residual(0);
            //for (int j = 0; j < residual_[bulkIdx][globalI].dimension; ++j)
            //    assert(std::isfinite(residual_[bulkIdx][globalI][j]));

            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

            int j = 0;
            for (auto&& globalJ : couplingStencil)
                (*matrix_)[bulkIdx][lowDimIdx][globalI][globalJ] = problem_().model().bulkLocalJacobian().mat(lowDimIdx, 0, j++);

            j = 0;
            for (auto&& globalJ : stencil)
                (*matrix_)[bulkIdx][bulkIdx][globalI][globalJ] = problem_().model().bulkLocalJacobian().mat(bulkIdx, 0, j++);
        }
        else
        {
            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

            for (unsigned int scvIdx = 0; scvIdx < element.subEntities(bulkDim); ++scvIdx)
            {
                auto globalI = problem_().model().bulkVertexMapper().subIndex(element, scvIdx, bulkDim);

                // update the right hand side
                residual_[bulkIdx][globalI] += problem_().model().bulkLocalJacobian().residual(scvIdx);
                for (int j = 0; j < residual_[bulkIdx][globalI].dimension; ++j)
                    assert(std::isfinite(residual_[bulkIdx][globalI][j]));

                int j = 0;
                for (auto&& globalJ : couplingStencil)
                    (*matrix_)[bulkIdx][lowDimIdx][globalI][globalJ] += problem_().model().bulkLocalJacobian().mat(lowDimIdx, scvIdx, j++);

                j = 0;
                for (auto&& globalJ : stencil)
                    (*matrix_)[bulkIdx][bulkIdx][globalI][globalJ] += problem_().model().bulkLocalJacobian().mat(bulkIdx, scvIdx, j++);
            }
        }
    }

    // assemble a bulk element
    void assembleElement_(const LowDimElement &element)
    {
        problem_().model().lowDimLocalJacobian().assemble(element);

        if (!lowDimIsBox)
        {
            int globalI = problem_().model().lowDimElementMapper().index(element);

            // update the right hand side
            residual_[lowDimIdx][globalI] =  problem_().model().lowDimLocalJacobian().residual(0);
            //for (int j = 0; j < residual_[lowDimIdx][globalI].dimension; ++j)
            //    assert(std::isfinite(residual_[lowDimIdx][globalI][j]));

            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

            int j = 0;
            for (auto&& globalJ : stencil)
                (*matrix_)[lowDimIdx][lowDimIdx][globalI][globalJ] = problem_().model().lowDimLocalJacobian().mat(lowDimIdx, 0, j++);

            j = 0;
            for (auto&& globalJ : couplingStencil)
                (*matrix_)[lowDimIdx][bulkIdx][globalI][globalJ] = problem_().model().lowDimLocalJacobian().mat(bulkIdx, 0, j++);
        }
        else
        {
            const auto& stencil = problem_().couplingManager().stencil(element);
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

            for (unsigned int scvIdx = 0; scvIdx < element.subEntities(lowDimDim); ++scvIdx)
            {
                auto globalI = problem_().model().lowDimVertexMapper().subIndex(element, scvIdx, lowDimDim);

                // update the right hand side
                residual_[lowDimIdx][globalI] +=  problem_().model().lowDimLocalJacobian().residual(scvIdx);
                for (int j = 0; j < residual_[lowDimIdx][globalI].dimension; ++j)
                    assert(std::isfinite(residual_[lowDimIdx][globalI][j]));

                int j = 0;
                for (auto&& globalJ : stencil)
                    (*matrix_)[lowDimIdx][lowDimIdx][globalI][globalJ] += problem_().model().lowDimLocalJacobian().mat(lowDimIdx, scvIdx, j++);

                j = 0;
                for (auto&& globalJ : couplingStencil)
                    (*matrix_)[lowDimIdx][bulkIdx][globalI][globalJ] += problem_().model().lowDimLocalJacobian().mat(bulkIdx, scvIdx, j++);
            }
        }
    }

    Problem &problem_()
    { return *problemPtr_; }
    const Problem &problem_() const
    { return *problemPtr_; }

    // the multidimension problem
    Problem *problemPtr_;

    // the jacobian matrix
    std::shared_ptr<JacobianMatrix> matrix_;
    // the right-hand side
    SolutionVector residual_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    bool printMatrix_;
};

} // namespace Dumux

#endif
