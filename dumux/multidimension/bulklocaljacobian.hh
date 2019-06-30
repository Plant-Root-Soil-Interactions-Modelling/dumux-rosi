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
 * \brief Caculates the Jacobian of the local residual for fully-implicit models
 */
#ifndef DUMUX_BULK_LOCAL_JACOBIAN_HH
#define DUMUX_BULK_LOCAL_JACOBIAN_HH

#include <dune/common/indices.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/multidimension/properties.hh>

namespace Dumux
{
/*!
 * \ingroup MultiDimension
 * \brief Calculates the Jacobian of the local residual for the bulk domain
 *
 * The default behavior is to use numeric differentiation, i.e.
 * forward or backward differences (2nd order), or central
 * differences (3rd order). The method used is determined by the
 * "NumericDifferenceMethod" property:
 *
 * - if the value of this property is smaller than 0, backward
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
 *   \f]
 *
 * - if the value of this property is 0, central
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
 *   \f]
 *
 * - if the value of this property is larger than 0, forward
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
 *   \f]
 *
 * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
 * is the value of a sub-control volume's primary variable at the
 * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
 */
template<class TypeTag>
class BulkLocalJacobian
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, BulkLocalJacobian) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) GlobalProblem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) GlobalModel;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP(TypeTag, SubProblemBlockIndices) SubProblemBlockIndices;

    typename SubProblemBlockIndices::BulkIdx bulkIdx;
    typename SubProblemBlockIndices::LowDimIdx lowDimIdx;

    // type of the local jacobian
    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixLittleBlockBulk MatrixLittleBlock;
    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixLittleBlockBulkCoupling CouplingMatrixLittleBlock;
    typedef typename Dune::Matrix<MatrixLittleBlock> MatrixBlock;
    typedef typename Dune::Matrix<CouplingMatrixLittleBlock> CouplingMatrixBlock;
    typedef Dune::MultiTypeBlockVector<MatrixBlock, CouplingMatrixBlock> LocalJacobian;

    // obtain the type tag of the bulk sub problems
    typedef typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag) BulkProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag) LowDimProblemTypeTag;

    // types of the bulk sub problem
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, LocalResidual) LocalResidual;

    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };
    enum { isBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox) };

    enum { bulkNumEq = GET_PROP_VALUE(BulkProblemTypeTag, NumEq) };
    enum { lowDimNumEq = GET_PROP_VALUE(LowDimProblemTypeTag, NumEq) };

    // copying a local jacobian is not a good idea
    BulkLocalJacobian(const BulkLocalJacobian &);

public:
    BulkLocalJacobian()
    {
        numericDifferenceMethod_ = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod);
        Valgrind::SetUndefined(globalProblemPtr_);
    }


    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(GlobalProblem &globalProblem)
    {
        globalProblemPtr_ = &globalProblem;
        localResidual_.init(problem_());
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * \param element The DUNE Codim<0> entity which we look at.
     */
    void assemble(const Element &element)
    {
        // set the current grid element and update the element's
        // finite volume geometry
        elemPtr_ = &element;
        fvElemGeom_.update(gridView_(), element, problem_());
        reset_();

        bcTypes_.update(problem_(), element_(), fvElemGeom_);

        // update the secondary variables for the element at the last
        // and the current time levels
        prevVolVars_.update(problem_(),
                            element_(),
                            fvElemGeom_,
                            true /* isOldSol? */);

        curVolVars_.update(problem_(),
                           element_(),
                           fvElemGeom_,
                           false /* isOldSol? */);

        // calculate the local residual
        localResidual().eval(element_(),
                             fvElemGeom_,
                             prevVolVars_,
                             curVolVars_,
                             bcTypes_);

        residual_ = localResidual().residual();

        // get stencil informations
        const auto& stencil = globalProblem_().couplingManager().stencil(element_());
        const auto& couplingStencil = globalProblem_().couplingManager().couplingStencil(element_());

        // set size of local jacobian matrix
        const std::size_t numCols = stencil.size();
        const std::size_t numColsCoupling = couplingStencil.size();
        std::size_t numRows;

        if (isBox)
            numRows = numCols;
        else
            numRows = 1;

        // set sizes of local jacobians
        A_[bulkIdx].setSize(numRows, numCols);
        A_[lowDimIdx].setSize(numRows, numColsCoupling);

        ElementSolutionVector partialDeriv(numRows);
        for (int col = 0; col < numCols; col++)
        {
            for (int pvIdx = 0; pvIdx < bulkNumEq; pvIdx++)
            {
                evalPartialDerivative_(partialDeriv,
                                       col,
                                       pvIdx);

                // update the local stiffness matrix with the current partial
                // derivatives
                updateLocalJacobian_(bulkIdx,
                                     col,
                                     pvIdx,
                                     partialDeriv);
            }
        }

        for (int col = 0; col < numColsCoupling; col++)
        {
            for (int pvIdx = 0; pvIdx < lowDimNumEq; pvIdx++)
            {
                evalPartialDerivativeCoupling_(partialDeriv,
                                               couplingStencil,
                                               col,
                                               pvIdx);

                // update the local stiffness matrix with the current partial
                // derivatives
                updateLocalJacobian_(lowDimIdx,
                                     col,
                                     pvIdx,
                                     partialDeriv);
            }
        }
    }

    /*!
     * \brief Returns a reference to the object which calculates the
     *        local residual.
     */
    const LocalResidual &localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns a reference to the object which calculates the
     *        local residual.
     */
    LocalResidual &localResidual()
    { return localResidual_; }

    /*!
     * \brief Returns the Jacobian of the equations at subcontrolvolume i
     * to the primary variables at subcontrolvolume j.
     *
     * \param i The local subcontrolvolume index on which
     *          the equations are defined
     * \param j The local subcontrolvolume index which holds
     *          primary variables
     */
    const MatrixLittleBlock &mat(typename SubProblemBlockIndices::BulkIdx idx,
                                 const int i, const int j) const
    { return A_[idx][i][j]; }

    const CouplingMatrixLittleBlock &mat(typename SubProblemBlockIndices::LowDimIdx idx,
                                         const int i, const int j) const
    { return A_[idx][i][j]; }

    /*!
     * \brief Returns the residual of the equations at subcontrolvolume i.
     *
     * \param i The local subcontrolvolume index on which
     *          the equations are defined
     */
    const PrimaryVariables &residual(const int i) const
    { return residual_[i]; }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     * \param domainIdx  The index of the sub problem domain
     *                   of the priVar with respect to which we want to derive
     * \param dofIdxGlobal  The global index of the element's subcontrolvolume for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    template <std::size_t index>
    Scalar numericEpsilon(Dune::index_constant<index> domainIdx,
                          const int dofIdxGlobal,
                          const int pvIdx) const
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        static const Scalar baseEps = 1e-10;
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        Scalar priVar;
        // derivative with respect to the bulk domain
        if (domainIdx == bulkIdx)
            priVar = globalProblem_().bulkProblem().model().curSol()[dofIdxGlobal][pvIdx];

        // derivative with respect to the low dimensional domain
        else
            priVar = globalProblem_().lowDimProblem().model().curSol()[dofIdxGlobal][pvIdx];

        return baseEps*(std::abs(priVar) + 1.0);
    }

protected:

    /*!
     * \brief Returns a reference to the problem.
     */
    const GlobalProblem &globalProblem_() const
    {
        Valgrind::CheckDefined(globalProblemPtr_);
        return *globalProblemPtr_;
    }

    GlobalProblem &globalProblem_()
    {
        Valgrind::CheckDefined(globalProblemPtr_);
        return *globalProblemPtr_;
    }

    /*!
     * \brief Returns a reference to the sub problem.
     */
    const Problem &problem_() const
    {
        Valgrind::CheckDefined(globalProblemPtr_);
        return globalProblemPtr_->bulkProblem();
    }

    /*!
     * \brief Returns a reference to the sub problem.
     */
    Problem &problem_()
    {
        Valgrind::CheckDefined(globalProblemPtr_);
        return globalProblemPtr_->bulkProblem();
    }

    /*!
     * \brief Returns a reference to the bulk grid view.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns a reference to the actual bulk element.
     */
    const Element &element_() const
    {
        Valgrind::CheckDefined(elemPtr_);
        return *elemPtr_;
    }

    /*!
     * \brief Returns a reference to the model.
     */
    const GlobalModel &globalModel_() const
    { return globalProblem_().model(); }

     /*!
     * \brief Returns a reference to the sub model.
     */
    const Model &model_() const
    { return problem_().model(); }

    Model &model_()
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the jacobian assembler.
     */
    const JacobianAssembler &jacAsm_() const
    { return globalModel_().jacobianAssembler(); }

    /*!
     * \brief Reset the local jacobian matrix to 0
     */
    void reset_()
    { A_ = 0.0; }

    /*!
     * \brief Compute the partial derivatives to a primary variable at
     *        an degree of freedom.
     *
     * This method can be overwritten by the implementation if a
     * better scheme than numerical differentiation is available.
     *
     * The default implementation of this method uses numeric
     * differentiation, i.e. forward or backward differences (2nd
     * order), or central differences (3rd order). The method used is
     * determined by the "NumericDifferenceMethod" property:
     *
     * - if the value of this property is smaller than 0, backward
     *   differences are used, i.e.:
     *   \f[
         \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
     *   \f]
     *
     * - if the value of this property is 0, central
     *   differences are used, i.e.:
     *   \f[
           \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
     *   \f]
     *
     * - if the value of this property is larger than 0, forward
     *   differences are used, i.e.:
     *   \f[
           \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
     *   \f]
     *
     * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
     * is the value of a sub-control volume's primary variable at the
     * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
     *
     * \param partialDeriv The vector storing the partial derivatives of all
     *              equations
     * \param storageDeriv the mass matrix contributions
     * \param col The block column index of the degree of freedom
     *            for which the partial derivative is calculated.
     *            Box: a sub-control volume index.
     *            Cell centered: a neighbor index.
     * \param pvIdx The index of the primary variable
     *              for which the partial derivative is calculated
     */
    void evalPartialDerivative_(ElementSolutionVector &partialDeriv,
                                const int col,
                                const int pvIdx)
    {
        int dofIdxGlobal;
        FVElementGeometry neighborFVGeom;
        auto neighbor = element_();
        if (isBox)
        {
            dofIdxGlobal = problem_().vertexMapper().subIndex(element_(), col, dim);
        }
        else
        {
            neighbor = fvElemGeom_.neighbors[col];
            neighborFVGeom.updateInner(neighbor);
            dofIdxGlobal = problem_().elementMapper().index(neighbor);
        }

        PrimaryVariables &priVars = model_().curSol()[dofIdxGlobal];
        PrimaryVariables origPriVars(priVars);
        VolumeVariables origVolVars(curVolVars_[col]);

        curVolVars_[col].setEvalPoint(&origVolVars);
        const Scalar eps = numericEpsilon(bulkIdx, dofIdxGlobal, pvIdx);
        Scalar delta = 0;

        if (numericDifferenceMethod_ >= 0)
        {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            // calculate the residual
            if (isBox)
            {
                curVolVars_[col].update(priVars,
                                        problem_(),
                                        element_(),
                                        fvElemGeom_,
                                        col,
                                        false);
            }
            else
            {
                curVolVars_[col].update(priVars,
                                        problem_(),
                                        neighbor,
                                        neighborFVGeom,
                                        /*scvIdx=*/0,
                                        false);
            }

            localResidual().eval(element_(),
                                 fvElemGeom_,
                                 prevVolVars_,
                                 curVolVars_,
                                 bcTypes_);

            // store the residual and the storage term
            partialDeriv = localResidual().residual();
        }
        else
        {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv = residual_;
        }


        if (numericDifferenceMethod_ <= 0)
        {
            // we are not using forward differences, i.e. we
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= 2*eps;
            delta += eps;

            // calculate residual again
            if (isBox)
                curVolVars_[col].update(priVars,
                                        problem_(),
                                        element_(),
                                        fvElemGeom_,
                                        col,
                                        false);
            else
                curVolVars_[col].update(priVars,
                                        problem_(),
                                        neighbor,
                                        neighborFVGeom,
                                        /*scvIdx=*/0,
                                        false);

            localResidual().eval(element_(),
                                 fvElemGeom_,
                                 prevVolVars_,
                                 curVolVars_,
                                 bcTypes_);
            partialDeriv -= localResidual().residual();
        }
        else
        {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv -= residual_;
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;

        // restore the original state of the solution vector and the
        // element's volume variables
        priVars = origPriVars;
        curVolVars_[col] = origVolVars;


#if HAVE_VALGRIND
        for (unsigned i = 0; i < partialDeriv.size(); ++i)
            Valgrind::CheckDefined(partialDeriv[i]);
#endif
    }

    void evalPartialDerivativeCoupling_(ElementSolutionVector &partialDeriv,
                                        const std::vector<unsigned int>& couplingStencil,
                                        const int col,
                                        const int pvIdx)
    {
        const unsigned int dofIdxGlobal = couplingStencil[col];
        auto& priVars = globalProblem_().lowDimProblem().model().curSol()[dofIdxGlobal];
        auto origPriVars = priVars;

        const Scalar eps = numericEpsilon(lowDimIdx, dofIdxGlobal, pvIdx);
        Scalar delta = 0;

        if (numericDifferenceMethod_ >= 0)
        {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            localResidual().eval(element_(),
                                 fvElemGeom_,
                                 prevVolVars_,
                                 curVolVars_,
                                 bcTypes_);

            // store the residual and the storage term
            partialDeriv = localResidual().residual();
        }
        else
        {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv = residual_;
        }


        if (numericDifferenceMethod_ <= 0)
        {
            // we are not using forward differences, i.e. we
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= 2*eps;
            delta += eps;

            localResidual().eval(element_(),
                                 fvElemGeom_,
                                 prevVolVars_,
                                 curVolVars_,
                                 bcTypes_);
            partialDeriv -= localResidual().residual();
        }
        else
        {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv -= residual_;
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;

        // restore the original state of the solution vector and the
        // element's volume variables
        priVars = origPriVars;

#if HAVE_VALGRIND
        for (unsigned i = 0; i < partialDeriv.size(); ++i)
            Valgrind::CheckDefined(partialDeriv[i]);
#endif
    }

    /*!
     * \brief Updates the current local Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at dof 'col' .
     */
    template <std::size_t index>
    void updateLocalJacobian_(Dune::index_constant<index> domainIdx,
                              const int col,
                              const int pvIdx,
                              const ElementSolutionVector &partialDeriv)
    {
        for (int i = 0; i < partialDeriv.size(); i++)
        {
            for (int eqIdx = 0; eqIdx < bulkNumEq; eqIdx++)
            {
                // A[domainIdx][i][col][eqIdx][pvIdx] is the rate of change of
                // the residual of equation 'eqIdx' at dof 'i'
                // depending on the primary variable 'pvIdx' at dof
                // 'col' with respect to a domainIdx element.
                this->A_[domainIdx][i][col][eqIdx][pvIdx] = partialDeriv[i][eqIdx];
                Valgrind::CheckDefined(this->A_[domainIdx][i][col][eqIdx][pvIdx]);
            }
        }
    }

    const Element *elemPtr_;
    FVElementGeometry fvElemGeom_;

    ElementBoundaryTypes bcTypes_;

    // The problem we would like to solve
    GlobalProblem *globalProblemPtr_;

    // secondary variables at the previous
    // and the current time levels
    ElementVolumeVariables prevVolVars_;
    ElementVolumeVariables curVolVars_;

    LocalResidual localResidual_;

    LocalJacobian A_;

    ElementSolutionVector residual_;

    int numericDifferenceMethod_;
};

}

#endif
