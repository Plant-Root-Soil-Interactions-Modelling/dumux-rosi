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
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_MPFA_LOCAL_RESIDUAL_HH
#define DUMUX_MPFA_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/localresidual.hh>

#include "mpfaproperties.hh"
#include "mpfaimplicitlocalresidual.hh"

namespace Dumux
{
/*!
 * \ingroup CCModel
 * \ingroup CCLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit cell-centered scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class MpfaLocalResidual : public MpfaImplicitLocalResidual<TypeTag>
{
    typedef MpfaImplicitLocalResidual<TypeTag> ParentType;
    friend class MpfaImplicitLocalResidual<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    // copying the local residual class is not a good idea
    MpfaLocalResidual(const MpfaLocalResidual &);

public:
    MpfaLocalResidual() : ParentType()
    { }

protected:

    /*!
     * \brief Add all fluxes resulting from Neumann, outflow and pure Dirichlet
     *        boundary conditions to the local residual.
     */
    void evalBoundaryFluxes_()
    {

    }

    /*!
     * \brief Evaluate Dirichlet conditions on faces that have mixed
     *        Dirichlet/Neumann boundary conditions.
     */
    void evalDirichlet_()
    {

    }

    /*!
     * \brief Add Neumann boundary conditions for a single intersection
     */
    void evalNeumannSegment_(const IntersectionIterator &isIt,
                             const BoundaryTypes &bcTypes)
    {

    }

    /*!
     * \brief Add outflow boundary conditions for a single intersection
     */
    void evalOutflowSegment_(const IntersectionIterator &isIt,
                             const BoundaryTypes &bcTypes)
    {

    }

    /*!
     * \brief Treat Dirichlet boundary conditions in a weak sense for a single
     *        intersection that only has Dirichlet boundary conditions
     */
    void evalDirichletSegment_(const IntersectionIterator &isIt,
                                   const BoundaryTypes &bcTypes)
    {

    }

    /*!
     * \brief Treat Dirichlet boundary conditions in a strong sense for a
     *        single intersection that has mixed D/N boundary conditions
     */
    void evalDirichletSegmentMixed_(const IntersectionIterator &isIt,
                                    const BoundaryTypes &bcTypes)
    {

    }

    /*!
     * \brief Add the flux terms to the local residual of the current element
     */
    void evalFluxes_()
    {
        // calculate the mass flux over the faces and subtract
        // it from the local rates
        IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
        IntersectionIterator isEndIt = this->gridView_().iend(this->element_());
        for (; isIt != isEndIt; ++isIt) {
            // Edit1: boundary is treated within the interaction volumes
            // therefore commented out the following two lines
            //if (!isIt->neighbor())
            //    continue;

            int fIdx = isIt->indexInInside();
            PrimaryVariables flux;

            Valgrind::SetUndefined(flux);
            this->asImp_().computeFlux(flux, fIdx);
            Valgrind::CheckDefined(flux);

            flux *= this->curVolVars_(0).extrusionFactor();

            this->residual_[0] += flux;
        }
    }
};

}

#endif
