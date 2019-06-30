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
 * \brief Dumux multidimension direct solver backend
 */
#ifndef DUMUX_MULTIDIMENSION_DIRECT_SOLVER_BACKEND_HH
#define DUMUX_MULTIDIMENSION_DIRECT_SOLVER_BACKEND_HH

#include <dumux/linear/seqsolverbackend.hh>

namespace Dumux {

#if HAVE_UMFPACK
/*! \brief UMFPackBackend for MultiTypeMatrices
 *         Copies the coupled matrix into a single BCRSMatrix.
 *         Very slow!! Only wise to use for verification purposes.
 */
template <class TypeTag>
using MultiDimensionUMFPackBackend DUNE_DEPRECATED_MSG("Simply use UMFPackBackend as linear solver!")
    = UMFPackBackend<TypeTag>;

#endif // HAVE_UMFPACK

} // end namespace Dumux

#endif
