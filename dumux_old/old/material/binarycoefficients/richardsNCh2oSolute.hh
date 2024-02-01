// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and heavy oil.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_SOLUTE_HH
#define DUMUX_BINARY_COEFF_H2O_SOLUTE_HH

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/richardsNCsolute.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and solute as in SAGD processes
 */
template<class Scalar, class Component>
class H2O_SOLUTE
{
    H2O_SOLUTE()
    {
        DUNE_THROW(Dune::NotImplemented, "The binary coefficients for H2O and your "
                     << "component are not implemented! Please implement the needed specialization.");
    }
};

template<class Scalar, int id>
class H2O_SOLUTE<Scalar, Components::Constant<id, Scalar>>
{
public:
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        static const Scalar D = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidDiffusionCoefficient", 1.0);
        return D;
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
