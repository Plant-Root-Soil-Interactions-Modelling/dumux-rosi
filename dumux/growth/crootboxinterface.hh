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
 * \brief Some helper classes to interface crootbox
 */
#ifndef DUMUX_CROOTBOX_INTERFACE_HH
#define DUMUX_CROOTBOX_INTERFACE_HH

#include <mymath.h>
#include <dune/common/fvector.hh>

namespace Dumux {

namespace GrowthModule {

class CRootBoxInterface
{
    using GlobalPosition = Dune::FieldVector<double, 3>;
public:
    //! scaling: default assumes cm units in CRootBox and m units in Dumux
    static GlobalPosition convert(const ::Vector3d& v, double scale = 0.01)
    { return GlobalPosition({v.x*scale, v.y*scale, v.z*scale}); }

    static std::vector<unsigned int> convert(const ::Vector2i& v)
    { return {static_cast<unsigned int>(v.x), static_cast<unsigned int>(v.y)}; }

    //! conversion factor from seconds to days
    static constexpr double sToDays = 1.0/86400.0;
};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
