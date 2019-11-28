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
 * \author Timo Koch <timo.koch@iws.uni-stuttgart.de>
 * \author Timo Koch
 * \ingroup GridGrowth
 * \brief Helper class to transform a network into the bounds of a periodic domain cell
 */
#ifndef DUMUX_PERIODIC_NETWORK_TRANSFORM_HH
#define DUMUX_PERIODIC_NETWORK_TRANSFORM_HH

#include <bitset>
#include <string>
#include <cmath>

#include <dune/common/fvector.hh>

namespace Dumux {

template<class GlobalCoordinate>
class PeriodicNetworkTransform
{
    static constexpr int dimWorld = GlobalCoordinate::dimension;
    using VertexMarker = Dune::FieldVector<int, dimWorld>;
public:
    /*!
     * \brief A helper class to do a periodic transform of a network grid
     * \param lowerLeft the lower left corner of the bounding box the grid should be restricted to
     * \param upperRight the upper right corner of the bounding box the grid should be restricted to
     * \param periodic if the grid is periodic in x and/or y and/or z
     * \param paramGroup the parameter group in which to preferably look up parameters
     */
    PeriodicNetworkTransform(const GlobalCoordinate& lowerLeft,
                             const GlobalCoordinate& upperRight,
                             const std::bitset<dimWorld>& periodic,
                             const std::string& paramGroup = "")
    : lowerLeft_(lowerLeft)
    , upperRight_(upperRight)
    , periodic_(periodic)
    {
        // ignore the bounds in the non-periodic directions
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            if (!periodic_[dimIdx])
            {
                upperRight_[dimIdx] = 1e100;
                lowerLeft_[dimIdx] = -1e100;
            }
        }

        shift_ = upperRight_ - lowerLeft_;
        eps_ = shift_; eps_ *= 1e-7;
    }

    //! check if a point is on the periodic boundary
    bool onBoundary(const GlobalCoordinate& point)
    {
        using std::abs;
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            if (periodic_[dimIdx])
                if (abs(point[dimIdx] - lowerLeft_[dimIdx]) < eps_[dimIdx])
                    return true;
        return false;
    }

    //! create a vertex marker denoting the shift factor in each periodic direction
    VertexMarker getVertexMarker(const GlobalCoordinate& vertexPos)
    {
        auto pos = vertexPos - lowerLeft_; // make lower left the origin

        using std::floor;
        VertexMarker marker(0);
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            if (periodic_[dimIdx])
                marker[dimIdx] = static_cast<int>(floor(pos[dimIdx]/shift_[dimIdx]));

        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            if (!periodic_[dimIdx] && marker[dimIdx] != 0)
                DUNE_THROW(Dune::InvalidStateException, "Computed non-zero marker ("
                           << marker[dimIdx] << ") in non-periodic direction (" << dimIdx << ")");

        return marker;
    }

    //! create the shift vector from the vertex marker
    GlobalCoordinate getShift(const VertexMarker& v)
    {
        auto shift = shift_;
        for (int i = 0; i < dimWorld; ++i)
            shift[i] *= v[i];
        return shift;
    }

    //! get periodic box (the boundary box of the shifted periodic cells) given a vertex marker
    std::pair<GlobalCoordinate, GlobalCoordinate> getPeriodicBBox(const VertexMarker& v)
    {
        auto lowerLeft = lowerLeft_ + getShift(v);
        const auto upperRight = lowerLeft + shift_;
        return {lowerLeft, upperRight};
    }

    //! intersect a directional ray with a periodic aabb
    GlobalCoordinate intersectRayBox(const GlobalCoordinate& lowerLeft, const GlobalCoordinate& upperRight,
                                     const GlobalCoordinate& origin, const GlobalCoordinate& dir) const
    {
        using std::signbit;
        auto t = signbit(dir[0]) ? (lowerLeft[0]-origin[0])/dir[0] : (upperRight[0]-origin[0])/dir[0];
        const auto ty = signbit(dir[1]) ? (lowerLeft[1]-origin[1])/dir[1] : (upperRight[1]-origin[1])/dir[1];
        const auto tz = signbit(dir[2]) ? (lowerLeft[2]-origin[2])/dir[2] : (upperRight[2]-origin[2])/dir[2];

        using std::min;
        t = min({t, ty, tz});

        if (t <= 0)
        {
            const auto tx = signbit(dir[0]) ? (lowerLeft[0]-origin[0])/dir[0] : (upperRight[0]-origin[0])/dir[0];
            DUNE_THROW(Dune::InvalidStateException, "parameter has to be greater than zero! "
                       << "tx: " << tx << ", ty: " << ty << ", tz: " << tz);
        }

        auto iPos = origin;
        iPos.axpy(t, dir);
        return iPos;
    }

    //! return bitset showing which axis directions are periodic
    const std::bitset<dimWorld>& periodic() const
    { return periodic_; }

    GlobalCoordinate lowerLeft_;
    GlobalCoordinate upperRight_;
    GlobalCoordinate eps_, shift_;
    std::bitset<dimWorld> periodic_;
};

} // end namespace Dumux

#endif
