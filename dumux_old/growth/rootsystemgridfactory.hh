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
 * \brief A grid factory using a CRootBox root system
 */
#ifndef DUMUX_ROOTSYSTEM_GRIDFACTORY_HH
#define DUMUX_ROOTSYSTEM_GRIDFACTORY_HH

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <RootSystem.h>

namespace Dumux {

namespace GrowthModule {

/**
 * Builds a grid (Dune::FoamGrid<1, 3>) from a root system (CPlantBox::RootSystem)
 *
 * use static member function: RootSystemGridFactory::makeGrid(RootSystem)
 */
class RootSystemGridFactory
{
    static constexpr int dim = 1;
    static constexpr int dimworld = 3;
    using GlobalPosition = Dune::FieldVector<double, dimworld>;

public:
    //! export grid type
    using Grid = Dune::FoamGrid<1, 3>;

    //! make the grid from the initial root system
    static std::shared_ptr<Grid> makeGrid(const CPlantBox::RootSystem& rs, double shootZ = 0., bool verbose = false)
    {
        // the grid factory creates the grid
        if (verbose) std::cout << "RootSystemGridFactory: " << std::endl;
        Dune::GridFactory<Grid> factory;
        constexpr auto line = Dune::GeometryTypes::line;

        auto nodes = rs.getNodes();
        nodes[0].z = shootZ; // optionally, fix depth of first node
        int counter = 0;
        for (const auto& n : nodes) {
            if (verbose) std::cout << "-- add vertex " << counter++ <<  " at " << n.toString() << std::endl;
            factory.insertVertex(convert_(n));
        }

        const auto shootSegments = rs.getShootSegments();
        for (const auto& s : shootSegments) {
            if (verbose) std::cout << "-- add element with vertices " << s.toString() << std::endl;
            factory.insertElement(line, convert_(s));
        }

        const auto segments = rs.getSegments();
        for (const auto& s : segments) {
            if (verbose) std::cout << "-- add element with vertices " << s.toString() << std::endl;
            factory.insertElement(line, convert_(s));
        }

        return std::shared_ptr<Grid>(factory.createGrid());
    }

private:

    static GlobalPosition convert_(const CPlantBox::Vector3d& v, double scale = 0.01) {
        return GlobalPosition( { v.x * scale, v.y * scale, v.z * scale });
    }

    static std::vector<unsigned int> convert_(const CPlantBox::Vector2i& v) {
        return {static_cast<unsigned int>(v.x), static_cast<unsigned int>(v.y)};
    }

};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
