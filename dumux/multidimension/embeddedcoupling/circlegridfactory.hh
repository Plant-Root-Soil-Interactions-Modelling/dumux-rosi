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
 * \ingroup EmbeddedCoupling
 * \brief A factory circle shaped grids embedded in three-dimensional space
 */
#ifndef DUMUX_CIRCLE_GRID_FACTORY_HH
#define DUMUX_CIRCLE_GRID_FACTORY_HH

#include <cmath>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dumux/common/math.hh>

namespace Dumux {


/*!
 * \ingroup EmbeddedCoupling
 * \brief A factory circle shaped grids embedded in three-dimensional space
 */
template <class Grid>
class CircleGridFactory
{
    static const int dimworld = Grid::dimensionworld;
    typedef typename Grid::ctype Scalar;
    typedef Dune::FieldVector<Scalar, dimworld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimworld, dimworld> GlobalMatrix;
    static constexpr Scalar eps = 1.5e-7;

public:

    /* \brief Create a circle 1d in 3d grid
     * \param center the position of the center point
     * \param normal a vector normal to the circle plane with the length of the circle radius
     * \param radius the radius of the circle
     * \param num number of vertices/elements on the circle
     * \param offset the offset from 0 radians for the first vertex
     */
    static Grid * createGrid(const GlobalPosition& center,
                             const GlobalPosition& normal,
                             const Scalar radius,
                             const unsigned int num = 20,
                             const Scalar offset = 0)
    {
        static_assert(Grid::dimensionworld == 3, "Only implemented for world dimension 3");

        // the grid factory
        Dune::GridFactory<Grid> gridFactory;

        // make sure n is a unit vector
        auto n = normal;
        n /= n.two_norm();

        // caculate a vector u perpendicular to n
        GlobalPosition u;
        if (std::abs(n[0]) < eps && std::abs(n[1]) < eps)
            if (std::abs(n[2]) < eps)
                DUNE_THROW(Dune::MathError, "The normal vector has to be non-zero!");
            else
                u = {0, 1, 0};
        else
            u = {-n[1], n[0], 0};

        u *= radius/u.two_norm();

        // the circle parameterization is p(t) = r*cos(t)*u + r*sin(t)*(n x u) + c
        auto tangent = crossProduct(u, n);
        tangent *= radius/tangent.two_norm();

        // the parameter with an offset
        Scalar t = offset;
        // insert the vertices
        for (unsigned int i = 0; i < num; ++i)
        {
            const GlobalPosition v = {u[0]*std::cos(t) + tangent[0]*std::sin(t) + center[0],
                                      u[1]*std::cos(t) + tangent[1]*std::sin(t) + center[1],
                                      u[2]*std::cos(t) + tangent[2]*std::sin(t) + center[2]};
            gridFactory.insertVertex(v);

            t += 2*M_PI/num;

            // periodic t
            if(t > 2*M_PI)
                t -= 2*M_PI;
        }
        // insert the elements
        for (unsigned int i = 1; i < num; ++i)
            gridFactory.insertElement(Dune::GeometryType(1), {i-1, i});
        // close the circle
        gridFactory.insertElement(Dune::GeometryType(1), {num-1, 0});

        return gridFactory.createGrid();
    }
};

} // end namespace Dumux

#endif
