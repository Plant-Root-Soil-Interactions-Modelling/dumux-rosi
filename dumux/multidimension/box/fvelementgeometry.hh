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
 * \brief Represents the finite volume geometry of a single element in
 *        the box scheme.
 *        \note As Dumux will be is in a transition regarding the structure
 *        of the FvElementGeometry this is temporary until a new interface
 *        is agreed on.
 */
#ifndef DUMUX_MULTIDIMENSION_BOX_FV_ELEMENTGEOMETRY_HH
#define DUMUX_MULTIDIMENSION_BOX_FV_ELEMENTGEOMETRY_HH

#include <dumux/implicit/box/fvelementgeometry.hh>

namespace Dumux
{

template<class TypeTag>
class MultiDimensionBoxFVElementGeometry : public BoxFVElementGeometry<TypeTag>
{
    typedef Dumux::BoxFVElementGeometry<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
public:
    void update(const GridView& gridView, const Element& element, const Problem& problem)
    {
        ParentType::update(gridView, element);
    }

    void update(const GridView& gridView, const Element& element)
    {
        ParentType::update(gridView, element);
    }
};

}

#endif
