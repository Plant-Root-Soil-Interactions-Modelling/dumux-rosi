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
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_ROOTSYSTEM_TEST_SPATIALPARAMS_HH
#define DUMUX_ROOTSYSTEM_TEST_SPATIALPARAMS_HH

//#include "rootspatialparams.hh"
#include <dumux/porousmediumflow/1d/rootsystem/spatialparamsDGF.hh>

namespace Dumux
{

/*!
 * \ingroup RootsystemBoxModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        rootsystem box model
 */
template<class TypeTag>
class RootsystemTestSpatialParams : public RootsystemSpatialParamsDGF<TypeTag>
{
    typedef RootsystemSpatialParamsDGF<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

public:
    RootsystemTestSpatialParams(const GridView& gridView)
        : ParentType(gridView)
    {
        Vmax_ = tree.template get<Scalar>("RootSystem.SpatialParams.Vmax", 0);
        Km_ = tree.template get<Scalar>("RootSystem.SpatialParams.Km", 0);
        Cmin_ = tree.template get<Scalar>("RootSystem.SpatialParams.Cmin", 0);
        PartitionCoeff_ = tree.template get<Scalar>("RootSystem.SpatialParams.PartitionCoeff");
        //rootRadius_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.rootRadius);
    }

    //  Scalar rootRadius(const Element &element,
    //                               const FVElementGeometry &fvGeometry,
    //                               const int scvIdx) const
    //  {
    //      return rootRadius_;
    //  }
//
    //  Scalar radius(int idx) const
    //  {
    //      return rootRadius_;
    //}

    Scalar PartitionCoeff() const
    {
        return PartitionCoeff_;
    }

    Scalar Vmax() const
    {
        return Vmax_;
    }

    Scalar Km() const
    {
        return Km_;
    }

    Scalar Cmin() const
    {
        return Cmin_;
    }

    Scalar porosity(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 const int scvIdx) const
    {
        return 1;
    }
    Scalar dispersivity(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 const int scvIdx) const
    {
        return 0;
    }
private:

    Scalar Vmax_, Km_, Cmin_, PartitionCoeff_, rootRadius_;
};

} // end namespace

#endif
