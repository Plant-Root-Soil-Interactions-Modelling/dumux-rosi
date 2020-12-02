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
 * \brief spatial parameters for the RichardsTestProblem
 */
#ifndef DUMUX_RICHARDS_TEST_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_TEST_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

// forward declaration
template<class TypeTag>
class RichardsTestSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(RichardsTestSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(RichardsTestSpatialParams, SpatialParams, Dumux::RichardsTestSpatialParams<TypeTag>);

// Set the material law
SET_PROP(RichardsTestSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef VanGenuchten<Scalar> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the RichardsTestProblem
 */
template<class TypeTag>
class RichardsTestSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    //! The parameters of the material law to be used
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief Constructor
     *
     * \param gridView The DUNE GridView representing the spatial
     *                 domain of the problem.
     */
    RichardsTestSpatialParams(const GridView& gridView)
        : ParentType(gridView)
    {
        permeability_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Permeability);
        porosity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Porosity);
        buffer_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.BufferPower);
        Swr_ = GET_RUNTIME_PARAM(TypeTag, Scalar, materialParams.Swr);
        VgAlpha_ = GET_RUNTIME_PARAM(TypeTag, Scalar, materialParams.VgAlpha);
        Vgn_ = GET_RUNTIME_PARAM(TypeTag, Scalar, materialParams.Vgn);

       // residual saturations
        materialParams_.setSwr(Swr_);
        materialParams_.setSnr(0.0);

        // parameters for the Van Genuchten law
        // alpha and n
        materialParams_.setVgAlpha(VgAlpha_);
        materialParams_.setVgn(Vgn_);
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 int scvIdx) const
    {
        return permeability_;
    }

    /*!
     * \brief Returns the porosity [] at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    { return porosity_; }
    /*!
     * \brief Returns the porosity [] at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    Scalar buffer(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    { return buffer_; }
    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvGeometry,
                                                int scvIdx) const
    {
        return materialLawParams(fvGeometry.subContVol[scvIdx].global);
    }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsTestProblem
     *
     * \param globalPos A global coordinate vector
     */
    const MaterialLawParams& materialLawParams(const GlobalPosition &globalPos) const
    {
        return materialParams_;
    }

    /*!
     * \brief Define the dispersivity.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    double dispersivity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        return 0;
    }


private:
    MaterialLawParams materialParams_;
    Scalar permeability_, VgAlpha_, Vgn_, Swr_, porosity_, buffer_;

};

} // end namespace

#endif
