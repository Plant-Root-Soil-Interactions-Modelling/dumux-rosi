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
#ifndef DUMUX_RICHARDS_BUFFER_TEST_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_BUFFER_TEST_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
//#include "csvReader.hh"
namespace Dumux
{

// forward declaration
template<class TypeTag>
class RichardsBufferTestSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(RichardsBufferTestSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(RichardsBufferTestSpatialParams, SpatialParams, Dumux::RichardsBufferTestSpatialParams<TypeTag>);

// Set the material law
SET_PROP(RichardsBufferTestSpatialParams, MaterialLaw)
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
class RichardsBufferTestSpatialParams : public ImplicitSpatialParams<TypeTag>
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

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

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
    RichardsBufferTestSpatialParams(const GridView& gridView)
        : ParentType(gridView)
    {
        //useSoilData_ = GET_RUNTIME_PARAM(TypeTag, bool, SpatialParams.UseSoilData);
        useSoilData_ = tree.template get<bool>("Soil.SpatialParams.UseSoilData", false);
        if (useSoilData_)
        {
            mu_ = 0.001;//WP::viscosity(temp,pnRef); // h2o: Dynamic viscosity of water 0.001 Pa·s (independent of temp and p)
            //std::cout<< "SpatialParams"<<"\n";
            //soilData_.getDataFromFile(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, SoilDataFileName));
            soilData_.getDataFromFile(tree.template get<std::string>("Soil.SpatialParams.SoilDataFileName"));
            soilData_.setDepthHeader("Depth[m]");
            std::vector<Scalar> soilWaterSaturation = soilData_.getValueInHeader("thetaS[-]");
            std::vector<Scalar> soilWaterResidual = soilData_.getValueInHeader("thetaR[-]");
            std::vector<Scalar> n = soilData_.getValueInHeader("n[-]");
            std::vector<Scalar> alpha = soilData_.getValueInHeader("alpha[1/m]");

            for (int layerNo = 0; layerNo < soilData_.getTotalSteps(); layerNo++)
            {
                materialParams_.push_back(MaterialLawParams());

                materialParams_[layerNo].setSwr(soilWaterResidual[layerNo]/soilWaterSaturation[layerNo]);
                materialParams_[layerNo].setSnr(0.0);
                materialParams_[layerNo].setVgAlpha(alpha[layerNo] /(1000*9.81)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
                materialParams_[layerNo].setVgn(n[layerNo]);
            }
        }
        else
	   {
            permeability_ = tree.template get<Scalar>("Soil.SpatialParams.Permeability");
            porosity_ = tree.template get<Scalar>("Soil.SpatialParams.Porosity");
            Swr_ = tree.template get<Scalar>("Soil.VanGenuchten.Swr");
            VgAlpha_ = tree.template get<Scalar>("Soil.VanGenuchten.VgAlpha");
            Vgn_ = tree.template get<Scalar>("Soil.VanGenuchten.Vgn");
            buffer_ = tree.template get<Scalar>("Soil.SpatialParams.BufferPower",0);
            bulkDensity_ = tree.template get<Scalar>("Soil.SpatialParams.BulkDensity", 1.4e3);
            FreundlichK_ = tree.template get<Scalar>("Soil.SpatialParams.FreundlichK",0);
            FreundlichN_ = tree.template get<Scalar>("Soil.SpatialParams.FreundlichN",1);
            //std::cout<< "1d spatial params"<< permeability_<<" "<< porosity_ <<" "<<Swr_ <<" "<<VgAlpha_ <<" "<<Vgn_ <<" "<<buffer_<<"\n";
           // residual saturations
            materialParams_.push_back(MaterialLawParams());
            materialParams_[0].setSwr(Swr_);
            materialParams_[0].setSnr(0.0);

            // parameters for the Van Genuchten law
            // alpha and n
            materialParams_[0].setVgAlpha(VgAlpha_);
            materialParams_[0].setVgn(Vgn_);
        }
    }

    /*!
     * \brief Returns the buffer power at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    Scalar buffer(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 int scvIdx) const
    {   if (useSoilData_)
        {
            GlobalPosition globalPos = element.geometry().center();
            return soilData_.getValueAtDepth("bufferCapacity[-]",-globalPos[2]);
        }
        else
            return buffer_;
    }

    Scalar FreundlichK(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 int scvIdx) const
    {   if (useSoilData_)
        {
            GlobalPosition globalPos = element.geometry().center();
            return soilData_.getValueAtDepth("FreundlichK",-globalPos[2]); // should be carefull with the unit
        }
        else
            return FreundlichK_;
    }

    Scalar FreundlichN(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 int scvIdx) const
    {   if (useSoilData_)
        {
            GlobalPosition globalPos = element.geometry().center();
            return soilData_.getValueAtDepth("FreundlichN[-]",-globalPos[2]);
        }
        else
            return FreundlichN_;
    }

    Scalar bulkDensity(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 int scvIdx) const
    {   if (useSoilData_)
        {
            GlobalPosition globalPos = element.geometry().center();
            return soilData_.getValueAtDepth("BulkDensity[gcm-3]",-globalPos[2])*1e3;
        }
        else
            return bulkDensity_;
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
    {   if (useSoilData_)
        {
            //GlobalPosition globalPos = element.geometry().center();
            //return soilData_.getValueAtDepth("Ksat",-globalPos[2])*mu_/(1000*9.81);
            return soilData_.getValueAtDepth("Ksat[m/day]",-element.geometry().center()[2])*mu_/(1000*9.81);
        }
        else
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
    {   if (useSoilData_)
        {
            //GlobalPosition globalPos = element.geometry().center();
            //return soilData_.getValueAtDepth("soilWaterSaturation",-globalPos[2]);
            return soilData_.getValueAtDepth("thetaS[-]",-element.geometry().center()[2]);
        }
        else
            return porosity_;

    }
    Scalar porosity(const GlobalPosition &globalPos) const
    {   if (useSoilData_)
        {
            //GlobalPosition globalPos = element.geometry().center();
            //return soilData_.getValueAtDepth("soilWaterSaturation",-globalPos[2]);
            return soilData_.getValueAtDepth("thetaS[-]",-globalPos[2]);
        }
        else
            return porosity_;
    }
    Scalar VgAlpha(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {   if (useSoilData_)
		return soilData_.getValueAtDepth("alpha[1/m]",-element.geometry().center()[2])/(1000*9.81);
        else
            return VgAlpha_;

    }
        Scalar Vgn(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {   if (useSoilData_)
		return soilData_.getValueAtDepth("n[-]",-element.geometry().center()[2]);
        else
            return Vgn_;

    }
        Scalar Swr(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {   if (useSoilData_)
		return soilData_.getValueAtDepth("thetaR[-]",-element.geometry().center()[2])/
                    soilData_.getValueAtDepth("thetaS[-]",-element.geometry().center()[2]);
        else
            return Swr_;

    }
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
        if (useSoilData_)
        {
            return materialParams_[soilData_.getDepthIndex(-fvGeometry.subContVol[scvIdx].global[2])];
        }
        else
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
        if (useSoilData_)
        {
            return materialParams_[soilData_.getDepthIndex(-globalPos[2])];
        }
        else
            return materialParams_[0];
    }

    Scalar initialMassFraction(const GlobalPosition &globalPos) const
    {
        if (useSoilData_)
        {
            return soilData_.getValueAtDepth("NutrientConcentration[mg/l]",-globalPos[2])*1e-6;
            ///(1 + soilData_.getValueAtDepth("bufferCapacity",-globalPos[2]));
        }
        else
            return -999;
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
    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */

private:
    Scalar permeability_, VgAlpha_, Vgn_, Swr_, porosity_, buffer_, mu_, initialSoluteMassFracInSoil_;
    Scalar bulkDensity_, FreundlichN_, FreundlichK_;
std::vector<MaterialLawParams> materialParams_;
    CSVReader soilData_;
    bool useSoilData_;
};

} // end namespace

#endif
