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
 * \ingroup OneTests
 * \brief Definition of the spatial parameters for the tissue problem
 */
#ifndef DUMUX_SOIL_SPATIAL_PARAMS_HH
#define DUMUX_SOIL_SPATIAL_PARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/spatialparams/fv.hh>
//#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Definition of the spatial parameters for the tissue problem
 */
template<class FVGridGeometry, class Scalar>
class SoilSpatialParams
: public FVSpatialParams<FVGridGeometry,Scalar, SoilSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = SoilSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;
    // the material law type
    using MaterialLaw = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
    // the material law params type
    using MaterialLawParams = typename MaterialLaw::Params;

    SoilSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // residual saturations
        materialParams_.setSwr(0.148);
        materialParams_.setSnr(0.0);

        // parameters for the Van Genuchten law
        // alpha and n
        materialParams_.setVgAlpha(2e-4);
        materialParams_.setVgn(1.41);

        // perm and poro
        permeability_ = getParam<Scalar>("Soil.SpatialParams.Permeability");
        porosity_ = getParam<Scalar>("Soil.SpatialParams.Porosity");

        double eps = 1.e-4;
        materialParams_.setPcLowSw(eps);
        materialParams_.setPcHighSw(1. - eps);
        materialParams_.setKrnLowSw(eps);
        materialParams_.setKrwHighSw(1 - eps);


        // output pc-sw and relperm curves
        if (getParam<bool>("Output.PlotVanGenuchtenCurves", false))
        {
            GnuplotInterface<double> gnuplot;
            gnuplot.resetPlot();
            gnuplot.setXlabel("log10(h [cm])");
            gnuplot.setYlabel("water content [-]");
            gnuplot.setOption("set y2label \"log10(hydr. cond. [cm/day])\"");

            constexpr int size = 300;
            std::vector<double> log10hcm(size);
            std::vector<double> watercontent(size);
            std::vector<double> log10conductivity(size);

            for (int i = 0; i <= size-1; ++i)
            {
                log10hcm[i] = double(i)/double(size-1)*4.0;
                const auto pc = std::pow(10, log10hcm[i])/100*9.81*1000;
                const auto sw = MaterialLaw::sw(materialParams_, pc);
                watercontent[i] = sw*porosity_;
                const auto krw = MaterialLaw::krw(materialParams_, sw);
                log10conductivity[i] = std::log10(krw*permeability_*9.81*1000/1e-3*100*86400);
            }

            gnuplot.setXRange(0, 4);
            gnuplot.setYRange(0, 0.5);
            gnuplot.addDataSetToPlot(log10hcm, watercontent, "watercontent", "with lines axes x1y1 lw 3");
            gnuplot.addDataSetToPlot(log10hcm, log10conductivity, "conductivity", "with lines axes x1y2 lw 3");
            gnuplot.setOption("set ytics nomirror");
            gnuplot.setOption("set y2tics");
            gnuplot.setOption("set y2range [-10 : 2.5]");
            gnuplot.setOption("set title \"Water retention curve and hydraulic conductivity\"");
            gnuplot.plot("vangenuchten");
        }
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol the element solution vector
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param element The current finite element
     * \param scv The sub control volume
     * \param elemSol The current element solution vector
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param globalPos A global coordinate vector
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition &globalPos) const
    { return materialParams_; }

    using ParentType::materialLawParams;
    const MaterialLawParams& materialLawParams() const
    { return materialParams_; }

private:
    MaterialLawParams materialParams_;
    Scalar permeability_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
