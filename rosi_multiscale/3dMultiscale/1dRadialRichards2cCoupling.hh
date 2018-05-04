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
  * \brief This is a class to calculate the nutrient root active uptake
  *        based on aproximation analytical approaches (Roose 2001).
  */
#ifndef ROSI_1D_RICHARDS_2C_SINGLE_ROOT_HH
#define ROSI_1D_RICHARDS_2C_SINGLE_ROOT_HH

#include <iomanip> // setprecision
#include <sstream>
#include <cassert>
#include <boost/math/special_functions/expint.hpp>
#include <dumux/common/basicproperties.hh>
#include <dumux/common/parameters.hh>
#include "1dRichards2cTestProblem.hh"
#include "1dRadialStart.hh"

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 20)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

namespace Dumux
{

template <class TypeTag>
class OneDRichards2CRadiallySymmetric
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

/*!
  * \brief This is a class to calculate the nutrient root active uptake
  *        based on aproximation analytical approaches (Roose 2001).
  */
public:
    OneDRichards2CRadiallySymmetric()
    {}

    void setVmax(Scalar Vmax)
    { Vmax_ = Vmax; }

    void setKm(Scalar Km)
    { Km_ = Km; }

    void setRootRadius(Scalar rootRadius)
    { rootRadius_ = rootRadius; }
    //! \brief Sets the ... \f$[m]\f$.
    void setEffDiffCoeff(Scalar effDiffCoeff)
    { effDiffCoeff_ = effDiffCoeff; }

    void setSaturation(Scalar saturation)
    { saturation_ = saturation; }

    void setPorosity(Scalar porosity)
    { porosity_ = porosity; }

    void setDensity(Scalar density)
    { density_ = density; }

    void setBuffer(Scalar buffer)
    { buffer_ = buffer; }

    /*!
     * \brief Returns the ActiveUptake
     */
    const Scalar ActiveUptakeAnalyticalROOSE2001()
    {
        Scalar time=restartTime_;
        Scalar activeUptake = 0.0;
        Cinf_ = initialSoilNutrient_*density_;
        Cinf_dl  = Cinf_/Km_;
        Scalar lambda = Vmax_*rootRadius_/(effDiffCoeff_*Km_);
                lambda /=saturation_*porosity_+buffer_;
        Scalar L = lambda/2*log(4*exp(-0.577215)*effDiffCoeff_*pow(rootRadius_,(-2))*time+1);
        activeUptake += -2*M_PI*rootRadius_*2*Vmax_*Cinf_dl/(1+Cinf_dl+L+sqrt(4*Cinf_dl+pow((1-Cinf_dl+L),2)));
        //std::cout<< "activeUptake" << activeUptake <<"\n";
        return activeUptake;
    }


    /*!
     * \brief Returns the MassFractionAtRootSureface based on analytical solution in Roose 2001
     */
    const Scalar MassFractionAtRootSurfaceAnalyticalROOSE2001()
    {
        Scalar time=restartTime_;
        Scalar massFractionAtRootSurface = 0.0;
        Cinf_ = initialSoilNutrient_*density_;
        Cinf_dl  = Cinf_/Km_;
        Scalar lambda = Vmax_*rootRadius_/(effDiffCoeff_*Km_);
                lambda /=saturation_*porosity_+buffer_;
        Scalar L = lambda/2*log(4*exp(-0.577215)*effDiffCoeff_*pow(rootRadius_,(-2))*time_+1);
        massFractionAtRootSurface += (Cinf_-Cinf_*lambda/(1+Cinf_dl+L+sqrt(4*Cinf_dl+pow((1-Cinf_dl+L),2)))*
                                                boost::math::expint(1,pow(rootRadius_,2)/(4*effDiffCoeff_*time_)))
                                                /density_;
        return massFractionAtRootSurface;
    }


    void setName(std::string name)
    { name_ = name; }

    void setRestartTime(Scalar RestartTime)
    { restartTime_ = RestartTime; }

    void setTEnd(Scalar TEnd)
    { tEnd_ = TEnd; }

    void setDtInitial(Scalar dtInitial)
    { dtInitial_ = dtInitial; }

    void setOuterBoundaryGrid(Scalar outerBoundary)
    { outerBoundary_ = outerBoundary; }

    void setZ(Scalar Z)
    { Z_ = Z; }

    void setVgAlpha(Scalar vgAlpha)
    { vgAlpha_ = vgAlpha; }

    void setVgn(Scalar vgn)
    { vgn_ = vgn; }

    void setSwr(Scalar swr)
    { swr_ = swr; }

    void setPermeability(Scalar permeability)
    { permeability_ = permeability; }

    void setKr(Scalar kr)
    { kr_ = kr; }

    void setDiffCoeffRootMembrane(Scalar DiffCoeffRootMembrane)
    { DiffCoeffRootMembrane_ = DiffCoeffRootMembrane; }

    void setPartitionCoeff(Scalar partitionCoeff)
    { partitionCoeff_ = partitionCoeff; }

    void setInitialSoilPressure(Scalar initialSoilPressure)
    { initialSoilPressure_ = initialSoilPressure;}

    void setInitialSoilNutrient(Scalar initialSoilNutrient)
    { initialSoilNutrient_ = initialSoilNutrient;}

    void setProblemName(std::string outputName)
    { outputName_ = outputName;}

    void setOuterFluxBoundaryConditions(PrimaryVariables fluxBoundaryCondition)
    { fluxBoundaryCondition_ = fluxBoundaryCondition;}

    void setInsideRootConditions(PrimaryVariables rootPriVars)
    { rootPriVars_ = rootPriVars;}

    /*!
     * \brief Provides an interface for customizing error messages associated with
     *        reading in parameters.
     *
     * \param progName  The name of the program, that was tried to be started.
     * \param errorMsg  The error message that was issued by the start function.
     *                  Comprises the thing that went wrong and a general help message.
     */
    void usage(const char *progName, const std::string &errorMsg)
    {
        // TODO
    }
    /*!
     * \brief Returns the ActiveUptakeNumerical
     */
    const auto ActiveUptakeNumerical()
    {
        std::vector<std::string> argv;
        argv.push_back("./"+name_);
        argv.push_back("-ParameterFile");
        argv.push_back(name_);
        argv.push_back("-radialProblemName");
        argv.push_back(outputName_);//0
        argv.push_back("-radialTimeManager.Restart");//1
        argv.push_back(to_string_with_precision(restartTime_));//2
        argv.push_back("-radialTimeManager.TEnd");
        argv.push_back(to_string_with_precision(tEnd_));
        argv.push_back("-radialTimeManager.DtInitial");
        argv.push_back(to_string_with_precision(dtInitial_));
        argv.push_back("-radialTimeManager.EpisodeTime");
        argv.push_back(to_string_with_precision(tEnd_));
        argv.push_back("-radialSoilGrid.Positions0");
        argv.push_back(to_string_with_precision(rootRadius_)+" "+to_string_with_precision(outerBoundary_)); //14
        std::cout<<"!!!! Coupling "<<to_string_with_precision(rootRadius_)+" "+to_string_with_precision(outerBoundary_)<<"\n";
        argv.push_back("-BoundaryConditions.OuterFluxesWater");
        argv.push_back(to_string_with_precision(fluxBoundaryCondition_[0])); //16
        argv.push_back("-BoundaryConditions.OuterFluxesNutrient");
        argv.push_back(to_string_with_precision(fluxBoundaryCondition_[1]));
        argv.push_back("-BoundaryConditions.InitialRootPressure");
        argv.push_back(to_string_with_precision(rootPriVars_[0])); //20
        argv.push_back("-BoundaryConditions.InitialSoluteMassFracInRoot");
        argv.push_back(to_string_with_precision(rootPriVars_[1]));
        argv.push_back("-BoundaryConditions.RootPressure");
        argv.push_back(to_string_with_precision(rootPriVars_[0]));
        argv.push_back("-materialParams.VgAlpha");
        argv.push_back(to_string_with_precision(vgAlpha_));
        argv.push_back("-materialParams.Vgn");
        argv.push_back(to_string_with_precision(vgn_));
        argv.push_back("-materialParams.Swr");
        argv.push_back(to_string_with_precision(swr_));
        argv.push_back("-SpatialParams.Porosity");
        argv.push_back(to_string_with_precision(porosity_));
        argv.push_back("-SpatialParams.Permeability");
        argv.push_back(to_string_with_precision(permeability_));
        argv.push_back("-SpatialParams.BufferPower");
        argv.push_back(to_string_with_precision(buffer_));
        argv.push_back("-SpatialParams.Kr");
        argv.push_back(to_string_with_precision(kr_));
        argv.push_back("-SpatialParams.Vmax");
        argv.push_back(to_string_with_precision(Vmax_));
        argv.push_back("-SpatialParams.Km");
        argv.push_back(to_string_with_precision(Km_));
        argv.push_back("-SpatialParams.rootRadius");
        argv.push_back(to_string_with_precision(rootRadius_));
        argv.push_back("-SpatialParams.DiffCoeffRootMembrane");
        argv.push_back(to_string_with_precision(DiffCoeffRootMembrane_));
        argv.push_back("-SpatialParams.PartitionCoeff");
        argv.push_back(to_string_with_precision(partitionCoeff_));
        argv.push_back("-BoundaryConditions.InitialSoluteMassFracInSoil");
        argv.push_back(to_string_with_precision(initialSoilNutrient_));
        argv.push_back("-BoundaryConditions.InitialSoilPressure");
        argv.push_back(to_string_with_precision(initialSoilPressure_)); //52
        argv.push_back("-radialSoilGrid.Z");
        argv.push_back(to_string_with_precision(Z_));

        char * argvector[argv.size()+1];
        for (int i = 0; i < argv.size(); ++i)
          argvector[i] = &argv[i][0];

        typedef TTAG(SoilRichardsTwoCBufferRadiallySymmetricTestCCProblem) ProblemTypeTag;
        return Dumux::radialStart<ProblemTypeTag>(argv.size(), argvector);
    }

private:
    Scalar Km_, Vmax_, density_, rootRadius_, time_, Cinf_, Cinf_dl, effDiffCoeff_, saturation_, porosity_, buffer_, activeUptake_;
    Scalar restartTime_, tEnd_, dtInitial_, outerBoundary_, Z_, vgAlpha_, vgn_, swr_, permeability_, kr_, partitionCoeff_, initialSoilPressure_, initialSoilNutrient_;
    std::string name_,outputName_;
    PrimaryVariables fluxBoundaryCondition_, rootPriVars_;
    Scalar DiffCoeffRootMembrane_;

};
} //end namespace

#endif // ROSI_ANALYTICAL_ACTIVE_UPTAKE_SINGLE_ROOT_HH
