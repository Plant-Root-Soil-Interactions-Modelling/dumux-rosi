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
#ifndef ROSI_ANALYTICAL_ACTIVE_UPTAKE_SINGLE_ROOT_HH
#define ROSI_ANALYTICAL_ACTIVE_UPTAKE_SINGLE_ROOT_HH


#include <cassert>
#include <boost/math/special_functions/expint.hpp>
#include <dumux/common/basicproperties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{

template <class TypeTag>
class AnalyticalSolutionROOSE2001
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

/*!
  * \brief This is a class to calculate the nutrient root active uptake
  *        based on aproximation analytical approaches (Roose 2001).
  */
public:
    AnalyticalSolutionROOSE2001(const Scalar Km, const Scalar Vmax,
                       const Scalar rootRadius, const Scalar Frac0, const Scalar time)
        : Km_(Km), Vmax_(Vmax), rootRadius_(rootRadius),Frac0_(Frac0), time_(time)
    {
        // check input whether be positive
        assert (rootRadius_ > 0);
        assert (Km_ > 0);
        assert (Vmax_ > 0);
        assert (Frac0_ > 0);
        assert (time_ >= 0);
    }

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
    const Scalar ActiveUptake()
    {
        Scalar activeUptake = 0.0;
        Cinf_ = Frac0_*density_;
        Cinf_dl  = Cinf_/Km_;
        Scalar lambda = Vmax_*rootRadius_/(effDiffCoeff_*Km_);
                lambda /=saturation_*porosity_+buffer_;
        Scalar L = lambda/2*log(4*exp(-0.577215)*effDiffCoeff_*pow(rootRadius_,(-2))*time_+1);
        activeUptake += -2*M_PI*rootRadius_*2*Vmax_*Cinf_dl/(1+Cinf_dl+L+sqrt(4*Cinf_dl+pow((1-Cinf_dl+L),2)));
    //    std::cout << "Cinf_ " << Cinf_ <<" "<<std::endl;
    //    std::cout << "Cinf_dl " << Cinf_dl <<" "<<std::endl;
    //    std::cout << "Vmax_ " << Vmax_ <<" "<<std::endl;
    //    std::cout << "rootRadius_ " << rootRadius_ <<" "<<std::endl;
    //    std::cout << "effDiffCoeff_ " << effDiffCoeff_ <<" "<<std::endl;
    //    std::cout << "Km_ " << Km_ <<" "<<std::endl;
    //    std::cout << "time_ " << time_ <<" "<<std::endl;
    //    std::cout << "saturation_ " << saturation_ <<" "<<std::endl;
    //    std::cout << "L " << L <<" "<<std::endl;
    //    std::cout << "lambda " << lambda <<" "<<std::endl;
    //    std::cout << "activeValue " << activeUptake<< M_PI <<" "<<std::endl;
        return activeUptake;
    }

private:
    Scalar Km_, Vmax_, density_, rootRadius_, Frac0_, time_, Cinf_, Cinf_dl, effDiffCoeff_, saturation_, porosity_, buffer_, activeUptake_;
};

} //end namespace

#endif // ROSI_ANALYTICAL_ACTIVE_UPTAKE_SINGLE_ROOT_HH
