/*
 * RichardsNCsolute.hh
 *
 *  Created on: 23.10.2019
 *      Author: Deepanshu Khare
 */

/*!
 * \file
 *
 * \brief A class for the transported soluted properties
 */
#ifndef DUMUX_RICHARDS_NC_SOLUTE_HH
#define DUMUX_RICHARDS_NC_SOLUTE_HH

#include <dumux/common/exceptions.hh>
//#include <dumux/material/components/component.hh>

#include <cmath>
#include <iostream>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux {
namespace Components {

template <class Scalar>
class SOLUTE
: public Components::Base<Scalar, SOLUTE<Scalar> >
, public Components::Liquid<Scalar, SOLUTE<Scalar> >
, public Components::Gas<Scalar, SOLUTE<Scalar> >
{
public:
    /*!
     * \brief Name of the transported solute.
     */
    static std::string name()
    { return getParam<std::string>("Component.Name"); }

    /*!
     * \brief The mass in [kg] of one mole of solute.
     */
    constexpr static Scalar molarMass()
    { return getParam<Scalar>("Component.MolarMass"); } // kg/mol

    /*!
     * \brief The diffusion Coefficient in water.
     */
    /*static Scalar liquidDiffCoeff()
    { return getParam<Scalar>("Component.liquidDiffCoeff"); }*/
};

} //end namespace Components

} // end namespace Dumux

#endif

