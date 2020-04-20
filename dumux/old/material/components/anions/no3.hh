/*
 * HCO3.hh
 *
 *  Created on: 14.03.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the HCO3 fluid properties
 */
#ifndef DUMUX_NO3_HH
#define DUMUX_NO3_HH

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/component.hh>

#include <cmath>
#include <iostream>


namespace Dumux
{
/*!
 * \brief A class for the NO3 fluid properties
 */
template <class Scalar>
class NO3 : public Component<Scalar, NO3<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the NO3.
     */
    static const char *name()
    { return "NO3"; }

    /*!
     * \brief The mass in [kg] of one mole of NO3.
     */
    static Scalar molarMass()
    { return 62.0049e-3; } // kg/mol

    /*!
     * \brief The diffusion Coefficient of NO3 in water.
     */
    static Scalar liquidDiffCoeff()
    { return 1.7e-9; }

    static Scalar charge()
    {
        return -1.0;
    }

//    /*!
//     * \brief Return the constant ai Parkhurst (1990) for the modified Debye-Hückel equation
//     */
//
//    static Scalar ai()
//    {
//        return 0;
//    }
//
//    /*!
//     * \brief Return the constant bi Parkhurst (1990) for the modified Debye-Hückel equation
//     */
//    static Scalar bi()
//    {
//        return 0.0;
//    }
};

} // end namespace

#endif

