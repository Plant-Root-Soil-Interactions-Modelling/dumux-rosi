/*
 * k.hh
 *
 *  Created on: 11.07.2016
 *      Author: Trung Hieu Mai
 */

/*!
 * \file
 *
 * \brief A class for the potassium K nutrient properties
 */
#ifndef DUMUX_K_HH
#define DUMUX_K_HH

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
class K : public Component<Scalar, K<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the HPO4.
     *\ http://www.endmemo.com/chem/compound/hpo4.php
     */
    static const char *name()
    { return "K"; }

    /*!
     * \brief The mass in [kg] of one mole of NO3.
     */
    static Scalar molarMass()
    { return 39.0983e-3; } // kg/mol

    /*!
     * \brief The diffusion Coefficient of HPO4 in water.
     * \http://www.unige.ch/cabe/dynamic/ESTDynamicPartISI.pdf
     */
    static Scalar liquidDiffCoeff()
    { return 19.57e-10; }

    static Scalar charge()
    {
        return 1.0;
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

