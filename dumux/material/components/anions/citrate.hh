/*
 * HCO3.hh
 *
 *  Created on: 15.09.2016
 *      Author: Trung Hieu Mai
 */

/*!
 * \file
 *
 * \brief A class for the HCO3 fluid properties
 */
#ifndef DUMUX_C6H5O7_HH
#define DUMUX_C6H5O7_HH

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/component.hh>

#include <cmath>
#include <iostream>


namespace Dumux
{
/*!
 * \brief A class for the C6H5O7 fluid properties
 https://pubchem.ncbi.nlm.nih.gov/compound/citrate#section=Top
 http://www.chemspider.com/Chemical-Structure.29081.html
 */
template <class Scalar>
class C6H5O7 : public Component<Scalar, C6H5O7<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the C6H5O7.
     */
    static const char *name()
    { return "Citrate"; }

    /*!
     * \brief The mass in [kg] of one mole of C6H5O7.
     */
    static Scalar molarMass()
    { return 189.101e-3; } // kg/mol

    /*!
     * \brief The diffusion Coefficient of C6H5O7 in water.
     http://www.unige.ch/cabe/dynamic/ESTDynamicPartISI.pdf
     */
    static Scalar liquidDiffCoeff()
    { return 6.23e-10; }

    static Scalar charge()
    {
        return -3.0;
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

