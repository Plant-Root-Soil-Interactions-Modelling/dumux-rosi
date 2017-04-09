/*
 *
 *  Created on: 04.04.2017
 *      Author: Trung Hieu Mai
 */

/*!
 * \file
 *
 * \brief A class for the Benzo[a]pyren fluid properties
 *	https://en.wikipedia.org/wiki/Benzo(a)pyrene
 */
#ifndef DUMUX_C20H12_HH
#define DUMUX_C20H12_HH

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/component.hh>

#include <cmath>
#include <iostream>


namespace Dumux
{
/*!
 * \brief A class for the Benzo(a)pyren fluid properties
 */
template <class Scalar>
class C20H12 : public Component<Scalar, C20H12<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the C20H12.
     */
    static const char *name()
    { return "benzo[a]pyren"; }

    /*!
     * \brief The mass in [kg] of one mole of C20H12.
     */
    static Scalar molarMass()
    { return 252.32e-3; } // kg/mol

    /*!
     * \brief The diffusion Coefficient of C20H12 in water.
     */
    static Scalar liquidDiffCoeff()
    { return 4.48e-10 ; } //https://publikationen.uni-tuebingen.de/xmlui/bitstream/handle/10900/48827/pdf/Diss_Henzler_TGA76.pdf?sequence=1

    static Scalar charge()
    {
        return 0.0;
    }
};

} // end namespace

#endif

