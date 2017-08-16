/*
 * solute.hh
 *
 *  Created on: 29.05.2017
 *      Author: Trung Hieu Mai
 */

/*!
 * \file
 *
 * \brief A class for the transported soluted properties
 */
#ifndef DUMUX_SOLUTE_HH
#define DUMUX_SOLUTE_HH

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/component.hh>

#include <cmath>
#include <iostream>


namespace Dumux
{
/*!
 * \brief A class for the fluid properties
 */
template<class TypeTag, class Scalar>
class SOLUTE : public Component<Scalar, SOLUTE<TypeTag, Scalar> >
{
public:
    /*!
     * \brief Name of the transported solute.
     */
    static const char *name()
    { return GET_PARAM_FROM_GROUP(TypeTag, std::string, Solute, Name); }

    /*!
     * \brief The mass in [kg] of one mole of solute.
     */
    static Scalar molarMass()
    { return GET_PARAM_FROM_GROUP(TypeTag, Scalar, Solute, MolarMass); } // kg/mol

    /*!
     * \brief The diffusion Coefficient in water.
     */
    static Scalar liquidDiffCoeff()
    { return GET_PARAM_FROM_GROUP(TypeTag, Scalar, Solute, liquidDiffCoeff); }
};

} // end namespace

#endif

