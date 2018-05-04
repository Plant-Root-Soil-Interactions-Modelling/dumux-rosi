#ifndef DUMUX_MPFAHELPER_HH
#define DUMUX_MPFAHELPER_HH

#include <dumux/common/math.hh>
#include "mpfaproperties.hh"
#include "mpfapropertydefaults.hh"

// In the following lines all the implemented interactionvolumes and interactionvolumfillers
// should be included. If a new combination is created the respective headers have to be included
// here and a corresponding specialization of the helper class needs to be implemented

#include "mpfao/2d/mpfao2dinteractionvolume.hh"
#include "mpfao/2d/mpfao2dmanager.hh"
#include "mpfao/2d/mpfao2dfluxcalculator.hh"
#include "mpfao/mpfaointeractionvolumecontainer.hh"

#include "mpfal/2d/mpfal2dinteractionvolume.hh"
#include "mpfal/2d/mpfal2dmanager.hh"
#include "mpfal/2d/mpfal2dfluxcalculator.hh"
#include "mpfal/mpfalinteractionvolumecontainer.hh"

namespace Dumux
{

//!
/*! \brief helper class to extract the necessary types needed for a specific MPFA method and dimension used.
 *
 * Types for the interactionvolume and the interactionvolume-filler class are specified.
 *
 */

template<class TypeTag, int dimension, typename Method = void>
class MpfaHelper
{};

// specialization of class for mpfa-o method and dimension 2
template<class TypeTag>
class MpfaHelper<TypeTag, 2,
    typename std::enable_if<GET_PROP_VALUE(TypeTag, MpfaMethod) == GET_PROP(TypeTag, MpfaMethods)::oMethod>::type >
{
public:
    typedef MpfaO2DInteractionVolume<TypeTag> InteractionVolume;
    typedef MpfaO2DManager<TypeTag> InteractionVolumeManager;
    typedef MpfaO2DFluxCalculator<TypeTag> FluxCalculator;
    typedef MpfaOInteractionVolumeContainer<TypeTag> InteractionVolumeContainer;

};

// specialization of class for mpfa-l method and dimension 2
template<class TypeTag>
class MpfaHelper<TypeTag, 2,
    typename std::enable_if<GET_PROP_VALUE(TypeTag, MpfaMethod) == GET_PROP(TypeTag, MpfaMethods)::lMethod>::type >
{
public:
    typedef MpfaL2DInteractionVolume<TypeTag> InteractionVolume;
    typedef MpfaL2DManager<TypeTag> InteractionVolumeManager;
    typedef MpfaL2DFluxCalculator<TypeTag> FluxCalculator;
    typedef MpfaLInteractionVolumeContainer<TypeTag> InteractionVolumeContainer;
};

}
#endif
