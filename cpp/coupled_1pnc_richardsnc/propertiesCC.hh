// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_COUPLED_1P2C__PROPERTIES_HH
#define DUMUX_COUPLED_1P2C__PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/geometry/quadraturerules.hh>
//#include <dumux/common/reorderingdofmapper.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "../roots_1pnc/properties.hh" // TypeTag:RootsOnePTwoCC
#include "../roots_1pnc/properties_1p2c.hh" // TypeTag:RootsOnePTwoCC

#include "../soil_richardsnc/properties.hh" // TypeTag:Richards2C
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

namespace Dumux {
namespace Properties {

using CouplingTransport = Embedded1d3dCouplingManager<MultiDomainTraits<
    Properties::TTag::Richards2CCC, Properties::TTag::RootsOnePTwoCCCTpfa>,
    Embedded1d3dCouplingMode::Line
>;

// tell the tissue sub-model about the coupling
template<class TypeTag> struct CouplingManager<TypeTag, TTag::Richards2CCC> { using type = CouplingTransport; };
template<class TypeTag> struct PointSource<TypeTag, TTag::Richards2CCC> { using type = CouplingTransport::PointSourceTraits::template PointSource<0>; };
template<class TypeTag> struct PointSourceHelper<TypeTag, TTag::Richards2CCC> { using type = CouplingTransport::PointSourceTraits::template PointSourceHelper<0>; };

// tell the network sub-model about the coupling
template<class TypeTag> struct CouplingManager<TypeTag, TTag::RootsOnePTwoCCCTpfa> { using type = CouplingTransport; };
template<class TypeTag> struct PointSource<TypeTag, TTag::RootsOnePTwoCCCTpfa> { using type = CouplingTransport::PointSourceTraits::template PointSource<1>; };
template<class TypeTag> struct PointSourceHelper<TypeTag, TTag::RootsOnePTwoCCCTpfa> { using type = CouplingTransport::PointSourceTraits::template PointSourceHelper<1>; };


} // namespace Properties
} // namespace Dumux

#endif
