// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_COUPLED_1P2C_PROPERTIES_HH
#define DUMUX_COUPLED_1P2C_PROPERTIES_HH

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

#include "../roots_1pnc/properties_stomata.hh"
#include "../soil_richards/properties.hh"

namespace Dumux {
namespace Properties {

/*
 * Define coupling manager according to dumux-rootgrowth
 */

// Coupling Properties for the Soil
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RichardsCC>
{
    using Traits = MultiDomainTraits<TypeTag, Properties::TTag::RootsOnePTwoCCCTpfa>;
    using type = EmbeddedCouplingManager1d3d<Traits, EmbeddedCouplingMode::line>;
};
// the point source type
template<class TypeTag>
struct PointSource<TypeTag, TTag::RichardsCC> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };
// the point source locater helper class
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::RichardsCC> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>;  };



// Coupling Properties for Roots
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RootsOnePTwoCCCTpfa>
{
    using Traits = MultiDomainTraits<Properties::TTag::RichardsCC, TypeTag>;
    using type = EmbeddedCouplingManager1d3d<Traits, EmbeddedCouplingMode::line>;
};
// the point source type
template<class TypeTag>
struct PointSource<TypeTag, TTag::RootsOnePTwoCCCTpfa> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };
// the point source locater helper class
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::RootsOnePTwoCCCTpfa> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };



} // namespace Properties
} // namespace Dumux

#endif
