// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_ROOT_PROPERTIES_PERIOIDC_STOMATA_HH
#define DUMUX_ROOT_PROPERTIES_PERIODIC STOMATA_HH

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/h2oABA.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/periodic/tpfa/fvgridgeometry.hh>

#include "properties.hh"

namespace Dumux {
namespace Properties {

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RootsOnePTwoC> {
    using type = RootsStomataProblem<TypeTag>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RootsOnePTwoC> {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2OABA = FluidSystems::H2OABA<Scalar, FluidSystems::H2OABADefaultPolicy</*simplified=*/true>>;
    using type = FluidSystems::OnePAdapter<H2OABA, H2OABA::liquidPhaseIdx>;
};

// for CC
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::RootsOnePTwoCCCTpfa> {
private:
    using GridView = typename FVGridGeometry::GridView;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>; // ReorderingDofMapper
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = PeriodicCCTpfaFVGridGeometry<GridView, /*enableCache=*/true>;
};

/**
 * Pick either RootSpatialParamsDGF (for static dgf files),
 * or RootSpatialParamsRB (for dynamic root growth) as SpatialParams.type,
 */
#if DGF
template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RootsOnePTwoC> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsCaviationDGF<FVGridGeometry, Scalar>;
};
int simtype = dgf;
#endif
#if ROOTBOX
template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RootsOnePTwoC> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsRB<FVGridGeometry, Scalar>; // TODO
};
int simtype = rootbox;
#endif


} // namespace Properties
} // namespace Dumux

#endif


