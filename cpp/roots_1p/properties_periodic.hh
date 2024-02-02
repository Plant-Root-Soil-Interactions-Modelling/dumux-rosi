// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_ROOTS_PROPERTIES_PERIODIC_HH
#define DUMUX_ROOTS_PROPERTIES_PERIODIC_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/periodic/tpfa/fvgridgeometry.hh>

namespace Dumux {
namespace Properties {

namespace TTag { // Create new type tags
struct Roots { using InheritsFrom = std::tuple<OneP>; };
struct RootsCCTpfa { using InheritsFrom = std::tuple<Roots, CCTpfaModel>; };
struct RootsBox { using InheritsFrom = std::tuple<Roots, BoxModel>; };
}

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Roots> {
    using type = Dune::FoamGrid<1, 3>;
};

// for CC
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::RootsCCTpfa> {
private:
    using GridView = typename FVGridGeometry::GridView;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>; // ReorderingDofMapper
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
     using type = PeriodicCCTpfaFVGridGeometry<GridView, /*enableCache=*/true>;
};

// for Box
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::RootsBox> {
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using GridView = typename FVGridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>; //ReorderingDofMapper
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, BoxDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Roots> {
    using type = RootsProblem<TypeTag>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Roots> {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
};

/**
 * Compile definitions are either DGF or ROOTBOX defined in CMakeLists
 */
enum modelType { dgf = 0, rootbox = 1 };

/**
 * Pick either RootSpatialParamsDGF (for static dgf files),
 * or RootSpatialParamsRB (for dynamic root growth) as SpatialParams.type,
 */
#if DGF
template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::Roots> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsCaviationDGF<FVGridGeometry, Scalar>;
};
int simtype = dgf;
#endif
#if ROOTBOX
template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::Roots> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsRB<FVGridGeometry, Scalar>;
};
int simtype = rootbox;
#endif

/**
 * to wrap a raw pointer into a shared pointer:
 * for not deleting it twice, an empty deleter must be defined
 */
template <typename T>
struct empty_delete {
    empty_delete() /* noexcept */
    { }
    template <typename U>
    empty_delete(const empty_delete<U>&,
        typename std::enable_if<
            std::is_convertible<U*, T*>::value
        >::type* = nullptr) /* noexcept */
    { }
    void operator()(T* const) const /* noexcept */
    { }// do nothing
};

} // namespace Properties
} // namespace Dumux

#endif
