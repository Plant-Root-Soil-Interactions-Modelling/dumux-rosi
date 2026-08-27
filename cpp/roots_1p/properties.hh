// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_ROOTS_PROPERTIES_HH
#define DUMUX_ROOTS_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>
// #include <dumux/common/reorderingdofmapper.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/geometry/quadraturerules.hh>
//#include <dumux/common/reorderingdofmapper.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

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
struct GridGeometry<TypeTag, TTag::RootsCCTpfa>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
	//todo: Y not use the ReorderingDofMapper?
	//todo: add the extrusion factor here?
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;// Dune::MultipleCodimMultipleGeomTypeMapper<GridView>?
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, CCTpfaDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

// for Box
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::RootsBox> {
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
	//todo: Y not use the ReorderingDofMapper?
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>; //ReorderingDofMapper
    // using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
	
	
    // We take the default traits as basis and exchange the extrusion type
    // The first axis (x-axis) is the radial axis, hence the zero. That means we rotate about the second axis (y-axis).
    // struct GGTraits : public BoxDefaultGridGeometryTraits<GridView>
    // { using Extrusion = RotationalExtrusion<0>; };
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, BoxDefaultGridGeometryTraits<GridView>>;
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
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsCaviationDGF<GridGeometry, Scalar>;
};
int simtype = dgf;
#endif
#if ROOTBOX
template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::Roots> {
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsRB<GridGeometry, Scalar>;
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
