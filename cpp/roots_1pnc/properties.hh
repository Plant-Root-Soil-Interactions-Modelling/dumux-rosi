// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_ROOTS_1PNC_PROPERTIES_HH
#define DUMUX_ROOTS_1PNC_PROPERTIES_HH

/**
 * General properties, shared by all roots_1pnc models
 */
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include <dune/foamgrid/foamgrid.hh>

namespace Dumux {
namespace Properties {

namespace TTag { // Create new type tags
struct RootsOnePTwoC { using InheritsFrom = std::tuple<OnePNC>; };
struct RootsOnePTwoCCCTpfa { using InheritsFrom = std::tuple<RootsOnePTwoC, CCTpfaModel>; };
struct RootsOnePTwoCBox { using InheritsFrom = std::tuple<RootsOnePTwoC, BoxModel>; };
}

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RootsOnePTwoC> {
    using type = Dune::FoamGrid<1, 3>;
};

// for Box (fixed for box, since box cannot be periodic),
// the type BoxFVGridGeometry is defined in dumux/discretization/box/fvgridgeometry.hh and is included in dumux/discretization/box.hh
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::RootsOnePTwoCBox> {
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

/**
 * Compile definitions are either DGF or ROOTBOX defined in CMakeLists
 */
enum modelType { dgf = 0, rootbox = 1 };

/*
 * Define whether mole (true) or mass (false) fractions are used
 * TODO I only understand false...
 */
template<class TypeTag>
struct UseMoles<TypeTag, TTag::RootsOnePTwoC> { static constexpr bool value = false; };

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
