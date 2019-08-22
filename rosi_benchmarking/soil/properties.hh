#ifndef DUMUX_SOIL_PROPERTIES_HH
#define DUMUX_SOIL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards/model.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

#include <RootSystem.h>

namespace Dumux {
namespace Properties {

namespace TTag { // Create new type tags
struct RichardsTT { using InheritsFrom = std::tuple<Richards>; };
struct RichardsBox { using InheritsFrom = std::tuple<RichardsTT, BoxModel>; };
struct RichardsCC { using InheritsFrom = std::tuple<RichardsTT, CCTpfaModel>; };
}

// Set grid type
#ifndef GRIDTYPE
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsTT> { using type = Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsTT> { using type = GRIDTYPE; };  // Use GRIDTYPE from CMakeLists.txt
#endif

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsTT> { using type = RichardsProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsTT> {
    using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>;
};

/**
 * no functionality (but Dumux wants its bindings)
 * ugly, but I found no other option...
 */
// The point source type (not used)
template<class TypeTag>
struct PointSource<TypeTag, TTag::RichardsTT> {
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using type = IntegrationPointSource<GlobalPosition, NumEqVector>;
};
/// Dummy types
class DummyPointSourceData {
public:
    double lowDimElementIdx() { throw 1; };
};
class DummySpatial {
public:
    double kr(int i) const { throw 1; };
    double radius(int i) const { throw 1; };
};
class DummyProblem {
public:
    DummySpatial spatialParams() { throw 1; };
};
class DummyCouplingManager {
public:
    std::vector<double> bulkPriVars(int i) { throw 1; };
    std::vector<double>  lowDimPriVars(int i) { throw 1; };
    DummyProblem& problem(int i) { throw 1; };
    DummyPointSourceData& pointSourceData(int i) { throw 1; };
};
// For a dummy manager
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RichardsTT> {
    using type = DummyCouplingManager;
};

} // end namespace properties
} // end namespace DUMUX

#endif
