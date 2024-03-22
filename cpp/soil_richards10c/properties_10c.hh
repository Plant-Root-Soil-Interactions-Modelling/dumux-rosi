// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS10C_PROPERTIES_HH
#define RICHARDS10C_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/spgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards10c/model.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

//#include <RootSystem.h>

namespace Dumux {
namespace Properties {

namespace TTag { // Create new type tags
struct RichardsNCTT { using InheritsFrom = std::tuple<Richards10C>; };
struct Richards10CBox { using InheritsFrom = std::tuple<RichardsNCTT, BoxModel>; };
struct Richards10CCC { using InheritsFrom = std::tuple<RichardsNCTT, CCTpfaModel>; };
struct RichardsNCBox { using InheritsFrom = std::tuple<RichardsNCTT, BoxModel>; };
struct RichardsNCCC { using InheritsFrom = std::tuple<RichardsNCTT, CCTpfaModel>; };
}

// Set grid type
#ifndef GRIDTYPE
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsNCTT> { using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>; }; // using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>;
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsNCTT> { using type = GRIDTYPE; };  // Use GRIDTYPE from CMakeLists.txt
#endif

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsNCTT> { using type = Richards1P10CProblem<TypeTag>; }; 

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsNCTT> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsParams<GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>>;
};

//// Set the physical problem to be solved
//template<class TypeTag>
//struct PointSource<TypeTag, TTag::Richards1CTT> { using type = SolDependentPointSource<TypeTag>; };

/*
 * Define whether mole (true) or mass (false) fractions are used
 * TODO I only understand false...
 */
template<class TypeTag>
struct UseMoles<TypeTag, TTag::RichardsNCTT> { static constexpr bool value = true; };

} // end namespace properties
} // end namespace DUMUX

#endif
