// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS_PROPERTIES_CYL_HH
#define RICHARDS_PROPERTIES_CYL_HH

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

#include <dumux/porousmediumflow/richardsCylindrical1d/model.hh>

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
struct Grid<TypeTag, TTag::RichardsTT> { using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>; }; // using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>;
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
    using type = RichardsParams<GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>>;
};

} // end namespace properties
} // end namespace DUMUX

#endif
