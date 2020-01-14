// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_SOIL_PROPERTIES_HH
#define DUMUX_SOIL_PROPERTIES_HH

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

#include <dumux/porousmediumflow/richardsnc/model.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

#include <RootSystem.h>

namespace Dumux {
namespace Properties {

namespace TTag { // Create new type tags
struct Richards1CTT { using InheritsFrom = std::tuple<RichardsNC>; };
struct Richards1CBox { using InheritsFrom = std::tuple<Richards1CTT, BoxModel>; };
struct Richards1CCC { using InheritsFrom = std::tuple<Richards1CTT, CCTpfaModel>; };
}

// Set grid type
#ifndef GRIDTYPE
template<class TypeTag>
struct Grid<TypeTag, TTag::Richards1CTT> { using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>; }; // using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>;
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::Richards1CTT> { using type = GRIDTYPE; };  // Use GRIDTYPE from CMakeLists.txt
#endif

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::Richards1CTT> { using type = Richards1CProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Richards1CTT> {
    using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>;
};

} // end namespace properties
} // end namespace DUMUX

#endif
