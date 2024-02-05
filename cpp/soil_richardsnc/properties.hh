// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDSNC_PROPERTIES_HH
#define RICHARDSNC_PROPERTIES_HH

//#include <dune/grid/yaspgrid.hh>
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

#include "../soil_richards/spatialparams.hh"

namespace Dumux {
namespace Properties {

namespace TTag { // Create new type tags
struct Richards2CTT { using InheritsFrom = std::tuple<RichardsNC>; };
struct Richards2CBox { using InheritsFrom = std::tuple<Richards2CTT, BoxModel>; };
struct Richards2CCC { using InheritsFrom = std::tuple<Richards2CTT, CCTpfaModel>; };
}

// Set grid type
#ifndef GRIDTYPE
template<class TypeTag>
struct Grid<TypeTag, TTag::Richards2CTT> { using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>; }; // using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>;
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::Richards2CTT> { using type = GRIDTYPE; };  // Use GRIDTYPE from CMakeLists.txt
#endif

// // Use 2d YaspGrid
// template<class TypeTag>
// struct Grid<TypeTag, TTag::Richards2CTT> { using type = Dune::YaspGrid<3>; };

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::Richards2CTT> { using type = Richards1P2CProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Richards2CTT> {
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
struct UseMoles<TypeTag, TTag::Richards2CTT> { static constexpr bool value = false; };

} // end namespace properties
} // end namespace DUMUX

#endif
