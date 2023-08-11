// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS5C_PROPERTIES_CYL_K_HH
#define RICHARDS5C_PROPERTIES_CYL_K_HH

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

#include <dumux/porousmediumflow/richards5cCylindrical1d/model.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

//#include <RootSystem.h>

namespace Dumux {
namespace Properties {

namespace TTag { // Create new type tags
struct RichardsNCTT { using InheritsFrom = std::tuple<Richards5C>; };
struct Richards2CBox { using InheritsFrom = std::tuple<RichardsNCTT, BoxModel>; };
struct Richards2CCC { using InheritsFrom = std::tuple<RichardsNCTT, CCTpfaModel>; };
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
struct Problem<TypeTag, TTag::RichardsNCTT> { using type = Richards1P5CProblem<TypeTag>; }; 

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsNCTT> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>;
};

/*
 * Define whether mole (true) or mass (false) fractions are used
 */
template<class TypeTag>
struct UseMoles<TypeTag, TTag::RichardsNCTT> { static constexpr bool value = true; };

// Set the physical problem to be solved
template<class TypeTag>
struct PointSource<TypeTag, TTag::RichardsNCTT> { 
	
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
	using type = SolDependentPointSource<TypeTag>; 
};


} // end namespace properties
} // end namespace DUMUX

#endif
