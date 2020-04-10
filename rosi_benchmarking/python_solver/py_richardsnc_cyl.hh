#ifndef PYTHON_RICHARDS_CYL_H_
#define PYTHON_RICHARDS_CYL_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
#include <dune/pybindxi/functional.h>
namespace py = pybind11;

#include <config.h> // configuration file

#include "richards_cyl.hh" // includes solverbase

#include "../soil_richardsnc/richards1p2cproblem.hh" // the problem class

#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richardsncCylindrical1d/model.hh>

#include <dune/grid/spgrid.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

/**
 * create type tags
 */
namespace Dumux { namespace Properties {

namespace TTag { // Create new type tags

struct RichardsTT { using InheritsFrom = std::tuple<RichardsNC>; }; // defaults,  <dumux/porousmediumflow/richardsCylindrical1d/model.hh>
struct RichardsNCCylFoamTT { using InheritsFrom = std::tuple<RichardsTT>; }; // Foam grid
struct RichardsNCCylFoamCC { using InheritsFrom = std::tuple<RichardsNCCylFoamTT, CCTpfaModel>; };

};

template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsTT> { using type = Richards1P2CProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsNCCylFoamTT> { using type = Dune::FoamGrid<1,1>; }; //  Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 1>

template<class TypeTag>
struct UseMoles<TypeTag, TTag::RichardsTT> { static constexpr bool value = false; };

} }

#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
using RCFoamTT = Dumux::Properties::TTag::RichardsNCCylFoamCC;
using RichardsCylFoamAssembler = Dumux::FVAssembler<RCFoamTT, Dumux::DiffMethod::numeric>;
using RichardsCylFoamLinearSolver = Dumux::AMGBackend<RCFoamTT>;
using RichardsCylFoamProblem = Dumux::Richards1P2CProblem<RCFoamTT>;


PYBIND11_MODULE(rosi_richardsnc_cyl, m) {
    init_solverbase<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1 /*dimension*/>(m, "BaseRichardsNCCylFoam");
    init_richards_cyl<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1 /*dimension*/>(m, "RichardsNCCylFoam");
}

#endif

