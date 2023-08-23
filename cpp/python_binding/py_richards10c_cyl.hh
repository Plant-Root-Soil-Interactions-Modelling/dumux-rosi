#ifndef PYTHON_RICHARDS10C_CYL_H_
#define PYTHON_RICHARDS10C_CYL_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
#include <dune/pybindxi/functional.h>
namespace py = pybind11;

#include <config.h> // configuration file

#include "richards_10cyl.hh" // includes solverbase

#include "../soil_richardsnc/richards1p10cproblem_cyl.hh" // the problem class

#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards10cCylindrical1d/model.hh>

#include <dune/grid/spgrid.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>



/**
 * create type tags
 */
namespace Dumux { namespace Properties {

namespace TTag { // Create new type tags

struct RichardsTT { using InheritsFrom = std::tuple<Richards10C>; }; // defaults,  <dumux/porousmediumflow/richardsCylindrical1d/model.hh>
struct Richards10CCylFoamTT { using InheritsFrom = std::tuple<RichardsTT>; }; // Foam grid
struct Richards10CCylFoamCC { using InheritsFrom = std::tuple<Richards10CCylFoamTT, CCTpfaModel>; };

};

template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsTT> { using type = Richards1P10CProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::Richards10CCylFoamTT> { using type = Dune::FoamGrid<1,1>; }; //  Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 1>

template<class TypeTag>
struct UseMoles<TypeTag, TTag::RichardsTT> { static constexpr bool value = true; };

} }

#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
using RCFoamTT = Dumux::Properties::TTag::Richards10CCylFoamCC;
using RichardsCylFoamAssembler = Dumux::FVAssembler<RCFoamTT, Dumux::DiffMethod::numeric>;
using RichardsCylFoamLinearSolver = Dumux::AMGBackend<RCFoamTT>;
using RichardsCylFoamProblem = Dumux::Richards1P10CProblem<RCFoamTT>;


PYBIND11_MODULE(rosi_richards10c_cyl, m) {
    init_solverbase<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1 /*dimension*/>(m, "BaseRichards10CCylFoam");
    init_richards_10cyl<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1 /*dimension*/>(m, "Richards10CCylFoam");
}


#endif

