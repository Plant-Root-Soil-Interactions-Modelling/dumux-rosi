#ifndef PYTHON_RICHARDS_CYL_H_
#define PYTHON_RICHARDS_CYL_H_

#include "external/pybind11/include/pybind11/pybind11.h"
namespace py = pybind11;

#include <config.h> // configuration file

#include "richards_cyl.hh" // includes solverbase

#include "../soil_richards/richardsproblem.hh" // the problem class

#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richardsCylindrical1d/model.hh>

#include <dune/grid/spgrid.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

/**
 * create type tags
 */
namespace Dumux { namespace Properties {

namespace TTag { // Create new type tags

struct RichardsTT { using InheritsFrom = std::tuple<Richards>; }; // defaults,  <dumux/porousmediumflow/richardsCylindrical1d/model.hh>
struct RichardsCylFoamTT { using InheritsFrom = std::tuple<RichardsTT>; }; // Foam grid
struct RichardsCylFoamCC { using InheritsFrom = std::tuple<RichardsCylFoamTT, CCTpfaModel>; };

};

template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsTT> { using type = RichardsProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsCylFoamTT> { using type = Dune::FoamGrid<1,1>; }; //  Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 1>

} }

#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
using RCFoamTT = Dumux::Properties::TTag::RichardsCylFoamCC;
using RichardsCylFoamAssembler = Dumux::FVAssembler<RCFoamTT, Dumux::DiffMethod::numeric>;
using RichardsCylFoamLinearSolver = Dumux::AMGBackend<RCFoamTT>;
using RichardsCylFoamProblem = Dumux::RichardsProblem<RCFoamTT>;


PYBIND11_MODULE(rosi_richards_cyl, m) {
    init_solverbase<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1 /*dimension*/>(m, "BaseRichardsCylFoam");
    init_richards_cyl<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1 /*dimension*/>(m, "RichardsCylFoam");
}

#endif

