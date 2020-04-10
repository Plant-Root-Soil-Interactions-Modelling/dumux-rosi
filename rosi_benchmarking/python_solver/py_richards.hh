#ifndef PYTHON_RICHARDS_H_
#define PYTHON_RICHARDS_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
#include <dune/pybindxi/functional.h>
namespace py = pybind11;

#include <config.h> // configuration file

#include "richards.hh" // includes solverbase

#include "../soil_richards/richardsproblem.hh" // the problem class

#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards/model.hh>

#include <dune/grid/spgrid.hh>
//#if HAVE_DUNE_ALUGRID
//#include <dune/alugrid/grid.hh>
//#endif
//#if HAVE_UG
//#include <dune/grid/uggrid.hh>
//#endif

/**
 * create type tags
 */
namespace Dumux { namespace Properties {

namespace TTag { // Create new type tags

struct RichardsTT { using InheritsFrom = std::tuple<Richards>; }; // defaults, dumux/porousmediumflow/richards/model.hh
struct RichardsSPTT { using InheritsFrom = std::tuple<RichardsTT>; }; // sp grid
//struct RichardsUGTT { using InheritsFrom = std::tuple<RichardsTT>; }; // ug grid
struct RichardsSPCC { using InheritsFrom = std::tuple<RichardsSPTT, CCTpfaModel>; };
struct RichardsSPBox { using InheritsFrom = std::tuple<RichardsSPTT, BoxModel>; };
// struct RichardsUGCC { using InheritsFrom = std::tuple<RichardsUGTT, CCTpfaModel>; };
//struct RichardsUGBox { using InheritsFrom = std::tuple<RichardsUGTT, BoxModel>; };

};


template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsTT> { using type = RichardsProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsSPTT> { using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>; };

//template<class TypeTag> // Set grid type
//struct Grid<TypeTag, TTag::RichardsUGTT> { using type = Dune::UGGrid<3>; };

} }

#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
using RSPTT = Dumux::Properties::TTag::RichardsSPCC; // CC!
using RichardsSPAssembler = Dumux::FVAssembler<RSPTT, Dumux::DiffMethod::numeric>;
using RichardsSPLinearSolver = Dumux::AMGBackend<RSPTT>;
using RichardsSPProblem = Dumux::RichardsProblem<RSPTT>;

//using RUGTT = Dumux::Properties::TTag::RichardsUGBox;
//using RichardsUGAssembler = Dumux::FVAssembler<RUGTT, Dumux::DiffMethod::numeric>;
//using RichardsUGLinearSolver = Dumux::AMGBackend<RUGTT>;
//using RichardsUGProblem = Dumux::RichardsProblem<RUGTT>;

PYBIND11_MODULE(rosi_richards, m) {
    init_solverbase<RichardsSPProblem, RichardsSPAssembler, RichardsSPLinearSolver>(m, "BaseRichardsSP");
    init_richardssp<RichardsSPProblem, RichardsSPAssembler, RichardsSPLinearSolver>(m, "RichardsSP");
//    init_solverbase<RichardsUGProblem, RichardsUGAssembler, RichardsUGLinearSolver>(m, "BaseRichardsUG");
//    init_richardssp<RichardsUGProblem, RichardsUGAssembler, RichardsUGLinearSolver>(m, "RichardsUG");
}

#endif

