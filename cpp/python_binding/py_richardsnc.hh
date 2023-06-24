#ifndef PYTHON_RICHARDSNC_H_
#define PYTHON_RICHARDSNC_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
#include <dune/pybindxi/functional.h>
namespace py = pybind11;

#include <config.h> // configuration file

#include "richards.hh" // includes solverbase

#include "../soil_richardsnc/richards1p2cproblem.hh" // the problem class

#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richardsnc/model.hh> // the model

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

struct RichardsNCTT { using InheritsFrom = std::tuple<RichardsNC>; }; // defaults, dumux/porousmediumflow/richards/model.hh
struct RichardsNCSPCC { using InheritsFrom = std::tuple<RichardsNCTT, CCTpfaModel>; };
struct RichardsNCSPBox { using InheritsFrom = std::tuple<RichardsNCTT, BoxModel>; };
};


template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsNCTT> { using type = Richards1P2CProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsNCTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsNCTT> { using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::RichardsNCTT> { static constexpr bool value = false; };


} }

#include "../soil_richardsnc/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
using RNCSPTT = Dumux::Properties::TTag::RichardsNCSPCC; // CC!
using RichardsNCSPAssembler = Dumux::FVAssembler<RNCSPTT, Dumux::DiffMethod::numeric>;
using RichardsNCSPLinearSolver = Dumux::AMGBackend<RNCSPTT>;
using RichardsNCSPProblem = Dumux::Richards1P2CProblem<RNCSPTT>;

PYBIND11_MODULE(rosi_richardsnc, m) {
    init_solverbase<RichardsNCSPProblem, RichardsNCSPAssembler, RichardsNCSPLinearSolver>(m, "BaseRichardsNCSP");
    init_richards<RichardsNCSPProblem, RichardsNCSPAssembler, RichardsNCSPLinearSolver>(m, "RichardsNCSP");
}

#endif
