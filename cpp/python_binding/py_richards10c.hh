#ifndef PYTHON_RICHARDS10C_H_
#define PYTHON_RICHARDS10C_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
#include <dune/pybindxi/functional.h>
namespace py = pybind11;

#include <config.h> // configuration file

#include "richards_10.hh" // includes solverbase

#include "../soil_richardsnc/richards1p10cproblem_cyl.hh" // the problem class

#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards10c/model.hh> // the model

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

struct RichardsNCTT { using InheritsFrom = std::tuple<Richards10C>; }; // defaults, dumux/porousmediumflow/richards/model.hh
struct Richards10CSPCC { using InheritsFrom = std::tuple<RichardsNCTT, CCTpfaModel>; };
struct Richards10CSPBox { using InheritsFrom = std::tuple<RichardsNCTT, BoxModel>; };
};


template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsNCTT> { using type = Richards1P10CProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsNCTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::FVGridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsNCTT> { using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::RichardsNCTT> { static constexpr bool value = true; };


} }

#include "../soil_richardsnc/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
using R10CSPTT = Dumux::Properties::TTag::Richards10CSPCC; // CC!
using Richards10CSPAssembler = Dumux::FVAssembler<R10CSPTT, Dumux::DiffMethod::numeric>;
using Richards10CSPLinearSolver = Dumux::AMGBackend<R10CSPTT>;
using Richards10CSPProblem = Dumux::Richards1P10CProblem<R10CSPTT>;

PYBIND11_MODULE(rosi_richards10c, m) {
    init_solverbase<Richards10CSPProblem, Richards10CSPAssembler, Richards10CSPLinearSolver>(m, "BaseRichards10CSP");
    init_richards_10<Richards10CSPProblem, Richards10CSPAssembler, Richards10CSPLinearSolver>(m, "Richards10CSP");
}

#endif
