#ifndef PYTHON_RICHARDSNC_H_
#define PYTHON_RICHARDSNC_H_

#include "../../../CPlantBox/src/external/pybind11/include/pybind11/pybind11.h"
#include "../../../CPlantBox/src/external/pybind11/include/pybind11/stl.h"
namespace py = pybind11;

#include <config.h> // configuration file

#include "../soil_richards/richardsparams.hh"							   
#include "richards.hh" // includes solverbase

#include "../soil_richardsnc/richards1p2cproblem.hh" // the problem class

#include <dumux/common/properties.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richardsnc/model.hh> // the model

#include <dune/grid/spgrid.hh>	  
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#define PYBIND11_DETAILED_ERROR_MESSAGES		
/**
 * create type tags
 */
namespace Dumux { namespace Properties {

namespace TTag { // Create new type tags

struct RichardsNCTT { using InheritsFrom = std::tuple<RichardsNC>; }; // defaults, dumux/porousmediumflow/richards/model.hh
struct RichardsSPTT { using InheritsFrom = std::tuple<RichardsNCTT>; }; // sp grid
struct RichardsUGTT { using InheritsFrom = std::tuple<RichardsNCTT>; }; // ug grid
struct RichardsNCSPCC { using InheritsFrom = std::tuple<RichardsSPTT, CCTpfaModel>; };
struct RichardsNCSPBox { using InheritsFrom = std::tuple<RichardsSPTT, BoxModel>; };
struct RichardsUGCC { using InheritsFrom = std::tuple<RichardsUGTT, CCTpfaModel>; };
struct RichardsUGBox { using InheritsFrom = std::tuple<RichardsUGTT, BoxModel>; };
};


template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsNCTT> { using type = Richards1P2CProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsNCTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

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
using GridGeometryRSPTT = typename Dumux::GetPropType<RNCSPTT, Dumux::Properties::GridGeometry>;
using RichardsNCSPAssembler = Dumux::FVAssembler<RNCSPTT, Dumux::DiffMethod::analytic>;
using RichardsNCSPAssemblerNum = Dumux::FVAssembler<RNCSPTT, Dumux::DiffMethod::numeric>;
using RichardsNCSPLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,// 
	Dumux::LinearAlgebraTraitsFromAssembler<RichardsNCSPAssembler>>;
using RichardsNCSPLinearSolverNum = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
	 Dumux::LinearAlgebraTraitsFromAssembler<RichardsNCSPAssemblerNum>>;
using RichardsNCSPProblem = Dumux::Richards1P2CProblem<RNCSPTT>;

using RUGTT = Dumux::Properties::TTag::RichardsUGCC; // choose CC or Box
using GridGeometryRUGTT = Dumux::GetPropType<RUGTT, Dumux::Properties::GridGeometry>; // typename ????
using RichardsUGAssembler = Dumux::FVAssembler<RUGTT, Dumux::DiffMethod::numeric>;
using RichardsUGLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRUGTT>,
		Dumux::LinearAlgebraTraitsFromAssembler<RichardsUGAssembler>>;
using RichardsUGProblem = Dumux::Richards1P2CProblem<RUGTT>;


PYBIND11_MODULE(rosi_richardsnc, m) {
    init_solverbase<RichardsNCSPProblem, RichardsNCSPAssembler, RichardsNCSPLinearSolver>(m, "BaseRichardsNCSPanalytic");
    init_richards<RichardsNCSPProblem, RichardsNCSPAssembler, RichardsNCSPLinearSolver>(m, "RichardsNCSPanalytic");
	
	init_solverbase<RichardsNCSPProblem, RichardsNCSPAssemblerNum, RichardsNCSPLinearSolver>(m, "BaseRichardsNCSP");
    init_richards<RichardsNCSPProblem, RichardsNCSPAssemblerNum, RichardsNCSPLinearSolver>(m, "RichardsNCSP");
}

#endif
