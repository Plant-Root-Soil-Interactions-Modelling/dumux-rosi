#ifndef PYTHON_RICHARDS4C_H_
#define PYTHON_RICHARDS4C_H_

#include "../../../CPlantBox/src/external/pybind11/include/pybind11/pybind11.h"
#include "../../../CPlantBox/src/external/pybind11/include/pybind11/stl.h"
namespace py = pybind11;

#include <config.h> // configuration file

#include "../soil_richards/richardsparams.hh"							   
#include "richards4.hh" // includes in solverbase

#include "../soil_richards4c/richards1p4cproblem.hh" // the problem class

#include <dumux/common/properties.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards4c/model.hh> // the model

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

struct RichardsNCTT { using InheritsFrom = std::tuple<Richards4C>; }; // defaults, dumux/porousmediumflow/richards/model.hh
struct RichardsSPTT { using InheritsFrom = std::tuple<RichardsNCTT>; }; // sp grid
struct RichardsUGTT { using InheritsFrom = std::tuple<RichardsNCTT>; }; // ug grid
struct RichardsNCSPCC { using InheritsFrom = std::tuple<RichardsSPTT, CCTpfaModel>; };
struct RichardsNCSPBox { using InheritsFrom = std::tuple<RichardsSPTT, BoxModel>; };
struct RichardsUGCC { using InheritsFrom = std::tuple<RichardsUGTT, CCTpfaModel>; };
struct RichardsUGBox { using InheritsFrom = std::tuple<RichardsUGTT, BoxModel>; };
};


template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsNCTT> { using type = Richards1P4CProblem<TypeTag>; };

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
using RichardsSPAssembler = Dumux::FVAssembler<RNCSPTT, Dumux::DiffMethod::analytic>;
using RichardsSPAssemblerNum = Dumux::FVAssembler<RNCSPTT, Dumux::DiffMethod::numeric>;

using RichardsSPLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssembler>>; 
using RichardsSPSeqLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::SeqLinearSolverTraits,
	 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssemblerNum>>;
using RichardsSPLinearSolverNum = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
	 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssemblerNum>>;


using RichardsILUIstlLinearSolver = Dumux::ILUBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
	 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssemblerNum>>; 
using RichardsSPSSORCGIstlLinearSolver = Dumux::SSORCGIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
	 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssemblerNum>>;
using RichardsSPILURestartedGMResIstlLinearSolver = Dumux::ILURestartedGMResIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
	 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssemblerNum>>;

using RichardsNCSPProblem = Dumux::Richards1P4CProblem<RNCSPTT>;


PYBIND11_MODULE(rosi_richards4c, m) {
    init_solverbase<RichardsNCSPProblem, RichardsSPAssembler, RichardsSPLinearSolver>(m, "BaseRichards4CSP");
    init_richards_4<RichardsNCSPProblem, RichardsSPAssembler, RichardsSPLinearSolver>(m, "Richards4CSP");

	init_solverbase<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPLinearSolverNum>(m, "BaseRichards4CSPnum");
    init_richards_4<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPLinearSolverNum>(m, "Richards4CSPnum");
	
	init_solverbase<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPSeqLinearSolver>(m, "BaseRichards4CSPSeq");
    init_richards_4<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPSeqLinearSolver>(m, "Richards4CSPSeq");
	
	init_solverbase<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsILUIstlLinearSolver>(m, "BaseRichards4CSPILU");
    init_richards_4<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsILUIstlLinearSolver>(m, "Richards4CSPILU");
	
	init_solverbase<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPSSORCGIstlLinearSolver>(m, "BaseRichards4CSPSSORC");
    init_richards_4<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPSSORCGIstlLinearSolver>(m, "Richards4CSPSSORC");
	
	init_solverbase<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPILURestartedGMResIstlLinearSolver>(m, "BaseRichards4CSPILURes");
    init_richards_4<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPILURestartedGMResIstlLinearSolver>(m, "Richards4CSPILURes");
}

#endif
