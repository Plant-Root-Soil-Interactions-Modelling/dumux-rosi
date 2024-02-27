#ifndef PYTHON_RICHARDS_H_
#define PYTHON_RICHARDS_H_

#include "../../../CPlantBox/src/external/pybind11/include/pybind11/pybind11.h"
#include "../../../CPlantBox/src/external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

#include <config.h> // configuration file

#include "../soil_richards/richardsparams.hh"
#include "../soil_richards/richardsproblem.hh" // the problem class
#include "richards.hh" // includes solverbase

#include <dumux/common/properties.hh>
// #include <dumux/linear/istlsolverfactorybackend.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/pdesolver.hh>
#include <dune/istl/umfpack.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards/model.hh> // the model

#include <dune/common/parallel/mpicommunication.hh>

#include <dune/grid/spgrid.hh>
//#if HAVE_DUNE_ALUGRID
//#include <dune/alugrid/grid.hh>
//#endif
//#if HAVE_UG
//#include <dune/grid/uggrid.hh>
//#endif

#define PYBIND11_DETAILED_ERROR_MESSAGES

/**
 * create type tags
 */
namespace Dumux { namespace Properties {

namespace TTag { // Create new type tags

struct RichardsTT { using InheritsFrom = std::tuple<Richards>; }; // defaults, dumux/porousmediumflow/richards/model.hh
struct RichardsSPTT { using InheritsFrom = std::tuple<RichardsTT>; }; // sp grid
struct RichardsUGTT { using InheritsFrom = std::tuple<RichardsTT>; }; // ug grid
struct RichardsSPCC { using InheritsFrom = std::tuple<RichardsSPTT, CCTpfaModel>; };
struct RichardsSPBox { using InheritsFrom = std::tuple<RichardsSPTT, BoxModel>; };
struct RichardsUGCC { using InheritsFrom = std::tuple<RichardsUGTT, CCTpfaModel>; };
struct RichardsUGBox { using InheritsFrom = std::tuple<RichardsUGTT, BoxModel>; };

};
template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsSPTT> { using type = Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 3>; };

//template<class TypeTag> // Set grid type
//struct Grid<TypeTag, TTag::RichardsUGTT> { using type = Dune::UGGrid<3>; }; //  type = Dune::UGGrid<3>; Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>;

template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsTT> { using type = RichardsProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

// TODO: remove after release (3.6)
// Set the primary variables type
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::RichardsTT>
{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

} }

#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
 using RSPTT = Dumux::Properties::TTag::RichardsSPCC; // choose CC or Box
 using GridGeometryRSPTT = Dumux::GetPropType<RSPTT, Dumux::Properties::GridGeometry>;
 using RichardsSPAssembler = Dumux::FVAssembler<RSPTT, Dumux::DiffMethod::analytic>;
 using RichardsSPAssemblerNum = Dumux::FVAssembler<RSPTT, Dumux::DiffMethod::numeric>;
 
 using RichardsSPLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssembler>>; // IstlSolverFactoryBackend (dynamic from input files)
 using RichardsSPSeqLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::SeqLinearSolverTraits,
		 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssembler>>;
 using RichardsSPLinearSolverNum = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssemblerNum>>; // IstlSolverFactoryBackend (dynamic from input files)
 
 using RichardsSPProblem = Dumux::RichardsProblem<RSPTT>;

 
 using RichardsILUIstlLinearSolver = Dumux::ILUBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssembler>>; 
 // using RichardsSPAMGBiCGSTABIstlLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 // Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssembler>>;
 // using RichardsSPILUBiCGSTABIstlLinearSolver = Dumux::ILUBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 // Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssembler>>;
 using RichardsSPSSORCGIstlLinearSolver = Dumux::SSORCGIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssembler>>;
 using RichardsSPILURestartedGMResIstlLinearSolver = Dumux::ILURestartedGMResIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssembler>>;
		 
 // using RichardsSPUMFPackIstlLinearSolverNumeric = Dumux::UMFPackIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 // Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssemblerNumeric>>; // IstlSolverFactoryBackend (dynamic from input files)
 // using RichardsSPAMGBiCGSTABIstlLinearSolverNumeric = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRSPTT>,
		 // Dumux::LinearAlgebraTraitsFromAssembler<RichardsSPAssemblerNumeric>>;
		 
//using RUGTT = Dumux::Properties::TTag::RichardsUGCC; // choose CC or Box
//using GridGeometryRUGTT = Dumux::GetPropType<RUGTT, Dumux::Properties::GridGeometry>; // typename ????
//using RichardsUGAssembler = Dumux::FVAssembler<RUGTT, Dumux::DiffMethod::numeric>;
//using RichardsUGLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRUGTT>,
//		Dumux::LinearAlgebraTraitsFromAssembler<RichardsUGAssembler>>;
//using RichardsUGProblem = Dumux::RichardsProblem<RUGTT>;

PYBIND11_MODULE(rosi_richards, m) {
    init_solverbase<RichardsSPProblem, RichardsSPAssembler, RichardsSPLinearSolver>(m, "BaseRichardsSP");
    init_richards<RichardsSPProblem, RichardsSPAssembler, RichardsSPLinearSolver>(m, "RichardsSP");
	
	init_solverbase<RichardsSPProblem, RichardsSPAssemblerNum, RichardsSPLinearSolverNum>(m, "BaseRichardsSPnum");
    init_richards<RichardsSPProblem, RichardsSPAssemblerNum, RichardsSPLinearSolverNum>(m, "RichardsSPnum");
	
	init_solverbase<RichardsSPProblem, RichardsSPAssembler, RichardsSPSeqLinearSolver>(m, "BaseRichardsSPSeq");
    init_richards<RichardsSPProblem, RichardsSPAssembler, RichardsSPSeqLinearSolver>(m, "RichardsSPSeq");
	
	init_solverbase<RichardsSPProblem, RichardsSPAssembler, RichardsILUIstlLinearSolver>(m, "BaseRichardsSPILU");
    init_richards<RichardsSPProblem, RichardsSPAssembler, RichardsILUIstlLinearSolver>(m, "RichardsSPILU");
	
	init_solverbase<RichardsSPProblem, RichardsSPAssembler, RichardsSPSSORCGIstlLinearSolver>(m, "BaseRichardsSPSSORC");
    init_richards<RichardsSPProblem, RichardsSPAssembler, RichardsSPSSORCGIstlLinearSolver>(m, "RichardsSPSSORC");
	
	init_solverbase<RichardsSPProblem, RichardsSPAssembler, RichardsSPILURestartedGMResIstlLinearSolver>(m, "BaseRichardsSPILURes");
    init_richards<RichardsSPProblem, RichardsSPAssembler, RichardsSPILURestartedGMResIstlLinearSolver>(m, "RichardsSPILURes");
    
//    init_solverbase<RichardsUGProblem, RichardsUGAssembler, RichardsUGLinearSolver>(m, "BaseRichardsUG");
//    init_richards<RichardsUGProblem, RichardsUGAssembler, RichardsUGLinearSolver>(m, "RichardsUG");
}

#endif
