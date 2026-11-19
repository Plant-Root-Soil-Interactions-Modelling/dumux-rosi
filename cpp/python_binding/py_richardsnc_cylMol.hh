#ifndef PYTHON_RICHARDSNC_CYLMOL_H_
#define PYTHON_RICHARDSNC_CYLMOL_H_

#include "../../../CPlantBox/src/external/pybind11/include/pybind11/pybind11.h"
#include "../../../CPlantBox/src/external/pybind11/include/pybind11/stl.h"
namespace py = pybind11;

#include <config.h> // configuration file

#include <dune/foamgrid/foamgrid.hh>		  
#include "richards_cyl.hh" // includes solverbase

#include "../soil_richardsnc/richards1p2cproblem.hh" // the problem class

#include <dumux/common/properties.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richardsncCylindrical1d/modelMoles.hh>


//#include <dumux/multidomain/traits.hh>
//#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

/**
 * create type tags
 */
namespace Dumux { namespace Properties {

namespace TTag { // Create new type tags

struct RichardsTT { using InheritsFrom = std::tuple<RichardsNC>; }; // defaults,  <dumux/porousmediumflow/richardsCylindrical1d/model.hh>
struct RichardsNCCylMolFoamTT { using InheritsFrom = std::tuple<RichardsTT>; }; // Foam grid
struct RichardsNCCylMolFoamCC { using InheritsFrom = std::tuple<RichardsNCCylMolFoamTT, CCTpfaModel>; };

};

template<class TypeTag> // Set Problem
struct Problem<TypeTag, TTag::RichardsTT> { using type = Richards1P2CProblem<TypeTag>; };

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::RichardsTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsNCCylMolFoamTT> { using type = Dune::FoamGrid<1,1>; }; //  Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 1>

//template<class TypeTag>
//struct UseMoles<TypeTag, TTag::RichardsNCTT> { static constexpr bool value = true; };

} }

#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
using RNCCFoamMolTT = Dumux::Properties::TTag::RichardsNCCylMolFoamCC;
using GridGeometryRNCCFoamMolTT = Dumux::GetPropType<RNCCFoamMolTT, Dumux::Properties::GridGeometry>;						
using RichardsNCCylMolFoamAssembler = Dumux::FVAssembler<RNCCFoamMolTT, Dumux::DiffMethod::numeric>;				
using RichardsNCCylMolFoamAssemblerAna = Dumux::FVAssembler<RNCCFoamMolTT, Dumux::DiffMethod::analytic>;
using RichardsNCCylFoamLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRNCCFoamMolTT>,//RCFoamTT, 
	Dumux::LinearAlgebraTraitsFromAssembler<RichardsNCCylMolFoamAssembler>>;
using RichardsNCCylFoamLinearSolverAna = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRNCCFoamMolTT>,//RCFoamTT, 
	Dumux::LinearAlgebraTraitsFromAssembler<RichardsNCCylMolFoamAssemblerAna>>;
using RichardsNCCylFoamProblem = Dumux::Richards1P2CProblem<RNCCFoamMolTT>;


PYBIND11_MODULE(rosi_richardsnc_cylMol, m) {
    init_solverbase<RichardsNCCylFoamProblem, RichardsNCCylMolFoamAssembler, RichardsNCCylFoamLinearSolver, 1 /*dimension*/>(m, "BaseRichardsNCCylMolFoam");
    init_richards_cyl<RichardsNCCylFoamProblem, RichardsNCCylMolFoamAssembler, RichardsNCCylFoamLinearSolver, 1 /*dimension*/>(m, "RichardsNCCylMolFoam");
    init_solverbase<RichardsNCCylFoamProblem, RichardsNCCylMolFoamAssemblerAna, RichardsNCCylFoamLinearSolverAna, 1 /*dimension*/>(m, "BaseRichardsNCCylFoamMolAna");
    init_richards_cyl<RichardsNCCylFoamProblem, RichardsNCCylMolFoamAssemblerAna, RichardsNCCylFoamLinearSolverAna, 1 /*dimension*/>(m, "RichardsNCCylFoamMolAna");
}

#endif

