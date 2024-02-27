#ifndef PYTHON_RICHARDSNC_CYL_H_
#define PYTHON_RICHARDSNC_CYL_H_

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

#include <dumux/porousmediumflow/richardsncCylindrical1d/model.hh>


//#include <dumux/multidomain/traits.hh>
//#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

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
struct SpatialParams<TypeTag, TTag::RichardsTT> { using type = RichardsParams<GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsNCCylFoamTT> { using type = Dune::FoamGrid<1,1>; }; //  Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 1>


} }

#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types (for RichardsTT)

/*
 * pick assembler, linear solver and problem
 */
using RNCCFoamTT = Dumux::Properties::TTag::RichardsNCCylFoamCC;
using GridGeometryRNCCFoamTT = Dumux::GetPropType<RNCCFoamTT, Dumux::Properties::GridGeometry>;						
using RichardsNCCylFoamAssembler = Dumux::FVAssembler<RNCCFoamTT, Dumux::DiffMethod::numeric>;
using RichardsNCCylFoamLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRNCCFoamTT>,//RCFoamTT, 
	Dumux::LinearAlgebraTraitsFromAssembler<RichardsNCCylFoamAssembler>>;
using RichardsNCCylFoamProblem = Dumux::Richards1P2CProblem<RNCCFoamTT>;


PYBIND11_MODULE(rosi_richardsnc_cyl, m) {
    init_solverbase<RichardsNCCylFoamProblem, RichardsNCCylFoamAssembler, RichardsNCCylFoamLinearSolver, 1 /*dimension*/>(m, "BaseRichardsNCCylFoam");
    init_richards_cyl<RichardsNCCylFoamProblem, RichardsNCCylFoamAssembler, RichardsNCCylFoamLinearSolver, 1 /*dimension*/>(m, "RichardsNCCylFoam");
}

#endif

