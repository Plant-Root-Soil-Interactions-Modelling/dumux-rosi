#ifndef PYTHON_RICHARDS_CYLExtrusion_H_
#define PYTHON_RICHARDS_CYLExtrusion_H_

#include "../../../external/pybind11/include/pybind11/pybind11.h"
#include "../../../external/pybind11/include/pybind11/stl.h"
namespace py = pybind11;

#include <config.h> // configuration file

#include <dune/foamgrid/foamgrid.hh>
#include "../soil_richards/richardsproblem.hh" // the problem class
#include "richards_cyl.hh" // includes solverbase

#include <dumux/common/properties.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richardsCylindrical1dExtrusion/model.hh>

//

/**
 * create type tags
 */
namespace Dumux { namespace Properties {

namespace TTag { // Create new type tags

struct RichardsTT { using InheritsFrom = std::tuple<Richards>; }; // defaults,  <dumux/porousmediumflow/richardsCylindrical1d/model.hh>
struct RichardsCylFoamTT { using InheritsFrom = std::tuple<RichardsTT>; }; // Foam grid
struct RichardsCylFoamCC { using InheritsFrom = std::tuple<RichardsCylFoamTT, CCTpfaModel>; };

};

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsCylFoamTT> { using type = Dune::FoamGrid<1,1>; }; //  Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 1>

// // [[/codeblock]]

// // As mentioned before, we modify the areas and volumes used
// // for integrals over the control volume and the control volume faces by changing the extrusion type.
// // Here, we pass these traits to the grid geometry of the box scheme (the scheme
// // that we use here) and specialize the `GridGeometry` property accordingly.
// // [[codeblock]]
// template<class TypeTag>
// struct GridGeometry<TypeTag, TTag::RichardsCylFoamCC>
// {
// private:
    // static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    // using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    // using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    // // We take the default traits as basis and exchange the extrusion type
    // // The first axis (x-axis) is the radial axis, hence the zero. That means we rotate about the second axis (y-axis).
    // struct GGTraits : public CCTpfaDefaultGridGeometryTraits<GridView>
    // { using Extrusion = RotationalExtrusion<0>; };

// public:
    // // Pass the above traits to the box grid geometry such that it uses the
    // // rotation-symmetric sub-control volumes and faces.
    // using type = CCTpfaFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
// };

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
using RCFoamTT = Dumux::Properties::TTag::RichardsCylFoamCC;
using GridGeometryRCFoamTT = Dumux::GetPropType<RCFoamTT, Dumux::Properties::GridGeometry>;
using RichardsCylFoamAssembler = Dumux::FVAssembler<RCFoamTT, Dumux::DiffMethod::numeric>;
using RichardsCylFoamLinearSolver = Dumux::AMGBiCGSTABIstlSolver<Dumux::LinearSolverTraits<GridGeometryRCFoamTT>,//RCFoamTT,
	Dumux::LinearAlgebraTraitsFromAssembler<RichardsCylFoamAssembler>>;
using RichardsCylFoamProblem = Dumux::RichardsProblem<RCFoamTT>;


PYBIND11_MODULE(rosi_richards_cylExtrusion, m) {
    init_solverbase<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1 /*dimension*/>(m, "BaseRichardsCylFoamExtrusion2");
    init_richards_cyl<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1 /*dimension*/>(m, "RichardsCylFoamExtrusion");
}

#endif

