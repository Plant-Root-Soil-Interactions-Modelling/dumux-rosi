// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Richards equation realized with Richards box model.
 */
#include <config.h>

#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dune/common/timer.hh> // to compute wall times
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid/coordinates.hh>
#include <dune/istl/io.hh>
// #include <dumux/common/properties.hh> // creates an undefined TypeTag types, and includes the property system
// #include <dumux/common/properties/propertysystem.hh>
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file
#include <dumux/common/valgrind.hh> // for debugging
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/defaultusagemessage.hh> // for information (a function)

#include <dumux/linear/amgbackend.hh> // linear solver (currently the only parallel solver available(?))
#include <dumux/porousmediumflow/richards/newtonsolver.hh>
/**
 * Some small adaption to <dumux/nonlinear/newtonsolver.hh>, which is the only nonlinear solver available.
 * The adaption is disabled per default (parameter EnableChop = false)
 */
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/assembly/fvassembler.hh> // assembles residual and Jacobian of the nonlinear system

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
// #include <dumux/io/loadsolution.hh> // functions to resume a simulation

#include "richards1p10cproblem_cyl.hh" // the problem class. Defines some TypeTag types and includes its spatialparams.hh class
#include "properties_cyl_10c.hh" // the property system related stuff (to pass types, used instead of polymorphism)
#include "properties_nocoupling.hh" // dummy types for replacing the coupling types

#include <dune/grid/utility/globalindexset.hh>
/**
 * here we go
 */
int main(int argc, char** argv) //try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::Richards10CCC; 
    using Problem = GetPropType<TypeTag, Properties::Problem>;
     using Grid = typename Problem::Grid;
    //using Grid = GetPropType<TypeTag, Problem::Grid>;
    using SoilGridType = GetPropType<TypeTag, Properties::Grid>;
    // using FVGridGeometry = typename Problem::FVGridGeometry;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>; // defined in discretization/fvproperties.hh, as Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>
    // using SolutionVector = typename Problem::SolutionVector;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    // using FluxVariables = typename Problem::FluxVariables;

    using GridData = Dumux::GridData<SoilGridType>;
    // using GridView = typename Grid::Traits::LeafGridView;
    using GridView = Dumux::GridData< SoilGridType::Traits::LeafGridView>;
    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)

    // print dumux start message
    if (mpiHelper.rank() == 0) { // rank is the process number
        DumuxMessage::print(/*firstCall=*/true);
    }

    // parse command line arguments and input file
    Parameters::init(argc, argv);
	
	
    /*
     * Parses the parameters from the .input file into a static data member of the static class Parameters,
     * the parameter tree can then be accessed with Parameters::paramTree() .
     * Furthermore, global auxiliary functions are defined e.g. getParam, setParam, and haveParam
     * (in parameters.hh).
     *
     * All Dumux classes access the global parameter tree where considered appropriate.
     */

    // try to create a grid (from the given grid file or the input file)
    /**
     * Properties::Grid is defined in the problem class richardsproblem.hh to Dune::YaspGrid<3>,
     * or if available to the compile definition GRIDTYPE that is given in CMakeLists, and is
     * YaspGrid<3> (for richards3d), Dune::FoamGrid<1,1> (for richards1d),
     * or Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming> (for richardsUG)
     */
    GridManager<SoilGridType> gridManager;
    gridManager.init("Soil"); // "Soil" is the parameter group name
    /**
     * Opens the grid file, or constructs the grid form the .input file parameters
     */


    /////////////////////////////////////////////////////////////////
    // run steady state or dynamic non-linear problem on this grid
    /////////////////////////////////////////////////////////////////


    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    /**
     * we work only on the leafGridView and do not need the grid ever again,
     * leafGridView.grid() returns a const reference to the grid (in case we need data attached to the grid).
     *
     * Of dune type: GridFamily::Traits::LeafGridView
     *
     * Have not found where the type is set (grid dependent),
     * but probably implements DefaultLeafGridView (in dune/grid/common/defaultgridview.hh)
     */

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    /**
     * The type is dependent on the discretization (e.g. box, tpfa, ...),
     *
     * For Box method:
     * Properties::FVGridGeometry is defined in problem.hh -> discretization/box.hh
     * The type is BoxFVGridGeometry (in discretization/box/fvgridgeometry.hh)
     * specialization of BaseFVGridGeometry (discretization/basefvgridgeometry.hh)
     */
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    /**
     * holds the vertexMapper() and elementMapper().
     *
     * ??? copies the leafGridView (but, how does it know when its updated?)
     * it seems i am only allowed to use fvGridGeometry->gridView() in the following
     * rendering leafGridView defined above, pointless
     */
    fvGridGeometry->update(); // update all Mappers(do this again after grid adaption)



    const auto rank = fvGridGeometry->gridView().comm().rank();
	std::shared_ptr<Dune::GlobalIndexSet<GridView>> cellIdx; // global index mappers
    cellIdx = std::make_shared<Dune::GlobalIndexSet<GridView>>(gridManager.grid().leafGridView(), 0);
	  // Create global index map
	  //Dune::GlobalIndexSet<GridView> cellIdx(gridManager.grid().leafGridView(),0);


		std::cout<<"getCellCenters rank:"<<rank<<"\n";
        for (const auto& e : elements(fvGridGeometry->gridView())) {
            auto p = e.geometry().center();
			 //const auto dofIdxGlobal = fvGridGeometry().elementMapper().index(element);
			int gIdx = cellIdx->index(e);
			std::cout<<"	gIdx "<<gIdx;
            for (int i=0; i< 2 ; i++) { // found no better way
				std::cout<<" "<<p[i];
            }std::cout<<std::endl;
        }


    return 0;
}
