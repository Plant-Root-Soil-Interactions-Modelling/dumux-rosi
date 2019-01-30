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
 * \brief Coupling
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/periodic/tpfa/periodicnetworkgridmanager.hh>
#include <dumux/periodic/tpfa/fvgridgeometry.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh>

#include <RootSystem.h>
#include "../rootsystem/rootsproblem.hh"
#include "../soil/richardsproblem.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>


namespace Dumux {
namespace Properties {

template<class TypeTag> // Set the spatial parameters for the root problem
struct SpatialParams<TypeTag, TTag::Roots> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsDGF<FVGridGeometry, Scalar>;
};

} // end namespace Properties
}



int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::RootsBox;
    // RootsCCTpfa RootsBox

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init("RootSystem");
    const auto gridData = gridManager.getGridData();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    std::cout << "i have the view \n" << "\n" << std::flush;

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();
    std::cout << "i have the geometry \n" << "\n" << std::flush;

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(fvGridGeometry);
    problem->spatialParams().initParameters(*gridData);
    // problem->spatialParams().analyseRootSystem();
    std::cout << "and i have a problem \n" << "\n" << std::flush;

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(fvGridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;
    std::cout << "no solution yet \n" << "\n" << std::flush;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);
    std::cout << "... but variables \n" << "\n" << std::flush;

    // get some time loop parameters & instantiate time loop
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    std::shared_ptr<CheckPointTimeLoop<Scalar>> timeLoop;
    if (tEnd > 0) { // dynamic problem
        const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
        auto initialDt = getParam<Scalar>("TimeLoop.DtInitial"); // initial time step
        timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(/*start time*/0., initialDt, tEnd);
        timeLoop->setMaxTimeStepSize(maxDt);
        try { // CheckPoints defined
            std::vector<double> checkPoints = getParam<std::vector<double>>("TimeLoop.CheckTimes");
            // insert check points
            for (auto p : checkPoints) { // don't know how to use the setCheckPoint( initializer list )
                timeLoop->setCheckPoint(p);
            }
        } catch (std::exception& e) {
            std::cout << "rootsysstem.cc: no check times (TimeLoop.CheckTimes) defined in the input file\n";
        }
    } else { // static
    }
    std::cout << "time might be an issue \n" << "\n" << std::flush;

    // intialize the vtk output module

    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    problem->axialFlux(x); // prepare fields
    problem->radialFlux(x); // prepare fields
    vtkWriter.addField(problem->axialFlux(), "axial flux");
    vtkWriter.addField(problem->radialFlux(), "radial flux");
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);
    std::cout << "vtk writer module initialized" << "\n" << std::flush;

//    using GridView = GetPropType<TypeTag, Properties::GridView>;
//    auto numDofs = fvGridGeometry->vertexMapper().size();
//    auto numVert = fvGridGeometry->gridView().size(GridView::dimension);
//    std::cout << "gridview dimension  " << GridView::dimension << "\n";
//    std::cout << "numDofs " << numDofs << "\n";
//    std::cout << "numVert  " << numVert << "\n";
//    return 0; // help

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    std::shared_ptr<Assembler> assembler;
    if (tEnd > 0) {
        assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop); // dynamic
    } else {
        assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables); // static
    }

    // the linear solver
    using LinearSolver = AMGBackend<TypeTag>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, fvGridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    std::cout << "i plan to actually start \n" << "\n" << std::flush;

    if (tEnd > 0) // dynamic
        {
        std::cout << "a time dependent model" << "\n" << std::flush;
        timeLoop->start();
        do {
            // set previous solution for storage evaluations
            assembler->setPreviousSolution(xOld);
            // solve the non-linear system with time step control
            nonLinearSolver.solve(x, *timeLoop);
            // make the new solution the old solution
            xOld = x;
            gridVariables->advanceTimeStep();
            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();
            // write vtk output (only at check points)
            if ((timeLoop->isCheckPoint()) || (timeLoop->finished())) {
                problem->axialFlux(x); // prepare fields
                problem->radialFlux(x); // prepare fields
                vtkWriter.write(timeLoop->time());
            }
            if (mpiHelper.rank() == 0) {
                problem->writeTranspirationRate(x);
            }
            // report statistics of this time step
            timeLoop->reportTimeStep();
            // set new dt as suggested by the newton solver
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
            // pass current time to the problem
            problem->setTime(timeLoop->time());
        } while (!timeLoop->finished());
        timeLoop->finalize(leafGridView.comm());
    } else // static
    {
        std::cout << "a static model" << "\n" << std::flush;
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);
        // solve the non-linear system
        nonLinearSolver.solve(x);
        // write vtk output
        problem->axialFlux(x); // prepare fields
        problem->radialFlux(x); // prepare fields
        problem->writeTranspirationRate(x);
        vtkWriter.write(1);
    }

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
        "). Most likely, the DGF file name is wrong "
        "or the DGF file is corrupted, "
        "e.g. missing hash at end of file or wrong number (dimensions) of entries." << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (std::exception &e)
{
    std::cerr << "Unknown exception thrown: " << e.what() << " ---> Abort!" << std::endl;
    return 4;
}

