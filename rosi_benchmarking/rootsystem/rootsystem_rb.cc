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
 * \brief test for the one-phase CC model
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <RootSystem.h>
#include <dumux/growth/rootsystemgridfactory.hh>
#include <dumux/growth/growthinterface.hh>
#include <dumux/growth/crootboxadapter.hh>
#include <dumux/growth/gridgrowth.hh>

#include "rootsproblem.hh"

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

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>



namespace Dumux {
namespace Properties {

template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::Roots> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsRB<FVGridGeometry, Scalar>;
};

} // end namespace Properties
}



int main(int argc, char** argv) try
{
    using namespace Dumux;
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // initialize MPI, finalize is done automatically on exit
    if (mpiHelper.rank() == 0) { // print dumux start message
        DumuxMessage::print(/*firstCall=*/true);
    } else {
        throw Dumux::ParameterException("Care! Foamgrid does not support parallel computation");
    }

    // define the type tag for this problem
    using TypeTag = Properties::TTag::RootsBox; // RootsCCTpfa RootsBox

    Parameters::init(argc, argv); // parse command line arguments and input file

    // initialize the CRootBox root system
    auto rootSystem = std::make_shared<CRootBox::RootSystem>();
    rootSystem->openFile(getParam<std::string>("RootSystem.Grid.File"), "modelparameter/");
    rootSystem->initialize();
    rootSystem->simulate(getParam<double>("RootSystem.Grid.InitialT"));

    // create a grid from a CRootBox root system
    using Grid = std::shared_ptr<Dune::FoamGrid<1, 3>>;
    Grid grid = GrowthModule::RootSystemGridFactory::makeGrid(*rootSystem); // in dumux/growth/rootsystemgridfactory.hh
    //    auto soilLookup = SoilLookUpBBoxTree<GrowthModule::Grid> (soilGridView, soilGridGeoemtry->boundingBoxTree(), saturation);
    //    rootSystem->setSoil(&soilLookup); todo

    //    using GridView = GetPropType<TypeTag, Properties::GridView>;
    //    using Element = typename GridView::template Codim<0>::Entity;
    //    using GlobalPosition = typename Element::Geometry::GlobalCoordinate; // the beauty

    using GlobalPosition = Dune::FieldVector<double, 3>;
    auto dumuxRootSystem = GrowthModule::CRootBoxAdapter<GlobalPosition>(*rootSystem);

    ////////////////////////////////////////////////////////////
    // run stationary or dynamic problem on this grid
    ////////////////////////////////////////////////////////////

    const auto& leafGridView = grid->leafGridView(); // we compute on the leaf grid view
    std::cout << "i have the view \n" << "\n" << std::flush;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>; // I assume BoxFVElementGeometry from fvelementgeometry.hh in dumux/discretization/box
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView); // but where is the property defined?
    fvGridGeometry->update();
    std::cout << "i have the geometry \n" << "\n" << std::flush;

    auto problem = std::make_shared<RootsProblem<TypeTag>>(fvGridGeometry);
    problem->spatialParams().updateParameters(dumuxRootSystem);
    // problem->spatialParams().analyseRootSystem(); // use for debugging
    // problem->spatialParams().initParameters(*gridData); // in case of grid
    std::cout << "... and, i have a problem \n" << "\n" << std::flush;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>; // in dumux/discretization/fvproperties.hh
    SolutionVector x(fvGridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;
    std::cout << "no solution, yet \n" << "\n" << std::flush;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>; // FVGridVariables, from dumux/discretization/fvproperties.hh
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);
    std::cout << "... but variables \n" << "\n" << std::flush;
    // grid variable knows the problem, the geometry and the solutionv vector

    bool grow = false;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    std::shared_ptr<CheckPointTimeLoop<Scalar>> timeLoop;
    if (tEnd > 0) { // dynamic problem
        grow = getParam<bool>("RootSystem.Grid.Grow", false); // use grid growth
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
            std::cout << "rootsystem.cc: no check times (TimeLoop.CheckTimes) defined in the input file\n";
        }
    } else { // static
    }
    std::cout << "time might be an issue \n" << "\n" << std::flush;

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
    std::cout << "vtk writer module initialized (how convenient)" << "\n" << std::flush;

    // class controlling the root growth
    // using Growth = typename GrowthModule::GridGrowth<TypeTag>;
    GrowthModule::GridGrowth<TypeTag> gridGrowth = GrowthModule::GridGrowth<TypeTag>(grid, fvGridGeometry, &dumuxRootSystem, x);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>; // dumux/assembly/fvassembler.hh
    std::shared_ptr<Assembler> assembler;
    if (tEnd > 0) {
        assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop); // dynamic
    } else {
        assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables); // static
    } // ? gridVariables already knows fvGridGeometry, and problem ?

    // the linear solver
    using LinearSolver = AMGBackend<TypeTag>; // dumux/linear/amgbackend.hh
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, fvGridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>; // dumux/nonlinear/newtonsolver.hh
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    std::cout << "planning to actually start \n" << "\n" << std::flush;
    if (tEnd > 0) // dynamic
    {
        std::cout << "a time dependent model" << "\n" << std::flush;

        timeLoop->start();
        do {

            if (grow) {
                std::cout << "grow() \n"<< std::flush;
                double dt = timeLoop->timeStepSize();
                gridGrowth.grow(dt);
                problem->spatialParams().updateParameters(dumuxRootSystem);
                std::cout << "grew \n"<< std::flush;

                gridVariables->updateAfterGridAdaption(x); // update the secondary variables
                assembler->setResidualSize(); // resize residual vector
                assembler->setJacobianPattern(); // resize and set Jacobian pattern
                assembler->setPreviousSolution(x);
                assembler->assembleJacobianAndResidual(x);
                std::cout << "hopefully modified assembler" << std::endl<< std::flush;;

                xOld = x;
            }

            assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations

            nonLinearSolver.solve(x); // solve the non-linear system with time step control // , *timeLoop

            gridVariables->advanceTimeStep();

            timeLoop->advanceTimeStep(); // advance to the time loop to the next step

            if ((timeLoop->isCheckPoint()) || (timeLoop->finished())) { // write vtk output (only at check points)
                problem->axialFlux(x); // prepare fields
                problem->radialFlux(x); // prepare fields
                vtkWriter.write(timeLoop->time());
            }
            problem->writeTranspirationRate(x); // always add transpiration data in the text file

            timeLoop->reportTimeStep(); // report statistics of this time step

            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize())); // set new dt as suggested by the newton solver

            problem->setTime(timeLoop->time()); // pass current time to the problem

        } while (!timeLoop->finished());

        timeLoop->finalize(leafGridView.comm());

    } else // static
    {
        std::cout << "a static model" << "\n" << std::flush;

        assembler->setPreviousSolution(xOld);  // set previous solution for storage evaluations

        nonLinearSolver.solve(x); // solve the non-linear system

        // write outputs
        problem->writeTranspirationRate(x);
        problem->axialFlux(x); // prepare fields
        problem->radialFlux(x); // prepare fields
        vtkWriter.write(1); // write vtk output
    }

    Parameters::print();
    DumuxMessage::print(/*firstCall=*/false);
    return 0;
} // end main
catch (Dumux::ParameterException &e) {
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e) {
    std::cerr << "DGF exception thrown (" << e <<
        "). Most likely, the DGF file name is wrong "
        "or the DGF file is corrupted, "
        "e.g. missing hash at end of file or wrong number (dimensions) of entries." << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (std::exception &e) {
    std::cerr << "Unknown exception thrown: " << e.what() << " ---> Abort!" << std::endl;
    return 4;
}

