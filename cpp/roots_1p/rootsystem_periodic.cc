// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:

/*!
 * \file
 *
 * \brief Doussan model for xylem flux (using dumux/porousmediumflow/1p/model.hh)
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <memory>

// Dune
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/istl/io.hh> // debug vector/matrix output
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>

// Dumux
#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/geometry/boundingboxtree.hh>
#include <dumux/common/geometry/geometricentityset.hh>
#include <dumux/common/geometry/intersectingentities.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh> // functions to resume a simulation

#include <dumux/periodic/tpfa/periodicnetworkgridmanager.hh>
#include <dumux/periodic/tpfa/fvgridgeometry.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "rootsproblem.hh"

#include "properties_periodic.hh" // the property system related stuff (to pass types, used instead of polymorphism)
#include "properties_nocoupling.hh" // dummy types for replacing the coupling types

/**
 * and so it begins...
 */
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::RootsCCTpfa; // RootsCCTpfa, RootsBox (TypeTag is defined in the problem class richardsproblem.hh)

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)
    if (mpiHelper.rank() == 0) { // print dumux start message
        DumuxMessage::print(/*firstCall=*/true);
    } else {
        throw Dumux::ParameterException("Care! Foamgrid does not support parallel computation");
    }

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    using GlobalPosition = Dune::FieldVector<double, 3>;
    std::bitset<3> periodic("110");
    if (hasParam("Grid.Periodic"))  {
        periodic = getParam<std::bitset<3>>("Grid.Periodic");
    }
    GlobalPosition lower = { -1.e9, -1.e9, -1.e9 };
    GlobalPosition upper = { 1.e9, 1.e9, 1.e9 };

    // Create the gridmanager and grid
    using Grid = Dune::FoamGrid<1, 3>;
    std::shared_ptr<Grid> grid;
    PeriodicNetworkGridManager<3> gridManager(lower, upper, periodic); // only for dgf
    std::cout << "\nSimulation type is dgf \n\n" << std::flush;
    gridManager.init("RootSystem");
    grid = std::shared_ptr<Grid>(&gridManager.grid(), Properties::empty_delete<Grid>());

    // we compute on the leaf grid view
    const auto& leafGridView = grid->leafGridView();
    std::cout << "i have the view \n"<< std::flush;

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(gridManager.grid().leafGridView());
    const auto periodicConnectivity = gridManager.getGridData()->createPeriodicConnectivity(fvGridGeometry->elementMapper(), fvGridGeometry->vertexMapper());
    fvGridGeometry->setExtraConnectivity(periodicConnectivity);
    fvGridGeometry->update();
    std::cout << "i have the geometry \n" << std::flush;

    ////////////////////////////////////////////////////////////
    // run stationary or dynamic problem on this grid
    ////////////////////////////////////////////////////////////

    // the problem (initial and boundary conditions)
    auto problem = std::make_shared<RootsProblem<TypeTag>>(fvGridGeometry);
    problem->spatialParams().initParameters(*gridManager.getGridData());

    // check if we are about to restart a previously interrupted simulation
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>; // defined in discretization/fvproperties.hh, as Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>
    SolutionVector x(fvGridGeometry->numDofs()); // degrees of freedoms
    if (restartTime > 0)
    {
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        const auto fileName = getParam<std::string>("Restart.RootFile");
        const auto pvName = createPVNameFunction<IOFields, PrimaryVariables, ModelTraits, FluidSystem>();
        loadSolution(x, fileName, pvName, *fvGridGeometry);
    }
    else
        problem->applyInitialSolution(x); // Dumux way of saying x = problem->applyInitialSolution()
    auto xOld = x;

    std::cout << "i have a problem \n" << std::flush;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);
    std::cout << "with variables \n" << std::flush;

    // get some time loop parameters & instantiate time loop
    bool grow = false;
    const auto tEnd = getParam<double>("TimeLoop.TEnd");
    std::shared_ptr<CheckPointTimeLoop<double>> timeLoop;
    if (tEnd > 0) { // dynamic problem
        grow = getParam<bool>("RootSystem.Grid.Grow", false); // use grid growth
        auto initialDt = getParam<double>("TimeLoop.DtInitial"); // initial time step
        timeLoop = std::make_shared<CheckPointTimeLoop<double>>(restartTime, initialDt, tEnd);
        timeLoop->setMaxTimeStepSize(getParam<double>("TimeLoop.MaxTimeStepSize"));
        if (hasParam("TimeLoop.CheckTimes")) {
            std::vector<double> checkPoints = getParam<std::vector<double>>("TimeLoop.CheckTimes");
            for (auto p : checkPoints) {
                timeLoop->setCheckPoint(p);
            }
        }
        if (hasParam("TimeLoop.PeriodicCheckTimes")) {
            std::cout << "using periodic check times \n";
            timeLoop->setPeriodicCheckPoint(getParam<double>("TimeLoop.PeriodicCheckTimes"));
        }
    } else { // static
    }
    std::cout << "time might be an issue \n" << std::flush;

    // intialize the vtk output module
    std::cout << "vtk writer module... \n" << std::flush;
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    problem->userData("pSoil", x);
    problem->userData("radius", x);
    problem->userData("order", x);
    problem->userData("id", x);
    problem->userData("axialFlux", x);// todo wrong (coarse approximation)
    problem->userData("radialFlux", x);
    problem->userData("age", x);
    problem->userData("initialPressure",x);
    problem->userData("kr", x);
    problem->userData("kx", x);
    vtkWriter.addField(problem->p(), "p [cm]"); // cm pressure head
    vtkWriter.addField(problem->radius(), "radius [m]"); // not in cm, because of tube plot
    vtkWriter.addField(problem->order(), "order [1]");
    vtkWriter.addField(problem->id(), "id [1]");
    vtkWriter.addField(problem->axialFlux(), "axial flux [cm3/d]");
    vtkWriter.addField(problem->radialFlux(), "radial flux [cm3/d]");
    vtkWriter.addField(problem->age(), "age [d]");
    vtkWriter.addField(problem->initialPressure(), "initial pressure [cm]");
    vtkWriter.addField(problem->kr(), "kr [cm/hPa/d]");
    vtkWriter.addField(problem->kx(), "kx [cm4/hPa/day]");
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(restartTime);
    std::cout << "vtk writer module initialized \n" << std::flush;

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    std::shared_ptr<Assembler> assembler;
    if (tEnd > 0) {
        assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop); // dynamic
    } else {
        assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables); // static
    }

    // the linear solver
    using LinearSolver = AMGBackend<TypeTag>; // how do i choose umfpack
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, fvGridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    std::cout << "\ni plan to actually start \n" << std::flush;

    if (tEnd > 0) // dynamic
    {
        std::cout << "a time dependent model\n\n" << std::flush;
        timeLoop->start();
        do {

            double t = timeLoop->time(); // dumux time
            double dt = timeLoop->timeStepSize(); // dumux time step
            problem->setTime(t, dt); // pass current time to the problem
            problem->postTimeStep(x, *gridVariables);
            problem->writeTranspirationRate(); // add transpiration data into the text file

            assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations
            nonLinearSolver.solve(x, *timeLoop); // solve the non-linear system with time step control
            xOld = x; // make the new solution the old solution

            gridVariables->advanceTimeStep();
            timeLoop->advanceTimeStep(); // advance to the time loop to the next step

            if ((timeLoop->isCheckPoint()) || (timeLoop->finished())) { // write vtk output (only at check points)
                if (grow) { // prepare static fields also
                    problem->userData("radius", x);
                    problem->userData("order", x);
                    problem->userData("id", x);
                    problem->userData("initialPressure", x);
                }
                problem->userData("p", x);
                problem->userData("axialFlux", x);
                problem->userData("radialFlux", x);
                problem->userData("age", x); // age changes with time
                problem->userData("kr", x);  // conductivities change with age
                problem->userData("kx", x);
                vtkWriter.write(timeLoop->time());
            }

            timeLoop->reportTimeStep();  // report statistics of this time step
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize())); // set new dt as suggested by the newton solver

        } while (!timeLoop->finished());

        timeLoop->finalize(leafGridView.comm());

    } else { // static

        std::cout << "a static model \n\n" << std::flush;

        assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations
        nonLinearSolver.solve(x); // solve the non-linear system

        // write outputs
        problem->userData("p", x);
        problem->userData("axialFlux", x);
        problem->userData("radialFlux", x);
        problem->userData("age", x); // prepare fields
        problem->userData("kr", x);  // conductivities change with age
        problem->userData("kx", x);
        vtkWriter.write(1); // write vtk output
        problem->postTimeStep(x, *gridVariables);
        problem->writeTranspirationRate();
    }

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0) {
        Parameters::print();
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e) {
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
} catch (Dune::DGFException & e) {
    std::cerr << "DGF exception thrown (" << e <<
        "). Most likely, the DGF file name is wrong "
        "or the DGF file is corrupted, "
        "e.g. missing hash at end of file or wrong number (dimensions) of entries." << " ---> Abort!" << std::endl;
    return 2;
} catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
} catch (std::exception &e) {
    std::cerr << "Unknown exception thrown: " << e.what() << " ---> Abort!" << std::endl;
    return 4;
}
