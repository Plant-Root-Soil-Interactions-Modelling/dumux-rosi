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

#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dune/common/timer.hh> // to compute wall times
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

// #include <dumux/common/properties.hh> // creates an undefined TypeTag types, and includes the property system
// #include <dumux/common/properties/propertysystem.hh>
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file
#include <dumux/common/valgrind.hh> // for debugging
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/defaultusagemessage.hh> // for information (a global function)

#include <dumux/linear/amgbackend.hh> // linear solver (currently the only solver available)
#include <dumux/nonlinear/newtonsolver.hh> // the only nonlinear solver available

#include <dumux/common/timeloop.hh>
#include <dumux/assembly/fvassembler.hh> // assembles residual and Jacobian of the nonlinear system

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
// #include <dumux/io/loadsolution.hh> // global functions to resume a simulation

#include <RootSystem.h>

#include <dumux/growth/rootsystemgridfactory.hh> // dumux-rosi growth ideas (modified from dumux-rootgrowth)
#include <dumux/growth/growthinterface.hh>
#include <dumux/growth/cplantboxadapter.hh>
#include <dumux/growth/gridgrowth.hh>

#include "rootsproblem.hh"

#include "properties.hh"
#include "properties_nocoupling.hh" // dummy types for replacing the coupling types

/**
 * and so it begins...
 */
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::RootsCCTpfa; // RootsCCTpfa, RootsBox (TypeTag is defined in the problem class richardsproblem.hh)
    int simtype = Properties::simtype;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)
    if (mpiHelper.rank() == 0) { // print dumux start message
        DumuxMessage::print(/*firstCall=*/true);
    } else {
        throw Dumux::ParameterException("Care! Foamgrid does not support parallel computation");
    }

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    //    using GridView = GetPropType<TypeTag, Properties::GridView>;
    //    using Element = typename GridView::template Codim<0>::Entity;
    //    using GlobalPosition = typename Element::Geometry::GlobalCoordinate; // the beauty (is there a correct & quicker way?)
    using GlobalPosition = Dune::FieldVector<double, 3>;

    // Create the gridmanager and grid
    using Grid = Dune::FoamGrid<1, 3>;
    std::shared_ptr<Grid> grid;
    GridManager<Grid> gridManager; // only for dgf
    std::shared_ptr<CPlantBox::RootSystem> rootSystem; // only for rootbox
    GrowthModule::GrowthInterface<GlobalPosition>* growth = nullptr; // in case of RootBox (or in future PlantBox)
    if (simtype==Properties::dgf) { // for a static dgf grid
        std::cout << "\nSimulation type is dgf \n\n" << std::flush;
        gridManager.init("RootSystem");
        grid = std::shared_ptr<Grid>(&gridManager.grid(), Properties::empty_delete<Grid>());
    } else if (simtype==Properties::rootbox) { // for a root model (static or dynamic)
        std::cout << "\nSimulation type is RootBox \n\n" << std::flush;
        rootSystem = std::make_shared<CPlantBox::RootSystem>();
        rootSystem->openFile(getParam<std::string>("RootSystem.Grid.File"), "../modelparameter/");
        if (hasParam("RootSystem.Grid.Confined")) {
            auto box = getParam<std::vector<double>>("RootSystem.Grid.Confined");
            rootSystem->setGeometry(std::make_shared<CPlantBox::SDF_PlantBox>(box.at(0)*100, box.at(1)*100, box.at(2)*100));
        } else { // half plane
            rootSystem->setGeometry(std::make_shared<CPlantBox::SDF_HalfPlane>(CPlantBox::Vector3d(0.,0.,0.5), CPlantBox::Vector3d(0.,0.,1.))); // care, collar needs to be top, make sure plant seed is located below -1 cm
        }
        rootSystem->initialize();
        double shootZ = getParam<double>("RootSystem.Grid.ShootZ", 0.); // root system initial time
        grid = GrowthModule::RootSystemGridFactory::makeGrid(*rootSystem, shootZ, true); // in dumux/growth/rootsystemgridfactory.hh
        //  todo static soil for hydrotropsim ...
        //    auto soilLookup = SoilLookUpBBoxTree<GrowthModule::Grid> (soilGridView, soilGridGeoemtry->boundingBoxTree(), saturation);
        //    rootSystem->setSoil(&soilLookup);
        growth = new GrowthModule::CPlantBoxAdapter<GlobalPosition>(rootSystem);
    }

    // we compute on the leaf grid view
    const auto& leafGridView = grid->leafGridView();
    std::cout << "i have the view \n"<< std::flush;

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();
    std::cout << "i have the geometry \n" << std::flush;

    ////////////////////////////////////////////////////////////
    // run stationary or dynamic problem on this grid
    ////////////////////////////////////////////////////////////

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>; // defined in discretization/fvproperties.hh, as Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>
    SolutionVector x(fvGridGeometry->numDofs()); // degrees of freedoms

    // root growth
    GrowthModule::GridGrowth<TypeTag>* gridGrowth = nullptr;
    double initialTime = 0.; // s
    if (simtype==Properties::rootbox) {
        gridGrowth = new GrowthModule::GridGrowth<TypeTag>(grid, fvGridGeometry, growth, x); // in growth/gridgrowth.hh
        std::cout << "...grid grower initialized \n" << std::flush;
        initialTime = getParam<double>("RootSystem.Grid.InitialT")*24*3600;
        gridGrowth->grow(initialTime);
        std::cout << "\ninitial growth performed... \n" << std::flush;
    }

    // the problem (initial and boundary conditions)
    auto problem = std::make_shared<RootsProblem<TypeTag>>(fvGridGeometry);
    if (simtype==Properties::dgf) {
        problem->spatialParams().initParameters(*gridManager.getGridData());
    } else if (simtype==Properties::rootbox){
        problem->spatialParams().updateParameters(*growth);
    }
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
        timeLoop = std::make_shared<CheckPointTimeLoop<double>>(/*start time*/0., initialDt, tEnd);
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
    problem->userData("age", x);
    problem->userData("initialPressure",x);
    problem->userData("kr", x);
    problem->userData("kx", x);
    vtkWriter.addField(problem->p(), "p [cm]"); // cm pressure head
    vtkWriter.addField(problem->radius(), "radius [m]"); // not in cm, because of tube plot
    vtkWriter.addField(problem->order(), "order [1]");
    vtkWriter.addField(problem->id(), "id [1]");
    vtkWriter.addField(problem->age(), "age [d]");
    vtkWriter.addField(problem->initialPressure(), "initial pressure [cm]");
    vtkWriter.addField(problem->kr(), "kr [cm/hPa/d]");
    vtkWriter.addField(problem->kx(), "kx [cm4/hPa/day]");
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);
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

            if (simtype==Properties::rootbox) {
                if (grow) {

                    // std::cout << "time " << growth->simTime()/24/3600 << " < " << (t+initialTime)/24/3600 << "\n";
                    while (growth->simTime()+dt<t+initialTime) {

                        std::cout << "\n grow ..."<< std::flush;
                        gridGrowth->grow(dt);
                        problem->spatialParams().updateParameters(*growth);
                        problem->applyInitialSolution(x); // reset todo (? does this make sense)?
                        std::cout << "grew \n"<< std::flush;

                        // what shall I update?
                        fvGridGeometry->update();
                        //gridVariables->update();
                        gridVariables->updateAfterGridAdaption(x); // update the secondary variables

                        // todo? what is necessary? no clue what i am doing ...
                        assembler->setResidualSize(); // resize residual vector
                        assembler->setJacobianPattern(); // resize and set Jacobian pattern
                        assembler->setPreviousSolution(x);
                        assembler->assembleJacobianAndResidual(x);

                        xOld = x;
                    }

                }
            }

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
                problem->userData("pSoil", x);
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
        problem->userData("pSoil", x);
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
