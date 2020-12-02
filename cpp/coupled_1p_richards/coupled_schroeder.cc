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
 * Monolythic coupling
 */
#include <config.h>

#include <ctime>
#include <iostream>

// Dune
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>
#include <dune/grid/common/rangegenerators.hh>

// Dumux
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh> // functions to resume a simulation

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

// growth model
#include <RootSystem.h>

// dumux-rosi
#include <dumux/growth/rootsystemgridfactory.hh>
#include <dumux/growth/growthinterface.hh>
#include <dumux/growth/cplantboxadapter.hh>
#include <dumux/growth/gridgrowth.hh>

#include "../roots_1p/rootsproblem_schroeder.hh"
#include "../soil_richards/richardsproblem_schroeder.hh"

#include "propertiesCC.hh" // includes root properties, soil properties, redefines coupling manager
// for Box                  properties.hh // <- not working for UG
// for CCTpfa               propertiesCC.hh // <- working, but bad results for UG
// for box soil, CC roots,  propertiesMix.hh (CC roots needed for periodicity)
// cahnge L70 & L71 accordingly

namespace Dumux {

using SoilTypeTag = Properties::TTag::RichardsCC; // RichardsCC //RichardsBox
using RootTypeTag = Properties::TTag::RootsCCTpfa; // RootsBox // RootsCCTpfa

/**
 * debugging
 */
using SoilFVGridGeometry = GetPropType<SoilTypeTag, Properties::FVGridGeometry>;
template<class SoilGridVariables, class SoilSolution>
void soilControl(const SoilFVGridGeometry& gridGeometry, const SoilGridVariables& gridVariables,
    const SoilSolution& sol, const SoilSolution& oldSol, double t, double dt) {
    double cVol = 0.;
    double oldVol = 0.;
    const auto& gridView = gridGeometry.gridView();  // soil
    for (const auto& element : elements(gridView)) { // soil elements
        auto fvGeometry = localView(gridGeometry); // soil solution -> volume variable
        fvGeometry.bindElement(element);
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        elemVolVars.bindElement(element, fvGeometry, sol);
        for (const auto& scv : scvs(fvGeometry)) {
            cVol += elemVolVars[scv].saturation(0)*scv.volume();
        }
        elemVolVars.bindElement(element, fvGeometry, oldSol);
        for (const auto& scv : scvs(fvGeometry)) {
            oldVol += elemVolVars[scv].saturation(0)*scv.volume();
        }
    }
    std::cout << "\nWater in domain: " << cVol*1.e6 << " g at day " << t/24/3600 << " \n";
    std::cout << "...   a change of: " << (oldVol-cVol)*1.e3 << " kg = " << (oldVol-cVol)*1.e6*24*3600/dt << " g/day \n" ;
}
} // namespace Dumux




/**
 * and so it begins...
 */
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)
    if (mpiHelper.rank() == 0) { // print dumux start message
        DumuxMessage::print(/*firstCall=*/true);
    } else {
        throw Dumux::ParameterException("Care! Foamgrid does not support parallel computation");
    }

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    std::string rootName = getParam<std::string>("Problem.RootName");
    Parameters::init(0, argv, rootName);
    std::string soilName = getParam<std::string>("Problem.SoilName");
    Parameters::init(0, argv, soilName);
    Parameters::init(argc, argv);

    // Define the sub problem type tags (see properties.hh)
    int simtype = Properties::simtype;

    // soil grid
    GridManager<GetPropType<SoilTypeTag, Properties::Grid>> soilGridManager;
    soilGridManager.init("Soil"); // pass parameter group (see input file)

    // soil grid geometry
    const auto& soilGridView = soilGridManager.grid().leafGridView();
    using SoilFVGridGeometry = GetPropType<SoilTypeTag, Properties::FVGridGeometry>;
    auto soilGridGeometry = std::make_shared<SoilFVGridGeometry>(soilGridView);
    soilGridGeometry->update();

    // root gridmanager and grid
    using GlobalPosition = Dune::FieldVector<double, 3>;
    using Grid = Dune::FoamGrid<1, 3>;
    std::shared_ptr<Grid> rootGrid;
    GridManager<Grid> rootGridManager; // only for dgf
    std::shared_ptr<CPlantBox::RootSystem> rootSystem; // only for rootbox
    GrowthModule::GrowthInterface<GlobalPosition>* growth = nullptr; // in case of RootBox (or in future PlantBox)
    if (simtype==Properties::dgf) { // for a static dgf grid
        std::cout << "\nSimulation type is dgf \n\n" << std::flush;
        rootGridManager.init("RootSystem");
        rootGrid = std::shared_ptr<Grid>(&rootGridManager.grid(), Properties::empty_delete<Grid>());
    } else if (simtype==Properties::rootbox) { // for a root model (static or dynamic)
        std::cout << "\nSimulation type is RootBox \n\n" << std::flush;
        rootSystem = std::make_shared<CPlantBox::RootSystem>();
        rootSystem->readParameters("modelparameter/rootsystem/"+getParam<std::string>("RootSystem.Grid.File")+".xml");
        // make sure we don't grow above the soil, but allow to grow in x and y because we will do the periodic mapping TODO
        // rootSystem->setGeometry(new CPlantBox::SDF_HalfPlane(CPlantBox::Vector3d(0.,0.,0.5), CPlantBox::Vector3d(0.,0.,1.))); // care, collar needs to be top, make sure plant seed is located below -1 cm
        const auto size = soilGridGeometry->bBoxMax() - soilGridGeometry->bBoxMin();
        // rootSystem->setGeometry(new CPlantBox::SDF_PlantBox(size[0]*100, size[1]*100, size[2]*100));
        rootSystem->initialize();
        double shootZ = getParam<double>("RootSystem.Grid.ShootZ", 0.); // root system initial time
        rootGrid = GrowthModule::RootSystemGridFactory::makeGrid(*rootSystem, shootZ, true); // in dumux/growth/rootsystemgridfactory.hh
        //  todo static soil for hydrotropsim ...
        //    auto soilLookup = SoilLookUpBBoxTree<GrowthModule::Grid> (soilGridView, soilGridGeoemtry->boundingBoxTree(), saturation);
        //    rootSystem->setSoil(&soilLookup);
        growth = new GrowthModule::CPlantBoxAdapter<GlobalPosition>(rootSystem);
    }

    // root grid geometry
    const auto& rootGridView = rootGrid->leafGridView();
    using FVGridGeometry = GetPropType<RootTypeTag, Properties::FVGridGeometry>;
    auto rootGridGeometry = std::make_shared<FVGridGeometry>(rootGridView);
    rootGridGeometry->update();

    ////////////////////////////////////////////////////////////
    // run stationary or dynamic problem on this grid
    ////////////////////////////////////////////////////////////

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<SoilTypeTag, RootTypeTag>;
    constexpr auto soilDomainIdx = Traits::template SubDomain<0>::Index();
    constexpr auto rootDomainIdx = Traits::template SubDomain<1>::Index();

    // the coupling manager
    using CouplingManager = GetPropType<SoilTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(soilGridGeometry, rootGridGeometry);

    // the problems
    using SoilProblem = GetPropType<SoilTypeTag, Properties::Problem>;
    auto soilProblem = std::make_shared<SoilProblem>(soilGridGeometry);
    soilProblem->setCouplingManager(&(*couplingManager));
    using RootProblem = GetPropType<RootTypeTag, Properties::Problem>;
    auto rootProblem = std::make_shared<RootProblem>(rootGridGeometry);
    rootProblem->setCouplingManager(&(*couplingManager));

    // check if we are about to restart a previously interrupted simulation
    double restartTime = getParam<double>("Restart.Time", 0);

    // the solution vector
    Traits::SolutionVector sol;
    sol[soilDomainIdx].resize(soilGridGeometry->numDofs());
    sol[rootDomainIdx].resize(rootGridGeometry->numDofs());
    soilProblem->applyInitialSolution(sol[soilDomainIdx]);
    rootProblem->applyInitialSolution(sol[rootDomainIdx]);
    auto oldSol = sol;

    // initial root growth
    GrowthModule::GridGrowth<RootTypeTag>* gridGrowth = nullptr;
    double initialTime = 0.; // s
    if (simtype==Properties::rootbox) {
        gridGrowth = new GrowthModule::GridGrowth<RootTypeTag>(rootGrid, rootGridGeometry, growth, sol[rootDomainIdx]); // in growth/gridgrowth.hh
        std::cout << "...grid grower initialized \n" << std::flush;
        initialTime = getParam<double>("RootSystem.Grid.InitialT")*24*3600;
        std::cout << "intialT " << initialTime/24/3600 << "\n" << std::flush;
        gridGrowth->grow(initialTime); ////////////////////////////////////////////////////////////////////
        std::cout << "\ninitial growth performed... \n" << std::flush;
    }

    // obtain parameters from the CPlantBox or dgf
    if (simtype==Properties::dgf) {
        rootProblem->spatialParams().initParameters(*rootGridManager.getGridData());
    } else if (simtype==Properties::rootbox){
        rootProblem->spatialParams().updateParameters(*growth);
    }

    // the solution vector
    sol[soilDomainIdx].resize(soilGridGeometry->numDofs());
    sol[rootDomainIdx].resize(rootGridGeometry->numDofs());
    if (restartTime > 0)
        {
        // soil
        using soilIOFields = GetPropType<SoilTypeTag, Properties::IOFields>;
        using soilPrimaryVariables = GetPropType<SoilTypeTag, Properties::PrimaryVariables>;
        using soilModelTraits = GetPropType<SoilTypeTag, Properties::ModelTraits>;
        using soilFluidSystem = GetPropType<SoilTypeTag, Properties::FluidSystem>;
        const auto soilfileName = getParam<std::string>("Restart.SoilFile");
        const auto soilpvName = createPVNameFunction<soilIOFields, soilPrimaryVariables, soilModelTraits, soilFluidSystem>();
        loadSolution(sol[soilDomainIdx], soilfileName, soilpvName, *soilGridGeometry);

        // root
        /*using rootIOFields = GetPropType<RootTypeTag, Properties::IOFields>;
        using rootPrimaryVariables = GetPropType<RootTypeTag, Properties::PrimaryVariables>;
        using rootModelTraits = GetPropType<RootTypeTag, Properties::ModelTraits>;
        using rootFluidSystem = GetPropType<RootTypeTag, Properties::FluidSystem>;
        const auto rootfileName = getParam<std::string>("Restart.RootFile");
        const auto rootpvName = createPVNameFunction<rootIOFields, rootPrimaryVariables, rootModelTraits, rootFluidSystem>();
        loadSolution(sol[rootDomainIdx], rootfileName, rootpvName, *rootGridGeometry);*/
        rootProblem->applyInitialSolution(sol[rootDomainIdx]);
        }
    else
        {
        soilProblem->applyInitialSolution(sol[soilDomainIdx]);
        rootProblem->applyInitialSolution(sol[rootDomainIdx]);
        }
    oldSol = sol;

    // coupling manager
    couplingManager->init(soilProblem, rootProblem, sol);
    soilProblem->computePointSourceMap();
    rootProblem->computePointSourceMap();

    // the grid variables
    using SoilGridVariables = GetPropType<SoilTypeTag, Properties::GridVariables>;
    auto soilGridVariables = std::make_shared<SoilGridVariables>(soilProblem, soilGridGeometry);
    soilGridVariables->init(sol[soilDomainIdx]);
    using RootGridVariables = GetPropType<RootTypeTag, Properties::GridVariables>;
    auto rootGridVariables = std::make_shared<RootGridVariables>(rootProblem, rootGridGeometry);
    rootGridVariables->init(sol[rootDomainIdx]);

    // update the saturation vector
    // RootSoil::updateSaturation(saturation, *soilGridGeoemtry, *soilGridVariables, sol[soilDomainIdx]);

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
            std::cout << "using "<< checkPoints.size() << "check times \n";
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

    // intialize the vtk output module
    using SoilSolution = std::decay_t<decltype(sol[soilDomainIdx])>;
    VtkOutputModule<SoilGridVariables, SoilSolution> soilVtkWriter(*soilGridVariables, sol[soilDomainIdx], soilProblem->name());
    GetPropType<SoilTypeTag, Properties::VtkOutputFields>::initOutputModule(soilVtkWriter);
    soilVtkWriter.write(0.0);

    using RootSolution = std::decay_t<decltype(sol[rootDomainIdx])>;
    VtkOutputModule<RootGridVariables, RootSolution> rootVtkWriter(*rootGridVariables, sol[rootDomainIdx], rootProblem->name());
    GetPropType<RootTypeTag, Properties::VtkOutputFields>::initOutputModule(rootVtkWriter);

    rootProblem->userData("pSoil", sol[rootDomainIdx]);
    rootProblem->userData("radius", sol[rootDomainIdx]);
    rootProblem->userData("order", sol[rootDomainIdx]);
    rootProblem->userData("id", sol[rootDomainIdx]);
    rootProblem->userData("axialFlux", sol[rootDomainIdx]); // todo wrong (coarse approximation)
    rootProblem->userData("radialFlux", sol[rootDomainIdx]);
    rootProblem->userData("age", sol[rootDomainIdx]);
    rootProblem->userData("initialPressure",sol[rootDomainIdx]);
    rootProblem->userData("kr", sol[rootDomainIdx]);
    rootProblem->userData("kx", sol[rootDomainIdx]);
    rootVtkWriter.addField(rootProblem->p(), "p [cm]");
    rootVtkWriter.addField(rootProblem->radius(), "radius [m]"); // not in cm, because of tube plot
    rootVtkWriter.addField(rootProblem->order(), "order [1]");
    rootVtkWriter.addField(rootProblem->id(), "id [1]");
    rootVtkWriter.addField(rootProblem->axialFlux(), "axial flux [cm3/d]");
    rootVtkWriter.addField(rootProblem->radialFlux(), "radial flux [cm3/d]");
    rootVtkWriter.addField(rootProblem->age(), "age [d]");
    rootVtkWriter.addField(rootProblem->initialPressure(), "initial pressure [cm]");
    rootVtkWriter.addField(rootProblem->kr(), "kr [cm/hPa/d]");
    rootVtkWriter.addField(rootProblem->kx(), "kx [cm4/hPa/day]");
    rootVtkWriter.write(restartTime);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    std::shared_ptr<Assembler> assembler;
    if (tEnd > 0) {
        assembler = std::make_shared<Assembler>(std::make_tuple(soilProblem, rootProblem),
            std::make_tuple(soilGridGeometry, rootGridGeometry),
            std::make_tuple(soilGridVariables, rootGridVariables),
            couplingManager, timeLoop); // dynamic
    } else {
        assembler = std::make_shared<Assembler>(std::make_tuple(soilProblem, rootProblem),
            std::make_tuple(soilGridGeometry, rootGridGeometry),
            std::make_tuple(soilGridVariables, rootGridVariables),
            couplingManager); // static
    }

    // the linear solver
    using LinearSolver = BlockDiagILU0BiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    std::cout << "\ni plan to actually start \n" << std::flush;
    if (tEnd > 0) // dynamic
    {
        std::cout << "a time dependent model\n\n" << std::flush;
        timeLoop->start();
        do {

            double t = timeLoop->time(); // dumux time
            double dt = timeLoop->timeStepSize(); // dumux time step
            rootProblem->setTime(t, dt); // pass current time to the root problem
            rootProblem->postTimeStep(sol[rootDomainIdx], *rootGridVariables);
            rootProblem->writeTranspirationRate(); // add transpiration data into the text file
            soilProblem->setTime(t, dt);
            soilProblem->postTimeStep(sol[soilDomainIdx], *soilGridVariables);
            soilProblem->writeBoundaryFluxes();

            if (simtype==Properties::rootbox) {
                if (grow) {

                    // std::cout << "time " << growth->simTime()/24/3600 << " < " << (t+initialTime)/24/3600 << "\n";
                    while (growth->simTime()+dt<t+initialTime) {

                        std::cout << "grow \n"<< std::flush;

                        gridGrowth->grow(dt);
                        rootProblem->spatialParams().updateParameters(*growth);
                        rootProblem->applyInitialSolution(sol[rootDomainIdx]);
                        rootGridVariables->updateAfterGridAdaption(sol[rootDomainIdx]); // update the secondary variables

                        couplingManager->updateAfterGridAdaption(soilGridGeometry, rootGridGeometry);
                        couplingManager->init(soilProblem, rootProblem, sol); // recompute coupling maps
                        couplingManager->updateSolution(sol); // update the solution vector for the coupling manager

                        soilProblem->computePointSourceMap(); // recompute the coupling sources
                        rootProblem->computePointSourceMap(); // recompute the coupling sources

                        assembler->setJacobianPattern(assembler->jacobian()); // resize and set Jacobian pattern
                        assembler->setResidualSize(assembler->residual()); // resize residual vector

                        oldSol[rootDomainIdx] = sol[rootDomainIdx]; // // update old solution to new grid

                        std::cout << "grew \n"<< std::flush;
                    }

                }
            }

            // set previous solution for storage evaluations
            assembler->setPreviousSolution(oldSol);

            nonLinearSolver.solve(sol, *timeLoop);

            soilControl(*soilGridGeometry, *soilGridVariables, sol[soilDomainIdx], oldSol[soilDomainIdx], t, dt); //debugging soil water content

            // make the new solution the old solution
            oldSol = sol;
            soilGridVariables->advanceTimeStep();
            rootGridVariables->advanceTimeStep();

            // update the saturation vector
            // updateSaturation(saturation, *soilGridGeoemtry, *soilGridVariables, sol[soilDomainIdx]);

            timeLoop->advanceTimeStep(); // advance to the time loop to the next step

            if ((timeLoop->isCheckPoint()) || (timeLoop->finished())) { // write vtk output (only at check points)
                if (grow) { // prepare static fields also
                    rootProblem->userData("radius", sol[rootDomainIdx]);
                    rootProblem->userData("order", sol[rootDomainIdx]);
                    rootProblem->userData("id", sol[rootDomainIdx]);
                    rootProblem->userData("initialPressure", sol[rootDomainIdx]);
                }
                rootProblem->userData("p", sol[rootDomainIdx]);
                rootProblem->userData("axialFlux", sol[rootDomainIdx]);
                rootProblem->userData("radialFlux", sol[rootDomainIdx]);
                rootProblem->userData("age", sol[rootDomainIdx]); // age changes with time
                rootProblem->userData("kr", sol[rootDomainIdx]);  // conductivities might change with age
                rootProblem->userData("kx", sol[rootDomainIdx]);
                rootVtkWriter.write(timeLoop->time());
                soilVtkWriter.write(timeLoop->time());
            }
            soilProblem->computeSourceIntegral(sol[soilDomainIdx], *soilGridVariables);
            rootProblem->computeSourceIntegral(sol[rootDomainIdx], *rootGridVariables);
            std::cout << "\n";

            timeLoop->reportTimeStep();  // report statistics of this time step
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize())); // set new dt as suggested by the newton solver

        } while (!timeLoop->finished());

        timeLoop->finalize();

    } else { // static

        std::cout << "a static model \n\n" << std::flush;

        assembler->setPreviousSolution(oldSol); // set previous solution for storage evaluations
        nonLinearSolver.solve(sol); // solve the non-linear system

        // write outputs
        rootProblem->userData("p", sol[rootDomainIdx]);
        rootProblem->userData("axialFlux", sol[rootDomainIdx]);
        rootProblem->userData("radialFlux", sol[rootDomainIdx]);
        rootProblem->userData("age", sol[rootDomainIdx]); // prepare fields
        rootProblem->userData("kr", sol[rootDomainIdx]);  // conductivities change with age
        rootProblem->userData("kx", sol[rootDomainIdx]);
        rootVtkWriter.write(1); // write vtk output
        soilVtkWriter.write(1);
        rootProblem->postTimeStep(sol[rootDomainIdx], *rootGridVariables);
        rootProblem->writeTranspirationRate();
        soilProblem->postTimeStep(sol[soilDomainIdx], *soilGridVariables);
        soilProblem->writeBoundaryFluxes();
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

