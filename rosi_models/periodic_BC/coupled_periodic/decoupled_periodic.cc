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
 * \brief Coupling (not)
 *
 * A soil and a root model are run in one simulation loop,
 * completely independent of each other (as first step) very boring...
 *
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <memory>

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

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/common/timeloop.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh> // solver
#include <dumux/porousmediumflow/richards/newtonsolver.hh> // solver
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>

#include <RootSystem.h>
#include "../rootsystem_periodic/rootsproblem_periodic.hh"
#include "../soil_periodic/richardsproblem_periodic.hh"
#include "../rootsystem_periodic/properties_periodic.hh" // TypeTag:Roots
#include "../rootsystem_periodic/properties_nocoupling_periodic.hh" // dummy types for replacing the coupling types
#include "../soil_periodic/properties_periodic.hh" // TypeTag:RichardsTT
#include "../soil_periodic/properties_nocoupling_periodic.hh" // dummy types for replacing the coupling types

/**
 * decoupled (for debugging, proof of concept)
 */
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using RootTypeTag = Properties::TTag::RootsBox;
    using SoilTypeTag = Properties::TTag::RichardsBox;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0) {
        DumuxMessage::print(/*firstCall=*/true);
    }

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    std::string soilName = getParam<std::string>("Problem.SoilName");
    Parameters::init(0, argv, soilName);
    std::string rootName = getParam<std::string>("Problem.RootName");
    Parameters::init(0, argv, rootName);

    // try to create a grid (from the given grid file or the input file)

    GridManager<GetPropType<RootTypeTag, Properties::Grid>> rootGridManager;
    rootGridManager.init("RootSystem");
    const auto rootGridData = rootGridManager.getGridData();

    using SoilGridType = Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>; // pick soil grid here
    GridManager<SoilGridType> soilGridManager;
    soilGridManager.init("Soil");

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& rootLGV = rootGridManager.grid().leafGridView();
    const auto& soilLGV = soilGridManager.grid().leafGridView();
    std::cout << "i have two views \n" << "\n" << std::flush;

    // create the finite volume grid geometry
    using RootFVGridGeometry = GetPropType<RootTypeTag, Properties::FVGridGeometry>;
    auto rootFVGridGeometry = std::make_shared<RootFVGridGeometry>(rootLGV);
    rootFVGridGeometry->update();
    using SoilFVGridGeometry = GetPropType<SoilTypeTag, Properties::FVGridGeometry>;
    auto soilFVGridGeometry = std::make_shared<SoilFVGridGeometry>(soilLGV);
    soilFVGridGeometry->update();
    std::cout << "i have two geometries built from the views \n" << "\n" << std::flush;

    // the problem (initial and boundary conditions)
    using RootProblem = GetPropType<RootTypeTag, Properties::Problem>;
    auto rootProblem = std::make_shared<RootProblem>(rootFVGridGeometry);
    rootProblem->spatialParams().initParameters(*rootGridData);
    // rootProblem->spatialParams().analyseRootSystem();
    using SoilProblem = GetPropType<SoilTypeTag, Properties::Problem>;
    auto soilProblem = std::make_shared<SoilProblem>(soilFVGridGeometry);
    std::cout << "and i have two problems \n" << "\n" << std::flush;

    // the solution vector
    using RootSolutionVector = GetPropType<RootTypeTag, Properties::SolutionVector>;
    RootSolutionVector r(rootFVGridGeometry->numDofs());
    rootProblem->applyInitialSolution(r);
    auto rOld = r;
    using SoilSolutionVector = GetPropType<SoilTypeTag, Properties::SolutionVector>;
    SoilSolutionVector s(soilFVGridGeometry->numDofs());
    soilProblem->applyInitialSolution(s);
    auto sOld = s;
    std::cout << "no solution yet \n" << "\n" << std::flush;

    // the grid variables
    using RootGridVariables = GetPropType<RootTypeTag, Properties::GridVariables>;
    auto rootGridVariables = std::make_shared<RootGridVariables>(rootProblem, rootFVGridGeometry);
    rootGridVariables->init(r);
    using SoilGridVariables = GetPropType<SoilTypeTag, Properties::GridVariables>;
    auto soilGridVariables = std::make_shared<SoilGridVariables>(soilProblem, soilFVGridGeometry);
    soilGridVariables->init(s);
    std::cout << "... but variables \n" << "\n" << std::flush;

    // get some time loop parameters & instantiate time loop
    using Scalar = GetPropType<RootTypeTag, Properties::Scalar>;
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
            std::cout << "decoupled.cc: no check times (TimeLoop.CheckTimes) defined in the input file\n";
        }
    } else { // static
    }
    std::cout << "time might be an issue \n" << "\n" << std::flush;

    // intialize the vtk output module
    using RootIOFields = GetPropType<RootTypeTag, Properties::IOFields>;
    VtkOutputModule<RootGridVariables, RootSolutionVector> rootVTKWriter(*rootGridVariables, r, rootProblem->name()+"R");
    using RootVelocityOutput = GetPropType<RootTypeTag, Properties::VelocityOutput>;
    rootVTKWriter.addVelocityOutput(std::make_shared<RootVelocityOutput>(*rootGridVariables));
    rootProblem->userData("axialFlux", r); // prepare fields // todo wrong (coarse approximation)
    rootProblem->userData("radialFlux", r); // prepare fields
    rootVTKWriter.addField(rootProblem->axialFlux(), "axial flux [cm3/d]");
    rootVTKWriter.addField(rootProblem->radialFlux(), "radial flux [cm3/d]");
    RootIOFields::initOutputModule(rootVTKWriter); //!< Add model specific output fields
    rootVTKWriter.write(0.0);
    using SoilIOFields = GetPropType<SoilTypeTag, Properties::IOFields>;
    VtkOutputModule<SoilGridVariables, SoilSolutionVector> soilVTKWriter(*soilGridVariables, s, soilProblem->name()+"S");
    using SoilVelocityOutput = GetPropType<SoilTypeTag, Properties::VelocityOutput>;
    soilVTKWriter.addVelocityOutput(std::make_shared<SoilVelocityOutput>(*soilGridVariables));
    SoilIOFields::initOutputModule(soilVTKWriter); //!< Add model specific output fields // TODO not workin...
    soilVTKWriter.write(0.0);
    std::cout << "vtk writer module initialized (in less than 20 lines)" << "\n" << std::flush;

    // the assembler with time loop for instationary problem
    using RootAssembler = FVAssembler<RootTypeTag, DiffMethod::numeric>;
    std::shared_ptr<RootAssembler> rootAssembler;
    if (tEnd > 0) {
        rootAssembler = std::make_shared<RootAssembler>(rootProblem, rootFVGridGeometry, rootGridVariables, timeLoop); // dynamic
    } else {
        rootAssembler = std::make_shared<RootAssembler>(rootProblem, rootFVGridGeometry, rootGridVariables); // static
    }
    using SoilAssembler = FVAssembler<SoilTypeTag, DiffMethod::numeric>;
    std::shared_ptr<SoilAssembler> soilAssembler;
    if (tEnd>0) {
        soilAssembler = std::make_shared<SoilAssembler>(soilProblem, soilFVGridGeometry, soilGridVariables, timeLoop); // dynamic
    } else {
        soilAssembler = std::make_shared<SoilAssembler>(soilProblem, soilFVGridGeometry, soilGridVariables); // static
    }

    // the linear solver
    using RootLinearSolver = AMGBackend<RootTypeTag>;
    auto rootLinearSolver = std::make_shared<RootLinearSolver>(rootLGV, rootFVGridGeometry->dofMapper());
    using SoilLinearSolver = Dumux::AMGBackend<SoilTypeTag>; // the only parallel linear solver available
    auto soilLinearSolver = std::make_shared<SoilLinearSolver>(soilLGV, soilFVGridGeometry->dofMapper());

    // the non-linear solver
    using RootNewtonSolver = Dumux::NewtonSolver<RootAssembler, RootLinearSolver>;
    RootNewtonSolver rootNonlinearSolver(rootAssembler, rootLinearSolver);
    using SoilNonlinearSolver = Dumux::RichardsNewtonSolver<SoilAssembler, SoilLinearSolver>; //Dumux::RichardsNewtonSolver<SoilAssembler, SoilLinearSolver>;
    SoilNonlinearSolver soilNonlinearSolver = SoilNonlinearSolver(soilAssembler, soilLinearSolver);

    std::cout << "i plan to actually start \n" << "\n" << std::flush;
    if (tEnd > 0) // dynamic
    {
        std::cout << "a time dependent model" << "\n" << std::flush;
        timeLoop->start();
        do {
            // set previous solution for storage evaluations
            rootAssembler->setPreviousSolution(rOld);
            soilAssembler->setPreviousSolution(sOld);
            // solve the non-linear system with time step control
            rootNonlinearSolver.solve(r, *timeLoop);
            soilNonlinearSolver.solve(s, *timeLoop);
            // make the new solution the old solution
            rOld = r;
            // sOld = s;
            rootGridVariables->advanceTimeStep();
            soilGridVariables->advanceTimeStep();
            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();
            // write vtk output (only at check points)
            if ((timeLoop->isCheckPoint()) || (timeLoop->finished())) {
                rootProblem->userData("axialFlux", r); // prepare fields // todo wrong (coarse approximation)
                rootProblem->userData("radialFlux", r); // prepare fields
                rootVTKWriter.write(timeLoop->time());
                soilVTKWriter.write(timeLoop->time());
            }
            if (mpiHelper.rank() == 0) {
                rootProblem->writeTranspirationRate(r);
            }

            timeLoop->reportTimeStep(); // report statistics of this time step

            double t = timeLoop->time();
            double dt = rootNonlinearSolver.suggestTimeStepSize(timeLoop->timeStepSize());
            timeLoop->setTimeStepSize(dt); // set new dt as suggested by the newton solver
            rootProblem->setTime(t, dt); // pass current time to the problem
            soilProblem->setTime(t);

        } while (!timeLoop->finished());
        timeLoop->finalize(rootLGV.comm());
    } else // static
    {
        std::cout << "a static model" << "\n" << std::flush;
        // set previous solution for storage evaluations
        rootAssembler->setPreviousSolution(rOld);
        soilAssembler->setPreviousSolution(sOld);
        // solve the non-linear system
        rootNonlinearSolver.solve(r);
        soilNonlinearSolver.solve(s);
        // write vtk output
        rootProblem->userData("axialFlux", r); // prepare fields // todo wrong (coarse approximation)
        rootProblem->userData("radialFlux", r); // prepare fields
        rootProblem->writeTranspirationRate(r);
        rootVTKWriter.write(1);
        soilVTKWriter.write(1);
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

