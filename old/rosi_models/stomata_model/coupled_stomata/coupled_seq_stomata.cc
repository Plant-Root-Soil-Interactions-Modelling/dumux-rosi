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
 *
 * We couple root and soil model sequentially (naive approach).
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
#include <dune/grid/yaspgrid/coordinates.hh>
#include <dune/grid/common/rangegenerators.hh> // elements, vertices, scv, scvfs, ...
#include <dune/istl/io.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/common/timeloop.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh> // solver
#include <dumux/porousmediumflow/richards/newtonsolver.hh> // solver

#include <dumux/assembly/fvassembler.hh>
//#include <dumux/assembly/diffmethod.hh>
//#include <dumux/discretization/method.hh>

// dumux-rosi
#include <dumux/growth/soillookup.hh> // for coupling
#include <dumux/growth/rootsystemgridfactory.hh>
#include <dumux/growth/growthinterface.hh>
#include <dumux/growth/crootboxadapter.hh>
#include <dumux/growth/gridgrowth.hh>

#include <RootSystem.h>
#include "../rootsystem_stomata/rootsproblem_stomata.hh"
#include "../soil_stomata/richardsproblem_stomata.hh"
#include "../rootsystem_stomata/properties_stomata.hh" // TypeTag:Roots
#include "../rootsystem_stomata/properties_nocoupling_stomata.hh" // dummy types for replacing the coupling types
#include "../soil_stomata/properties_stomata.hh" // TypeTag:RichardsTT
#include "../soil_stomata/properties_nocoupling_stomata.hh" // dummy types for replacing the coupling types

namespace Dumux {

using RootTypeTag = Properties::TTag::RootsBox;
using SoilTypeTag = Properties::TTag::RichardsBox;

//using SoilGridType = Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>; // pick soil grid here (its in compile definition in the soil model)
//using RootGridType = GetPropType<RootTypeTag, Properties::Grid>;

using RootFVGridGeometry = GetPropType<RootTypeTag, Properties::FVGridGeometry>;
using SoilFVGridGeometry = GetPropType<SoilTypeTag, Properties::FVGridGeometry>;

using SoilLookUp = GrowthModule::SoilLookUpBBoxTree<SoilFVGridGeometry>;

/**
 * debugging
 */
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
    std::cout << "Water in domain: " << cVol*1.e3 << " kg at time " << t/24/3600 << " days\n";
    std::cout << "change of " << (oldVol-cVol)*1.e3 << " kg = " << (oldVol-cVol)*1.e3/dt << " kg/s \n" ;
}


/*!
 * calculates saturation from the solution vector
 */
template<class SoilGridVariables, class SoilSolution>
void updateSaturation(std::vector<double>& saturation, const SoilFVGridGeometry& gridGeometry, const SoilGridVariables& gridVariables,
    const SoilSolution& sol) {
    double vol = 0;

    const auto& gridView = gridGeometry.gridView();  // soil
    for (const auto& element : elements(gridView)) { // soil elements

        auto fvGeometry = localView(gridGeometry); // soil solution -> volume variable
        fvGeometry.bindElement(element);
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (const auto& scv : scvs(fvGeometry)) { // i dont quite get that..
            saturation[scv.dofIndex()] = elemVolVars[scv].saturation(0);
            vol += elemVolVars[scv].saturation(0)*scv.volume();
        }
    }
    std::cout << "updateSaturation: Water volume: " << vol << " m3 \n";
}

/*!
 * calculates the radial fluxes of the segments
 */
template< class RootGridVariables, class RootSolution, class RootProlem>
void radialFlux2soilSink(std::vector<double>& source, const RootFVGridGeometry& gridGeometry, const RootGridVariables& gridVariables,
    const RootSolution& r, const RootProlem& rootProblem, SoilLookUp* soilLookUp) {

    std::fill(source.begin(), source.end(), 0);
    const auto& gridView = gridGeometry.gridView(); // root

    for (const auto& element : elements(gridView)) { // root elementsrootbox

        // root solution -> volume variable
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        elemVolVars.bindElement(element, fvGeometry, r);

        for (const auto& scv : scvs(fvGeometry)) {  //root sub control volumes
            double s = -rootProblem.source(element, gridGeometry, elemVolVars, scv); // pass to problem class [kg/s/m^3]
            s *=  scv.volume()*elemVolVars[scv].extrusionFactor();  // [kg / s]
            auto pos = scv.center();
            int eIdx = soilLookUp->pick(pos); // find element index of soil, for root each root element
            if (eIdx>=0) { // just to be sure...
                source.at(eIdx) += s; //  [kg/s] accumulate source term
            } else {
                std::cout << "root at position " << pos << " not within soil";
            }
        }
    }
    double sum = std::accumulate(source.begin(), source.end(), 0.);
    std::cout << "radialFlux2soilSink: summed source: " << sum << " [kg / s] \n";

}

} // end namespace Dumux



int main(int argc, char** argv) try
{
    using namespace Dumux;

    int simtype = Properties::simtype; // definition in ../rootsystem/properties.hh, use DGF or ROOTBOX in CMakeLists.txt

    // initialize MPI, finalize is done automatically#include <dumux/growth/soillookup.hh> on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0) {
        DumuxMessage::print(/*firstCall=*/true);
    }

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    std::string rootName = getParam<std::string>("Problem.RootName");
    Parameters::init(0, argv, rootName);
    std::string soilName = getParam<std::string>("Problem.SoilName");
    Parameters::init(0, argv, soilName);
    Parameters::init(argc, argv); // reread to overwrite parameters

    // soil grid
    GridManager<GetPropType<SoilTypeTag, Properties::Grid>> soilGridManager;
    soilGridManager.init("Soil"); // pass parameter group (see input file)

    // soil grid geometry
    const auto& soilGridView = soilGridManager.grid().leafGridView();
    using SoilFVGridGeometry = GetPropType<SoilTypeTag, Properties::FVGridGeometry>;
    auto soilGridGeometry = std::make_shared<SoilFVGridGeometry>(soilGridView);
    soilGridGeometry->update();

    // Create the gridmanager and grid
    using GlobalPosition = Dune::FieldVector<double, 3>;
    using Grid = Dune::FoamGrid<1, 3>;
    std::shared_ptr<Grid> rootGrid;
    GridManager<Grid> gridManager; // only for dgf
    std::shared_ptr<CRootBox::RootSystem> rootSystem; // only for rootbox
    GrowthModule::GrowthInterface<GlobalPosition>* growth = nullptr; // in case of RootBox (or in future PlantBox)
    if (simtype==Properties::dgf) { // for a static dgf grid
        std::cout << "\nSimulation type is dgf \n\n" << std::flush;
        gridManager.init("RootSystem");
        rootGrid = std::shared_ptr<Grid>(&gridManager.grid(), Properties::empty_delete<Grid>());
    } else if (simtype==Properties::rootbox) { // for a root model (static or dynamic)
        std::cout << "\nSimulation type is RootBox \n\n" << std::flush;
        rootSystem = std::make_shared<CRootBox::RootSystem>();
        rootSystem->openFile(getParam<std::string>("RootSystem.Grid.File"), "modelparameter/");
        if (hasParam("RootSystem.Grid.Confined")) {
            auto box = getParam<std::vector<double>>("RootSystem.Grid.Confined");
            rootSystem->setGeometry(new CRootBox::SDF_PlantBox(box.at(0)*100, box.at(1)*100, box.at(2)*100));
        } else { // half plane
            rootSystem->setGeometry(new CRootBox::SDF_HalfPlane(CRootBox::Vector3d(0.,0.,0.5), CRootBox::Vector3d(0.,0.,1.))); // care, collar needs to be top, make sure plant seed is located below -1 cm
        }
        rootSystem->initialize();
        double shootZ = getParam<double>("RootSystem.Grid.ShootZ", 0.); // root system initial time
        rootGrid = GrowthModule::RootSystemGridFactory::makeGrid(*rootSystem, shootZ, true); // in dumux/growth/rootsystemgridfactory.hh
        //  todo static soil for hydrotropsim ...
        //    auto soilLookup = SoilLookUpBBoxTree<GrowthModule::Grid> (soilGridView, soilGridGeoemtry->boundingBoxTree(), saturation);
        //    rootSystem->setSoil(&soilLookup);
        growth = new GrowthModule::CRootBoxAdapter<GlobalPosition>(*rootSystem);
    }

    // root grid geometry
    const auto& rootGridView = rootGrid->leafGridView();
    using FVGridGeometry = GetPropType<RootTypeTag, Properties::FVGridGeometry>;
    auto rootGridGeometry = std::make_shared<FVGridGeometry>(rootGridView);
    rootGridGeometry->update();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // the problems (initial and boundary conditions)
    using RootProblem = GetPropType<RootTypeTag, Properties::Problem>;
    auto rootProblem = std::make_shared<RootProblem>(rootGridGeometry);

    using SoilProblem = GetPropType<SoilTypeTag, Properties::Problem>;
    auto soilProblem = std::make_shared<SoilProblem>(soilGridGeometry);
    std::cout << "\nI have two problems " << "\n" << std::flush;

    // the solution vector
    using RootSolutionVector = GetPropType<RootTypeTag, Properties::SolutionVector>;
    RootSolutionVector r(rootGridGeometry->numDofs());
    rootProblem->applyInitialSolution(r);
    auto rOld = r;
    using SoilSolutionVector = GetPropType<SoilTypeTag, Properties::SolutionVector>;
    SoilSolutionVector s(soilGridGeometry->numDofs());
    soilProblem->applyInitialSolution(s);
    auto sOld = s;
    std::cout << "no solution yet \n" << std::flush;

    // initial root growth
    GrowthModule::GridGrowth<RootTypeTag>* gridGrowth = nullptr;
    double initialTime = 0.; // s
    if (simtype==Properties::rootbox) {
        gridGrowth = new GrowthModule::GridGrowth<RootTypeTag>(rootGrid, rootGridGeometry, growth, r); // in growth/gridgrowth.hh
        std::cout << "...grid grower initialized \n" << std::flush;
        initialTime = getParam<double>("RootSystem.Grid.InitialT")*24*3600;
        gridGrowth->grow(initialTime);
        std::cout << "\ninitial growth performed... \n" << std::flush;
    }

    // obtain parameters from the crootbox or dgf
    if (simtype==Properties::dgf) {
        rootProblem->spatialParams().initParameters(*gridManager.getGridData());
    } else if (simtype==Properties::rootbox){
        rootProblem->spatialParams().updateParameters(*growth);
    }

    // the grid variables
    using RootGridVariables = GetPropType<RootTypeTag, Properties::GridVariables>;
    auto rootGridVariables = std::make_shared<RootGridVariables>(rootProblem, rootGridGeometry);
    rootGridVariables->init(r);

    using SoilGridVariables = GetPropType<SoilTypeTag, Properties::GridVariables>;
    auto soilGridVariables = std::make_shared<SoilGridVariables>(soilProblem, soilGridGeometry);
    soilGridVariables->init(s);
    std::cout << "... but grid variables \n" << std::flush;

    // sequentially coupling the problems ...
    static constexpr int soilDim = SoilFVGridGeometry::GridView::dimension;
    std::vector<double> saturation(soilGridView.size(soilDim), 1.0);
    updateSaturation(saturation, *soilGridGeometry, *soilGridVariables, s);
    auto soilLookUp = SoilLookUp(*soilGridGeometry, saturation);

    rootProblem->setSoil(&soilLookUp);
    std::cout << "roots know the soil \n" << "\n" << std::flush;

    std::vector<double> soilSink = std::vector<double>(soilGridGeometry->gridView().size(0));
    radialFlux2soilSink(soilSink, *rootGridGeometry, *rootGridVariables, r, *rootProblem, &soilLookUp); // precomputes the sink for the soil problem
    soilProblem->setSource(&soilSink);
    std::cout << "and the soil knows the roots \n" << "\n" << std::flush;
    std::cout << "value = " << soilLookUp.getValue(CRootBox::Vector3d()) << "\n";

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

    std::cout << "\npress button \n";
    std::string sss = "";     std::getline(std::cin, sss); // for debugging

    // intialize the vtk output module
    using RootIOFields = GetPropType<RootTypeTag, Properties::IOFields>;
    VtkOutputModule<RootGridVariables, RootSolutionVector> rootVTKWriter(*rootGridVariables, r, rootProblem->name()+"R");
    using RootVelocityOutput = GetPropType<RootTypeTag, Properties::VelocityOutput>;
    rootProblem->userData("p", r);
    rootProblem->userData("radius", r);
    rootProblem->userData("order", r);
    rootProblem->userData("id", r);
    rootProblem->userData("axialFlux", r); // todo wrong (coarse approximation)
    rootProblem->userData("radialFlux", r);
    rootProblem->userData("age", r);
    rootProblem->userData("initialPressure", r);
    rootProblem->userData("kr", r);
    rootProblem->userData("kx", r);
    rootVTKWriter.addField(rootProblem->p(), "p [cm]");
    rootVTKWriter.addField(rootProblem->radius(), "radius [m]"); // not in cm, because of tube plot
    rootVTKWriter.addField(rootProblem->order(), "order [1]");
    rootVTKWriter.addField(rootProblem->id(), "id [1]");
    rootVTKWriter.addField(rootProblem->axialFlux(), "axial flux [cm3/d]");
    rootVTKWriter.addField(rootProblem->radialFlux(), "radial flux [cm3/d]");
    rootVTKWriter.addField(rootProblem->age(), "age [d]");
    rootVTKWriter.addField(rootProblem->initialPressure(), "initial pressure [cm]");
    rootVTKWriter.addField(rootProblem->kr(), "kr [cm/hPa/d]");
    rootVTKWriter.addField(rootProblem->kx(), "kx [cm4/hPa/day]");
    rootVTKWriter.write(0.0);
    RootIOFields::initOutputModule(rootVTKWriter); //!< Add model specific output fields
    rootVTKWriter.write(0.0);

    using SoilIOFields = GetPropType<SoilTypeTag, Properties::IOFields>;
    VtkOutputModule<SoilGridVariables, SoilSolutionVector> soilVTKWriter(*soilGridVariables, s, soilProblem->name()+"S");
    using SoilVelocityOutput = GetPropType<SoilTypeTag, Properties::VelocityOutput>;
//    soilVTKWriter.addVelocityOutput(std::make_shared<SoilVelocityOutput>(*soilGridVariables));
    SoilIOFields::initOutputModule(soilVTKWriter); //!< Add model specific output fields
    soilVTKWriter.write(0.0);
    std::cout << "vtk writer module initialized (in less than 20 lines)" << "\n" << std::flush;

    // the assembler with time loop for instationary problem
    using RootAssembler = FVAssembler<RootTypeTag, DiffMethod::numeric>;
    using SoilAssembler = FVAssembler<SoilTypeTag, DiffMethod::numeric>;
    std::shared_ptr<RootAssembler> rootAssembler;
    std::shared_ptr<SoilAssembler> soilAssembler;
    if (tEnd > 0) {
        rootAssembler = std::make_shared<RootAssembler>(rootProblem, rootGridGeometry, rootGridVariables, timeLoop); // dynamic
        soilAssembler = std::make_shared<SoilAssembler>(soilProblem, soilGridGeometry, soilGridVariables, timeLoop); // dynamic
    } else {
        rootAssembler = std::make_shared<RootAssembler>(rootProblem, rootGridGeometry, rootGridVariables); // static
        soilAssembler = std::make_shared<SoilAssembler>(soilProblem, soilGridGeometry, soilGridVariables); // static
    }

    // the linear solver
    using RootLinearSolver = AMGBackend<RootTypeTag>;
    using SoilLinearSolver = AMGBackend<SoilTypeTag>; // the only parallel linear solver available
    auto rootLinearSolver = std::make_shared<RootLinearSolver>(rootGridView, rootGridGeometry->dofMapper());
    auto soilLinearSolver = std::make_shared<SoilLinearSolver>(soilGridView, soilGridGeometry->dofMapper());

    // the non-linear solver
    using RootNewtonSolver = NewtonSolver<RootAssembler, RootLinearSolver>;
    using SoilNonlinearSolver = RichardsNewtonSolver<SoilAssembler, SoilLinearSolver>; //Dumux::RichardsNewtonSolver<SoilAssembler, SoilLinearSolver>;
    RootNewtonSolver rootNonlinearSolver(rootAssembler, rootLinearSolver);
    SoilNonlinearSolver soilNonlinearSolver = SoilNonlinearSolver(soilAssembler, soilLinearSolver);

    std::cout << "i plan to actually start \n" << "\n" << std::flush;
    if (tEnd > 0) // dynamic
    {
        std::cout << "a time dependent model" << "\n" << std::flush;
        timeLoop->start();
        do {

            double t = timeLoop->time(); // dumux time
            double dt = timeLoop->timeStepSize(); // dumux time step
            rootProblem->setTime(t, dt); // pass current time to the root problem
            soilProblem->setTime(t);

            if (grow) {

                // std::cout << "time " << growth->simTime()/24/3600 << " < " << (t+initialTime)/24/3600 << "\n";
                while (growth->simTime()+dt<t+initialTime) {

                    std::cout << "grow \n"<< std::flush;
                    gridGrowth->grow(dt);
                    rootProblem->spatialParams().updateParameters(*growth);
                    rootProblem->applyInitialSolution(r); // reset todo (? does this make sense)?
                    std::cout << "grew \n"<< std::flush;

                    // what shall I update?
                    rootGridGeometry->update();
                    //gridVariables->update();
                    rootGridVariables->updateAfterGridAdaption(r); // update the secondary variables

                    // todo? what is necessary? no clue what i am doing ...
                    rootAssembler->setResidualSize(); // resize residual vector
                    rootAssembler->setJacobianPattern(); // resize and set Jacobian pattern
                    rootAssembler->setPreviousSolution(r);
                    rootAssembler->assembleJacobianAndResidual(r);

                    rOld = r;

                }

            }

            // set previous solution for storage evaluations
            soilAssembler->setPreviousSolution(sOld);
            rootAssembler->setPreviousSolution(rOld);

            // solves the soil problem
            radialFlux2soilSink(soilSink, *rootGridGeometry, *rootGridVariables, rOld, *rootProblem, &soilLookUp); // precomputes the sink for the soil problem
            std::cout << "solve soil\n";
            soilNonlinearSolver.solve(s, *timeLoop);

            // solves the root problem
            updateSaturation(saturation, *soilGridGeometry, *soilGridVariables, sOld); // updates soil look up for the root problem#
            std::cout << "solve roots\n";
            rootNonlinearSolver.solve(r, *timeLoop);

            soilControl(*soilGridGeometry, *soilGridVariables, s, sOld, t, dt); //debugging soil water content


            // make the new solution the old solution
            rOld = r;
            sOld = s;

            rootGridVariables->advanceTimeStep();
            soilGridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();
            // write vtk output (only at check points)
            if ((timeLoop->isCheckPoint()) || (timeLoop->finished())) {
                if (grow) { // prepare static fields also
                    rootProblem->userData("radius", r);
                    rootProblem->userData("order", r);
                    rootProblem->userData("id", r);
                    rootProblem->userData("initialPressure", r);
                }
                rootProblem->userData("p", r);
                rootProblem->userData("axialFlux", r);
                rootProblem->userData("radialFlux", r);
                rootProblem->userData("age", r); // age changes with time
                rootProblem->userData("kr", r);  // conductivities change with age
                rootProblem->userData("kx", r);
                rootVTKWriter.write(timeLoop->time());
                soilVTKWriter.write(timeLoop->time());
            }
            if (mpiHelper.rank() == 0) {
                rootProblem->writeTranspirationRate(r);
            }

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by the newton solver
            timeLoop->setTimeStepSize(rootNonlinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        } while (!timeLoop->finished());

        timeLoop->finalize(rootGridView.comm());

    } else { // static

        std::cout << "a static model" << "\n" << std::flush;
        // set previous solution for storage evaluations
        rootAssembler->setPreviousSolution(rOld);
        soilAssembler->setPreviousSolution(sOld);
        // solve the non-linear system
        rootNonlinearSolver.solve(r);
        soilNonlinearSolver.solve(s);
        // write vtk output
        rootProblem->userData("p", r);
        rootProblem->userData("axialFlux", r);
        rootProblem->userData("radialFlux", r);
        rootProblem->userData("age", r); // prepare fields
        rootProblem->userData("kr", r);  // conductivities change with age
        rootProblem->userData("kx", r);
        rootVTKWriter.write(1); // write vtk output
        soilVTKWriter.write(1);
        rootProblem->writeTranspirationRate(r);

    }

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0) {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

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

