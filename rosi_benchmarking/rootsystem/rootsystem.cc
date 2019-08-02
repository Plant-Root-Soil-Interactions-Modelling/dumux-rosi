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
#include <dumux/growth/crootboxadapter.hh>
#include <dumux/growth/gridgrowth.hh>

#include "rootsproblem.hh"



/**
 * Compile definitions are either DGF or ROOTBOX defined in CMakeLists
 */
enum modelType { dgf=0, rootbox=1 };

/**
 * Pick either RootSpatialParamsDGF (for static dgf files),
 * or RootSpatialParamsRB (for dynamic root growth) as SpatialParams.type,
 */
namespace Dumux {
namespace Properties {
#if DGF
template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::Roots> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsDGF<FVGridGeometry, Scalar>;
};
int simtype = dgf;
#endif
#if ROOTBOX
template<class TypeTag> // Set the spatial parameters
struct SpatialParams<TypeTag, TTag::Roots> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootSpatialParamsRB<FVGridGeometry, Scalar>;
};
int simtype = rootbox;
#endif
} // end namespace Properties
} // end namespace Dumux



/**
 * to wrap a raw pointer into a shared pointer:
 * for not deleting it twice, an empty deleter must be defined
 */
template <typename T>
struct empty_delete {
    empty_delete() /* noexcept */
    { }
    template <typename U>
    empty_delete(const empty_delete<U>&,
        typename std::enable_if<
            std::is_convertible<U*, T*>::value
        >::type* = nullptr) /* noexcept */
    { }
    void operator()(T* const) const /* noexcept */
    { }// do nothing
};

/**
 *
 */

//template <class Assembler, class LinearSolver>
//class MyNewton :public Dumux::NewtonSolver<Assembler,LinearSolver> {
//
//    using GlobalPosition = Dune::FieldVector<double, 3>;
//
//public:
//
//    virtual ~MyNewton() { }
//
//    virtual void newtonFail(SolutionVector& u) {
//        std::cout << "i failed \n";
//        grow->restore();
//    }
//
//    GrowthModule::GrowthInterface<GlobalPosition>* grow;
//
//
//};



/**
 *
 */
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::RootsBox; // RootsCC, RootsBox (TypeTag is defined in the problem class richardsproblem.hh)
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
    std::shared_ptr<CRootBox::RootSystem> rootSystem; // only for rootbox
    GrowthModule::GrowthInterface<GlobalPosition>* growth = nullptr; // in case of RootBox (or in future PlantBox)
    if (simtype==dgf) { // for a static dgf grid
        std::cout << "\nSimulation type is dgf \n\n" << std::flush;
        gridManager.init("RootSystem");
        grid = std::shared_ptr<Grid>(&gridManager.grid(), empty_delete<Grid>());
    } else if (simtype==rootbox) { // for a root model (static or dynamic)
        std::cout << "\nSimulation type is RootBox \n\n" << std::flush;
        rootSystem = std::make_shared<CRootBox::RootSystem>();
        rootSystem->openFile(getParam<std::string>("RootSystem.Grid.File"), "modelparameter/");
        rootSystem->setGeometry(new CRootBox::SDF_HalfPlane(CRootBox::Vector3d(0.,0.,0.5), CRootBox::Vector3d(0.,0.,1.))); // care, collar needs to be top, make sure plant seed is located below -1 cm
        rootSystem->initialize();
        double shootZ = getParam<double>("RootSystem.Grid.ShootZ", 0.); // root system initial time
        grid = GrowthModule::RootSystemGridFactory::makeGrid(*rootSystem, shootZ, true); // in dumux/growth/rootsystemgridfactory.hh
        //  todo static soil for hydrotropsim ...
        //    auto soilLookup = SoilLookUpBBoxTree<GrowthModule::Grid> (soilGridView, soilGridGeoemtry->boundingBoxTree(), saturation);
        //    rootSystem->setSoil(&soilLookup);
        growth = new GrowthModule::CRootBoxAdapter<GlobalPosition>(*rootSystem);
    }

    ////////////////////////////////////////////////////////////
    // run stationary or dynamic problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = grid->leafGridView();
    std::cout << "i have the view \n"<< std::flush;

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();
    std::cout << "i have the geometry \n" << std::flush;

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>; // defined in discretization/fvproperties.hh, as Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>
    SolutionVector x(fvGridGeometry->numDofs()); // degrees of freedoms

    // root growth
    GrowthModule::GridGrowth<TypeTag>* gridGrowth = nullptr;
    double initialTime = 0.; // s
    if (simtype==rootbox) {
        gridGrowth = new GrowthModule::GridGrowth<TypeTag>(grid, fvGridGeometry, growth, x); // in growth/gridgrowth.hh
        std::cout << "...grid grower initialized \n" << std::flush;
        initialTime = getParam<double>("RootSystem.Grid.InitialT")*24*3600;
        gridGrowth->grow(initialTime);
        std::cout << "initial growth performed... \n" << std::flush;
    }

    // the problem (initial and boundary conditions)
    auto problem = std::make_shared<RootsProblem<TypeTag>>(fvGridGeometry);
    if (simtype==dgf) {
        problem->spatialParams().initParameters(*gridManager.getGridData());
    } else if (simtype==rootbox){
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
    std::cout << "time might be an issue \n" << std::flush;

    // intialize the vtk output module
    std::cout << "vtk writer module... \n" << std::flush;
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    problem->userData("radius", x); // prepare fields
    problem->userData("order", x); // prepare fields
    problem->userData("id", x); // prepare fields
    problem->userData("axialFlux", x); // prepare fields // todo wrong (coarse approximation)
    problem->userData("radialFlux", x); // prepare fields
    problem->userData("age", x); // prepare fields
    problem->userData("initialPressure",x); //prepare fields
    problem->userData("kr", x); // prepare fields
    problem->userData("kx", x); //prepare fields
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

    double dt2 = 3600; // root box time step todo this is NEW

    if (tEnd > 0) // dynamic
    {
        std::cout << "a time dependent model\n\n" << std::flush;
        timeLoop->start();
        do {

            double t = timeLoop->time(); // dumux time
            double dt = timeLoop->timeStepSize(); // dumux time step
            problem->setTime(t, dt); // pass current time to the problem

            if (grow) {

                std::cout << "time " << growth->simTime()/24/3600 << " < " << (t+initialTime)/24/3600 << "\n";
                while (growth->simTime()<t+initialTime) {

                    std::cout << "grow \n"<< std::flush;
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

            assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations
            nonLinearSolver.solve(x); // solve the non-linear system with time step control
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
                problem->userData("axialFlux", x);
                problem->userData("radialFlux", x);
                problem->userData("age", x); // age changes with time
                problem->userData("kr", x);  // conductivities change with age
                problem->userData("kx", x);
                vtkWriter.write(timeLoop->time());
            }
            problem->writeTranspirationRate(x); // always add transpiration data into the text file
            timeLoop->reportTimeStep();  // report statistics of this time step

            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize())); // set new dt as suggested by the newton solver

        } while (!timeLoop->finished());

        timeLoop->finalize(leafGridView.comm());

    } else { // static

        std::cout << "a static model \n\n" << std::flush;

        assembler->setPreviousSolution(xOld); // set previous solution for storage evaluations

        nonLinearSolver.solve(x); // solve the non-linear system

        // write outputs
        problem->userData("axialFlux", x);
        problem->userData("radialFlux", x);
        problem->userData("age", x); // prepare fields
        problem->userData("kr", x);  // conductivities change with age
        problem->userData("kx", x);
        vtkWriter.write(1); // write vtk output
        problem->writeTranspirationRate(x);
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

