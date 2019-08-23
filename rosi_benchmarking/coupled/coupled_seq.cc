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
#include <dune/grid/common/rangegenerators.hh> // elements, vertices, scv, scvfs, ...
#include <dune/istl/io.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/periodic/tpfa/periodicnetworkgridmanager.hh>
#include <dumux/periodic/tpfa/fvgridgeometry.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh>

#include <RootSystem.h>

#include <dumux/common/timeloop.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh> // solver
#include <dumux/porousmediumflow/richards/newtonsolver.hh> // solver
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>

#include <dumux/growth/soillookup.hh> // for coupling

#include "../rootsystem/rootsproblem.hh"
#include "../soil/richardsproblem.hh"
#include "../rootsystem/properties.hh" // TypeTag:Roots
#include "../soil/properties.hh" // TypeTag:RichardsTT

namespace Dumux { namespace Properties {


using RootsTag = Properties::TTag::RootsBox;
using SoilTag = Properties::TTag::RichardsBox;

using SoilGridType = Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>; // pick soil grid here (its in compile definition in the soil model)
using RootGridType = GetPropType<RootsTag, Properties::Grid>;

using RootFVGridGeometry = GetPropType<RootsTag, Properties::FVGridGeometry>;
using SoilFVGridGeometry = GetPropType<SoilTag, Properties::FVGridGeometry>;

using SoilLookUp = GrowthModule::SoilLookUpBBoxTree<SoilFVGridGeometry>;



/*!
 * calculates saturation from the solution vector (todo which dof is associated to scv?)
 * todo move to soil
 */
template<class SoilGridVariables, class SoilSolution>
void updateSaturation(std::vector<double>& saturation, const SoilFVGridGeometry& gridGeometry, const SoilGridVariables& gridVariables,
    const SoilSolution& sol) {

    const auto& gridView = gridGeometry.gridView();  // soil
    for (const auto& element : elements(gridView)) { // soil elements

        // soil solution -> volume variable
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (const auto& scv : scvs(fvGeometry)) {
            saturation[scv.dofIndex()] = elemVolVars[scv].saturation(0); // saturation(0);
        }
    }
}

/*!
 * calculates the radial fluxes of the segments
 * todo: move to RootProlbem...
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

            double s = rootProblem.source(element, gridGeometry, elemVolVars, scv); // pass to problem class

            // define the type tag for this problem
            auto pos = scv.center();
            int eIdx = soilLookUp->pick(pos); // find element index of soil, for root each root element
            if (eIdx>=0) {
                //std::cout << "eIdx " << eIdx << ", " << s <<"\n"<<std::flush;
                source.at(eIdx) += s; // accumulate source term
            } else {
                std::cout << "root at position " << pos << " not within soil";
            }

        }

    }

}

} // end namespace Dumux





int main(int argc, char** argv) try
{
    using namespace Dumux;
    int simtype = Properties::simtype;

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
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<RootGridType> rootGridManager;
    rootGridManager.init("RootSystem");
    const auto rootGridData = rootGridManager.getGridData();

    GridManager<SoilGridType> soilGridManager;
    soilGridManager.init("Soil");

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // compute on the leaf grid views
    const RootGridType::LeafGridView& rootLGV = rootGridManager.grid().leafGridView();
    const SoilGridType::LeafGridView& soilLGV = soilGridManager.grid().leafGridView();
    std::cout << "i have two views \n" << std::flush;

    // create the finite volume grid geometry
    auto rootFVGridGeometry = std::make_shared<RootFVGridGeometry>(rootLGV);
    rootFVGridGeometry->update();

    auto soilFVGridGeometry = std::make_shared<SoilFVGridGeometry>(soilLGV);
    soilFVGridGeometry->update();
    std::cout << "i have two geometries built from the views \n" << std::flush;

    // the problem (initial and boundary conditions)
    using RootProblem = GetPropType<RootsTag, Properties::Problem>;
    auto rootProblem = std::make_shared<RootProblem>(rootFVGridGeometry);
    rootProblem->spatialParams().initParameters(*rootGridData);

    using SoilProblem = GetPropType<SoilTag, Properties::Problem>;
    auto soilProblem = std::make_shared<SoilProblem>(soilFVGridGeometry);
    std::cout << "and i have two problems \n" << "\n" << std::flush;

    // the solution vector
    using RootSolutionVector = GetPropType<RootsTag, Properties::SolutionVector>;
    RootSolutionVector r(rootFVGridGeometry->numDofs());
    rootProblem->applyInitialSolution(r);
    auto rOld = r;

    using SoilSolutionVector = GetPropType<SoilTag, Properties::SolutionVector>;
    SoilSolutionVector s(soilFVGridGeometry->numDofs());
    soilProblem->applyInitialSolution(s);
    auto sOld = s;
    std::cout << "no solution yet \n" << std::flush;

    // the grid variables
    using RootGridVariables = GetPropType<RootsTag, Properties::GridVariables>;
    auto rootGridVariables = std::make_shared<RootGridVariables>(rootProblem, rootFVGridGeometry);
    rootGridVariables->init(r);

    using SoilGridVariables = GetPropType<SoilTag, Properties::GridVariables>;
    auto soilGridVariables = std::make_shared<SoilGridVariables>(soilProblem, soilFVGridGeometry);
    soilGridVariables->init(s);
    std::cout << "... but grid variables \n" << std::flush;

    /**
     * COUPLING
     */
    static constexpr int soilDim = SoilFVGridGeometry::GridView::dimension;
    std::vector<double> saturation(soilLGV.size(soilDim), 1.0);
    updateSaturation(saturation, *soilFVGridGeometry, *soilGridVariables, s);
    auto soilLookUp = SoilLookUp(*soilFVGridGeometry, saturation);

    rootProblem->setSoil(&soilLookUp);
    std::cout << "roots know the soil \n" << "\n" << std::flush;

    std::vector<double> soilSink = std::vector<double>(soilFVGridGeometry->gridView().size(0));
    radialFlux2soilSink(soilSink, *rootFVGridGeometry, *rootGridVariables, r, *rootProblem, &soilLookUp); // precomputes the sink for the soil problem
    soilProblem->setSource(&soilSink);
    std::cout << "and the soil knows the roots \n" << "\n" << std::flush;

    std::cout << "value = " << soilLookUp.getValue(CRootBox::Vector3d());

    std::string sss = "";
    std::getline(std::cin, sss);



    // get some time loop parameters & instantiate time loop
    using Scalar = GetPropType<RootsTag, Properties::Scalar>;
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
    using RootIOFields = GetPropType<RootsTag, Properties::IOFields>;
    VtkOutputModule<RootGridVariables, RootSolutionVector> rootVTKWriter(*rootGridVariables, r, rootProblem->name()+"R");
    using RootVelocityOutput = GetPropType<RootsTag, Properties::VelocityOutput>;
//    rootVTKWriter.addVelocityOutput(std::make_shared<RootVelocityOutput>(*rootGridVariables));
//    rootProblem->axialFlux(r); // prepare fields
//    rootProblem->radialFlux(r); // prepare fields
//    rootVTKWriter.addField(rootProblem->axialFlux(), "axial flux");
//    rootVTKWriter.addField(rootProblem->radialFlux(), "radial flux");
    RootIOFields::initOutputModule(rootVTKWriter); //!< Add model specific output fields
    rootVTKWriter.write(0.0);

    using SoilIOFields = GetPropType<SoilTag, Properties::IOFields>;
    VtkOutputModule<SoilGridVariables, SoilSolutionVector> soilVTKWriter(*soilGridVariables, s, soilProblem->name()+"S");
    using SoilVelocityOutput = GetPropType<SoilTag, Properties::VelocityOutput>;
//    soilVTKWriter.addVelocityOutput(std::make_shared<SoilVelocityOutput>(*soilGridVariables));
    SoilIOFields::initOutputModule(soilVTKWriter); //!< Add model specific output fields
    soilVTKWriter.write(0.0);
    std::cout << "vtk writer module initialized (in less than 20 lines)" << "\n" << std::flush;

    // the assembler with time loop for instationary problem
    using RootAssembler = FVAssembler<RootsTag, DiffMethod::numeric>;
    std::shared_ptr<RootAssembler> rootAssembler;
    if (tEnd > 0) {
        rootAssembler = std::make_shared<RootAssembler>(rootProblem, rootFVGridGeometry, rootGridVariables, timeLoop); // dynamic
    } else {
        rootAssembler = std::make_shared<RootAssembler>(rootProblem, rootFVGridGeometry, rootGridVariables); // static
    }
    using SoilAssembler = FVAssembler<SoilTag, DiffMethod::numeric>;
    std::shared_ptr<SoilAssembler> soilAssembler;
    if (tEnd>0) {
        soilAssembler = std::make_shared<SoilAssembler>(soilProblem, soilFVGridGeometry, soilGridVariables, timeLoop); // dynamic
    } else {
        soilAssembler = std::make_shared<SoilAssembler>(soilProblem, soilFVGridGeometry, soilGridVariables); // static
    }

    // the linear solver
    using RootLinearSolver = AMGBackend<RootsTag>;
    auto rootLinearSolver = std::make_shared<RootLinearSolver>(rootLGV, rootFVGridGeometry->dofMapper());
    using SoilLinearSolver = Dumux::AMGBackend<SoilTag>; // the only parallel linear solver available
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
            soilAssembler->setPreviousSolution(sOld);
            rootAssembler->setPreviousSolution(rOld);

            // solves the soil problem
            radialFlux2soilSink(soilSink, *rootFVGridGeometry, *rootGridVariables, r, *rootProblem, &soilLookUp); // precomputes the sink for the soil problem
            soilNonlinearSolver.solve(s, *timeLoop);

            // solves the root problem
            updateSaturation(saturation, *soilFVGridGeometry, *soilGridVariables, s); // updates soil look up for the root problem
            rootNonlinearSolver.solve(r, *timeLoop);

            // make the new solution the old solution
            rOld = r;
            sOld = s;

            rootGridVariables->advanceTimeStep();
            soilGridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();
            // write vtk output (only at check points)
            if ((timeLoop->isCheckPoint()) || (timeLoop->finished())) {
//                rootProblem->axialFlux(r); // prepare fields
//                rootProblem->radialFlux(r); // prepare fields
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

            // pass current time to the problem
            double t = timeLoop->time();
            double dt = timeLoop->timeStepSize();
            rootProblem->setTime(t,dt);
            soilProblem->setTime(t);

        } while (!timeLoop->finished());

        timeLoop->finalize(rootLGV.comm());

    } else { // static

        std::cout << "a static model" << "\n" << std::flush;
        // set previous solution for storage evaluations
        rootAssembler->setPreviousSolution(rOld);
        soilAssembler->setPreviousSolution(sOld);
        // solve the non-linear system
        rootNonlinearSolver.solve(r);
        soilNonlinearSolver.solve(s);
        // write vtk output
//        rootProblem->axialFlux(r); // prepare fields
//        rootProblem->radialFlux(r); // prepare fields
        rootProblem->writeTranspirationRate(r);
        rootVTKWriter.write(1);
        soilVTKWriter.write(1);

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

