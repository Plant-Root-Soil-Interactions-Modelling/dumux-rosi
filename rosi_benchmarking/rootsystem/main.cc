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
 * \brief Test for the 1d-3d embedded mixed-dimension model coupling two
 *        one-phase porous medium flow problems
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <memory>
#include <algorithm>

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

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

// CRootBox
#include <RootSystem.h>

// Dumux Growth
#include <dumux/growth/gridgrowth.hh>
#include <dumux/growth/crootboxalgorithm.hh>
#include <dumux/growth/rootsystemgridfactory.hh>
#include <dumux/growth/soillookup.hh>
#include <dumux/growth/localresidual.hh>
#include <dumux/growth/problem.hh>
#include <dumux/growth/couplingmanager.hh>

// problem specific includes
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// the problems to solve
#include "rootproblem_.hh"
#include "soilproblem.hh"

namespace Dumux {
namespace Properties {

////////////////////////////////////////////////////////////
// SubTypeTags /////////////////////////////////////////////
////////////////////////////////////////////////////////////
NEW_TYPE_TAG(RootTypeTag, INHERITS_FROM(CCTpfaModel, OneP));
NEW_TYPE_TAG(SoilTypeTag, INHERITS_FROM(CCTpfaModel, Richards));

////////////////////////////////////////////////////////////
// Root Properties /////////////////////////////////////////
////////////////////////////////////////////////////////////
// Set the grid type
SET_TYPE_PROP(RootTypeTag, Grid, Dune::FoamGrid<1, 3>);

SET_BOOL_PROP(RootTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(RootTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(RootTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(RootTypeTag, SolutionDependentAdvection, false);
SET_BOOL_PROP(RootTypeTag, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(RootTypeTag, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(RootTypeTag, Problem, GrowthModule::GrowthProblemAdapter<RootProblem<TypeTag>>);

// the fluid system
SET_PROP(RootTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
};

// Set the problem property
SET_TYPE_PROP(RootTypeTag, LocalResidual, GrowthModule::GrowthLocalResidualAdapter<OnePIncompressibleLocalResidual<TypeTag>>);

// Set the spatial parameters
SET_TYPE_PROP(RootTypeTag, SpatialParams, RootSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                                            typename GET_PROP_TYPE(TypeTag, Scalar)>);

// coupling properties
SET_PROP(RootTypeTag, CouplingManager){
using Traits = MultiDomainTraits<TTAG(SoilTypeTag), TypeTag>;
using type = SoilRootCouplingManager<Traits>;
};
// the point source type
SET_TYPE_PROP(RootTypeTag, PointSource, typename GET_PROP_TYPE(TypeTag, CouplingManager)::PointSourceTraits::template PointSource<1>);
// the point source locater helper class
SET_TYPE_PROP(RootTypeTag, PointSourceHelper, typename GET_PROP_TYPE(TypeTag, CouplingManager)::PointSourceTraits::template PointSourceHelper<1>);

////////////////////////////////////////////////////////////
// Soil Properties /////////////////////////////////////////
////////////////////////////////////////////////////////////
// Set the grid type
SET_TYPE_PROP(SoilTypeTag, Grid, Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);

SET_BOOL_PROP(SoilTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(SoilTypeTag, EnableGridVolumeVariablesCache, true);SET_BOOL_PROP(SoilTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(SoilTypeTag, SolutionDependentAdvection, false);SET_BOOL_PROP(SoilTypeTag, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(SoilTypeTag, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(SoilTypeTag, Problem, SoilProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(SoilTypeTag, SpatialParams, SoilSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry), typename GET_PROP_TYPE(TypeTag, Scalar)>);

// coupling properties
SET_PROP(SoilTypeTag, CouplingManager){
using Traits = MultiDomainTraits<TypeTag, TTAG(RootTypeTag)>;
using type = SoilRootCouplingManager<Traits>;
};
// the point source type
SET_TYPE_PROP(SoilTypeTag, PointSource, typename GET_PROP_TYPE(TypeTag, CouplingManager)::PointSourceTraits::template PointSource<0>);
// the point source locater helper class
SET_TYPE_PROP(SoilTypeTag, PointSourceHelper, typename GET_PROP_TYPE(TypeTag, CouplingManager)::PointSourceTraits::template PointSourceHelper<0>);

} // end namespace Properties

template<class SoilFVGridGeometry, class SoilGridVariables, class SoilSolution>
void updateSaturation(std::vector<double>& saturation, const SoilFVGridGeometry& gridGeometry, const SoilGridVariables& gridVariables,
    const SoilSolution& sol) {
    const auto& gridView = gridGeometry.gridView();
    for (const auto& element : elements(gridView)) {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);

        auto elemVolVars = localView(gridVariables.curGridVolVars());
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (const auto& scv : scvs(fvGeometry))
            saturation[scv.dofIndex()] = elemVolVars[scv].saturation(0);
    }
}

} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;
    using namespace GrowthModule;

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // start the simulation timer
    Dune::Timer timer;

    // Define the sub problem type tags
    using SoilTypeTag = TTAG(SoilTypeTag);
    using RootTypeTag = TTAG(RootTypeTag);

    // create the soil grid
    GridManager<typename GET_PROP_TYPE(SoilTypeTag, Grid)> soilGridManager;
    soilGridManager.init("Soil"); // pass parameter group (see input file)

    // create a croot box rootsystem
    auto rootSystem = std::make_shared<CRootBox::RootSystem>();

    // read in the plant and root parameter files and initialize root system
    rootSystem->openFile(getParam<std::string>("Parameters.File"), "params/");

    // soil grid lookup
    const auto& soilGridView = soilGridManager.grid().leafGridView();
    using SoilFVGridGeometry = typename GET_PROP_TYPE(SoilTypeTag, FVGridGeometry);
    auto soilGridGeoemtry = std::make_shared<SoilFVGridGeometry>(soilGridView);
    soilGridGeoemtry->update();

    // set the soil lookup
    static constexpr int soilDim = SoilFVGridGeometry::GridView::dimension;
    std::vector<double> saturation(soilGridView.size(soilDim), 1.0); // initialize with 1.0
    // using SoilGridType = typename GET_PROP_TYPE(SoilTypeTag, Grid);
    // auto soilLookup = SoilLookUpBBoxTree<SoilGridType>(soilGridView, soilGridGeoemtry->boundingBoxTree(), saturation);
    // rootSystem->setSoil(&soilLookup);

    // set the domain as confining boundary
    const auto size = soilGridGeoemtry->bBoxMax() - soilGridGeoemtry->bBoxMin();
    auto box = std::make_unique<CRootBox::SDF_PlantBox>(size[0] * 100, size[1] * 100, size[2] * 100); // cm (!)
    rootSystem->setGeometry(box.get());

    // intialize the root system
    rootSystem->initialize();

    // simulate one day so we have a small structure (our initial structure for the fluid simulation)
    const double initialGrowthTime = getParam<double>("RootSystem.InitialGrowthTime", 1.0);
    rootSystem->simulate(initialGrowthTime);

    // write initial rootsystem as vtp
    std::ofstream initialRootSystem; initialRootSystem.open("initial.vtp");
    rootSystem->writeVTP(initialRootSystem);
    initialRootSystem.close();

    // create the dumux grid (using dune-foamgrid) from the initial crootbox grid
    auto rootGrid = RootSystemGridFactory::makeGrid(*rootSystem);
    using RootFVGridGeometry = typename GET_PROP_TYPE(RootTypeTag, FVGridGeometry);
    auto rootGridGeometry = std::make_shared<RootFVGridGeometry>(rootGrid->leafGridView());
    rootGridGeometry->update();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<SoilTypeTag, RootTypeTag>;
    constexpr auto soilDomainIdx = Traits::template DomainIdx<0>();
    constexpr auto rootDomainIdx = Traits::template DomainIdx<1>();

    // the coupling manager
    using CouplingManager = typename GET_PROP_TYPE(SoilTypeTag, CouplingManager);
    auto couplingManager = std::make_shared<CouplingManager>(soilGridGeoemtry, rootGridGeometry);

    // the problem (initial and boundary conditions)
    using SoilProblem = typename GET_PROP_TYPE(SoilTypeTag, Problem);
    auto soilProblem = std::make_shared<SoilProblem>(soilGridGeoemtry, couplingManager);
    using RootProblem = typename GET_PROP_TYPE(RootTypeTag, Problem);
    const auto domainSize = soilGridGeoemtry->bBoxMax()-soilGridGeoemtry->bBoxMin();
    auto rootProblem = std::make_shared<RootProblem>(rootGridGeometry, couplingManager, domainSize);
    // obtain parameters from the crootbox root system
    rootProblem->spatialParams().initParameters(*rootSystem);
    rootProblem->spatialParams().analyseRootSystem(*rootSystem);
    rootProblem->updateRootVolume();

    // the solution vector
    Traits::SolutionVector sol;
    sol[soilDomainIdx].resize(soilGridGeoemtry->numDofs());
    sol[rootDomainIdx].resize(rootGridGeometry->numDofs());
    soilProblem->applyInitialSolution(sol[soilDomainIdx]);
    rootProblem->applyInitialSolution(sol[rootDomainIdx]);
    auto oldSol = sol;

    couplingManager->init(soilProblem, rootProblem, sol);
    soilProblem->computePointSourceMap();
    rootProblem->computePointSourceMap();

    // the grid variables
    using SoilGridVariables = typename GET_PROP_TYPE(SoilTypeTag, GridVariables);
    auto soilGridVariables = std::make_shared<SoilGridVariables>(soilProblem, soilGridGeoemtry);
    soilGridVariables->init(sol[soilDomainIdx], oldSol[soilDomainIdx]);
    using RootGridVariables = typename GET_PROP_TYPE(RootTypeTag, GridVariables);
    auto rootGridVariables = std::make_shared<RootGridVariables>(rootProblem, rootGridGeometry);
    rootGridVariables->init(sol[rootDomainIdx], oldSol[rootDomainIdx]);

    // update the saturation vector
    // updateSaturation(saturation, *soilGridGeoemtry, *soilGridVariables, sol[soilDomainIdx]);

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // intialize the vtk output module
    using SoilSolution = std::decay_t<decltype(sol[soilDomainIdx])>;
    VtkOutputModule<SoilGridVariables, SoilSolution> soilVtkWriter(*soilGridVariables, sol[soilDomainIdx], soilProblem->name());
    GET_PROP_TYPE(SoilTypeTag, VtkOutputFields)::init(soilVtkWriter);
    soilVtkWriter.write(0.0);

    using RootSolution = std::decay_t<decltype(sol[rootDomainIdx])>;
    VtkOutputModule<RootGridVariables, RootSolution> rootVtkWriter(*rootGridVariables, sol[rootDomainIdx], rootProblem->name());
    GET_PROP_TYPE(RootTypeTag, VtkOutputFields)::init(rootVtkWriter);
    rootVtkWriter.addField(rootProblem->spatialParams().radii(), "radius"); //! Add radius field
    rootVtkWriter.addField(rootProblem->spatialParams().orders(), "order"); //! Add root order field
    rootVtkWriter.addField(rootProblem->spatialParams().ages(), "age"); //! Add age field
    rootVtkWriter.addField(rootProblem->spatialParams().radialConductivities(), "Kr"); //! Add Kr field
    rootVtkWriter.addField(rootProblem->spatialParams().axialConductivities(), "Kx"); //! Add Kx field
    rootVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(soilProblem, rootProblem),
                                                 std::make_tuple(soilGridGeoemtry, rootGridGeometry),
                                                 std::make_tuple(soilGridVariables, rootGridVariables),
                                                 couplingManager, timeLoop);

    // the linear solver
    using LinearSolver = BlockDiagILU0BiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // the growth algorithm
    auto growthAlgorithm = std::make_shared<CRootBoxGrowthAlgorithm<RootTypeTag>>
          (rootGrid, rootGridGeometry, rootSystem, *rootProblem, sol[rootDomainIdx]);
    GridGrowth growth(growthAlgorithm);
    const bool enableGrowth = getParam<bool>("RootSystem.EnableGrowth", false);

    // time loop
    // timeLoop->setPeriodicCheckPoint(episodeLength);
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // compute the segment lengths at the beginning of the time step
        for (const auto& element : elements(rootGrid->leafGridView()))
        {
            const auto eIdx = rootGridGeometry->elementMapper().index(element);
            rootProblem->spatialParams().setPreviousLength(eIdx, element.geometry().volume());
        }

        // set previous solution for storage evaluations
        assembler->setPreviousSolution(oldSol);

        // set time for time dependent boundary conditions
        soilProblem->setTime(timeLoop->time() + timeLoop->timeStepSize());
        rootProblem->setTime(timeLoop->time() + timeLoop->timeStepSize());

        // solve the non-linear system
        // TODO: implement time step control
        // for this it has to be possible to reset the crootbox rootsystem to the last step
        nonLinearSolver.solve(sol);

        // make the new solution the old solution
        oldSol = sol;
        soilGridVariables->advanceTimeStep();
        rootGridVariables->advanceTimeStep();

        // update the saturation vector
        // updateSaturation(saturation, *soilGridGeoemtry, *soilGridVariables, sol[soilDomainIdx]);

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // output the source terms
        soilProblem->computeSourceIntegral(sol[soilDomainIdx], *soilGridVariables);
        rootProblem->computeSourceIntegral(sol[rootDomainIdx], *rootGridVariables);
        rootProblem->plotTranspirationRate(timeLoop->timeStepSize(), sol[rootDomainIdx], *rootGridVariables);

        // write vtk output
        soilVtkWriter.write(timeLoop->time());
        rootVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }

    timeLoop->finalize();

    // report used and unused parameters
    Parameters::print();

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
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (std::exception &e)
{
    std::cerr << "Standard library exception thrown: " << e.what() << " ---> Abort!" << std::endl;
    return 4;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 5;
}
