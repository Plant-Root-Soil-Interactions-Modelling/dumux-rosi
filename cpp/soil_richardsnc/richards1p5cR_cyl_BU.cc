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
 * \brief Richards equation realized with Richards box model.
 */
#include <config.h>

#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
#include <dune/common/timer.hh> // to compute wall times
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid/coordinates.hh>
#include <dune/istl/io.hh>
// #include <dumux/common/properties.hh> // creates an undefined TypeTag types, and includes the property system
// #include <dumux/common/properties/propertysystem.hh>
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file
#include <dumux/common/valgrind.hh> // for debugging
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/defaultusagemessage.hh> // for information (a function)

#include <dumux/linear/amgbackend.hh> // linear solver (currently the only parallel solver available(?))
#include <dumux/porousmediumflow/richards/newtonsolver.hh>
/**
 * Some small adaption to <dumux/nonlinear/newtonsolver.hh>, which is the only nonlinear solver available.
 * The adaption is disabled per default (parameter EnableChop = false)
 */
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/assembly/fvassembler.hh> // assembles residual and Jacobian of the nonlinear system

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
// #include <dumux/io/loadsolution.hh> // functions to resume a simulation

#include "richards1p5cproblemReaction_cyl.hh" // the problem class. Defines some TypeTag types and includes its spatialparams.hh class
#include "properties_cyl_5c.hh" // the property system related stuff (to pass types, used instead of polymorphism)
#include "properties_nocoupling.hh" // dummy types for replacing the coupling types

/**
 * here we go
 */
int main(int argc, char** argv) //try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::RichardsNCBox; // Richards2CCC, Richards2CBox, (TypeTag is defined in the problem class richardsproblem.hh)

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)

    // print dumux start message
    if (mpiHelper.rank() == 0) { // rank is the process number
        DumuxMessage::print(/*firstCall=*/true);
    }

    // parse command line arguments and input file
    Parameters::init(argc, argv);
	// const auto& myParams = Parameters::paramTree();
	//myParams.report();
	//const auto& paramKeys = myParams.getValueKeys();
	// for (int i = 1; i < paramKeys.size(); ++i)
    // {
		// std::cout<<paramKeys[i]<<" "<<myParams[paramKeys[i]]<<std::endl;
	// }
    /*
     * Parses the parameters from the .input file into a static data member of the static class Parameters,
     * the parameter tree can then be accessed with Parameters::paramTree() .
     * Furthermore, global auxiliary functions are defined e.g. getParam, setParam, and haveParam
     * (in parameters.hh).
     *
     * All Dumux classes access the global parameter tree where considered appropriate.
     */

    // try to create a grid (from the given grid file or the input file)
    using SoilGridType = GetPropType<TypeTag, Properties::Grid>;
    /**
     * Properties::Grid is defined in the problem class richardsproblem.hh to Dune::YaspGrid<3>,
     * or if available to the compile definition GRIDTYPE that is given in CMakeLists, and is
     * YaspGrid<3> (for richards3d), Dune::FoamGrid<1,1> (for richards1d),
     * or Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming> (for richardsUG)
     */
    GridManager<SoilGridType> gridManager;
    gridManager.init("Soil"); // "Soil" is the parameter group name
    /**
     * Opens the grid file, or constructs the grid form the .input file parameters
     */


    /////////////////////////////////////////////////////////////////
    // run steady state or dynamic non-linear problem on this grid
    /////////////////////////////////////////////////////////////////


    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    /**
     * we work only on the leafGridView and do not need the grid ever again,
     * leafGridView.grid() returns a const reference to the grid (in case we need data attached to the grid).
     *
     * Of dune type: GridFamily::Traits::LeafGridView
     *
     * Have not found where the type is set (grid dependent),
     * but probably implements DefaultLeafGridView (in dune/grid/common/defaultgridview.hh)
     */

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    /**
     * The type is dependent on the discretization (e.g. box, tpfa, ...),
     *
     * For Box method:
     * Properties::FVGridGeometry is defined in problem.hh -> discretization/box.hh
     * The type is BoxFVGridGeometry (in discretization/box/fvgridgeometry.hh)
     * specialization of BaseFVGridGeometry (discretization/basefvgridgeometry.hh)
     */
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    /**
     * holds the vertexMapper() and elementMapper().
     *
     * ??? copies the leafGridView (but, how does it know when its updated?)
     * it seems i am only allowed to use fvGridGeometry->gridView() in the following
     * rendering leafGridView defined above, pointless
     */
    fvGridGeometry->update(); // update all Mappers(do this again after grid adaption)

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
	auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
	using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>; // defined in discretization/fvproperties.hh, as Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>
    SolutionVector x(fvGridGeometry->numDofs()); // degrees of freedoms
	problem->applyInitialSolution(x); // Dumux way of saying x = problem->applyInitialSolution()
    auto xOld = x;
	std::cout<<"shape of solution matrix "<<x.size()<<" "<<x[0].size()<<std::endl;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    /**
     * The type is defined GridVariables in discretization/fvproperties.hh, as
     * FVGridVariables<FVGridGeometry, GVV, GFVC>
     * where GVV = grid volume variables, and GFVC = grid flux variables cache. I have not found where these are set,
     * but I assume for box method types are BoxGridFluxVariablesCache.
     *
     * FVGridVariables is defined in fvgridvariables.hh.
     *
     * Manages the grid volume variables, and the flux variable cache, stores a pointer to problem (why?)
     */
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x); // initialize all variables , updates volume variables to the current solution, and updates the flux variable cache

    // get some time loop parameters  & instantiate time loop
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
	std::string extenction = getParam<std::string>("Problem.extenction", "csv");
	
    std::shared_ptr<CheckPointTimeLoop<Scalar>> timeLoop; // defined in common/timeloop.hh, everything clear, easy to use.
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
        } catch(std::exception& e) {
            std::cout<< "richards1p5c.cc: no check times (TimeLoop.CheckTimes) defined in the input file\n";
        }
    } else { // static
    }

    // instantiate time loop & intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);//somehow throws error
    /**
     * home grown vkt output, rather uninteresting
     */

    // the assembler with time loop for instationary or stationary problem (assembles resdiual, and Jacobian for Newton)
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>; //  FVAssembler (in fvassembler.hh)
    std::shared_ptr<Assembler> assembler;
    if (tEnd>0) {
        assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop); // dynamic
    } else {
        assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables); // static
    }

    // the linear solver
	//std::cout<<"the linear solver"<<std::endl;
    using LinearSolver = Dumux::AMGBackend<TypeTag>; // the only linear solver available, files located in located in dumux/linear/*
    auto linearSolver = std::make_shared<LinearSolver>(fvGridGeometry->gridView(), fvGridGeometry->dofMapper());

    // the non-linear solver
	//std::cout<<"the NON-linear solver"<<std::endl;
    using NonLinearSolver = Dumux::RichardsNewtonSolver<Assembler, LinearSolver>;
    NonLinearSolver nonLinearSolver = NonLinearSolver(assembler, linearSolver);

    // std::cin.ignore();  // wait for key (debugging)
	//std::cout<<"tEnd>0? "<<tEnd<<std::endl;
    // print dumux end message
	double pRefPa = 101300; 
	double rho__ = 1.e3;
	double g__ = 9.80665;
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
		std::ofstream myfile_;
		std::string filestr = problem->name() + "_1p3cR_cyl_T0." + extenction; // output file
		myfile_.open(filestr.c_str());
		for (const auto& vertex : Dune::vertices(fvGridGeometry->gridView()))
		{
			const auto dofIdxGlobal = fvGridGeometry->vertexMapper().index(vertex);
			auto rVal = vertex.geometry().center()[0];
			myfile_ <<rVal<<", ";
			
			myfile_ <<x[dofIdxGlobal][0]<<", "<<(( x[dofIdxGlobal][0] - pRefPa)*100/rho__/g__) << ", ";
				for(int j = 1; j < 3; j++)
				{
					myfile_ << x[dofIdxGlobal][j]/18.<<", ";
				} 
				
				for(int j = 3; j < x[dofIdxGlobal].size(); j++)
				{
					myfile_ << x[dofIdxGlobal][j]*25615.8<<", ";
				} myfile_ << "\n";
			
		}
		myfile_.close();
        DumuxMessage::print(/*firstCall=*/false);
    }
	
	double bulkSoilDensity = (2700. / 60.08e-3) * (1 - 0.43) ;
	
			std::cout<<"\n\nwe get: "<< timeLoop->time()<<" "
			<<timeLoop->timeStepSize()   
			<<"\n\t"<<x[0][3]*bulkSoilDensity
			<<" "<< problem->getDecay() <<" "
			<<x[0][3]*bulkSoilDensity + problem->getDecay() * timeLoop->timeStepSize() 
			<<"\n\t"<<x[0][1]*problem->getTheta()* (1e6/18.) //mol /mol* g/m3 / g/mol
			<<" "<< problem->getDecayCs() <<" "
			<<x[0][1]*problem->getTheta() * (1e6/18.)+ problem->getDecayCs() * timeLoop->timeStepSize() 
			<<std::endl<<std::endl;
	
    if (tEnd>0)  { // dynamic
        timeLoop->start();

        do {
            // set previous solution for storage evaluations
			//std::cout<<"set previous solution for storage evaluations "<<xOld.size()<<" "<<xOld[0].size()<<std::endl;
            assembler->setPreviousSolution(xOld);
            // solve the non-linear system with time step control
			//std::cout<<"solve the non-linear system with time step control"<<std::endl;
            nonLinearSolver.solve(x, *timeLoop);
            // make the new solution the old solution
			//std::cout<<"make the new solution the old solution"<<std::endl;
            xOld = x;
			//std::cout<<"advance time step"<<std::endl;
            gridVariables->advanceTimeStep();
            // advance to the time loop to the next step
			//std::cout<<"advance to the time loop to the next step"<<std::endl;
            timeLoop->advanceTimeStep();
            // write vtk output (only at check points)
            //if ((timeLoop->isCheckPoint()) || (timeLoop->finished())) {
            //    vtkWriter.write(timeLoop->time());
            //}
            // report statistics of this time step
			std::cout<<"\n\nwe get: "<< timeLoop->time()<<" "
			<<timeLoop->timeStepSize()   
			<<"\n\t"<<x[0][3]*bulkSoilDensity
			<<" "<< problem->getDecay() <<" "
			<<x[0][3]*bulkSoilDensity + problem->getDecay() * timeLoop->timeStepSize() 
			<<"\n\t"<<x[0][1]*problem->getTheta()* (1e6/18.) //mol /mol* g/m3 / g/mol
			<<" "<< problem->getDecayCs() <<" "
			<<x[0][1]*problem->getTheta() * (1e6/18.)+ problem->getDecayCs() * timeLoop->timeStepSize() 
			<<std::endl<<std::endl;
			std::cout<<"report statistics of this time step"<<std::endl;
            timeLoop->reportTimeStep();
            // set new dt as suggested by the newton solver
			std::cout<<"set new dt as suggested by the newton solve"<<std::endl;
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
            // pass current time to the problem
			std::cout<<"pass current time to the problem"<<std::endl;
            problem->setTime( timeLoop->time() , timeLoop->timeStepSize()  );
            problem->postTimeStep(x, *gridVariables);
			//DUNE_THROW(Dune::InvalidStateException, "problem->postTimeStep(x, *gridVariables);");
            problem->writeBoundaryFluxes();

        } while (!timeLoop->finished());

        timeLoop->finalize(fvGridGeometry->gridView().comm());

    } else { // static
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);
        // solve the non-linear system
        nonLinearSolver.solve(x);
        vtkWriter.write(1);
        problem->postTimeStep(x, *gridVariables);
        problem->writeBoundaryFluxes();
    }

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
		//[kg soil / m3 soil] / [kg soil / mol soil] * [m3space /m3 space - m3 pores / m3 space] * [m3/cm3]
		// == [mol soil/ m3 soil] * [m3 soil/m3 space] * [m3/cm3] = mol soil /cm3 space
		double bulkSoilDensity = (2700. / 60.08e-3) * (1 - 0.43) /1e6; 
        Parameters::print();
		std::ofstream myfile_;
		std::string filestr = problem->name() + "_1p5cR_cyl_end."+extenction; // output file
		// myfile_.open(filestr.c_str());
		// //auto allnames = std::vector<std::string>{"P", "Cs", "Cl", "r"}
		// for(int i = 0; i < x.size(); i++)
		// {
			// //myfile_<<allnames[i]<<",\n";
			// myfile_ <<x[i][0]<<", "<<(( x[i][0] - pRefPa)*100/rho__/g__) << ", ";
			// for(int j = 1; j < 3; j++)
			// {
				// myfile_ << "("<<x[i][j] << ", "<<x[i][j]/18.<<"); ";//mol/cm3 water
			// } for(int j = 3; j < x[i].size(); j++)
			// {
				// myfile_ << "("<<x[i][j] << ", "<<x[i][j]*25615.8<<"); ";
			// }myfile_ << "\n";
		// }
		// myfile_.close();
		
		// filestr = problem->name() + "_1p5cR_cyl_end2."+extenction; // output file
		myfile_.open(filestr.c_str());
		
		for (const auto& vertex : Dune::vertices(fvGridGeometry->gridView()))
		{
			const auto dofIdxGlobal = fvGridGeometry->vertexMapper().index(vertex);
			auto rVal = vertex.geometry().center()[0];
			myfile_ <<rVal<<", ";
			
			myfile_ <<x[dofIdxGlobal][0]<<", "<<(( x[dofIdxGlobal][0] - pRefPa)*100/rho__/g__) << ", ";
				for(int j = 1; j < 3; j++)
				{
					myfile_ << x[dofIdxGlobal][j] * (1/18.) <<", ";//[mol / mol] * (g/cm3) / (g/mol) = [mol solute / cm3 water]
				} 
				
				for(int j = 3; j < x[dofIdxGlobal].size(); j++)
				{
					myfile_ << x[dofIdxGlobal][j]*bulkSoilDensity<<", "; // [mol C / mol soil] * [mol / m3 space] 
				} myfile_ << "\n";
			
		}
		
		myfile_.close();

        DumuxMessage::print(/*firstCall=*/false);
    }
	
    return 0;
}