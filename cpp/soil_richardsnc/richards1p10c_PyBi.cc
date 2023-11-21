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
// #include <config.h>

// #include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
// #include <dune/common/timer.hh> // to compute wall times
// #include <dune/grid/io/file/dgfparser/dgfexception.hh>
// #include <dune/grid/io/file/vtk.hh>
// #include <dune/grid/yaspgrid/coordinates.hh>
// #include <dune/istl/io.hh>
// // #include <dumux/common/properties.hh> // creates an undefined TypeTag types, and includes the property system
// // #include <dumux/common/properties/propertysystem.hh>
// #include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file
 #include <dumux/common/valgrind.hh> // for debugging
// #include <dumux/common/dumuxmessage.hh> // for fun (a static class)
// #include <dumux/common/defaultusagemessage.hh> // for information (a function)

// #include <dumux/linear/amgbackend.hh> // linear solver (currently the only parallel solver available(?))
// #include <dumux/porousmediumflow/richards/newtonsolver.hh>
// /**
 // * Some small adaption to <dumux/nonlinear/newtonsolver.hh>, which is the only nonlinear solver available.
 // * The adaption is disabled per default (parameter EnableChop = false)
 // */
// #include <dumux/common/parameters.hh>
// #include <dumux/common/timeloop.hh>
// #include <dumux/assembly/fvassembler.hh> // assembles residual and Jacobian of the nonlinear system

// #include <dumux/io/vtkoutputmodule.hh>
// #include <dumux/io/grid/gridmanager.hh>
// // #include <dumux/io/loadsolution.hh> // functions to resume a simulation

// #include "richards1p10cproblem_cyl.hh" // the problem class. Defines some TypeTag types and includes its spatialparams.hh class
// #include "properties_10c.hh" // the property system related stuff (to pass types, used instead of polymorphism)
// #include "properties_nocoupling.hh" // dummy types for replacing the coupling types

// // getDofIndices, getPointIndices, getCellIndices
// #include <dune/grid/utility/globalindexset.hh>



#include "py_richards.hh"
/**
 * here we go
 */
int main(int argc, char** argv) //try
{
    //using namespace Dumux;
    // parse command line arguments and input file
	auto s = Richards<RichardsSPProblem, RichardsSPAssembler, RichardsSPLinearSolver, 3>();
	std::vector<std::string> args_{""};
	s.initialize(args_);
	std::string icz = "-100";
    s.setParameter("Soil.IC.P", icz) ; // cm pressure head
	s.setParameter("Soil.BC.Top.Type", "2");
	s.setParameter("Soil.BC.Top.Value", "0.");
	s.setParameter("Soil.BC.Bot.Type", "2");
	s.setParameter("Soil.BC.Bot.Value", "0.");
    
	std::cout<<"to create grid"<<std::endl;
	//double l = 0.53/100.;
	//int N = 41;
	std::array<double, 3> boundsMin{-0.53, -0.53, -1./100.};//std::vector<double> 
	std::array<double, 3> boundsMax{0.53, 0.53, 0.};// std::vector<double> 
	std::array<int, 3> numberOfCells{41, 41, 1};
    s.createGrid(boundsMin,boundsMax, numberOfCells)  ;// [m]
    
	std::cout<<"grid created"<<std::endl;
	
	s.setParameter( "Soil.VanGenuchten.Qr", "0.08");
	s.setParameter( "Soil.VanGenuchten.Qs", "0.43");
	s.setParameter( "Soil.VanGenuchten.Alpha", "0.04");
	s.setParameter("Soil.VanGenuchten.N", "1.6");
	s.setParameter("Soil.VanGenuchten.Ks", "50");

    
	std::cout<<"to initializeProblem"<<std::endl;
    s.initializeProblem();
    
	std::cout<<"initializeProblem created"<<std::endl;
    //s.setCriticalPressure(-15000);
	
    return 0;
}
