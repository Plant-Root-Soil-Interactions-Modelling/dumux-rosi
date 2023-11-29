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

#include "../python_binding/py_richards10c.hh"
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI


int main(int argc, char** argv) //try
{
    using namespace Dumux;
    // parse command line arguments and input file
	auto s = Richards<Richards10CSPProblem, Richards10CSPAssembler, Richards10CSPLinearSolver, 3>();
    
    // already done in s.initialize
    //const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    for(int i = 0; i < argc; ++i)
	{
        std::cout << argv[i] << '\n';
	}
	std::vector<std::string> args_{""};
	s.initialize(args_);
	std::string icz = "-1000";
    s.setParameter("Soil.IC.P", icz) ; // cm pressure head
	s.setParameter("Soil.BC.Top.Type", "3");
	s.setParameter("Soil.BC.Top.Value", "1.");
	s.setParameter("Soil.BC.Bot.Type", "3");
	s.setParameter("Soil.BC.Bot.Value", "0.");
    
    
	//s.setParameter("Problem.verboseIndexSet", "0");
	std::cout<<"to create grid"<<std::endl;
	
	std::array<double, 3> boundsMin{-0.05, -0.05, -0.1};//std::vector<double> 
	std::array<double, 3> boundsMax{0.05, 0.05, 0.};// std::vector<double> 
	std::array<int, 3> numberOfCells{5, 5, 20};
	s.setParameter("Soil.Grid.Cells", "5 5 20");
    s.createGrid(boundsMin,boundsMax, numberOfCells)  ;// [m]
    
	std::cout<<"grid created"<<std::endl;
	
	s.setParameter( "Soil.VanGenuchten.Qr", "0.08");
	s.setParameter( "Soil.VanGenuchten.Qs", "0.43");
	s.setParameter( "Soil.VanGenuchten.Alpha", "0.04");
	s.setParameter("Soil.VanGenuchten.N", "1.6");
	s.setParameter("Soil.VanGenuchten.Ks", "50");
	
	s.setParameter("Problem.EnableGravity", "true");
	s.setParameter("Problem.reactionExclusive", "1");
	s.setParameter("Soil.MolarMass", "0.06008");
	s.setParameter("Soil.solidDensity", "2700");
    
    
	std::cout<<"to initializeProblem"<<std::endl;
    s.initializeProblem();
    
	std::cout<<"initializeProblem created"<<std::endl;
    s.solve(60.);
	std::cout<<"s.solve finished"<<std::endl;
	
    return 0;
}
