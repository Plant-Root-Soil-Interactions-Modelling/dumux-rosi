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

#include "../../cpp/python_binding/py_richards10c.hh"
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file


static constexpr double g_ = 9.81; // cm / s^2 (for type conversions)
static constexpr double rho_ = 1.e3; // kg / m^3 (for type conversions)
static constexpr double pRef_ = 1.e5; // Pa

//! cm pressure head -> Pascal
static double toPa(double ph) {
	return pRef_ + ph / 100. * rho_ * g_;
}
//! Pascal -> cm pressure head
static double toHead(double p) {
	return (p - pRef_) * 100. / rho_ / g_;
}
int main(int argc, char** argv) //try
{
    using namespace Dumux;
	// parse command line arguments and input file
	auto s = Richards10<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPLinearSolverNum, 3>();
    
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
	
	std::cout<<"to create grid"<<std::endl;
    
	double bulkSoilDensity = (2650. / 60.08e-3) * (1 - getParam<double>("Soil.VanGenuchten.Qs")) ;
    
	std::array<double, 3> boundsMin{-0.015, -0.06, -0.4}; 
	std::array<double, 3> boundsMax{0.015, 0.06, 0}; 
	int nCellX = getParam<int>("Soil.nCellX",2);
	int nCellY = getParam<int>("Soil.nCellY",2);
	int nCellZ = getParam<int>("Soil.nCellZ",2);
	std::array<int, 3> numberOfCells{nCellX,nCellY,nCellZ};//{1,1,1};
    s.createGrid(boundsMin,boundsMax, numberOfCells)  ;// [m]
	// Soil.Grid.Cells
	std::cout<<"grid created"<<std::endl;
	std::cout<<"to initializeProblem"<<std::endl;
	
	double maxDt = getParam<double>("Problem.maxDt",250);
    s.initializeProblem(maxDt);
	
	for(int cs = 0; cs < C_S_fr.size(); cs ++)
	{
		WatVol.at(cs) = cellVol.at(cs) * watCont.at(cs);// m3 wat
		C_S_W.at(cs) = C_S_fr.at(cs) * rho_;// kg solute/m3 wat
		C_S.at(cs) = C_S_W.at(cs) * WatVol.at(cs);// kg solute
		C_all_init +=  C_tot_init.at(cs)  ;// kg solute
	}
	
	std::vector<double> watCont = s.getWaterContent(); // m3 wat/m3 scv
	std::cout<<"watCont "<<watCont.at(0)<<std::endl;
	
	std::vector<double> phead = s.getSolutionHead(); // m3 wat/m3 scv
	std::cout<<"phead "<<phead.at(0)<<std::endl;
	double dt = getParam<double>("Problem.dt",60);
	s.solve(dt);
    return 0;
}
