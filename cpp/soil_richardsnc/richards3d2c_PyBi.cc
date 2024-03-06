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

#include "../python_binding/py_richardsnc.hh"


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
    // initialize (e.g. multi-threading backend)
    initialize(argc, argv);
	// parse command line arguments and input file
	auto s = Richards<RichardsNCSPProblem, 
						RichardsNCSPAssembler, 
						RichardsNCSPLinearSolver, 3>();
    
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
	
	std::cout<<"create grid"<<std::endl;
	
	bool periodic = false;
	
	
	// std::string icz = "-1000";
    // s.setParameter("Soil.IC.P", icz) ; // cm pressure head
	// s.setParameter("Soil.BC.Top.Type", "3");
	// s.setParameter("Soil.BC.Top.Value", "1.");
	// s.setParameter("Soil.BC.Bot.Type", "3");
	// s.setParameter("Soil.BC.Bot.Value", "0.");
	
	int const dim = 3;
	std::array<double, dim> min_b= {-4./100, -4./100, -25./100};
	std::array<double, dim> max_b={4./100, 4./100, 0./100};
	std::array<int, dim> cell_number={8, 8, 25} ;
    s.createGrid(min_b, max_b, cell_number, periodic);// [m]

	s.setParameter("Soil.Grid.Cells", 
			std::to_string(cell_number[0]+cell_number[1]+cell_number[2]));
	std::cout<<"grid created"<<std::endl;
	std::cout<<"to initializeProblem"<<std::endl;
    s.initializeProblem();
	// std::map<int, double> sourceMap;// mol/s
	// double flux = getParam<double>("Soil.sourceCs");
	// sourceMap[0] = flux; //2.2486750867416232e-07; mol/m3/s
	// s.setSource( sourceMap, 1);
	std::cout<<"initializeProblem created "<<std::endl;
		 
	std::vector<double> C_S_fr = s.getSolution(1);//kg/kg
	std::vector<double> C_S_W(C_S_fr.size());//kg/m3 wat
	std::vector<double> C_S(C_S_fr.size());//kg
	std::vector<double> C_tot_init(C_S_fr.size());//kg
	std::vector<double> watCont = s.getWaterContent(); // m3 wat/m3 scv
	std::vector<double> WatVol(C_S_fr.size());//m3 wat
	std::vector<double> cellVol = s.getCellVolumes();//m3 scv
	double C_all_init = 0.;
	
	for(int cs = 0; cs < C_S_fr.size(); cs ++)
	{
		WatVol.at(cs) = cellVol.at(cs) * watCont.at(cs);// m3 wat
		C_S_W.at(cs) = C_S_fr.at(cs) * rho_;// kg solute/m3 wat
		C_S.at(cs) = C_S_W.at(cs) * WatVol.at(cs);// kg solute
		C_all_init +=  C_tot_init.at(cs)  ;// kg solute
	}
	std::cout<<"Cs all cells "<< C_all_init << std::endl;
	double simTime = getParam<double>("Problem.simTime",1200); //s
    s.solve(simTime);
	std::cout<<"s.solve finished"<<std::endl;
	C_S_fr = s.getSolution(1);
	std::vector<double> C_S_tot_end(C_S_fr.size());//mol
	watCont = s.getWaterContent();//m3 wat
	std::vector<double> C_tot_end(C_S_fr.size());//mol
	double C_all_end = 0.;
	double Error_all_end = 0.;
	std::vector<double> errorMass(C_S_fr.size());//mol
	for(int cs = 0; cs < C_S_fr.size(); cs ++)
	{
		WatVol.at(cs) = cellVol.at(cs) * watCont.at(cs);
		C_S_W.at(cs) = C_S_fr.at(cs) * rho_;
		C_S.at(cs) = C_S_W.at(cs) * WatVol.at(cs);
		
		C_all_end +=  C_tot_end.at(cs) ;
		// errorMass.at(cs) = C_tot_end.at(cs) 
		// -( C_tot_init.at(cs)  + bulkSoil_sources.at(cs)[1] - flux10c.at(cs)[1] );
		// std::cout<<"error "<<errorMass.at(cs) << std::endl;
		// Error_all_end += errorMass.at(cs);
		
	}
	std::cout<<"C_all_end "<<C_all_end<<std::endl;
	// std::cout<<"Cs all cells "<< C_all_end << std::endl;
	// std::cout<<"error all cells "<< Error_all_end << std::endl;
    return 0;
}
