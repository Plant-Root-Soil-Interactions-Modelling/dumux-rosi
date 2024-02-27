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

#include "../python_binding/py_richardsnc_cyl.hh"

std::vector<std::array<double, 1>> setShape(
	//Richards10Cyl<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1> s, 
	int nVertices, double r_in, double r_out) 
{
    //int nVertices = 10;
    double logbase = Dumux::getParam<double>("Grid.logbase", 0.5);// 0.5;
	bool doLogarithm = Dumux::getParam<bool>("Grid.doLogarithm", true);

	std::vector<std::array<double, 1>> points(nVertices);
	for (int i = 0; i < nVertices; ++i) {
		double point;
		if(doLogarithm)
		{
			point = std::pow(logbase, (std::log(r_in) / std::log(logbase)) + i * ((std::log(r_out) / std::log(logbase)) - (std::log(r_in) / std::log(logbase))) / nVertices);
		}else{
			point = r_in + (r_out - r_in) * (double(i)/double(nVertices - 1));
		}
		points.at(i).at(0) = (point);
	}
	return points;

}

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
	auto s = RichardsCyl<RichardsNCCylFoamProblem, 
						RichardsNCCylFoamAssembler, 
						RichardsNCCylFoamLinearSolver, 1>();
    
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
	
    double l = getParam<double>("Problem.segLength",1.); // length in m
	s.setParameter("Problem.segLength", std::to_string(l));// m
	double r_in = getParam<double>("Grid.LowerLeft", 0.02/100.);// m
	double r_out = getParam<double>("Grid.UpperRight", 0.2104870824716315/100.);// m
	int nVertices = getParam<int>("Grid.Cells", 9) + 1;
	int nCells = getParam<int>("Grid.Cells", 9);
	auto points = setShape(nVertices, r_in, r_out);
	s.createGrid1d(points);

    s.setParameter("Soil.Grid.Cells", std::to_string(nCells));
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
	double simTime = getParam<double>("Problem.segLength",1200); //s
    s.solveNoMPI(simTime);
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
	std::array<double, 1> innerCell;
	innerCell[0] = r_in;
	int idx = s.pickCell(innerCell);
	double f_inW = s.getNeumann(idx, 0);
	double f_inC = s.getNeumann(idx, 1);
	std::cout<<"C_all_end "<<C_all_end<<"f_inW_end "<<f_inW<<" f_inC "<<f_inC<<std::endl;
	// std::cout<<"Cs all cells "<< C_all_end << std::endl;
	// std::cout<<"error all cells "<< Error_all_end << std::endl;
    return 0;
}
