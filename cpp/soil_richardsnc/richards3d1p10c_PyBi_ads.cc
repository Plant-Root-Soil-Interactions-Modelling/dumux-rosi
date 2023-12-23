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
	
	std::cout<<"to create grid"<<std::endl;
	
	std::array<double, 3> boundsMin{-0.05, -0.05, -0.1}; 
	std::array<double, 3> boundsMax{0.05, 0.05, 0.}; 
	int nCellX = getParam<int>("Soil.nCellX",1);
	int nCellY = getParam<int>("Soil.nCellY",1);
	int nCellZ = getParam<int>("Soil.nCellZ",1);
	std::array<int, 3> numberOfCells{nCellX,nCellY,nCellZ};//{1,1,1};
    s.createGrid(boundsMin,boundsMax, numberOfCells)  ;// [m]
	
	std::cout<<"grid created"<<std::endl;
	std::cout<<"to initializeProblem"<<std::endl;
    s.initializeProblem();
	std::map<int, double> sourceMap;// mol/s
	double flux = getParam<double>("Soil.sourceCs");
	sourceMap[0] = flux; //2.2486750867416232e-07;
	s.setSource( sourceMap, 1);
	std::cout<<"initializeProblem created"<<std::endl;
	double m3_2_cm3 = 1e6;
	//double f_sorp = getParam<double>("Soil.f_sorp");//[-]
	//double k_sorp = getParam<double>("Soil.k_sorp")* m3_2_cm3;// mol/m3 water
	//double CSSmax = getParam<double>("Soil.CSSmax")* m3_2_cm3;//mol/m3
		 
    double molarMassWat = 18.; // [g/mol]
    double densityWat_m3 = 1e6 ;//[g/m3]
    //[mol/m3] = [g/m3] /  [g/mol] 
    double molarDensityWat_m3 =  densityWat_m3 / molarMassWat;
	std::vector<double> C_S_fr = s.getSolution(1);//mol/mol
	std::vector<double> C_S_W(C_S_fr.size());//mol/m3 wat
	std::vector<double> CSS1_init = s.computeCSS1s();//(C_S_fr.size());//mol/m3 scv
	std::vector<double> C_S_tot(C_S_fr.size());//mol
	std::vector<double> C_SS1_tot(C_S_fr.size());//mol
	std::vector<double> watCont = s.getWaterContent(); // m3 wat/m3 scv
	std::vector<double> WatVol(C_S_fr.size());//m3 wat
	std::vector<double> cellVol = s.getCellVolumes();//m3 scv
	double C_tot_init = 0.;
	for(int cs = 0; cs < C_S_fr.size(); cs ++)
	{
		WatVol.at(cs) = cellVol.at(cs) * watCont.at(cs);
		C_S_W.at(cs) = C_S_fr.at(cs) * molarDensityWat_m3;
		C_S_tot.at(cs) = C_S_W.at(cs) * WatVol.at(cs);
		// CSS1_init.at(cs) = CSSmax*(C_S_W.at(cs) /(C_S_W.at(cs) +k_sorp))*f_sorp;
		C_SS1_tot.at(cs) = CSS1_init.at(cs) * cellVol.at(cs);
		std::cout<<"C_S_tot "<<C_S_tot.at(cs) <<" CSS1_tot "<<C_SS1_tot.at(cs) 
		<<" tot "<< C_S_tot.at(cs) + C_SS1_tot.at(cs) <<std::endl;
		C_tot_init = C_tot_init + C_S_tot.at(cs) + C_SS1_tot.at(cs) ;
	}
    s.solve(1200);
	std::cout<<"s.solve finished"<<std::endl;
	
	C_S_fr = s.getSolution(1);
	watCont = s.getWaterContent();//m3 wat
	double C_tot_end = 0.;
	std::vector<double> CSS1_end = s.computeCSS1s();
	for(int cs = 0; cs < C_S_fr.size(); cs ++)
	{
		WatVol.at(cs) = cellVol.at(cs) * watCont.at(cs);
		C_S_W.at(cs) = C_S_fr.at(cs) * molarDensityWat_m3;
		C_S_tot.at(cs) = C_S_W.at(cs) * WatVol.at(cs);
		// CSS1_end.at(cs) = CSSmax*(C_S_W.at(cs) /(C_S_W.at(cs) +k_sorp))*f_sorp;
		C_SS1_tot.at(cs) = CSS1_end.at(cs) * cellVol.at(cs);
		std::cout<<"C_S_tot "<<C_S_tot.at(cs) <<" CSS1_tot "<<C_SS1_tot.at(cs) 
		<<" tot "<< C_S_tot.at(cs) + C_SS1_tot.at(cs) <<std::endl;
		C_tot_end = C_tot_end + C_S_tot.at(cs) + C_SS1_tot.at(cs) ;
	}
	std::cout<<"balance: init "<<C_tot_init<<" C_tot_end "<<C_tot_end
	<<" input "<<flux*1200.<<" diff "<<C_tot_end - (C_tot_init + flux*1200.)<<std::endl;
    return 0;
}
