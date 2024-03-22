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



int main(int argc, char** argv) //try
{
    using namespace Dumux;
	
	//typedef FakeMPIHelper MPIHelper;
	// parse command line arguments and input file
	auto s = Richards10<RichardsNCSPProblem, RichardsSPAssemblerNum, RichardsSPSSORCGIstlLinearSolver>();
    
    
    // parse command line arguments and input file
    Parameters::init(argc, argv);
    for(int i = 0; i < argc; ++i)
	{
        std::cout << argv[i] << '\n';
	}
	std::vector<std::string> args_{""};
	s.initialize(args_);
	
	std::cout<<"to create grid"<<std::endl;
    
	double MolarMass = getParam<double>("Soil.MolarMass");
	double solidDensity = getParam<double>("Soil.solidDensity");
	double css1Function = getParam<double>("Soil.css1Function");
	
	double maxDt = getParam<double>("Problem.maxDt",  250);//s
	double simTime =getParam<int>("Problem.simTime",1200);//s
	double wilting_point =getParam<int>("Problem.wilting_point",-15000);//Pa?
	bool doMPIsolve =getParam<int>("Problem.doMPIsolve",true);//
	bool saveInnerDumuxValues =getParam<int>("Problem.saveInnerDumuxValues",false);//
	
    s.createGrid()  ;// [m]
	
	std::cout<<"grid created"<<std::endl;
	std::cout<<"to initializeProblem"<<std::endl;
    
    
    s.initializeProblem(maxDt);
	
    s.setCriticalPressure(wilting_point)  ;
	
    s.solve(simTime, doMPIsolve, saveInnerDumuxValues);
	
	
	std::vector<double> C_S_fr = s.getSolution(1);
	std::vector<double> watCont = s.getWaterContent();//m3 wat/m3 scv
	std::vector<double> cellVol = s.getCellVolumes();//m3 scv
    double molarDensityWat_m3 =  s.molarDensityWat_m3;
	std::vector<double> WatVol(C_S_fr.size());//m3 wat
	std::vector<double> C_S_W(C_S_fr.size());//mol/m3 wat
	std::vector<double> C_S(C_S_fr.size());//mol
	double C_all_end = 0.;
	
	
	for(int cs = 0; cs < C_S_fr.size(); cs ++)
	{
		WatVol.at(cs) = cellVol.at(cs) * watCont.at(cs);
		C_S_W.at(cs) = C_S_fr.at(cs) * molarDensityWat_m3;
		C_S.at(cs) = C_S_W.at(cs) * WatVol.at(cs);
		
		C_all_end +=  C_S.at(cs) ;
		
	}
	std::cout<<"Cs all cells "<< C_all_end << std::endl;
    return 0;
}
