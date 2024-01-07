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

#include "../python_binding/py_richards10c_cyl.hh"
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI


std::vector<std::vector<double>> getFluxOrSource_10c(
	Richards10Cyl<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1> s, 
	bool doGetSource) 
{
        
        std::vector<std::vector<std::vector<double>>> inSources;
        std::vector<double> vols; // m3 scv or [-]
		if(doGetSource)
		{
			inSources = s.inSources;
			vols = s.getCellVolumes(); // m3 scv,  att for MPI
		}else{
			inSources = s.inFluxes;
			vols = std::vector<double> (s.getCellVolumes().size(),1.); // dummy,  att for MPI
		}
        std::vector<double> inFluxes_ddt = s.inFluxes_ddt;

		int numberOfCellsTot = vols.size();
		
        try {
            assert(inSources.size() == inFluxes_ddt.size());
            assert(inSources[0].size() == numberOfCellsTot);
            assert(inSources[0][0].size() == s.numComp());
        } catch (std::exception& e) {
            std::cerr << "Assertion failed: " << e.what() << std::endl;
            throw e;
        }

        std::vector<std::vector<double>> inSources_tot(numberOfCellsTot, std::vector<double>(s.numComp() , 0.0));

        for (size_t isrcs = 0; isrcs < inFluxes_ddt.size(); ++isrcs) {
            for (size_t idc = 0; idc < numberOfCellsTot; ++idc) {
                for (size_t comp = 0; comp < s.numComp(); ++comp) {
                    inSources_tot[idc][comp] += inSources[isrcs][idc][comp] * vols[idc] * inFluxes_ddt[isrcs];
                }
            }
        }


        return inSources_tot;
}
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

int main(int argc, char** argv) //try
{
    using namespace Dumux;
	// parse command line arguments and input file
	auto s = Richards10Cyl<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1>();
    
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
	std::map<int, double> sourceMap;// mol/s
	double flux = getParam<double>("Soil.sourceCs");
	sourceMap[0] = flux; //2.2486750867416232e-07; mol/m3/s
	s.setSource( sourceMap, 1);
	std::cout<<"initializeProblem created "<<s.numComp()<<std::endl;
	double f_sorp = getParam<double>("Soil.f_sorp");//[-]
		 
    double molarDensityWat_m3 =  s.molarDensityWat_m3;//densityWat_m3 / molarMassWat;
	std::vector<double> C_S_fr = s.getSolution(1);//mol/mol
	std::vector<double> C_S_W(C_S_fr.size());//mol/m3 wat
	std::vector<double> CSS1_init = s.computeCSS1s();//(C_S_fr.size());//mol/m3 scv
	std::vector<double> C_S(C_S_fr.size());//mol
	std::vector<double> C_tot_init(C_S_fr.size());//mol
	std::vector<double> C_SS1(C_S_fr.size());//mol
	std::vector<double> watCont = s.getWaterContent(); // m3 wat/m3 scv
	std::vector<double> WatVol(C_S_fr.size());//m3 wat
	std::vector<double> cellVol = s.getCellVolumes();//m3 scv
	double C_all_init = 0.;
	
	for(int cs = 0; cs < C_S_fr.size(); cs ++)
	{
		WatVol.at(cs) = cellVol.at(cs) * watCont.at(cs);
		C_S_W.at(cs) = C_S_fr.at(cs) * molarDensityWat_m3;
		C_S.at(cs) = C_S_W.at(cs) * WatVol.at(cs);
		// CSS1_init.at(cs) = CSSmax*(C_S_W.at(cs) /(C_S_W.at(cs) +k_sorp))*f_sorp;
		C_SS1.at(cs) = f_sorp * CSS1_init.at(cs) * cellVol.at(cs);
		C_tot_init.at(cs) = C_SS1.at(cs) + C_S.at(cs);
		std::cout<<"cellId "<<cs<<", C_S "<<C_S.at(cs) <<" CSS1 "<<C_SS1.at(cs) 
		<<" tot "<< C_tot_init.at(cs) <<std::endl;
		C_all_init +=  C_tot_init.at(cs)  ;
	}
	std::cout<<"Cs all cells "<< C_all_init << std::endl;
    s.solveNoMPI(1200, -1, false, true, true);
	std::cout<<"s.solve finished"<<std::endl;
	std::vector<std::vector<double> > flux10c = getFluxOrSource_10c(s,false);
    std::vector<std::vector<double> > bulkSoil_sources = getFluxOrSource_10c(s, true);
	C_S_fr = s.getSolution(1);
	std::vector<double> C_S_tot_end(C_S_fr.size());//mol
	watCont = s.getWaterContent();//m3 wat
	std::vector<double> C_tot_end(C_S_fr.size());//mol
	double C_all_end = 0.;
	double Error_all_end = 0.;
	std::vector<double> CSS1_end = s.computeCSS1s();
	std::vector<double> errorMass(C_S_fr.size());//mol
	for(int cs = 0; cs < C_S_fr.size(); cs ++)
	{
		WatVol.at(cs) = cellVol.at(cs) * watCont.at(cs);
		C_S_W.at(cs) = C_S_fr.at(cs) * molarDensityWat_m3;
		C_S.at(cs) = C_S_W.at(cs) * WatVol.at(cs);
		
		C_SS1.at(cs) = f_sorp * CSS1_end.at(cs) * cellVol.at(cs);
		C_tot_end.at(cs) = C_SS1.at(cs) + C_S.at(cs);
		std::cout<<"cellId "<<cs<<", C_S "<<C_S.at(cs) <<" CSS1 "<<C_SS1.at(cs) 
		<<" tot "<< C_tot_end.at(cs) <<std::endl;
		std::cout<<"source "<<bulkSoil_sources.at(cs)[1]<<" flux cs "<<flux10c.at(cs)[1]<<std::endl;
		C_all_end +=  C_tot_end.at(cs) ;
		errorMass.at(cs) = C_tot_end.at(cs) 
		-( C_tot_init.at(cs)  + bulkSoil_sources.at(cs)[1] - flux10c.at(cs)[1] );
		std::cout<<"error "<<errorMass.at(cs) << std::endl;
		Error_all_end += errorMass.at(cs);
		
	}
	std::cout<<"Cs all cells "<< C_all_end << std::endl;
	std::cout<<"error all cells "<< Error_all_end << std::endl;
    return 0;
}
