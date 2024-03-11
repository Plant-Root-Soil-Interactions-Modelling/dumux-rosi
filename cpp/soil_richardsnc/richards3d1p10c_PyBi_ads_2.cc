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


std::vector<std::vector<double>> getFluxOrSource_10c(
	Richards10<Richards10CSPProblem, Richards10CSPAssembler, Richards10CSPLinearSolver, 3> s, 
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
            assert(inSources[0][0].size() == s.numComp());//+ 1
        } catch (std::exception& e) {
            std::cerr << "Assertion failed: " << e.what() << std::endl;
            throw e;
        }

        std::vector<std::vector<double>> inSources_tot(numberOfCellsTot, std::vector<double>(s.numComp() , 0.0));//+ 2

        for (size_t isrcs = 0; isrcs < inFluxes_ddt.size(); ++isrcs) {
            for (size_t idc = 0; idc < numberOfCellsTot; ++idc) {
                for (size_t comp = 0; comp < s.numComp() ; ++comp) {//+ 1
                    inSources_tot[idc][comp] += inSources[isrcs][idc][comp] * vols[idc] * inFluxes_ddt[isrcs];
                }
            }
        }

        // Sum along the first dimension
        // for (size_t idc = 0; idc < numberOfCellsTot; ++idc) {
            // inSources_tot[idc][s.numComp() ] = 0.0;
            // for (size_t comp = 0; comp < s.numComp() + 1; ++comp) {
                // inSources_tot[idc][s.numComp() + 1] += inSources_tot[idc][comp];
            // }
        // }

        //std::vector<std::vector<double>> d_css1(numberOfCellsTot, std::vector<double>(1, 0.0));

        

        //for (size_t idc = 0; idc < numberOfCellsTot; ++idc) {
        //    inSources_tot[idc].insert(inSources_tot[idc].end(), d_css1[idc].begin(), d_css1[idc].end());
        //}

        return inSources_tot;
}

double computeCSS1(double C_S_W, double theta, double svc_volume) 
	{// [ mol / m^3]
		
    using namespace Dumux;
	double m3_2_cm3 = 1e6;
	double k_sorp =  getParam<double>("Soil.k_sorp")* m3_2_cm3;// mol/m3 water
	double CSSmax =  getParam<double>("Soil.CSSmax")* m3_2_cm3;//mol/m3
	int css1Function = getParam<double>("Soil.css1Function");
		switch(css1Function) {
		  case 0:
			return CSSmax*(C_S_W/(C_S_W+k_sorp));
		  case 1:
			return 0.;//none
		  case 2:
		  // [mol C/m3 scv zone 1] = [m3 soil water/m3 scv zone 1] * [mol C/m3 soil water]
			return CSSmax*C_S_W/k_sorp;//linear, with CSSmax in [m3 wat/m3 scv zone 1]
		  case 3:
			return 0.;//only pde
		  case 4:
			return CSSmax*(C_S_W/(C_S_W+k_sorp));
		  case 5:
			return CSSmax*(svc_volume*C_S_W*theta/(svc_volume*C_S_W*theta+k_sorp));
		  case 6://CSSmax is content
			return CSSmax*(svc_volume*C_S_W*theta/(svc_volume*C_S_W*theta+k_sorp))/svc_volume;
		  case 7://linear with depends on content
			return CSSmax*(svc_volume*C_S_W*theta/k_sorp);
		  case 8://linear with CSSmax is content
			return CSSmax*C_S_W*theta/k_sorp;
		  default:
			DUNE_THROW(Dune::InvalidStateException, "css1Function not recognised (0, 1, or 2)"+ std::to_string(css1Function));
		}
		return 1.;
	}

int main(int argc, char** argv) //try
{
    using namespace Dumux;
	// parse command line arguments and input file
	auto s = Richards10<Richards10CSPProblem, Richards10CSPAssembler, Richards10CSPLinearSolver, 3>();
    
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
    
	double bulkSoilDensity = (2650. / 60.08e-3) * (1 - getParam<int>("Soil.VanGenuchten.Qs")) ;
    double CS_init_mFr =  getParam<double>("Soil.IC.C1");//mol/mol
	double CS_init = CS_init_mFr * bulkSoilDensity;//mol/m3
    double CSS2Init = s.computeInitCSS2(0., CS_init);// mol/m3
	
	std::array<double, 3> boundsMin{-0.05, -0.05, -0.1}; 
	std::array<double, 3> boundsMax{0.05, 0.05, 0.}; 
	int nCellX = getParam<int>("Soil.nCellX",2);
	int nCellY = getParam<int>("Soil.nCellY",2);
	int nCellZ = getParam<int>("Soil.nCellZ",2);
	std::array<int, 3> numberOfCells{nCellX,nCellY,nCellZ};//{1,1,1};
    s.createGrid(boundsMin,boundsMax, numberOfCells)  ;// [m]
	// Soil.Grid.Cells
	std::cout<<"grid created"<<std::endl;
	std::cout<<"to initializeProblem"<<std::endl;
    double f_sorp = 0;
    
    s.initializeProblem();
	
    //double molarMassWat = 18.; // [g/mol]
    //double densityWat_m3 = 1e6 ;//[g/m3]
    //[mol/m3] = [g/m3] /  [g/mol] 
    double molarDensityWat_m3 =  s.molarDensityWat_m3;//densityWat_m3 / molarMassWat;
	std::vector<double> C_S_fr = s.getSolution(1);//mol/mol
	std::vector<double> C_SS2_fr = s.getSolution(7);//mol/mol
    
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
    s.solve(1200, -1, false, true);
	std::cout<<"s.solve finished"<<std::endl;
	std::vector<std::vector<double> > flux10c = getFluxOrSource_10c(s,false);
    std::vector<std::vector<double> > bulkSoil_sources = getFluxOrSource_10c(s, true);
	C_S_fr = s.getSolution(1);
	C_SS2_fr = s.getSolution(7);//mol/mol
    
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
		double flux_css1 = 0;
		if(false)//getFluxcss1)
		{
			flux_css1 = f_sorp *computeCSS1(flux10c.at(cs)[1], 1., cellVol.at(cs));
		}
		std::cout<<"source "<<bulkSoil_sources.at(cs)[1]<<" flux cs "<<flux10c.at(cs)[1]<<" flux css1 "<<flux_css1<<std::endl;
		C_all_end +=  C_tot_end.at(cs) ;
		errorMass.at(cs) = C_tot_end.at(cs) 
		-( C_tot_init.at(cs)  + bulkSoil_sources.at(cs)[1] - flux10c.at(cs)[1] + flux_css1);
		std::cout<<"error "<<errorMass.at(cs) << std::endl;
		Error_all_end += errorMass.at(cs);
		
	}
	std::cout<<"Cs all cells "<< C_all_end << std::endl;
	std::cout<<"error all cells "<< Error_all_end << std::endl;
    return 0;
}
