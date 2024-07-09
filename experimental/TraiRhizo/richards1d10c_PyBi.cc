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

#include "../../cpp/python_binding/py_richards10c_cyl.hh"

std::vector<std::array<double, 1>> setShape(
	//Richards10Cyl<RichardsCylFoamProblem, RichardsCylFoamAssembler, RichardsCylFoamLinearSolver, 1> s, 
	int nVertices, double r_in, double r_out) 
{
    //int nVertices = 10;
    double logbase = Dumux::getParam<double>("Grid.logbase", 0.5);// 0.5;
	bool doLogarithm = Dumux::getParam<bool>("Grid.doLogarithm", true);

	std::vector<std::array<double, 1>> points(nVertices);
	
	
    double start = std::log(r_in) / std::log(logbase);
    double end = std::log(r_out) / std::log(logbase);
    double step = (end - start) / (nVertices - 1);
	
	for (int i = 0; i < nVertices; ++i) {
        double point = std::pow(logbase, start + i * step);
		
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
	
	double maxDt = getParam<double>("Problem.maxDt",250);
    s.initializeProblem(maxDt);
	
	std::cout<<"initializeProblem created "<<std::endl;
		 
	double simTime = getParam<double>("Problem.simTime",1200); //s
	s.ddt= getParam<double>("Problem.ddt",1); //s
    s.solveNoMPI(simTime);
	std::cout<<"s.solve finished"<<std::endl;
	for (int cs = 1; cs < s.numComp(); cs ++)
	{
		std::vector<double> C_S_fr = s.getSolution(cs);
		
		double mincfr = *min_element(C_S_fr.begin(), C_S_fr.end());
		std::cout << cs << " " << mincfr << std::endl;
	}
	
    return 0;
}
