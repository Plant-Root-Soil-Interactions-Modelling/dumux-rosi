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
 * \brief test for the rootsystem coupled model
 */
#include "config.h"
//#include "rosiRichards2cbuffertestproblem.hh"
#include "2soilRichards2cbuffertestproblem.hh"
//#include <dumux/multidimension/embeddedcoupling/start.hh>
#include <dumux/common/start.hh>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    // TODO
}

int main(int argc, char** argv)
{
	//for (int i=0; i<argc; ++i)
	//	std::cout<<argv[i]<<std::endl;
	//std::cout << "Press Enter to Continue" << std::to_string(1.234);
	//std::cin.ignore();

#if HAVE_UMFPACK
    typedef TTAG(SoilRichardsTwoCBufferTestCCProblem) ProblemTypeTag;
    //typedef TTAG(SoilRichardsTwoCBufferTestBoxProblem) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv, usage);
#else
#warning "You need to have the UMFPack solver library installed to run this test!"
    std::cout << "This test needs the UMFPack solver library to run!" << std::endl;
#endif
}
