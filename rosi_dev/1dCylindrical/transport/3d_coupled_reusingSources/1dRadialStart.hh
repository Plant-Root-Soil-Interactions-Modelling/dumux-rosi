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
 * \brief Provides a few default main functions for convenience.
 */
#ifndef DUMUX_RADIAL_START_HH
#define DUMUX_RADIAL_START_HH

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/parameterparser.hh>

namespace Dumux
{
// forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridCreator);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(TimeManager);
}

/*!
 * \ingroup Start
 * \ingroup EmbeddedCoupling
 * \brief Read the command line arguments and write them into the parameter tree.
 *        Do some syntax checks.
 *
 * \param   argc      The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
 * \param   argv      The 'argv' argument of the main function: array of pointers to the argument strings
 * \param   paramTree The parameterTree. It can be filled from an input file or the command line.
 * \return            Empty string if everything worked out. Otherwise the thing that could not be read.
 */
std::string radialreadOptions_(int argc, char **argv, Dune::ParameterTree &paramTree)
{
    // All command line options need to start with '-'
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            std::ostringstream oss;
            oss << "\n -> Command line argument " << i << " (='" << argv[i] << "') is invalid. <- \n\n\n\n";
            return oss.str();
        }

        std::string paramName, paramValue;

        // read a --my-opt=VALUE option. This gets transformed
        // into the parameter "MyOpt" with the value being "VALUE"
        if (argv[i][1] == '-') {
            std::string s(argv[i] + 2);
            // There is nothing after the '='
            if (s.size() == 0 || s[0] == '=')
            {
                std::ostringstream oss;
                oss << "\n -> Parameter name of argument " << i << " (='" << argv[i] << "')"
                    << " is empty. <- \n\n\n\n";
                return oss.str();
            }

            // capitalize the first character
            s[0] = toupper(s[0]);

            // parse argument
            unsigned int j = 0;
            while (true) {
                if (j >= s.size()) {
                    // encountered the end of the string, i.e. we
                    // have a parameter where the argument is empty
                    paramName = s;
                    paramValue = "";
                    break;
                }
                else if (s[j] == '=') {
                    // we encountered a '=' character. everything
                    // before is the name of the parameter,
                    // everything after is the value.
                    paramName = s.substr(0, j);
                    paramValue = s.substr(j+1);
                    break;
                }
                else if (s[j] == '.') {
                    // we encountered a '.' character, indicating
                    // the end of a group name, and need
                    // to captitalize the following character
                    if (s.size() == j)
                    {
                        std::ostringstream oss;
                        oss << "\n -> Parameter name of argument " << i << " ('" << argv[i] << "')"
                            << " is invalid (ends with a '.' character). <- \n\n\n\n";
                        return oss.str();
                    }
                    s[j+1] = toupper(s[j+1]);
                }
                else if (s[j] == '-') {
                    // remove all "-" characters and capitalize the
                    // character after them
                    s.erase(j, 1);
                    if (s.size() == j)
                    {
                        std::ostringstream oss;
                        oss << "\n -> Parameter name of argument " << i << " ('" << argv[i] << "')"
                            << " is invalid (ends with a '-' character). <- \n\n\n\n";
                        return oss.str();
                    }
                    else if (s[j] == '-')
                    {
                        std::ostringstream oss;
                        oss << "\n -> Malformed parameter name name in argument " << i << " ('" << argv[i] << "'): "
                            << "'--' in parameter name. <- \n\n\n\n";
                        return oss.str();
                    }
                    s[j] = toupper(s[j]);
                }

                ++j;
            }
        }
        else {
            // read a -MyOpt VALUE option
            paramName = argv[i] + 1;

            if (argc == i + 1 || argv[i+1][0] == '-') {
                std::ostringstream oss;
                oss << "\n -> No argument given for parameter '" << argv[i] << "'! <- \n\n\n\n";
                return oss.str();
            }

            paramValue = argv[i+1];
            //std::cout<<"paramName "<<paramName<<" paramValue "<<paramValue<<std::endl;
            ++i; // In the case of '-MyOpt VALUE' each pair counts as two arguments
        }

        // Put the key=value pair into the parameter tree
        paramTree[paramName] = paramValue;
    }
    //paramTree["Restart"] = argv[2];
    //paramTree["TEnd"] = argv[4];
    //paramTree["DtInitial"] = argv[6];
    //paramTree["EpisodeTime"] = argv[8];
    //paramTree["Positions0"] = argv[10];
    //paramTree["VgAlpha"] = argv[12];
    //paramTree["Vgn"] = argv[14];
    //paramTree["Swr"] = argv[16];
    //paramTree["Porosity"] = argv[18];
    //paramTree["Permeability"] = argv[20];
    //paramTree["BufferPower"] = argv[22];
    //paramTree["Kr"] = argv[24];
    //paramTree["Vmax"] = argv[26];
    //paramTree["Km"] = argv[28];
    //paramTree["PartitionCoeff"] = argv[30];
    //paramTree["InitialSoilFracK"] = argv[32];
    //paramTree["InitialSoilPressure"] = argv[34];;

    return "";
}

/*!
 * \ingroup Start
 *
 * \brief Provides a main function which reads in parameters from the
 *        command line and a parameter file.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param   argc    The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
 * \param   argv    The 'argv' argument of the main function: array of pointers to the argument strings
 * \param   usage   Callback function for printing the usage message
 */
template <class TypeTag>
double radialStart_(int argc,
           char **argv)
{
    //for(int i = 0; i < argc; ++i)
    //{
    //    //    argvector[i] = new char[argv[i].size() + 1];
    //    //    std::strcpy(argvector[i], argv[i].c_str());
    //    std::cout<<argv[i]<<" ";
    //}
    // some aliases for better readability
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using radialGridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using radialProblem = typename GET_PROP_TYPE(TypeTag, Problem);
    using radialTimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using radialParameterTree = typename GET_PROP(TypeTag, ParameterTree);

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////
    // fill the parameter tree with the options from the command line
    std::string s = radialreadOptions_(argc, argv, radialParameterTree::tree());
    //std::cout<<" S: "<<s<<std::endl;
    // parse the input file into the parameter tree
    // check first if the user provided an input file through the command line, if not use the default
    const auto radialParameterFileName = radialParameterTree::tree().hasKey("ParameterFile") ? radialParameterTree::tree().template get<std::string>("ParameterFile") : "";
    ParameterParser::parseInputFile(argc, argv, radialParameterTree::tree(), radialParameterFileName);
    // open and check whether the parameter file exists.
    std::ifstream radialParameterFile(radialParameterFileName.c_str());
    if (not radialParameterFile.is_open()) {
        // if the name of the file has been specified,
        // this must be an error. Otherwise proceed.
        if (radialParameterTree::tree().hasKey("ParameterFile"))
        {
            std::cout << "\n\t -> Could not open file '"
                      << radialParameterFileName
                      << "'. <- \n\n\n\n";
            return 1;
        }
    }
    else
    {
        // read parameters from the file without overwriting
        Dune::ParameterTreeParser::readINITree(radialParameterFileName,
                                               radialParameterTree::tree(),
                                               /*overwrite=*/false);
    }
    radialParameterFile.close();
    ////////////////////////////////////////////////////////////
    // check for some user debugging parameters
    ////////////////////////////////////////////////////////////


    bool printProps = false; // per default don't print all properties
    if (radialParameterTree::tree().hasKey("PrintProperties") || radialParameterTree::tree().hasKey("radialTimeManager.PrintProperties"))
        //printProps = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, radialTimeManager, PrintProperties);
        printProps = radialParameterTree::tree().template get<bool>("radialTimeManager.PrintProperties");

    if (printProps && mpiHelper.rank() == 0)
        Properties::print<TypeTag>();

    bool printParams = false; // per default print all properties
    if (radialParameterTree::tree().hasKey("PrintParameters") || radialParameterTree::tree().hasKey("radialTimeManager.PrintParameters"))
        //printParams = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, radialTimeManager, PrintParameters);
        printParams = radialParameterTree::tree().template get<bool>("radialTimeManager.PrintParameters");
    //bool printParams = false;
    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    try { radialGridCreator::makeGrid(); }
    catch (...) {
        std::string usageMessage = "\n\t -> Creation of the grid failed! <- \n\n";
        usageMessage += defaultUsageMessage(argv[0]);
        throw;
    }
    radialGridCreator::loadBalance();

    //////////////////////////////////////////////////////////////////////
    // run the simulation
    /////////////////////////////////////////////////////////////////////

    // check if we are about to restart a previously interrupted simulation
    bool restart = false;
    Scalar restartTime = 0;
    if (radialParameterTree::tree().hasKey("Restart") || radialParameterTree::tree().hasKey("radialTimeManager.Restart"))
    {
        //restartTime = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, radialTimeManager, Restart);
        restartTime = radialParameterTree::tree().template get<Scalar>("radialTimeManager.Restart");
        if (restartTime != 0) restart = true;
    }

    // read the initial time step and the end time (mandatory parameters)
    Scalar tEnd = radialParameterTree::tree().template get<Scalar>("radialTimeManager.TEnd"); //= GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, radialTimeManager, TEnd);
    Scalar dt = radialParameterTree::tree().template get<Scalar>("radialTimeManager.DtInitial");//= GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, radialTimeManager, DtInitial);

    //std::cout<<"tEnd: "<<tEnd<<std::endl;
    //std::cout<<"restartTime: "<<restartTime<<std::endl;
    //std::cout<<"dt: "<<dt<<std::endl;
    //restartTime= std::stod(radialParameterTree::tree()["radialTimeManager.Restart"]);
    //tEnd= std::stod(radialParameterTree::tree()["radialTimeManager.TEnd"]);
    //dt= std::stod(radialParameterTree::tree()["radialTimeManager.DtInitial"]);

    //std::string problemName(argv[4]);
    //std::cout<<"problemName: "<<problemName<<std::endl;
    //radialParameterTree::tree()["Name"] = problemName;
    //std::cout<<"problemName: "<<radialParameterTree::tree()["Name"]<<std::endl;

    //if (restartTime != 0) restart = true;
    //std::cout<<"tEnd: "<<tEnd<<std::endl;
    //std::cout<<"restartTime: "<<restartTime<<std::endl;
    //std::cout<<"dt: "<<dt<<std::endl;

    std::cout << "\033[1;33m" << "The radially symmetric simulation is called in time range: " << restartTime <<" - "
                << tEnd << " with inital timestep = " << dt << "\033[0m" << '\n';

    // instantiate and run the problem
    radialTimeManager radialtimeManager;
    //radialtimeManager.setTimeStepSize(6);
    //radialtimeManager.setEndTime(36);
    radialProblem radialproblem(radialtimeManager, radialGridCreator::grid().leafGridView());
    radialtimeManager.init(radialproblem, restartTime, dt, tEnd, restart);
    radialtimeManager.run();

    // print dumux end message and maybe the parameters for debugging
    if (mpiHelper.rank() == 0)
    {
        DumuxMessage::print(/*firstCall=*/false);

        if (printParams)
            Parameters::print<TypeTag>();
    }
    //std::cout<<"problem.uptake() "<< problemRadial.uptake() <<std::cin.ignore();
    std::cout << "\033[1;33m" << "The radially symmetric simulation is done !" << "\033[0m" << '\n';
    return radialproblem.uptakeRate();
    //return radialproblem.cumulativeUptake()/(tEnd-restartTime);
}

/*!
 * \ingroup Start
 *
 * \brief Provides a main function with error handling
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc  The number of command line arguments of the program
 * \param argv  The contents of the command line arguments of the program
 * \param usage Callback function for printing the usage message
 */
template <class TypeTag>
double radialStart(int argc,
          char **argv)
{
    try {
        return radialStart_<TypeTag>(argc, argv);
    }
    catch (ParameterException &e) {
        Parameters::print<TypeTag>();
        std::cerr << std::endl << e << ". Abort!" << std::endl;
        return 1;
    }
    catch (Dune::DGFException & e) {
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << std::endl;
    return 2;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
        return 3;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        return 4;
    }
}

} // end namespace Dumux

#endif
