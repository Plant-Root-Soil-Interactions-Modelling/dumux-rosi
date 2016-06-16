
/*!
 * \file
 *
 * \brief Test for the Richards box model.
 */
#include <config.h> // from cmake

#include "myrichardsproblem.hh"  // all relevant is in here

#include <dumux/common/start.hh> // start function and processing of input file

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
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe list of mandatory options for this program is:\n"
                                        "\t-TimeManager.TEnd      End of the simulation [s] \n"
                                        "\t-TimeManager.DtInitial Initial timestep size [s] \n"
                                        "\t-Grid.File             Name of the file containing the grid \n"
                                        "\t                       definition in DGF format\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
    typedef TTAG(RichardsBoxProblem) PTT;
    return Dumux::start<PTT>(argc, argv, usage);
}
