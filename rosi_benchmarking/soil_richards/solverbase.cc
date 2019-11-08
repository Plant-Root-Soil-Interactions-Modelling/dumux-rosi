#include "solverbase.hh"

// initialize
#include <dune/common/parallel/mpihelper.hh> // in dune parallelization is realized with MPI
int DUMUX_VERSION = 100; // the future (needed in dumuxmessage.hh)
#include <dumux/common/dumuxmessage.hh> // for fun (a static class)
#include <dumux/common/parameters.hh> // global parameter tree with defaults and parsed from args and .input file

// createGrid
#include <dumux/io/grid/gridmanager.hh>



/**
 * Writes the Dumux welcome message, and creates the global Dumux parameter tree
 *
 * Normally you state an input file, that contains all parameters that are needed for the simulation.
 * SolverBase will try to set some of them dynamically.
 */
void SolverBase::initialize(std::vector<std::string> args)
{
    std::vector<char*> cargs;
    cargs.reserve(args.size());
    for(size_t i = 0; i < args.size(); ++i) {
        cargs.push_back(const_cast<char*>(args[i].c_str()));
    } // its a beautiful language

	auto& mpiHelper = Dune::MPIHelper::instance(cargs.size(), &cargs[0]); // of type MPIHelper, or FakeMPIHelper (in mpihelper.hh)

    // print dumux start message
    if (mpiHelper.rank() == 0) { // rank is the process number
    	Dumux::DumuxMessage::print(/*firstCall=*/true);
    }

    // parse command line arguments and input file
    Dumux::Parameters::init(cargs.size(), &cargs[0]);
}

/**
 * Creates a rectangular grid with a given resolution
 *
 * depending on the Grid you choose at compile time this will work, or not
 * e.g. it does not make a lot of sense for unstructured grids
 */
void SolverBase::createGrid(VectorType boundsMin, VectorType boundsMax, VectorType numberOfCells)
{
	// pass via (hasParamInGroup(modelParamGroup, "??? ", call init

	// makeStructuredGrid
}

/**
 * Creates a grid from a file
 *
 * depending on the Grid you choose at compile time it will accept the file type, or not.
 */
void SolverBase::createGrid(std::string file)
{
	GridManagerFix gridManager;
	// pass via (hasParamInGroup(modelParamGroup, "Grid.File", call init
	gridManager.init();
	grid = gridManager.gridPtr();
}

void SolverBase::simulate()
{
// todo
}

std::vector<VectorType> SolverBase::getPoints()
{
// todo
	return { VectorType(0) };
}

std::vector<VectorType> SolverBase::getCells()
{
// todo
	return { VectorType(0) };
}

int SolverBase::pickIndex(VectorType pos)
{
	// todo
	return 0;
}


double SolverBase::solutionAt(VectorType pos)
{
	return solution[pickIndex(pos)];
}
