#ifndef RICHARDS_YASP_SOLVER_H_
#define RICHARDS_YASP_SOLVER_H_

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>
namespace py = pybind11;

#include <config.h> // configuration file

// #include <dumux/discretization/box/fvelementgeometry.hh> // free function scvs(...)

// pick assembler and linear solver
#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>


// most includes are in solverbase
#include "solverbase.hh"

#include "../soil_richards/richardsproblem.hh" // the problem class. Defines some TypeTag types and includes its spatialparams.hh class
#include "../soil_richards/propertiesYasp.hh" // the property system related stuff (to pass types, used instead of polymorphism)
#include "../soil_richards/properties_nocoupling.hh" // dummy types for replacing the coupling types

/*
 * Define the type tag for this problem (in propertiesYasp.hh)
 */
using TypeTag = Dumux::Properties::TTag::RichardsBox; // RichardsCC, RichardsBox

/*
 * The problem
 */
using Problem = Dumux::RichardsProblem<TypeTag>;

/*
 * Pick assembler and solver
 */
using Assembler = Dumux::FVAssembler<TypeTag, Dumux::DiffMethod::numeric>;
using LinearSolver = Dumux::AMGBackend<TypeTag>;

/*
 * Name of the configuration of Problem, Assembler, LinearSolver
 */
std::string name = "RichardsYaspSolver";

/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
class RichardsYaspSolver : public SolverBase<Problem, Assembler, LinearSolver> {
public:

    std::vector<double> ls; // local sink

    virtual ~RichardsYaspSolver()
    { }

    /**
     * Sets the source term of the problem.
     * The source is given per cell (element),
     * as a map with element index as key, and source as value
     */
    virtual void setSource(std::map<int, double> source)
    {
        checkInitialized();
        int n = gridGeometry->gridView().size(0);
        std::vector<int> indices;
        ls.resize(n);
        indices.reserve(n);
        for (const auto& e : elements(gridGeometry->gridView())) { // local elements
            int gIdx = cellIdx->index(e); // global index
            auto eIdx = gridGeometry->elementMapper().index(e);
            if (source.count(gIdx)) {
                ls[eIdx] = source[gIdx]/24.*3600./1.e3; // g/day -> kg/s
            } else {
                ls[eIdx] = 0.;
            }
        }
        problem->setSource(&ls); // why a pointer? todo
    }

    /**
     * Returns the current solution for a single mpi process.
     * Gathering and mapping is done in Python
     */
    virtual std::vector<double> getSaturation()
    {
        checkInitialized();
        std::vector<double> s;
        int n;
        if (isBox) {
            n = gridGeometry->gridView().size(dim);
        } else {
            n = gridGeometry->gridView().size(0);
        }
        s.resize(n);
        for (const auto& element : Dune::elements(gridGeometry->gridView())) { // soil elements
            auto eIdx = gridGeometry->elementMapper().index(element);
            s[eIdx] = 0;
            auto fvGeometry = Dumux::localView(*gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(element);
            auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, x);
            for (const auto& scv : scvs(fvGeometry)) {
                s[eIdx] += elemVolVars[scv].saturation(0)*scv.volume();
            }
        }
        return s;
    }

    /**
     * Returns the total water volume [cm3] in domain of the current solution
     */
    virtual double getWaterVolume()
    {
        checkInitialized();
        double cVol = 0.;
        for (const auto& element : Dune::elements(gridGeometry->gridView())) { // soil elements
            auto fvGeometry = Dumux::localView(*gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(element);
            auto elemVolVars = Dumux::localView(gridVariables->curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, x);
            for (const auto& scv : scvs(fvGeometry)) {
                cVol += elemVolVars[scv].saturation(0)*scv.volume();
            }

        }
        return gridGeometry->gridView().comm().sum(cVol)*1.e6; // m3 -> cm3
    }

    /**
     * Uses the Dumux VTK Writer
     */
    virtual void writeDumuxVTK(std::string name)
    {
        Dumux::VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, name);
        // using VelocityOutput = PorousMediumFlowVelocityOutput<GridVariables> // <- can't get this type without TTAG :-(
        // vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
        vtkWriter.addVolumeVariable([](const auto& volVars){ return volVars.pressure(); }, "p (Pa)");
        vtkWriter.addVolumeVariable([](const auto& volVars){ return volVars.saturation(); }, "saturation");
        vtkWriter.write(0.0);
    }

};

using Solver = RichardsYaspSolver;

/**
 * Python binding of the Dumux solver base class
 */
PYBIND11_MODULE(richards_yasp_solver, m) {

    py::class_<Solver>(m, name.c_str())
    // initialization
        .def(py::init<>())
        .def("initialize", &Solver::initialize)
        .def("createGrid", (void (Solver::*)(std::string)) &Solver::createGrid, py::arg("modelParamGroup") = "") // overloads, defaults
        .def("createGrid", (void (Solver::*)(VectorType, VectorType, VectorType, std::string)) &Solver::createGrid,
            py::arg("boundsMin"), py::arg("boundsMax"), py::arg("numberOfCells"), py::arg("periodic") = "false false false") // overloads, defaults
        .def("readGrid", &Solver::readGrid)
        .def("getGridBounds", &Solver::getGridBounds)
        .def("setParameter", &Solver::setParameter)
        .def("getParameter", &Solver::getParameter)
        .def("initializeProblem", &Solver::initializeProblem)
    // simulation
        .def("simulate", &Solver::simulate, py::arg("dt"), py::arg("maxDt") = -1)
    // post processing (vtk naming)
        .def("getPoints", &Solver::getPoints) //
        .def("getCellCenters", &Solver::getCellCenters)
        .def("getDofCoordinates", &Solver::getDofCoordinates)
        .def("getPointIndices", &Solver::getPointIndices)
        .def("getCellIndices", &Solver::getCellIndices)
        .def("getDofIndices", &Solver::getDofIndices)
        .def("getSolution", &Solver::getSolution)
        .def("pickCell", &Solver::pickCell)
    // members
        .def_readonly("simTime", &Solver::simTime) // read only
        .def_readonly("rank", &Solver::rank) // read only
        .def_readonly("maxRank", &Solver::maxRank) // read only
        .def_readwrite("ddt", &Solver::ddt) // initial internal time step
    // useful
        .def("__str__",&Solver::toString)
        .def("checkInitialized", &Solver::checkInitialized)
    // added by class specialization
        .def("setSource", &Solver::setSource)
        .def("getSaturation",&Solver::getSaturation)
        .def("getWaterVolume",&Solver::getWaterVolume)
        .def("writeDumuxVTK",&Solver::writeDumuxVTK);
}

#endif

// #include <dumux/discretization/localview.hh> // free fucntion localView(...)
