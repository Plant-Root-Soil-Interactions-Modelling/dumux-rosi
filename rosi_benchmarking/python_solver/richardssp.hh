#ifndef RICHARDS_SP_SOLVER_H_
#define RICHARDS_SP_SOLVER_H_

// most includes are in solverbase
#include "solverbase.hh"

// writeDumuxVTK
#include <dumux/io/vtkoutputmodule.hh>


/**
 * Adds solver functionality, that specifically makes sense for Richards equation
 */
template<class Problem, class Assembler, class LinearSolver>
class RichardsSP : public SolverBase<Problem, Assembler, LinearSolver> {
public:

    std::vector<double> ls; // local sink [kg/s]

    virtual ~RichardsSP() { }

    /**
     * Sets the source term of the problem [kg/s].
     *
     * The source is given per cell (Dumux element),
     * as a map with global element index as key, and source as value
     *
     * for simplicity avoiding mpi broadcasting or scattering
     */
    virtual void setSource(std::map<int, double> source) {
        this->checkInitialized();
        int n = this->gridGeometry->gridView().size(0);
        std::vector<int> indices;
        ls.resize(n);
        indices.reserve(n);
        for (const auto& e : elements(this->gridGeometry->gridView())) { // local elements
            int gIdx = this->cellIdx->index(e); // global index
            auto eIdx = this->gridGeometry->elementMapper().index(e);
            if (source.count(gIdx)>0) {
                ls[eIdx] = source[gIdx];
            } else {
                ls[eIdx] = 0.;
            }
        }
        this->problem->setSource(&ls); // why a pointer? todo
    }

    /**
     * TODO
     * applies source (operator splitting)
     * limits with wilting point from below, and with full saturation
     */
    virtual void applySource(std::map<int, double> init) {
        if (isBox) {
            throw std::invalid_argument("SolverBase::setInitialCondition: Not implemented yet (sorry)");
        } else {
            for (const auto& e : Dune::elements(gridGeometry->gridView())) {
                int eIdx = gridGeometry->elementMapper().index(e);
                int gIdx = cellIdx->index(e);
                /* equations... */

                x[eIdx] = init[gIdx];
            }
        }
    }

    /**
     * Returns the current solution for a single mpi process in cm pressure head. @see SolverBase::getSolution
     * Gathering and mapping is done in Python TODO pass equ id, return vector<double>
     */
    virtual std::vector<double> getSolutionHead(int eqIdx = 0) {
        std::vector<double> sol = this->getSolution(eqIdx);
        std::transform(sol.begin(), sol.end(), sol.begin(), std::bind1st(std::plus<double>(), -1.e5));
        std::transform(sol.begin(), sol.end(), sol.begin(), std::bind1st(std::multiplies<double>(), 100. / 1.e3 / 9.81));
        return sol;
    }

    /*
     * TODO setInitialConditions, std::map<int, double> ic,
     * TODO setLayers(std::map<int, int> l)
     */

    /**
     * Returns the current solution for a single mpi process.
     * Gathering and mapping is done in Python
     *
     * TODO does this work for box ??? should i loop over vertices?
     */
    virtual std::vector<double> getWaterContent() {
        this->checkInitialized();
        std::vector<double> s;
        int n;
        if (this->isBox) {
            n = this->gridGeometry->gridView().size(this->dim);
        } else {
            n = this->gridGeometry->gridView().size(0);
        }
        s.resize(n);
        for (const auto& element : Dune::elements(this->gridGeometry->gridView())) { // soil elements

            auto eIdx = this->gridGeometry->elementMapper().index(element);
            s[eIdx] = 0;
            auto fvGeometry = Dumux::localView(*this->gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(element);
            auto elemVolVars = Dumux::localView(this->gridVariables->curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, this->x);

            int c = 0;
            for (const auto& scv : scvs(fvGeometry)) {
                c++;
                s[eIdx] += elemVolVars[scv].waterContent();
            }
            s[eIdx] /= c; // mean value
        }
        return s;
    }

    /**
     * Returns the total water volume [m3] within the domain
     */
    virtual double getWaterVolume()
    {
        this->checkInitialized();
        double cVol = 0.;
        for (const auto& element : Dune::elements(this->gridGeometry->gridView())) { // soil elements
            auto fvGeometry = Dumux::localView(*this->gridGeometry); // soil solution -> volume variable
            fvGeometry.bindElement(element);
            auto elemVolVars = Dumux::localView(this->gridVariables->curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, this->x);
            for (const auto& scv : scvs(fvGeometry)) {
                cVol += elemVolVars[scv].saturation(0)*scv.volume();
            }
        }
        return this->gridGeometry->gridView().comm().sum(cVol);
    }

    /**
     * Uses the Dumux VTK Writer
     */
    virtual void writeDumuxVTK(std::string name)
    {
        Dumux::VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*this->gridVariables, this->x, name);
        // using VelocityOutput = PorousMediumFlowVelocityOutput<GridVariables> // <- can't get this type without TTAG :-(
        // vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
        vtkWriter.addVolumeVariable([](const auto& volVars){ return volVars.pressure(); }, "p (Pa)");
        vtkWriter.addVolumeVariable([](const auto& volVars){ return volVars.saturation(); }, "saturation");
        vtkWriter.write(0.0);
    }

protected:
    using SolutionVector = typename Problem::SolutionVector;
    using GridVariables = typename Problem::GridVariables;

};

/**
 * pybind11
 */
template<class Problem, class Assembler, class LinearSolver>
void init_richardssp(py::module &m, std::string name) {
    using RichardsSP = RichardsSP<Problem, Assembler, LinearSolver>;
	py::class_<RichardsSP, SolverBase<Problem, Assembler, LinearSolver>>(m, name.c_str())
   .def(py::init<>())
   .def("setSource", &RichardsSP::setSource)
   .def("getSolutionHead", &RichardsSP::getSolutionHead, py::arg("eqIdx") = 0)
   .def("getWaterContent",&RichardsSP::getWaterContent)
   .def("getWaterVolume",&RichardsSP::getWaterVolume)
   .def("writeDumuxVTK",&RichardsSP::writeDumuxVTK);
}


#endif
