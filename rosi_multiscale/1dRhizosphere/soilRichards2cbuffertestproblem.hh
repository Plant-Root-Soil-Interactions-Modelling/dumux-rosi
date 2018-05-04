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
 * \brief A sub problem for the richards problem
 */
#ifndef DUMUX_SOIL_RICHARD2C_BUFFER_TEST_PROBLEM_HH
#define DUMUX_SOIL_RICHARD2C_BUFFER_TEST_PROBLEM_HH

#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <boost/math/special_functions/expint.hpp>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/richards2cCylindrical1d/richards2cbuffermodel.hh>
#include <dumux/material/fluidsystems/h2osolute.hh>
#include "richardsbuffertestspatialparams.hh"

#define NONISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class SoilRichardsTwoCBufferTestProblem;

namespace Properties
{
NEW_TYPE_TAG(SoilRichardsTwoCBufferTestProblem, INHERITS_FROM(RichardsTwoCBufferRadiallySymmetric, RichardsBufferTestSpatialParams)); //RichardsTwoC from properties file
NEW_TYPE_TAG(SoilRichardsTwoCBufferTestBoxProblem, INHERITS_FROM(BoxModel, SoilRichardsTwoCBufferTestProblem));
NEW_TYPE_TAG(SoilRichardsTwoCBufferTestCCProblem, INHERITS_FROM(CCModel, SoilRichardsTwoCBufferTestProblem));

// Set the wetting phase
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OSOLUTE<TypeTag>);

// Set the grid type
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::YaspGrid<1, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 1> >);

// Set the problem property
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Problem, Dumux::SoilRichardsTwoCBufferTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, SpatialParams, Dumux::RichardsBufferTestSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(SoilRichardsTwoCBufferTestProblem, ProblemEnableGravity, false);

// Enable velocity output
SET_BOOL_PROP(SoilRichardsTwoCBufferTestProblem, VtkAddVelocity, true);

// Set the grid parameter group
SET_STRING_PROP(SoilRichardsTwoCBufferTestProblem, GridParameterGroup, "SoilGrid");

//SET_BOOL_PROP(SoilRichardsTwoCTestProblem, UseHead, false);
//SET_BOOL_PROP(SoilOnePTwoCTestProblem, UseMoles, true);
SET_BOOL_PROP(SoilRichardsTwoCBufferTestProblem, UseMoles, false);
}

/*!
 * \ingroup OnePBoxModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * In the middle of the domain, a lens with low permeability (\f$K=10e-12\f$)
 * compared to the surrounding material (\f$ K=10e-10\f$) is defined.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p -parameterFile test_box1p.input</tt> or
 * <tt>./test_cc1p -parameterFile test_cc1p.input</tt>
 *
 * The same parameter file can be also used for 3d simulation but you need to change line
 * <tt>typedef Dune::SGrid<2,2> type;</tt> to
 * <tt>typedef Dune::SGrid<3,3> type;</tt> in the problem file
 * and use <tt>1p_3d.dgf</tt> in the parameter file.
 */
template <class TypeTag>
class SoilRichardsTwoCBufferTestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx
#endif
    };
     enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
        transportCompIdx = Indices::transportCompIdx,
        phaseIdx = Indices::phaseIdx,
#if NONISOTHERMAL
        energyEqIdx = Indices::energyEqIdx
#endif
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes; //added
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables; //added

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    SoilRichardsTwoCBufferTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "-soil";
        nameUptakeRate_ = this->name() + "_uptakeRate.log";
        if(useMoles)
        {
            std::cout<<"problem uses mole fractions"<<std::endl;
        }
        else
        {
            std::cout<<"problem uses mass fractions"<<std::endl;
        }
        pnRef_ = 1e5;
        uptakeC_ana_= 0;
        uptakeC_num_= 0;
        lastTime_=0;
        if (std::ifstream(nameUptakeRate_))
            std::remove(nameUptakeRate_.c_str());
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    {
        return name_;
    }

    /*!
     * \brief Return the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
//#if !NONISOTHERMAL
//
    Scalar temperature() const
    {   Scalar soilTemp_;
        soilTemp_ = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.SoilTemperature);
    return soilTemp_ ; } // 10C
//
//#endif
    /*
      * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar referencePressure(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    { return pnRef_; }


    // TODO identify source term
    // for example by reading from a vector that contains the source term for
    // each element
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx) const
    {
        values = 0.0;
    }

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of Dumux::PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { }


    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
            const GlobalPosition &globalPos) const
    {
        bool ZeroNeumann = GET_RUNTIME_PARAM(TypeTag,
                                        bool,
                                        BoundaryConditions.ZeroNeumannSoil);
        if (ZeroNeumann)
            values.setAllNeumann();
        else
            if (globalPos[0] < this->bBoxMin()[0] + eps_)
                values.setAllNeumann();
            else
                values.setAllDirichlet();

    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &priVars,
                        const GlobalPosition &globalPos) const
    {
        initial_(priVars, globalPos);
    }


    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
     void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
        auto posInnerRadius = std::min(element.geometry().corner(0)[0],element.geometry().corner(1)[0]);

        if (posInnerRadius < this->bBoxMin()[0] + eps_)
        {
            Scalar TranspirationRate= GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.TranspirationRate);
            if (GET_RUNTIME_PARAM(TypeTag, bool, BoundaryConditions.DiurnalFluctuation))
                TranspirationRate *= std::abs(std::cos(this->timeManager().time()/(24*3600)*(2*M_PI))+1);
            values[conti0EqIdx] = TranspirationRate*posInnerRadius;

            Scalar Vmax = GET_RUNTIME_PARAM(TypeTag,
                                            Scalar,
                                            SpatialParams.Vmax);
            Scalar Km = GET_RUNTIME_PARAM(TypeTag,
                                            Scalar,
                                            SpatialParams.Km);

            Scalar compFracRoot = GET_RUNTIME_PARAM(TypeTag,
                                            Scalar,
                                            BoundaryConditions.InitialSoluteMassFracInRoot);
            Scalar diffCoeffRootMembrane = GET_RUNTIME_PARAM(TypeTag,
                                            Scalar,
                                            SpatialParams.DiffCoeffRootMembrane);
            Scalar PartitionCoeff = GET_RUNTIME_PARAM(TypeTag,
                                            Scalar,
                                            SpatialParams.PartitionCoeff);

            Scalar c3D;
            if(useMoles)
                c3D = elemVolVars[scvIdx].moleFraction(1);
            else
                c3D = elemVolVars[scvIdx].massFraction(1);

            Scalar passive_uptake = 0;
            if (values[conti0EqIdx]>0)
            {
            	passive_uptake += (values[conti0EqIdx]*c3D + diffCoeffRootMembrane*(c3D - compFracRoot))
            						*elemVolVars[scvIdx].density()*posInnerRadius;
            }

            Scalar active_uptake = 0;
            active_uptake =  Vmax*c3D*elemVolVars[scvIdx].density()
                                    /(Km+c3D*elemVolVars[scvIdx].density())
                                    *posInnerRadius;
            values[transportEqIdx] = PartitionCoeff*passive_uptake + (1-PartitionCoeff)*active_uptake;
        }
        else
        {
            values[transportEqIdx] = 0;
            values[conti0EqIdx] = 0;
        }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param values Storage for all primary variables of the initial condition
     * \param globalPos The position for which the boundary type is set
     */
    void initialAtPos(PrimaryVariables &values,
                 const GlobalPosition &globalPos) const
    {
      initial_(values, globalPos);
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    void addOutputVtkFields()
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        unsigned numDofs = this->model().numDofs();

        // create required scalar fields for the vtk output
        ScalarField& diffPressureHead = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& diffWaterContent = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ratioMassFraction = *(this->resultWriter().allocateManagedBuffer(numDofs));

        diffPressureHead = 0.0;
        diffWaterContent = 0.0;
        ratioMassFraction = 0.0;
        Scalar pressureMax=0;

        // iterate over all elements
        for (const auto& element : elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

            auto dofGlobalIdx = this->model().dofMapper().subIndex(element, 0, dofCodim);

            // only call the source function if we are not initializing
            if (!(this->timeManager().time() < 0))
            {
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
    		        auto dofGlobalIdx = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);

    	            PrimaryVariables initialPriVars;
    	            GlobalPosition globalPos = element.geometry().center();
    	            initial_(initialPriVars, globalPos);
    	            diffPressureHead[dofGlobalIdx] = (elemVolVars[scvIdx].pressure()-initialPriVars[conti0EqIdx])/9.81/elemVolVars[scvIdx].density();

    	            ratioMassFraction[dofGlobalIdx] = elemVolVars[scvIdx].massFraction(transportCompIdx)/initialPriVars[transportEqIdx];

    	            Scalar initialPc_ = pnRef_-initialPriVars[conti0EqIdx];
    	            Scalar initialSw_ = MaterialLaw::sw(this->spatialParams().materialLawParams(globalPos),initialPc_);
    	            diffWaterContent[dofGlobalIdx] = (elemVolVars[scvIdx].saturation(phaseIdx)-initialSw_)*elemVolVars[scvIdx].porosity();
                }
            }
        }
        // attach data to the vtk output
        this->resultWriter().attachDofData(diffPressureHead, "diffPressureHead", isBox);
        this->resultWriter().attachDofData(diffWaterContent, "diffWaterContent", isBox);
        this->resultWriter().attachDofData(ratioMassFraction, "w/w0", isBox);
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        ElementVolumeVariables elemVolVars;
        FVElementGeometry fvGeometry;
        ElementBoundaryTypes bcTypes;

        PrimaryVariables flux(0.0);
        bool LoopOverElements = true;

        Scalar rootRadius = GET_RUNTIME_PARAM(TypeTag,
                                            Scalar,
                                            SpatialParams.rootRadius);

        // Loop over elements
        for (const auto& element : elements(this->gridView()))
        {
        	if (LoopOverElements == false) break ;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (!intersection.boundary())
                        continue;
		        Scalar posInnerRadius = std::min(element.geometry().corner(0)[0],element.geometry().corner(1)[0]);
		        if (posInnerRadius < this->bBoxMin()[0] + eps_)
		        {
		        	fvGeometry.update(this->gridView(), element);
	                elemVolVars.update(*this, element, fvGeometry,false);
	                unsigned bfIdx = intersection.indexInInside();
	                PrimaryVariables values(0.0);
		            this->solDependentNeumann(values,
	                                                  element,
	                                                  fvGeometry,
	                                                  intersection,
	                                                  /*scvIdx=*/0,
	                                                  bfIdx,
	                                                  elemVolVars);
	                //values *= intersection.geometry().volume();
	                values /= rootRadius;
	                flux += values;
	                LoopOverElements = false;
	            }
	        	break;
            }
        }
        logFile_.open(name_ + "_uptakeRate.log", std::ios::app);
        logFile_ << this->timeManager().time()+this->timeManager().timeStepSize()
                << " " <<flux[0]
                << " " <<flux[1]
                << std::endl;
        logFile_.close();
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behavior is to write out the solution for
     * every time step. This function is intended to be overwritten by the
     * implementation.
     */
    //bool shouldWriteOutput() const
    //{ return (this->timeManager().willBeFinished()) or (this->timeManager().time()<0); }

private:
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        Scalar sw = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilWaterContent)/
                    GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        SpatialParams.Porosity);
        Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParams(globalPos),sw);

        priVars[pressureIdx] = 1e5 - pc;


        priVars[massOrMoleFracIdx]= GET_RUNTIME_PARAM(TypeTag,
                                       Scalar,
                                       BoundaryConditions.InitialSoluteMassFracInSoil);
    };
    std::ofstream logFile_;
    const Scalar eps_ = 1e-9;
    Scalar episodeTime, temperature_, pnRef_, pc_, sw_, pw_, uptakeC_ana_, uptakeC_num_, lastTime_;
    std::string name_, nameUptakeRate_;
};

} //end namespace

#endif
