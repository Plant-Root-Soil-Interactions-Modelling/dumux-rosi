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
//#include "math_expint.hh"
#include <boost/math/special_functions/expint.hpp>

// explicitly set the fvgeometry as we currently use another one as in stable
#include <dumux/multidimension/cellcentered/fvelementgeometry.hh>
#include <dumux/multidimension/box/fvelementgeometry.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/richards2cbuffer/richards2cbuffermodel.hh>

#include <dumux/material/fluidsystems/h2osolute.hh>
//! get the properties needed for subproblems
#include <dumux/multidimension/subproblemproperties.hh>

//#include "soiltestspatialparams1p2c.hh"
 #include "richardsbuffertestspatialparams.hh"

#define NONISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class SoilRichardsTwoCBufferTestProblem;

namespace Properties
{
NEW_TYPE_TAG(SoilRichardsTwoCBufferTestProblem, INHERITS_FROM(RichardsTwoCBuffer, RichardsBufferTestSpatialParams)); //RichardsTwoC from properties file
NEW_TYPE_TAG(SoilRichardsTwoCBufferTestBoxProblem, INHERITS_FROM(BoxModel, SoilRichardsTwoCBufferTestProblem));
NEW_TYPE_TAG(SoilRichardsTwoCBufferTestCCProblem, INHERITS_FROM(CCModel, SoilRichardsTwoCBufferTestProblem));

// Set the wetting phase
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OSOLUTE<TypeTag>);

// Set the grid type
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::YaspGrid<3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);

// explicitly set the fvgeometry as we currently use another one as in stable
SET_TYPE_PROP(SoilRichardsTwoCBufferTestBoxProblem, FVElementGeometry, Dumux::MultiDimensionBoxFVElementGeometry<TypeTag>);
SET_TYPE_PROP(SoilRichardsTwoCBufferTestCCProblem, FVElementGeometry, Dumux::MultiDimensionCCFVElementGeometry<TypeTag>);

// Set the problem property
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Problem, Dumux::SoilRichardsTwoCBufferTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, SpatialParams, Dumux::RichardsBufferTestSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(SoilRichardsTwoCBufferTestProblem, ProblemEnableGravity, true);

// Enable velocity output
SET_BOOL_PROP(SoilRichardsTwoCBufferTestProblem, VtkAddVelocity, true);

// Set the grid parameter group
SET_STRING_PROP(SoilRichardsTwoCBufferTestProblem, GridParameterGroup, "SoilGrid");

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
    //typedef typename GET_PROP_TYPE(TypeTag, PointSource) IntegrationPointSourve<TypeTag, unsigned int>

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables; //added

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag) GlobalProblemTypeTag;
    typedef typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager) CouplingManager;

public:
    SoilRichardsTwoCBufferTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "-soil";
        nameUptakeRate_ = this->name() + "_uptakeRate.log";
        //stating in the console whether mole or mass fractions are used

        if(useMoles)
        {
            std::cout<<"problem uses mole fractions"<<std::endl;
        }
        else
        {
            std::cout<<"problem uses mass fractions"<<std::endl;
        }
        pnRef_ = 1e5;
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
    { pointSources = this->couplingManager().bulkPointSources(); }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param source A single point sources for the scv
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void solDependentPointSource( PointSource& source,//PrimaryVariables &source, ///
                                  const Element &element,
                                  const FVElementGeometry &fvGeometry,
                                  const int scvIdx,
                                  const ElementVolumeVariables &elemVolVars) const
    {
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[conti0EqIdx];
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[conti0EqIdx] ;

        const auto& spatialParams = this->couplingManager().lowDimProblem().spatialParams();
        unsigned int rootEIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = spatialParams.Kr(rootEIdx);
        const Scalar rootRadius = spatialParams.radius(rootEIdx);
        const Scalar Vmax = spatialParams.Vmax();
        const Scalar Km = spatialParams.Km();
        const Scalar sigma = spatialParams.PartitionCoeff();

        PrimaryVariables sourceValues;
        sourceValues = 0.0;
        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3] = [m-1 kg s-1]
        if (GET_RUNTIME_PARAM(TypeTag, bool,
                                        BoundaryConditions.WithWaterUptake))
        	sourceValues[conti0EqIdx] = 2* M_PI *rootRadius * Kr *(pressure1D - pressure3D)
                                   *elemVolVars[scvIdx].density();
        //std::cout<<std::setprecision(17)<<pressure1D-pressure3D<<"\n";
        // sourceValues positive mean flow from root to soil
        // needs concentrations in soil and root
        Scalar c1D;
        if(useMoles)
            c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];
        else
            c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];

        Scalar c3D;
        if(useMoles)
            c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];
        else
            c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];

        const Scalar DiffValue = 0.0;
        //2* M_PI *rootRadius *DiffCoef_*(c1D - c3D)*elemVolVars[scvIdx].density(/*phaseIdx=*/0);

        //Advective flux term of transport
        Scalar AdvValue;
        if (sourceValues[conti0EqIdx]>0) // flow from root to soil
            //AdvValue = 0;
            AdvValue = sourceValues[conti0EqIdx]*c1D; // [m-1 kg s-1]* [kg m-3] =[m-4 kg2 s-1]
        else // flow from soil to root
            AdvValue = sourceValues[conti0EqIdx]*c3D;

        //Active flux - active uptake based on Michaeles Menten
        Scalar ActiveValue;
        ActiveValue = 0.;
        ActiveValue = -2*M_PI*rootRadius*Vmax*c3D*elemVolVars[scvIdx].density()
                    /(Km+c3D*elemVolVars[scvIdx].density()); //[m][kg m-2 s-1] [kg m-3] / [kg m-3] =[kg m-1 s-1]

        sourceValues[transportEqIdx] = (sigma*(AdvValue + DiffValue) + (1-sigma)*ActiveValue)*source.quadratureWeight()*source.integrationElement();
        sourceValues[conti0EqIdx] *= source.quadratureWeight()*source.integrationElement();
        source =  sourceValues;
    }

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
        {
            if((globalPos[1] > this->bBoxMax()[1] - eps_)||(globalPos[1] < this->bBoxMin()[1] + eps_)||
              (globalPos[0] > this->bBoxMax()[0] - eps_)||(globalPos[0] < this->bBoxMin()[0] + eps_))
            {
                values.setDirichlet(transportEqIdx);
                values.setDirichlet(conti0EqIdx);
            }
            else
                values.setAllNeumann();
        }
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
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    using ParentType::neumann;
    void neumann(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        priVars[conti0EqIdx] = 0.0;
        priVars[transportEqIdx] = 0.0;
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
        ScalarField& ratioPressure = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ratioWaterContent = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ratioMassFraction = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceP = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& PeclectNo = *(this->resultWriter().allocateManagedBuffer(numDofs));

	PeclectNo = 0.0;
        sourceP = 0.0;
        sourceC = 0.0;
        Scalar totalSourceC, totalSourceP;
        totalSourceC = 0.0;
        totalSourceP = 0.0;

        ratioPressure = 0.0;
        ratioWaterContent = 0.0;
        ratioMassFraction = 0.0;


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
                PrimaryVariables priVars;
                GlobalPosition globalPos = element.geometry().center();
                initial_(priVars, globalPos);
                ratioPressure[dofGlobalIdx] = elemVolVars[0].pressure()/priVars[conti0EqIdx];

                ratioMassFraction[dofGlobalIdx] = elemVolVars[0].massFraction(transportCompIdx)/priVars[transportEqIdx];

                Scalar initialPc_ = pnRef_-priVars[conti0EqIdx];
                Scalar initialSw_ = MaterialLaw::sw(this->spatialParams().materialLawParams(globalPos),initialPc_);
                ratioWaterContent[dofGlobalIdx] = elemVolVars[0].saturation(phaseIdx)/initialSw_;

                PrimaryVariables values;
                this->scvPointSources(values, element, fvGeometry, 0, elemVolVars);
                sourceP[dofGlobalIdx] += values[conti0EqIdx] * fvGeometry.subContVol[0].volume
                                             * this->boxExtrusionFactor(element, fvGeometry, 0);
                totalSourceP += sourceP[dofGlobalIdx];
                sourceC[dofGlobalIdx] += values[transportEqIdx] * fvGeometry.subContVol[0].volume
                                             * this->boxExtrusionFactor(element, fvGeometry, 0);
                totalSourceC += sourceC[dofGlobalIdx];
            }
            //////
            FluxVariables fluxVars;
            Scalar advFluxes = 0;
            Scalar diffFluxes = 0;

            for (int scvfIdx = 0; scvfIdx < fvGeometry.numScvf; scvfIdx++)
            {
                fluxVars.update(*this, element, fvGeometry, scvfIdx, elemVolVars, false /* oldSol? */);

                advFluxes +=
                        -fluxVars.faceVolumeFlux(0) *elemVolVars[0].density()*elemVolVars[0].massFraction(transportCompIdx);

                Scalar tmp = (fluxVars.massFractionGrad(transportCompIdx)*fluxVars.face().normal);
                tmp *= fluxVars.porousDiffCoeff() * fluxVars.density();
                diffFluxes += tmp;
            }

            PeclectNo[dofGlobalIdx] = advFluxes/diffFluxes;
            //////
        }
        // attach data to the vtk output
        std::cout << "Integrated mass source (3D): " << totalSourceP << std::endl;
        std::cout << "Integrated concentration source (3D): " << totalSourceC << std::endl;

        // attach data to the vtk output
        this->resultWriter().attachDofData(sourceP, "mass source (kg/s)", isBox);
        this->resultWriter().attachDofData(sourceC, "concentration source (kg/s)", isBox);
        this->resultWriter().attachDofData(ratioPressure, "pw/pw0", isBox);
        this->resultWriter().attachDofData(ratioWaterContent, "theta/theta0", isBox);
        this->resultWriter().attachDofData(ratioMassFraction, "w/w0", isBox);
        this->resultWriter().attachDofData(PeclectNo, "PeclectNo", isBox);
    }


    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        unsigned numDofs = this->model().numDofs();

        // create required scalar fields for the vtk output
        ScalarField& sourceP = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceC = *(this->resultWriter().allocateManagedBuffer(numDofs));

        sourceP = 0.0;
        sourceC = 0.0;

        Scalar totalSourceC, totalSourceP;
        totalSourceC = 0.0;
        totalSourceP = 0.0;

        // iterate over all elements
        for (const auto& element : elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            ElementBoundaryTypes bcTypes;
            bcTypes.update(*this, element, fvGeometry);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

            // output pressure
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                auto dofGlobalIdx = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);
                // only call the source function if we are not initializing
                if (!(this->timeManager().time() < 0.0))
                {
                    PrimaryVariables values;
                    this->scvPointSources(values, element, fvGeometry, scvIdx, elemVolVars);
                    sourceP[dofGlobalIdx] += values[conti0EqIdx] * fvGeometry.subContVol[scvIdx].volume
                                            * this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                    totalSourceP += sourceP[dofGlobalIdx];
                    sourceC[dofGlobalIdx] += values[transportEqIdx] * fvGeometry.subContVol[scvIdx].volume
                                            * this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                    totalSourceC += sourceC[dofGlobalIdx];
                }
            }
        }
        logFile_.open(this->name() + "_uptakeRate.log", std::ios::app);
        logFile_ << this->timeManager().time() + this->timeManager().timeStepSize()
                << " " << -totalSourceP
                << " " << -totalSourceC
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
    bool shouldWriteOutput() const
    { return (this->timeManager().willBeFinished()) or (this->timeManager().time()<0); }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

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

        priVars[pressureIdx] = 1e5 - pc - 1000*9.81*globalPos[2];

        priVars[massOrMoleFracIdx]= GET_RUNTIME_PARAM(TypeTag,
                                       Scalar,
                                       BoundaryConditions.InitialSoluteMassFracInSoil);
    };

    std::ofstream logFile_;
    const Scalar eps_ = 1e-9;
    Scalar episodeTime, temperature_, pnRef_, pw_;
    std::string name_, nameUptakeRate_;
    std::shared_ptr<CouplingManager> couplingManager_;
    //Scalar DiffCoef_;
};

} //end namespace

#endif
