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
#include "1dRadialRichards2cCoupling.hh"

// explicitly set the fvgeometry as we currently use another one as in stable
#include <dumux/multidimension/cellcentered/fvelementgeometry.hh>
#include <dumux/multidimension/box/fvelementgeometry.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/richards2cbuffer/richards2cbuffermodel.hh>

#include <dumux/material/fluidsystems/h2ok.hh>
//#include "h2oX.hh"
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
              Dumux::FluidSystems::H2OK<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

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
        uptakeC_num_ = 0.0;
        uptakeC_ana_ = 0.0;
    }

    void preTimeStep()
    {
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
        unsigned int rootEIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        unsigned int soilEIdx = this->couplingManager().pointSourceData(source.id()).bulkElementIdx();
        std::pair<unsigned int, unsigned int> Intersection = std::make_pair(rootEIdx, soilEIdx);

        PrimaryVariables sourceValues;
        sourceValues=0.0;
        if (this->couplingManager().reuseCouplingSources_[Intersection])
        {
            sourceValues[conti0EqIdx] = 0;
            sourceValues[transportEqIdx] = this->couplingManager().couplingSources_[Intersection];
            //std::cout<<" Reuse sources in SOIL system : "<<this->couplingManager().couplingSources_[Intersection]<<"\n";
        }
        else
        {
            const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[conti0EqIdx];
            const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[conti0EqIdx] ;

            const auto& spatialParams = this->couplingManager().lowDimProblem().spatialParams();

            //std::cout << "\033[1;33m" << "SOIL - rootEIdx: " << rootEIdx <<" - soilEIdx: " << soilEIdx << "- Position: " << source.position() << "\033[0m" << '\n';
            const Scalar Kr = spatialParams.Kr(rootEIdx);
            const Scalar rootRadius = spatialParams.radius(rootEIdx);
            const Scalar Vmax = spatialParams.Vmax();
            const Scalar Km = spatialParams.Km();
            const Scalar sigma = spatialParams.PartitionCoeff();


            // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3] = [m-1 kg s-1]
            sourceValues[conti0EqIdx] = 2* M_PI *rootRadius * Kr *(pressure1D - pressure3D)
                                       *elemVolVars[scvIdx].density();

            // sourceValues positive mean flow from root to soil
            // needs concentrations in soil and root
            Scalar c1D;
            if(useMoles)
                //c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx]*elemVolVars[0].molarDensity();
                c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];
            else
                //c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx]*1000;//elemVolVars[0].density(); //unit: kg/m3
                c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];
            //std::cout << "concentrations c1D " <<c1D<< std::endl;
            Scalar c3D;
            if(useMoles)
                //c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx]*elemVolVars[0].molarDensity();
                c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];
            else
                //c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx]*1000;//elemVolVars[0].density(); //unit: kg/m3
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
            //ActiveValue = -2*M_PI*rootRadius*Vmax*c3D*elemVolVars[scvIdx].density()
            //            /(Km+c3D*elemVolVars[scvIdx].density()); //[m][kg m-2 s-1] [kg m-3] / [kg m-3] =[kg m-1 s-1]

            Scalar Frac0=GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.InitialSoilFracK);
            OneDRichards2CRadiallySymmetric<TypeTag> OnedRadialSol(Km, Vmax, rootRadius, Frac0, this->timeManager().timeStepSize());
            OnedRadialSol.setEffDiffCoeff(elemVolVars[scvIdx].effDiffCoeff());
            OnedRadialSol.setSaturation(elemVolVars[scvIdx].saturation(phaseIdx));
            OnedRadialSol.setPorosity(elemVolVars[scvIdx].porosity());
            OnedRadialSol.setDensity(elemVolVars[scvIdx].density());
            OnedRadialSol.setBuffer(elemVolVars[scvIdx].buffer());

            //ActiveValue = OnedRadialSol.ActiveUptakeAnalyticalROOSE2001();
            //OnedRadialSol.setName("test_rosi_RadialCoupled");
            OnedRadialSol.setName(GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile));
            Scalar restartTimeRhizosphere = this->timeManager().time();
            Scalar endTimeRhizosphere;
            if (this->timeManager().time()+this->timeManager().timeStepSize()>this->timeManager().endTime())
                endTimeRhizosphere = this->timeManager().time();
            else
                endTimeRhizosphere = this->timeManager().time()+this->timeManager().timeStepSize();
            Scalar dtRhizosphere = (endTimeRhizosphere-this->timeManager().time())/GET_RUNTIME_PARAM(TypeTag, Scalar, radialTimeManager.TimeStepRatio);
            OnedRadialSol.setRestartTime(restartTimeRhizosphere);
            OnedRadialSol.setTEnd(endTimeRhizosphere);
            OnedRadialSol.setDtInitial(dtRhizosphere);

            //estimate and set the rhizosphere Radius
            Scalar rootLength = this->couplingManager().lowDimVolume(element)/(M_PI*rootRadius*rootRadius);
            Scalar rhizosphereRadius = pow(element.geometry().volume()/rootLength/M_PI,0.5);
            //if (rhizosphereRadius> 10*rootRadius)
            //    rhizosphereRadius = 10*rootRadius;
            OnedRadialSol.setOuterBoundary(rhizosphereRadius);

            OnedRadialSol.setVgAlpha(GET_RUNTIME_PARAM(TypeTag, Scalar, materialParams.VgAlpha));
            OnedRadialSol.setVgn(GET_RUNTIME_PARAM(TypeTag, Scalar, materialParams.Vgn));
            OnedRadialSol.setSwr(GET_RUNTIME_PARAM(TypeTag, Scalar, materialParams.Swr));
            OnedRadialSol.setPermeability(this->spatialParams().intrinsicPermeability(element,fvGeometry,scvIdx));
            OnedRadialSol.setKr(spatialParams.Kr(rootEIdx));
            OnedRadialSol.setPartitionCoeff(spatialParams.PartitionCoeff());
            OnedRadialSol.setInitialSoilPressure(GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.InitialSoilPressure));
            OnedRadialSol.setInitialSoilNutrient(Frac0);
            OnedRadialSol.setProblemName(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name)
                                        +"_rhizosphere_"+std::to_string(rootEIdx)+"_"+std::to_string(soilEIdx));

            ActiveValue = -2*M_PI*rootRadius*OnedRadialSol.ActiveUptakeNumerical();
            Scalar TransportValue = (sigma*(AdvValue + DiffValue) + (1-sigma)*ActiveValue)*source.quadratureWeight()*source.integrationElement();
            sourceValues[transportEqIdx] = TransportValue;
            sourceValues[conti0EqIdx] *= source.quadratureWeight()*source.integrationElement();
            source =  sourceValues;
            this->couplingManager().setCouplingSources(Intersection, TransportValue);
            this->couplingManager().setReuseCouplingSources(Intersection);
        }
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
        if((globalPos[1] > this->bBoxMax()[1] - eps_)||(globalPos[1] < this->bBoxMin()[1] + eps_)||
        	(globalPos[0] > this->bBoxMax()[0] - eps_)||(globalPos[0] < this->bBoxMin()[0] + eps_))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
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
        Scalar pw_ = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilPressure);
        priVars[pressureIdx] = pw_;

        priVars[massOrMoleFracIdx] = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilFracK);
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
        ScalarField& sourceP = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& relErrC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& absErrC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& AnalyticalC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& AnalyticalS = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& Distance = *(this->resultWriter().allocateManagedBuffer(numDofs));

        sourceP = 0.0;
        sourceC = 0.0;
        relErrC = 0.0;
        absErrC = 0.0;
        AnalyticalC = 0.0;
        AnalyticalS = 0.0;
        Distance = 0.0;

        Scalar Km_=GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Km);
        Scalar Vmax_=GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Vmax);
        Scalar rootRadius_=GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.rootRadius);
        Scalar Frac0=GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.InitialSoilFracK);

        Scalar totalSourceC, totalSourceP;
        totalSourceC = 0.0;
        totalSourceP = 0.0;

        Scalar anaUptakeRate;
        anaUptakeRate = 0;

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
                Scalar Cinf_ = Frac0*elemVolVars[scvIdx].density();
                Scalar Cinf_dl  = Cinf_/Km_;

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

                    Scalar x = pow((pow(element.geometry().center()[0],2) + pow(element.geometry().center()[1],2)),0.5);
                    if (x<rootRadius_)
                        x=rootRadius_;
                    Distance[dofGlobalIdx] = x;

                    Scalar lambda = Vmax_*rootRadius_/(elemVolVars[scvIdx].effDiffCoeff()*Km_);
                            lambda /=(elemVolVars[scvIdx].saturation(phaseIdx)* elemVolVars[scvIdx].porosity()+elemVolVars[scvIdx].buffer());
                    Scalar L = lambda/2*log(4*exp(-0.577215)*elemVolVars[scvIdx].effDiffCoeff()*pow(rootRadius_,(-2))*(this->timeManager().time()+this->timeManager().timeStepSize())+1);
                    AnalyticalC[dofGlobalIdx] += (Cinf_-Cinf_*lambda/(1+Cinf_dl+L+sqrt(4*Cinf_dl+pow((1-Cinf_dl+L),2)))*
                                                boost::math::expint(1,pow(x,2)/(4*elemVolVars[scvIdx].effDiffCoeff()*(this->timeManager().time()+this->timeManager().timeStepSize()))))
                                                /elemVolVars[scvIdx].density();

                    if (AnalyticalC[dofGlobalIdx] != 0)
                    {
                        relErrC=std::abs((elemVolVars[scvIdx].massFraction(transportCompIdx)-AnalyticalC[dofGlobalIdx])/AnalyticalC[dofGlobalIdx]);
                        absErrC=std::abs((elemVolVars[scvIdx].massFraction(transportCompIdx))/AnalyticalC[dofGlobalIdx]);
                    }

                    AnalyticalS[dofGlobalIdx] += 2*Vmax_*Cinf_dl/(1+Cinf_dl+L+sqrt(4*Cinf_dl+pow((1-Cinf_dl+L),2)));
                    if (((AnalyticalS[dofGlobalIdx]) < anaUptakeRate)or(anaUptakeRate==0.0))
                        anaUptakeRate = AnalyticalS[dofGlobalIdx];
                }
            }
        }

        std::cout << "Integrated mass source (3D): " << totalSourceP << std::endl;
        std::cout << "Integrated concentration source (3D): " << totalSourceC << std::endl;

        this->resultWriter().attachDofData(AnalyticalC, "AnalyticalC", isBox);
        this->resultWriter().attachDofData(Distance, "distance", isBox);
        this->resultWriter().attachDofData(relErrC, "relErrC", isBox);
        this->resultWriter().attachDofData(absErrC, "absErrC", isBox);

        uptakeC_num_ += -(totalSourceC/(rootRadius_*2*M_PI)/0.02)*(this->timeManager().timeStepSize());
        uptakeC_ana_ += anaUptakeRate*(this->timeManager().timeStepSize());

        logFile_.open(this->name() + ".log", std::ios::app);
        logFile_ << "time = " << this->timeManager().time()+this->timeManager().timeStepSize()
                << " uptakeRate_num = " << -totalSourceC/(rootRadius_*2*M_PI)/0.02
                << " uptakeRate_ana = "<< anaUptakeRate
                << " cummulativeUptake_num = "<< uptakeC_num_
                << " cummulativeUptake_ana = "<< uptakeC_ana_
                << " ratioRate = " << -totalSourceC/(rootRadius_*2*M_PI)/0.02/anaUptakeRate
                << " ratioUptake = " << uptakeC_num_/uptakeC_ana_
                << std::endl;
        logFile_.close();

        logFile_.open(this->name() + "_value.log", std::ios::app);
        logFile_ << this->timeManager().time()+this->timeManager().timeStepSize()
                << " " << -totalSourceC/(rootRadius_*2*M_PI)/0.02
                << " " << anaUptakeRate
                << " " << uptakeC_num_
                << " " << uptakeC_ana_
                << " " << -totalSourceC/(rootRadius_*2*M_PI)/0.02/anaUptakeRate
                << " " << uptakeC_num_/uptakeC_ana_
                << std::endl;
        logFile_.close();
    }


    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        this->couplingManager().resetCouplingSources();
        if (!(this->timeManager().time() < 0.0))
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
            logFile_.open(this->name() + "_post.log", std::ios::app);
            logFile_ << "time = " << this->timeManager().time()
                    << " uptakeRate_num = " << -totalSourceC
                    << std::endl;
            logFile_.close();

            logFile_.open(this->name() + "_valuepost.log", std::ios::app);
            logFile_ << this->timeManager().time()
                    << " " << -totalSourceC
                    << std::endl;
            logFile_.close();
        }
    }

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
        Scalar pw_ = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilPressure);
        priVars[pressureIdx] = pw_;

        priVars[massOrMoleFracIdx]= GET_RUNTIME_PARAM(TypeTag,
                                       Scalar,
                                       BoundaryConditions.InitialSoilFracK);

};

    std::ofstream logFile_;
    const Scalar eps_ = 1e-9;
    Scalar episodeTime, temperature_, pnRef_, pc_, sw_, pw_, uptakeC_num_, uptakeC_ana_;
    std::string name_;
    std::shared_ptr<CouplingManager> couplingManager_;
    //Scalar DiffCoef_;
};

} //end namespace

#endif
