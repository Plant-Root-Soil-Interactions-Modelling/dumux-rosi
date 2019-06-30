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
 * \brief A sub problem for the rootsystem
 */
#ifndef DUMUX_ROOTSYSTEM_1P2C_TEST_PROBLEM_HH
#define DUMUX_ROOTSYSTEM_1P2C_TEST_PROBLEM_HH

#include <cmath>

#include <dumux/porousmediumflow/rootmodel1p2c/model1p2c.hh>
#include <dumux/porousmediumflow/rootmodel1p2c/problem1p2c.hh>
#include <dumux/material/fluidsystems/h2osolute.hh>
//#include "h2oX.hh"
//! get the properties needed for subproblems
#include <dumux/multidimension/subproblemproperties.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include "1dRadialRichards2cCoupling.hh"

#include "rootsystemtestspatialparams.hh"
#include "dumux/io/csvReader.hh"
namespace Dumux
{
template <class TypeTag>
class RootsystemOnePTwoCTestProblem;

namespace Properties
{
NEW_TYPE_TAG(RootsystemOnePTwoCTestProblem, INHERITS_FROM(RootsystemOnePTwoC));
NEW_TYPE_TAG(RootsystemOnePTwoCTestBoxProblem, INHERITS_FROM(BoxModel, RootsystemOnePTwoCTestProblem));
NEW_TYPE_TAG(RootsystemOnePTwoCTestCCProblem, INHERITS_FROM(CCModel, RootsystemOnePTwoCTestProblem));

SET_TYPE_PROP(RootsystemOnePTwoCTestProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OSOLUTE<TypeTag>);
              //Dumux::H2OXFluidSystem<TypeTag>);

// Set the grid type
SET_TYPE_PROP(RootsystemOnePTwoCTestProblem, Grid, Dune::FoamGrid<1, 3>);

// Set the problem property
SET_TYPE_PROP(RootsystemOnePTwoCTestProblem, Problem, Dumux::RootsystemOnePTwoCTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(RootsystemOnePTwoCTestProblem, SpatialParams, Dumux::RootsystemTestSpatialParams<TypeTag>);

// Linear solver settings
#if HAVE_UMFPACK
SET_TYPE_PROP(RootsystemOnePTwoCTestProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif

// Enable gravity
SET_BOOL_PROP(RootsystemOnePTwoCTestProblem, ProblemEnableGravity, true);

// write newton convergence to vtk
SET_BOOL_PROP(RootsystemOnePTwoCTestProblem, NewtonWriteConvergence, false);

// Enable velocity output
SET_BOOL_PROP(RootsystemOnePTwoCTestProblem, VtkAddVelocity, true);

//SET_BOOL_PROP(RootsystemTestProblem, GrowingGrid, false);
//SET_BOOL_PROP(RootsystemOnePTwoCTestProblem, UseMoles, true);
SET_BOOL_PROP(RootsystemOnePTwoCTestProblem, UseMoles, false);
}

/*!
 * \ingroup OneDRootSystem
 * \ingroup ImplicitTestProblems
 * \brief TODO
 */
template <class TypeTag>
class RootsystemOnePTwoCTestProblem : public RootsystemOnePTwoCProblem<TypeTag>
{
    typedef RootsystemOnePTwoCProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
        transportEqIdx = Indices::transportEqIdx,
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;
    //typedef typename GET_PROP_TYPE(TypeTag, IntegrationPointSource) PointSource;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    // Coupling
    typedef typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag) GlobalProblemTypeTag;
    typedef typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager) CouplingManager;

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();


public:
    RootsystemOnePTwoCTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = tree.template get<std::string>("Problem.Name") + "-root";
        mimicGrowth_ = tree.template get<bool>("Grid.MimicGrowth", false);
        mimicGrowthEndTime_ = tree.template get<Scalar>("Grid.MimicGrowthEndTime", 1e99);

        this->spatialParams().setParams();

        ////read weather data
        useWeatherData_= tree.template get<bool>("Soil.BoundaryConditions.UseWeatherData", false);
        if (useWeatherData_)
        {
            weatherData_.getDataFromFile(tree.template get<std::string>("Soil.BoundaryConditions.WeatherDataFileName"));
            weatherData_.setTimeHeader("Days");
            //max_ = weatherData_.getMaxValueOverTime("NumberOfLeave");
        }
        //read whether the rhizosphere models are use
        useRhizoScaleModel_ = tree.template get<bool>("Rhizosphere.Grid.EnableRhizosphereScaleModel", false);
        reuseCouplingTerms_  = tree.template get<bool>("Rhizosphere.Grid.ReuseCouplingTerms_", true);
        if (useRhizoScaleModel_)
        {
            std::string command = "mkdir ";
            std::string path = tree.template get<std::string>("Problem.Name")+"_rhizo";
            const int dir = system(command.append(path).c_str());
            if (dir< 0)
            {
                 DUNE_THROW(Dune::IOError,
                    " Cannot create folders for the results from the rhiziosphere simulation !");
            }
        }

        neglectRootCollarSegment_ = tree.template get<bool>("Grid.NeglectRootCollarSegment", false);
        wiltingPointPressure_ = tree.template get<Scalar>("RootSystem.BoundaryConditions.CriticalCollarPressure", -1.4e6);
        collarPressure_ = tree.template get<Scalar>("RootSystem.BoundaryConditions.CollarPressure");
        dirichletRootCollar_ = tree.template get<bool>("RootSystem.BoundaryConditions.DirichletRootCollar", false);
        startEpisodeWithTimeStepSize_ = tree.template get<Scalar>("TimeManager.StartEpisodeWithTimeStepSize", 0);
        nutrientUptakeRate_ = 0.0;
        startEpisode_ = false;
        outputEveryTimeStep_ = tree.template get<bool>("TimeManager.OutputEveryTimeStep", false);
        transpirationRateRatio_ = tree.template get<Scalar>("RootSystem.BoundaryConditions.TranspirationRateRatio", 1);
        diurnalFluctuation_ = tree.template get<bool>("RootSystem.BoundaryConditions.DiurnalFluctuation", false);
        rootHydraulicRedistribution_ = tree.template get<bool>("RootSystem.BoundaryConditions.RootHydraulicRedistribution", true);
        transpirationRate_ = tree.template get<Scalar>("RootSystem.BoundaryConditions.TranspirationRate", 0);
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
    Scalar temperature() const
    {   Scalar rootTemp_;
        rootTemp_ = tree.template get<Scalar>("RootSystem.InitialConditions.RootTemperature");
        return rootTemp_;
    } // 10C
//    Scalar temperature() const
//    { temperatureRoot_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.TemperatureRoot);
//    return temperatureRoot_; } // 10C

    /*!
     * \brief Called by the time manager after the time integration.
     */
    void preTimeStep()
    {
        ParentType::preTimeStep();
        // If gridGrowth is used, this method grows the grid.
        // Remeber to call the parent class function if this is overwritten
        // on a lower problem level when using an adaptive grid
        //if (growingGrid && this->timeManager().timeStepIndex() > 0){
        //    this->gridGrowth().growGrid();
        //    this->resultWriter().gridChanged();
        //}
        preSol_ = this->model().curSol();

        if ((startEpisode_) and (startEpisodeWithTimeStepSize_ > 0))
        {
            this->timeManager().setTimeStepSize(startEpisodeWithTimeStepSize_);
            std::cout<<" ROOT set time step "<< this->timeManager().timeStepSize()<<"\n";
            startEpisode_ = false;
        }
        if (useRhizoScaleModel_)
        {
            if (this->timeManager().time()==0)
            {
                rhizoStencils_.clear();
	            for (auto&& stencil : this->couplingManager().lowDimCouplingStencils())
	            {
	                for (auto&& second : stencil.second)
	                {
	                    int rootEIdx = stencil.first;
	                    int soilEIdx = second;
	                    std::string rhizoOutputName = tree.template get<std::string>("Problem.Name")
	                                                +"_rhizo/r"+std::to_string(rootEIdx)+"_s"+std::to_string(soilEIdx);
	                    std::string command = "mkdir ";
	                    std::string path = rhizoOutputName;
	                    const int dir = system(command.append(path).c_str());
	                    if (dir< 0)
	                    {
	                         DUNE_THROW(Dune::IOError," Cannot create folders for the results from the rhiziosphere simulation !");
	                    }
                        Scalar rootAge = this->timeManager().time()-this->spatialParams().rootBornTime(rootEIdx);
                        if ((!mimicGrowth_) or ((rootAge >= 0) and (this->spatialParams().rootBornTime(rootEIdx) < this->mimicGrowthEndTime())))
                            rhizoStencils_.push_back({rootEIdx, soilEIdx});
	                }
	            }
            }
            else
            {
                if ((mimicGrowth_ ) and (this->timeManager().time() < mimicGrowthEndTime_))
                this->updaterhizoStencils();
            }

            this->couplingManager().resetCouplingSources();
            std::cout<<" updateCouplingTerms \n";
            //std::cin.ignore();
            for (int i = 0; i < rhizoStencils_.size(); i++)
            {
                std::vector<int> rhizoStencil = rhizoStencils_[i];
                int rootEIdx = rhizoStencil[0];
                //int soilEIdx = rhizoStencil[1];
                Dune::FieldVector<double, 2> TransportValue (0.0);
                Dune::FieldVector<double, 2> rootPriVars (0.0);
                rootPriVars[pressureIdx] = this->model().curSol()[rootEIdx][pressureIdx];//this->couplingManager().lowDimPriVars(rhizoStencil[0])[pressureIdx] ;
                //rootPriVars[massOrMoleFracIdx] = 0;
                rootPriVars[massOrMoleFracIdx] = this->couplingManager().lowDimPriVars(rhizoStencil[0])[massOrMoleFracIdx];
                //std::cout<<" root-soil "<<rootEIdx<<" - "<< soilEIdx<<" "<<this->couplingManager().lowDimPriVars(rootEIdx)[massOrMoleFracIdx]
                //<<" - "<<this->couplingManager().bulkPriVars(soilEIdx)[massOrMoleFracIdx]<<"\n";
                //std::cin.ignore();
                std::string rhizoOutputName = tree.template get<std::string>("Problem.Name")
                                            +"_rhizo/r"+std::to_string(rhizoStencil[0])+"_s"+std::to_string(rhizoStencil[1])+"/r";

                TransportValue = this->couplingManager().bulkProblem().runRhizosphereModelFromRoot(rhizoStencil[0], rhizoStencil[1],
                                                                                                        rootPriVars, rhizoOutputName);

                //std::pair<unsigned int, unsigned int> Intersection = std::make_pair(rootEIdx, soilEIdx);
                //TransportValue *= this->couplingManager().lowDimStencilLength(Intersection);

                this->couplingManager().setCouplingSources(std::make_pair(rhizoStencil[0], rhizoStencil[1]), TransportValue);
                this->couplingManager().setReuseCouplingSources(std::make_pair(rhizoStencil[0], rhizoStencil[1]));
            }
        }
    }

   void updaterhizoStencils() //called at postTimeStep to update new root segment
    {
        //const SpatialParams &spatialParams = this->spatialParams();
        rhizoStencils_.clear();
        for (auto&& stencil : this->couplingManager().lowDimCouplingStencils())
        {
            for (auto&& second : stencil.second)
            {
                int rootEIdx = stencil.first;
                int soilEIdx = second;
                Scalar rootAge = this->timeManager().time()-this->spatialParams().rootBornTime(rootEIdx);
                if ((!mimicGrowth_) or ((rootAge >= 0) and (this->spatialParams().rootBornTime(rootEIdx) < this->mimicGrowthEndTime())))
                    rhizoStencils_.push_back({rootEIdx, soilEIdx});
            }
        }
    }


    void runRhizo(std::vector<std::vector<int>> rhizoStencils, int begin, int end) const
    {

    }
    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    void boundaryTypesAtPos (BoundaryTypes &values,
                             const GlobalPosition &globalPos ) const
    {
        if (globalPos[2] + eps_ >  this->bBoxMax()[2] )
        {
            //values.setOutflow(transportEqIdx);
            //values.setNeumann(transportEqIdx);
            if ((tree.template get<bool>("RootSystem.BoundaryConditions.DirichletRootCollar")))
                values.setDirichlet(conti0EqIdx);
            else
            {
                //get element index Eid of root segment at root colar
                int Eid=-1;
                for (const auto& element : elements(this->gridView()))
                {
                    Eid ++;
                    auto posZ = std::max(element.geometry().corner(0)[2],element.geometry().corner(1)[2]);
                    if (posZ + eps_ > this->bBoxMax()[2])
                    {
                        break;
                    }

                }
                if (this->timeManager().time()>=0)
                {
                    if ((preSol_[Eid][conti0EqIdx] < wiltingPointPressure_ ))
                    {
                        std::cout<<"Collar pressure: "<<preSol_[Eid][conti0EqIdx]<<" < " <<wiltingPointPressure_<<"\n";
                        std::cout<<"WATER STRESS !! SET BC at collar as Dirichlet !!"<<"\n";
                        values.setDirichlet(conti0EqIdx);
                    }
                    else
                    {
                        //std::cout<<"Collar pressure: "<<preSol_[Eid][conti0EqIdx]<<" > " <<criticalCollarPressure<<"\n";
                        //std::cout<<"NO water stress !! SET BC at collar as Neumann !!"<<"\n";
                        values.setNeumann(conti0EqIdx);
                    }
                }
                else
                {
                    std::cout<<"SET BC at collar as Neumann !!"<<"\n";
                    values.setNeumann(conti0EqIdx);
                }
            }
        }
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */

//    void dirichlet(PrimaryVariables &values,
//                   const Intersection &intersection) const
//    {
//        GlobalPosition globalPos = intersection.geometry().center();
//        if (globalPos[2] + eps_ >  this->bBoxMax()[2] )
//              values[pIdx] = GET_RUNTIME_PARAM(TypeTag,
//                                               Scalar,
//                                               RootSystem.BoundaryConditions.CriticalCollarPressure);
//
//          //values[pIdx] *= 1.01;
//    }

    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {   if (dirichletRootCollar_)
            values[conti0EqIdx] = tree.template get<Scalar>("RootSystem.BoundaryConditions.CollarPressure");
        else
            values[conti0EqIdx] = wiltingPointPressure_;
    }
    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        Scalar TranspirationRate = transpirationRate_;
        GlobalPosition globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;
        if ( globalPos[2] + eps_ >  this->bBoxMax()[2] )
        {
            if (useWeatherData_)
            {
                Scalar simTime = this->timeManager().time();
                //TranspirationRate *= weatherData_.getValueAtTime("NumberOfLeave",simTime/(24*3600))/maxLeafNumber_;
                TranspirationRate = weatherData_.getValueAtTime("TranspirationRate_[kg/day]",simTime/(24*3600))/(24*3600)*transpirationRateRatio_;
                //std::cout<<" !! TranspirationRate "<<TranspirationRate<<"\n";
            }
            else
            {   bool isMultiRoot = tree.template get<bool>("RootSystem.BoundaryConditions.MultiRoot");
                if (isMultiRoot)
                {
                    const Scalar rootSurface = this->spatialParams().rootSurface(element, fvGeometry, scvIdx);
                    TranspirationRate *= rootSurface;
                }
            }
            if (diurnalFluctuation_)
                TranspirationRate *= (std::sin(this->timeManager().time()/(24*3600)*(2*M_PI)-M_PI/2)+1);
            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);
            values[conti0EqIdx] = TranspirationRate;
            values[transportEqIdx] = values[conti0EqIdx]*elemVolVars[scvIdx].massFraction(massOrMoleFracIdx);
        }
        else
            values = 0.0;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        GlobalPosition globalPos = element.geometry().center();
        Dune::FieldVector<double, 2> priVars_;
        this->couplingManager().bulkProblem().initialAtPos(priVars_, globalPos);
        priVars[pressureIdx] = priVars_[pressureIdx];
	//priVars[pressureIdx] =  tree.template get<Scalar>("BoundaryConditions.InitialRootPressure");
        priVars[massOrMoleFracIdx] = 0.0;
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
    { pointSources = this->couplingManager().lowDimPointSources(); }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param pointSources The vector of point sources
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void solDependentPointSource(PointSource& source,
                                  const Element &element,
                                  const FVElementGeometry &fvGeometry,
                                  const int scvIdx,
                                  const ElementVolumeVariables &elemVolVars) const
    {
        source = 0.0;
        Scalar nutrientUptakeRate = 0;
        unsigned int rootEIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        if ((rootEIdx==0) and neglectRootCollarSegment_)
            return;

        // compute source at every integration point
        const SpatialParams &spatialParams = this->spatialParams();
        Scalar rootAge = this->timeManager().time() - spatialParams.rootBornTime(rootEIdx);
        //if (((rootAge >= 0) and (rootEIdx != 0)) or (!mimicGrowth_))
        if ((rootAge >= 0) and (spatialParams.rootBornTime(rootEIdx) < mimicGrowthEndTime_))
        {
	        unsigned int soilEIdx = this->couplingManager().pointSourceData(source.id()).bulkElementIdx();
	        std::pair<unsigned int, unsigned int> Intersection = std::make_pair(rootEIdx, soilEIdx);

	        const Scalar Kr = spatialParams.Kr(element, fvGeometry, scvIdx);
	        const Scalar rootRadius = spatialParams.rootRadius(element, fvGeometry, scvIdx);
	        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[pressureIdx];
	        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[pressureIdx];

			PrimaryVariables sourceValues;
			sourceValues = 0.0;
			// sink defined as radial flow Jr [m^3 s-1]*density
            if (((pressure1D < pressure3D)  and  (pressure3D > wiltingPointPressure_))
                     or ((pressure1D > pressure3D) and (rootHydraulicRedistribution_)))
                sourceValues[conti0EqIdx] = 2 * M_PI *rootRadius * Kr *(pressure3D - pressure1D)
			                               * elemVolVars[scvIdx].density();
            if (!useRhizoScaleModel_)
            {
                Scalar c3D;
                c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];//*1000;//elemVolVars[0].density();
                ////Difussive flux term of transport
                const Scalar DiffValue = 0.0;
                ////2* M_PI *rootRadius *DiffCoef_*(c1D - c3D)*elemVolVars[scvIdx].density(/*phaseIdx=*/0);
                ////Advective flux term of transport
                Scalar AdvValue;
                if (sourceValues[conti0EqIdx]>0)
                    AdvValue = sourceValues[conti0EqIdx]*c3D;
                else
                    AdvValue = 0;
                //    AdvValue = sourceValues[conti0EqIdx]*c1D;
                //Active flux - active uptake based on Michaeles Menten
                Scalar ActiveValue = 0.0;

                const Scalar Vmax = spatialParams.Vmax();
                const Scalar Km = spatialParams.Km();
                const Scalar Cmin = spatialParams.Cmin();
                if (c3D*elemVolVars[scvIdx].density()-Cmin >= 0)
                    ActiveValue = 2 * M_PI*rootRadius*Vmax*(c3D*elemVolVars[scvIdx].density()-Cmin)/(Km+c3D*elemVolVars[scvIdx].density()-Cmin);
                Scalar sigma;
                //sigma = spatialParams.PartitionCoeff();
                sigma = 0; //active uptake
                sourceValues[transportEqIdx] = (sigma*(AdvValue + DiffValue) + (1-sigma)*ActiveValue);
                //setNutrientUptakeRate(nutrientUptakeRate_, (sigma*(AdvValue + DiffValue) + (1-sigma)*ActiveValue));
                nutrientUptakeRate = (sigma*(AdvValue + DiffValue) + (1-sigma)*ActiveValue);
            }
            else
            {
    	        // REUSE approach
    	        if (this->couplingManager().reuseCouplingSources(Intersection))
    	        {
    	            //sourceValues[conti0EqIdx] = this->couplingManager().couplingSources(Intersection)[conti0EqIdx];
    	            sourceValues[transportEqIdx] = this->couplingManager().couplingSources(Intersection)[transportEqIdx];
                    nutrientUptakeRate = this->couplingManager().couplingSources(Intersection)[1];
    	        }
    	        else
    	        {
    	            DUNE_THROW(Dune::NotImplemented,"  The coupling sources need to be evaluated in soil problem !");

    	            Dune::FieldVector<double, 2> rootPriVars;
    	            rootPriVars[pressureIdx] = this->couplingManager().lowDimPriVars(source.id())[pressureIdx] ;
    	            //rootPriVars[massOrMoleFracIdx] = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];

                        std::string rhizoOutputName = tree.template get<std::string>("Problem.Name")
                                                    +"_rhizo/r"+std::to_string(rootEIdx)+"_s"+std::to_string(soilEIdx)+"/r";

	            Dune::FieldVector<double, 2> TransportValue = this->couplingManager().bulkProblem().runRhizosphereModelFromRoot(rootEIdx, soilEIdx, rootPriVars, rhizoOutputName);
                    nutrientUptakeRate = TransportValue[1];
	            TransportValue *= source.quadratureWeight()*source.integrationElement();
	            sourceValues[transportEqIdx] = TransportValue[transportEqIdx];
    	            this->couplingManager().setCouplingSources(Intersection, TransportValue);
    	            this->couplingManager().setReuseCouplingSources(Intersection);
    	        }
            }
        	sourceValues *= source.quadratureWeight()*source.integrationElement();
			nutrientUptakeRate *= source.quadratureWeight()*source.integrationElement();
	        source = sourceValues;
			nutrientUptakeRate_ += nutrientUptakeRate;
        }
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The source and sink values for the conservation equations in
     units of \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void solDependentSource(PrimaryVariables &values,
                     const Element &element,
                     const FVElementGeometry &fvGeometry,
                     const int scvIdx,
                     const ElementVolumeVariables &elemVolVars) const
    {
        values = 0.0;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    // add source terms to the output
    void addOutputVtkFields()
    {
        const SpatialParams &spatialParams = this->spatialParams();

        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        unsigned numDofs = this->model().numDofs();

        // create required scalar fields for the vtk output
        ScalarField& sourceP = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceCFlux = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourcePFlux = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& rootAge = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& distanceFromOrigin = *(this->resultWriter().allocateManagedBuffer(numDofs));
        sourceP = 0.0;
        sourceC = 0.0;
        sourceCFlux = 0.0;
        sourcePFlux = 0.0;

        Scalar totalSourceC, totalSourceP;
        totalSourceC = 0.0;
        totalSourceP = 0.0;

        // iterate over all elements
        for (const auto& element : elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

            int branchId = spatialParams.rootBranch(element, fvGeometry, 0);

            PrimaryVariables Sources;
            // output pressure
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                auto dofGlobalIdx = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);

                // only call the source function if we are not initializing
                if (!(this->timeManager().time() < 0.0))
                {
                    nutrientUptakeRate_ = 0;
                    this->scvPointSources(Sources, element, fvGeometry, scvIdx, elemVolVars);
                    Sources *= fvGeometry.subContVol[scvIdx].volume
                                             * this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                    sourceP[dofGlobalIdx] += Sources[0];
                    totalSourceP += sourceP[dofGlobalIdx];

                    sourceC[dofGlobalIdx] += nutrientUptakeRate_;// * fvGeometry.subContVol[scvIdx].volume
                                             //* this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                    totalSourceC += sourceC[dofGlobalIdx];

                    sourcePFlux[dofGlobalIdx] += sourceP[dofGlobalIdx]/fvGeometry.subContVol[scvIdx].volume;
                    sourceCFlux[dofGlobalIdx] += sourceC[dofGlobalIdx]/fvGeometry.subContVol[scvIdx].volume;
                    rootAge[dofGlobalIdx] = this->timeManager().time()-spatialParams.rootBornTime(element, fvGeometry, scvIdx);
                    distanceFromOrigin[dofGlobalIdx] = spatialParams.distanceFromOrigin(element, fvGeometry, scvIdx);
                }
            }
            //int eIdx = this->elementMapper().index(element);
            //elementIndicesInBranch_[branchId].push_back(eIdx);
        }

        std::cout << "Integrated mass source (1D): " << totalSourceP << std::endl;
        std::cout << "Integrated concentration source (1D): " << totalSourceC << std::endl;
        // attach data to the vtk output
        this->resultWriter().attachDofData(sourceP, "waterSource_[kg/s)]", isBox);
        this->resultWriter().attachDofData(sourceC, "nutrientSource_[kg/s]", isBox);
        this->resultWriter().attachDofData(sourcePFlux, "waterSourceFlux_[kg/s/m]", isBox);
        this->resultWriter().attachDofData(sourceCFlux, "nutrientSourceFlux_[kg/s/m]", isBox);
        this->resultWriter().attachDofData(rootAge, "rootAge_[s]", isBox);
        this->resultWriter().attachDofData(distanceFromOrigin, "distanceFromOrigin_[m]", isBox);
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        if (this->timeManager().timeStepIndex()==1)
        {
            logFile_.open(this->name() + "_pointsource.log", std::ofstream::out | std::ofstream::trunc);
            logFile_<< "rootEIdx"<<"\t" <<"soilEIdx"<<"\t"<<"Position"<<"\t"<<"Length"<<"\t"<<"rootRadius"<<"\t"<<"Vector"
                    <<"\t"<<"rootOrder"<<"\t"<<"rootBranch"<<"\t"<<"rootBornTime"<<"\t"<<"distanceFromOrigin"<<"\n";
            for (auto& pointSource : this->couplingManager().lowDimPointSources())
                {
                    unsigned int rootEIdx = this->couplingManager().pointSourceData(pointSource.id()).lowDimElementIdx();
                    unsigned int soilEIdx = this->couplingManager().pointSourceData(pointSource.id()).bulkElementIdx();
                    auto radius = this->spatialParams().radius(rootEIdx);
                    auto rootOrder = this->spatialParams().rootOrder(rootEIdx);
                    auto rootBranch = this->spatialParams().rootBranch(rootEIdx);
                    auto rootBornTime = this->spatialParams().rootBornTime(rootEIdx);
                    auto distanceFromOrigin = this->spatialParams().distanceFromOrigin(rootEIdx);
                    int i=0;
                    for (const auto& element : elements(this->gridView()))
                    {
                        if (i==rootEIdx)
                        {
                            auto rootGeometry = element.geometry();
                            auto vector= rootGeometry.corner(1)-rootGeometry.corner(0);
                            logFile_<< rootEIdx << "\t" <<  soilEIdx << "\t" <<  pointSource.position() << "\t"
                                    << pointSource.quadratureWeight()*pointSource.integrationElement() << "\t" << radius
                                    << "\t"<<vector<< "\t"<<rootOrder<< "\t"<<rootBranch<< "\t"<<rootBornTime<<"\t"<<distanceFromOrigin<<'\n';
                            break;
                        }
                        i++;
                    }
                }
            logFile_.close();

            std::string OutputFileName = this->name() +"_rootOrder.log";
            if (std::ifstream(OutputFileName))
                std::remove(OutputFileName.c_str());
            logFile_.open(OutputFileName, std::ios::app);
            for (auto&& order : this->spatialParams().orderOfBranches())
            {
                for (auto&& branchId : order.second)
                    logFile_ <<branchId<<"\t";
                logFile_<< std::endl;
            }
            logFile_.close();
        }
        std::string OutputFileName = this->name() +"_collarPressure.log";
        logFile_.open(OutputFileName, std::ios::app);
        if (std::ifstream(OutputFileName))
            std::remove(OutputFileName.c_str());
        logFile_.open(OutputFileName, std::ios::app);
        logFile_ << this->timeManager().time() + this->timeManager().timeStepSize()
            << " " << this->model().curSol()[0][conti0EqIdx]
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
    {
        if (this->timeManager().time()<0) return true;
        if (outputEveryTimeStep_)
            return true;
        else
            return (this->timeManager().episodeWillBeFinished() || this->timeManager().willBeFinished());
    }

    void episodeEnd()
    {
        if (startEpisodeWithTimeStepSize_ > 0)
        {
            //this->timeManager().setTimeStepSize(startEpisodeWithTimeStepSize_);
            std::cout<<" ROOT set next time step "<< this->timeManager().timeStepSize()<<"\n";
            startEpisode_ = true;
        }

    }
    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    Scalar mimicGrowthEndTime() const
    { return mimicGrowthEndTime_;}
private:
    std::string name_;
    std::ofstream logFile_;
    Scalar wiltingPointPressure_, collarPressure_, startEpisodeWithTimeStepSize_, mimicGrowthEndTime_, transpirationRate_;
    mutable Scalar nutrientUptakeRate_;
    const Scalar eps_ = 1e-9;
    std::shared_ptr<CouplingManager> couplingManager_;
    SolutionVector preSol_;
    CSVReader weatherData_;
    bool useWeatherData_, mimicGrowth_, useRhizoScaleModel_, neglectRootCollarSegment_, reuseCouplingTerms_, dirichletRootCollar_,
            startEpisode_, outputEveryTimeStep_, diurnalFluctuation_,
            rootHydraulicRedistribution_;
    Scalar maxLeafNumber_;
    std::vector<std::vector<int>> rhizoStencils_;
    //std::map<unsigned int, std::vector<unsigned int>> elementIndicesInBranch_;
    Scalar transpirationRateRatio_;
};

} //end namespace

#endif
