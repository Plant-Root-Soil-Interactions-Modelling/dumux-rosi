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

#include <dumux/material/fluidsystems/h2osolute.hh>
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
              Dumux::FluidSystems::H2OSOLUTE<TypeTag>);

// Set the grid type
//SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::YaspGrid<3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::UGGrid<3>);

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
SET_STRING_PROP(SoilRichardsTwoCBufferTestProblem, GridParameterGroup, "Soil.Grid");

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

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

public:
    SoilRichardsTwoCBufferTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();
        name_ = tree.template get<std::string>("Problem.Name") + "-soil";
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
        uptakeC_num_ = 0.0;
        uptakeC_ana_ = 0.0;
        if (std::ifstream(nameUptakeRate_))
            std::remove(nameUptakeRate_.c_str());

        ////read weather data
        useWeatherData_ = tree.template get<bool>("Soil.BoundaryConditions.UseWeatherData", false);
        if (useWeatherData_)
        {
            weatherData_.getDataFromFile(tree.template get<std::string>("Soil.BoundaryConditions.WeatherDataFileName"));
            weatherData_.setTimeHeader("Days");
        }
        useSoilData_ = tree.template get<bool>("Soil.SpatialParams.UseSoilData", false);
        mimicGrowth_ = tree.template get<bool>("Grid.MimicGrowth", false);
        mimicGrowthEndTime_ = tree.template get<Scalar>("Grid.MimicGrowthEndTime", 1e99);
        useRhizoScaleModel_ = tree.template get<bool>("Rhizosphere.Grid.EnableRhizosphereScaleModel", false);
        reuseCouplingTerms_  = tree.template get<bool>("Rhizosphere.Grid.ReuseCouplingTerms_", true);
        neglectRootCollarSegment_ = tree.template get<bool>("Grid.NeglectRootCollarSegment", false);
        radialTimeStepRatio_ =tree.template get<Scalar>("radialTimeManager.TimeStepRatio", 1);
        wiltingPointPressure_ = tree.template get<Scalar>("RootSystem.BoundaryConditions.CriticalCollarPressure", -1.4e6);
        periodicVerticalFlow_  = tree.template get<bool>("Soil.BoundaryConditions.PeriodicVerticalFlow", false);
        bottomInFlow_  = tree.template get<bool>("Soil.BoundaryConditions.BottomInFlow", true);
        bottomFreeDrainageFlow_  = tree.template get<bool>("Soil.BoundaryConditions.BottomFreeDrainageFlow", false);
   		//startEpisodeWithInitialTimeStep_ =tree.template get<Scalar>("TimeManager.StartEpisodeWithTimeStepSize", 0);
   		startEpisodeWithTimeStepSize_ = tree.template get<Scalar>("TimeManager.StartEpisodeWithTimeStepSize", 0);
   		startEpisode_ = false;
        outputEveryTimeStep_ = tree.template get<bool>("TimeManager.OutputEveryTimeStep",false);
        diurnalFluctuation_ = tree.template get<bool>("RootSystem.BoundaryConditions.DiurnalFluctuation", false);
        rootHydraulicRedistribution_ = tree.template get<bool>("RootSystem.BoundaryConditions.RootHydraulicRedistribution", true);
        topSoilSaturated_ = false;
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
        soilTemp_ = tree.template get<Scalar>("Soil.InitialConditions.SoilTemperature");
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
        source = 0.0;
        unsigned int rootEIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        if ((rootEIdx==0) and neglectRootCollarSegment_)
            return;
        // compute source at every integration point
        const auto& spatialParams = this->couplingManager().lowDimProblem().spatialParams();
        Scalar rootAge = this->timeManager().time() - spatialParams.rootBornTime(rootEIdx);
        //if (((rootAge >= 0) and (rootEIdx != 0)) or (!mimicGrowth_))
        if ((rootAge >= 0) and (spatialParams.rootBornTime(rootEIdx) < mimicGrowthEndTime_))
        {
            unsigned int soilEIdx = this->couplingManager().pointSourceData(source.id()).bulkElementIdx();
            std::pair<unsigned int, unsigned int> Intersection = std::make_pair(rootEIdx, soilEIdx);
            const Scalar Kr = spatialParams.Kr(rootEIdx);
            const Scalar rootRadius = spatialParams.radius(rootEIdx);
            const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[conti0EqIdx];
            const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[conti0EqIdx] ;

            PrimaryVariables sourceValues;
            sourceValues = 0.0;

            // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3] = [m-1 kg s-1]
            if (((pressure1D < pressure3D)  and  (pressure3D > wiltingPointPressure_))
                     or ((pressure1D > pressure3D) and (rootHydraulicRedistribution_)))
                sourceValues[conti0EqIdx] = 2* M_PI *rootRadius * Kr *(pressure1D - pressure3D)
                                       *elemVolVars[scvIdx].density();
            if (!useRhizoScaleModel_)
            {
                //Scalar c1D;
                //if(useMoles)
                //    c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];
                //else
                //    c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];

                Scalar c3D;
                if(useMoles)
                    c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];
                else
                    c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];

                const Scalar DiffValue = 0.0;
                //2* M_PI *rootRadius *DiffCoef_*(c1D - c3D)*elemVolVars[scvIdx].density(/*phaseIdx=*/0);

                //Advective flux term of transport
                Scalar AdvValue;
                //if (sourceValues[conti0EqIdx]>0) // flow from root to soil
                //    //AdvValue = 0;
                //    AdvValue = sourceValues[conti0EqIdx]*c1D; // [m-1 kg s-1]* [kg m-3] =[m-4 kg2 s-1]
                //else // flow from soil to root
                    AdvValue = sourceValues[conti0EqIdx]*c3D;

                //Active flux - active uptake based on Michaeles Menten
                Scalar ActiveValue;
                const Scalar Vmax = spatialParams.Vmax();
                const Scalar Km = spatialParams.Km();
                const Scalar Cmin = spatialParams.Cmin();
                ActiveValue = 0.;
                if (c3D*elemVolVars[scvIdx].density()-Cmin >= 0)
                    ActiveValue = -2*M_PI*rootRadius*Vmax*(c3D*elemVolVars[scvIdx].density()-Cmin)
                            /(Km+c3D*elemVolVars[scvIdx].density()-Cmin); //[m][kg m-2 s-1] [kg m-3] / [kg m-3] =[kg m-1 s-1]
                const Scalar sigma = 0;
                sourceValues[transportEqIdx] = (sigma*(AdvValue + DiffValue) + (1-sigma)*ActiveValue);
            }
            else
            {
                if (this->couplingManager().reuseCouplingSources(Intersection))
                {
                    //sourceValues[conti0EqIdx] = -this->couplingManager().couplingSources(Intersection)[conti0EqIdx];
                    sourceValues[transportEqIdx] = -this->couplingManager().couplingSources(Intersection)[transportEqIdx];
                }
                else
                {
                    //std::cout<<rootEIdx<<" "<< soilEIdx;
                    DUNE_THROW(Dune::NotImplemented,"  The coupling sources need to be evaluated in root problem !");

                    PrimaryVariables rootPriVars;
                    rootPriVars[pressureIdx] = this->couplingManager().lowDimPriVars(source.id())[pressureIdx] ;
                    rootPriVars[massOrMoleFracIdx] = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];

                    std::string rhizoOutputName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name)
                                                +"_rhizo_s"+std::to_string(soilEIdx)+"_r"+std::to_string(rootEIdx);

                    PrimaryVariables TransportValue = runRhizosphereModel(rootEIdx, soilEIdx, rootPriVars, element, fvGeometry, elemVolVars, rhizoOutputName);

                    sourceValues[transportEqIdx] = TransportValue[transportEqIdx];
                    //sourceValues[conti0EqIdx] = TransportValue[conti0EqIdx];

                    this->couplingManager().setCouplingSources(Intersection, TransportValue);
                    this->couplingManager().setReuseCouplingSources(Intersection);
                }
            }
    		sourceValues *= source.quadratureWeight()*source.integrationElement();
        	source = sourceValues;
	   }
    }

    PrimaryVariables runRhizosphereModel(const int rootEIdx, const int soilEIdx, const PrimaryVariables rootPriVars,
                                        const Element element, const FVElementGeometry fvGeometry, const ElementVolumeVariables &elemVolVars,
                                        const std::string outputName) const
    {
            std::pair<unsigned int, unsigned int> Intersection = std::make_pair(rootEIdx, soilEIdx);
            //std::cout<<"!!!! !!!! rootEIdx "<<rootEIdx<<" soilEIdx "<<soilEIdx<<"\n";
            // GETTING PARAMETERTERS FOR RHIZOSPHERE MODEL
            const auto& spatialParams = this->couplingManager().lowDimProblem().spatialParams();
            //const Scalar Kr = spatialParams.Kr(rootEIdx);
            const Scalar rootRadius = spatialParams.radius(rootEIdx);
            const Scalar Vmax = spatialParams.Vmax();
            const Scalar Km = spatialParams.Km();
            const Scalar Cmin = spatialParams.Cmin();
            const Scalar sigma = spatialParams.PartitionCoeff();
            auto restartSols = this->couplingManager().rhizoSolutions(Intersection);

            bool restart = false;
            if (restartSols.size() != 0)
                restart = true;

            Scalar restartTimeRhizosphere = this->timeManager().time();
            Scalar endTimeRhizosphere = this->timeManager().time() + this->timeManager().timeStepSize();;

            Scalar dtRhizosphere = this->timeManager().timeStepSize()
                                    /radialTimeStepRatio_;
            //if (this->timeManager().willBeFinished())
            //    dtRhizosphere = this->timeManager().timeStepSize();

            Scalar rootIntersectionLength = this->couplingManager().lowDimStencilLength(Intersection);
            //Scalar rhizosphereVolume = rootIntersectionLength*(M_PI*rootRadius*rootRadius)
            //                            /this->couplingManager().lowDimVolume(element)
            //                            *element.geometry().volume();
            //Scalar rhizosphereRadius = pow(rhizosphereVolume/rootIntersectionLength/M_PI,0.5)+rootRadius;

            Scalar rhizosphereRadius = this->couplingManager().rhizosphereRadius(rootEIdx, soilEIdx);

            //std::cout<<"!!!! !!!! rootEIdx "<<rootEIdx<<" rootRadius "<<rootRadius<<" rhizosphereRadius "<<rhizosphereRadius<<"\n";
            // Calculation of coupled fluxes for Boundary Conditions
            PrimaryVariables coupledFluxes;
            coupledFluxes = 0;
            //Scalar diffFluxTotal = 0;
            //Scalar avdFluxTotal = 0;
            FluxVariables fluxVars;

            for (int scvfIdx = 0; scvfIdx < fvGeometry.numScvf; scvfIdx++)
            {
                fluxVars.update(*this, element, fvGeometry, scvfIdx, elemVolVars, false /* oldSol? */);

                //const ElementVolumeVariables &up = this->curVolVars(fluxVars.upstreamIdx(phaseIdx));
                auto elementUp = this->boundingBoxTree().entity(fluxVars.upstreamIdx(0));
                FVElementGeometry fvGeometryUp;
                fvGeometryUp.update(this->gridView(), elementUp);

                ElementVolumeVariables up;
                up.update(*this, elementUp, fvGeometryUp, false /* oldSol? */);

                auto elementDn = this->boundingBoxTree().entity(fluxVars.downstreamIdx(0));
                FVElementGeometry fvGeometryDn;
                fvGeometryDn.update(this->gridView(), elementDn);

                ElementVolumeVariables dn;
                dn.update(*this, elementDn, fvGeometryDn, false /* oldSol? */);

                Scalar massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);

                Scalar advFlux = fluxVars.volumeFlux(0);
                advFlux *=
                            ((     massUpwindWeight_)*up[0].density()
                             +
                             ((1 - massUpwindWeight_)*dn[0].density()));
                coupledFluxes[conti0EqIdx] += advFlux;

                coupledFluxes[transportEqIdx] +=
                        fluxVars.volumeFlux(0) *
                        ((    massUpwindWeight_)*up[0].density() * up[0].massFraction(transportCompIdx)
                         +
                         (1 - massUpwindWeight_)*dn[0].density()*dn[0].massFraction(transportCompIdx));

                Scalar tmp = -(fluxVars.massFractionGrad(transportCompIdx)*fluxVars.face().normal);
                tmp *= fluxVars.porousDiffCoeff() * fluxVars.density();
                coupledFluxes[transportEqIdx] += tmp;
            }
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (!intersection.boundary())
                        continue;
                //GlobalPosition globalPos = intersection.geometry().center();
                unsigned bfIdx = intersection.indexInInside();

                PrimaryVariables BCFlux(0.0);
                this->neumann(BCFlux, element,fvGeometry, intersection, 0, bfIdx);
                BCFlux *= intersection.geometry().volume();

                coupledFluxes -= BCFlux;
            }

            Scalar surfaceWeight = (2*M_PI*rootRadius*rootIntersectionLength)/this->couplingManager().lowDimSurface(element);

            coupledFluxes *= surfaceWeight/(2*M_PI*rhizosphereRadius*rootIntersectionLength);

            // SETTING RHIZOSPHERE MODEL
            OneDRichards2CRadiallySymmetric<TypeTag> OnedRadialSol;

            OnedRadialSol.setName(tree.template get<std::string>("ParameterFile"));
            OnedRadialSol.setProblemName(outputName);

            // set ROOT spatial params
            OnedRadialSol.setVmax(Vmax);
            OnedRadialSol.setKm(Km);
            OnedRadialSol.setCmin(Cmin);
            OnedRadialSol.setRootRadius(rootRadius);
            OnedRadialSol.setKr(spatialParams.Kr(rootEIdx));
            OnedRadialSol.setPartitionCoeff(sigma);
            OnedRadialSol.setDiffCoeffRootMembrane(tree.template get<Scalar>("RootSystem.SpatialParams.DiffCoeffRootMembrane"));

            // set SOIL spatial params
            OnedRadialSol.setEffDiffCoeff(elemVolVars[0].effDiffCoeff());
            OnedRadialSol.setSaturation(elemVolVars[0].saturation(phaseIdx));
            OnedRadialSol.setPorosity(elemVolVars[0].porosity());
            OnedRadialSol.setDensity(elemVolVars[0].density());
            OnedRadialSol.setBuffer(elemVolVars[0].buffer());
            OnedRadialSol.setFreundlichK(this->spatialParams().FreundlichK(element, fvGeometry, 0));
            OnedRadialSol.setFreundlichN(this->spatialParams().FreundlichN(element, fvGeometry, 0));
            OnedRadialSol.setBulkDensity(this->spatialParams().bulkDensity(element, fvGeometry, 0));
            OnedRadialSol.setPermeability(this->spatialParams().intrinsicPermeability(element,fvGeometry,0));
            OnedRadialSol.setVgAlpha(this->spatialParams().VgAlpha(element,fvGeometry,0));
            OnedRadialSol.setVgn(this->spatialParams().Vgn(element,fvGeometry,0));
            OnedRadialSol.setSwr(this->spatialParams().Swr(element,fvGeometry,0));

            //set time and grid
            OnedRadialSol.setRestart(restart);
            OnedRadialSol.setRestartSolutions(restartSols);
            OnedRadialSol.setRestartTime(restartTimeRhizosphere);
            OnedRadialSol.setTEnd(endTimeRhizosphere);
            OnedRadialSol.setDtInitial(dtRhizosphere);
            OnedRadialSol.setOuterBoundaryGrid(rhizosphereRadius);
            OnedRadialSol.setFinalRhizoRadius(this->couplingManager().finalRhizoRadius(Intersection));
            OnedRadialSol.setZ(element.geometry().center()[2]);

            //set initial conditions
            //OnedRadialSol.setInitialSoilPressure(GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.InitialSoilPressure));
            OnedRadialSol.setInitialSoilPressure(elemVolVars[0].pressure());
           //Scalar initialTotalMassConcentrationInRhizosphere =
           //        (elemVolVars[0].totalMassConcentration()*element.geometry().volume()
           //                            - this->couplingManager().rhizoMassInBulkElement(soilEIdx))/
           //        (element.geometry().volume() - this->couplingManager().lowDimVolume(element));

            //OnedRadialSol.setInitialSoilNutrient(elemVolVars[0].dissolveMassFraction(initialTotalMassConcentrationInRhizosphere));
            OnedRadialSol.setInitialSoilNutrient(initialMassFractionInRhizospheres_[soilEIdx]);
            //OnedRadialSol.setInitialSoilNutrient(elemVolVars[0].massFraction(transportCompIdx));
            //set boundary conditions
            OnedRadialSol.setOuterFluxBoundaryConditions(coupledFluxes);
            OnedRadialSol.setInsideRootConditions(rootPriVars);

            std::vector<PrimaryVariables> OneDSolutions = OnedRadialSol.ActiveUptakeNumerical();

            this->couplingManager().setRhizoSolutions(Intersection, OneDSolutions);
            //std::cout<<"setRhizoSolutions uptake: \n";
            PrimaryVariables uptake;
            uptake[0] = OneDSolutions.back()[0];
            uptake[1] = OneDSolutions.back()[1];

            //std::cout<<"SOIL uptake: "<<uptake[0]<<" "<<uptake[1]<<"\n";
            uptake *= -2*M_PI*rootRadius;

            return uptake;
    }

    PrimaryVariables runRhizosphereModelFromRoot(const int rootEIdx, const int soilEIdx, const PrimaryVariables rootPriVars, const std::string outputName) const
    {
            //std::cout<< "in 3D pressureRoot_ : !!!!! "<< rootPriVars[conti0EqIdx]<< "\n";
            //std::pair<unsigned int, unsigned int> Intersection = std::make_pair(rootEIdx, soilEIdx);
            auto element = this->boundingBoxTree().entity(soilEIdx);
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

            PrimaryVariables uptake = runRhizosphereModel(rootEIdx, soilEIdx, rootPriVars, element, fvGeometry, elemVolVars, outputName);
            uptake *= -1;
            return uptake;
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
        bool ZeroNeumann = tree.template get<bool>("Soil.BoundaryConditions.ZeroNeumannSoil");
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
        if ((useWeatherData_))
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
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        values = 0;
        bool ZeroNeumann = tree.template get<bool>("Soil.BoundaryConditions.ZeroNeumannSoil");

        GlobalPosition globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        if (ZeroNeumann)
        {
            values[conti0EqIdx] = 0.0;
            values[transportEqIdx] = 0.0;
        }

        if (useWeatherData_)
        {
            // set precipitation at top BC
            if(globalPos[2] > this->bBoxMax()[2] - eps_)
            {
                if (periodicVerticalFlow_)
                    values[conti0EqIdx] = -mFluxOut_[0];
                else
                {
                    Scalar simTime = this->timeManager().time();
                    Scalar topInRate = weatherData_.getValueAtTime("TopIn_[kg/m2/day]",simTime/(24*60*60))/(24*60*60);
                    topInRate -=weatherData_.getValueAtTime("Evaporation_[kg/m2/day]",simTime/(24*60*60))/(24*60*60);
                    //precipitationValue /= (24*60*60); //mm/day -> kg/(s m2)
                    //values[conti0EqIdx] = -topInRate; //inflow is negative
                    if (topSoilSaturated_ and (topInRate > mFluxOut_[0]))
                    {
                        values[conti0EqIdx] = -mFluxOut_[0];
                    }
                    else
                        values[conti0EqIdx] = -topInRate;
                    //std::cout<<topSoilSaturated_<<" topInRate "<<topInRate<<" :-mFluxOut_[0] "<<-mFluxOut_[0]<<" values[conti0EqIdx] "<<values[conti0EqIdx]<<"\n";
                }
                //std::cout<<"set precipitation at top BC"<<values[conti0EqIdx]<<"\n";
                values[transportEqIdx] = 0.0;
            }
            //// set ZeroNeumann at the sides
            //if((globalPos[1] > this->bBoxMax()[1] - eps_)||(globalPos[1] < this->bBoxMin()[1] + eps_)||
            //  (globalPos[0] > this->bBoxMax()[0] - eps_)||(globalPos[0] < this->bBoxMin()[0] + eps_))
            //{
            //    values[conti0EqIdx] = 0.0;
            //    //std::cout<<"set  ZeroNeumann at the sides"<<values[conti0EqIdx]<<"\n";
            //    values[transportEqIdx] = 0.0;
            //}
            // set free drainage at bottom BC
            if(globalPos[2] < this->bBoxMin()[2] + eps_)
            {
                if (bottomFreeDrainageFlow_)
                {
                    ElementVolumeVariables elemVolVars;
                    elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

                    Scalar intrinsicPermeability = this->spatialParams().intrinsicPermeability(element, fvGeometry, 0);
                    Scalar mobility = elemVolVars[scvIdx].mobility(phaseIdx);

                    values[conti0EqIdx] = mobility*intrinsicPermeability*1000; // * 1 [m]
                    //std::cout<<"set free drainage at bottom BC "<<values[conti0EqIdx]<<"\n";
                    //values[conti0EqIdx] = -precipitationValue;
                    values[transportEqIdx] = values[conti0EqIdx]*elemVolVars[scvIdx].massFraction(massOrMoleFracIdx);
                }
                else if (bottomInFlow_)
                {
                    Scalar simTime = this->timeManager().time();
                    Scalar bottomInRate = weatherData_.getValueAtTime("BottomIn_[kg/m2/day]",simTime/(24*60*60))/(24*60*60);
                    //precipitationValue /= (24*60*60); //mm/day -> kg/(s m2)
                    values[conti0EqIdx] = -bottomInRate;
                }
            }

            if ((diurnalFluctuation_) and (!bottomFreeDrainageFlow_))
                values[conti0EqIdx] *= (std::sin(this->timeManager().time()/(24*3600)*(2*M_PI)-M_PI/2)+1);
        }
    }
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
        //ScalarField& ratioPressure = *(this->resultWriter().allocateManagedBuffer(numDofs));
        //ScalarField& ratioWaterContent = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ratioMassFraction = *(this->resultWriter().allocateManagedBuffer(numDofs));
        //ScalarField& diffWaterContent = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& diffMassCon = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceP = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& rootSurface = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& rootLength = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& rootLengthDensity = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& volElement =*(this->resultWriter().allocateManagedBuffer(numDofs));
        //ScalarField& CoNumber = *(this->resultWriter().allocateManagedBuffer(numDofs));

        //ratioPressure = 0.0;
        //ratioWaterContent = 0.0;
        ratioMassFraction = 0.0;
        //diffWaterContent = 0.0;
        diffMassCon = 0.0;
        sourceP = 0.0;
        sourceC = 0.0;
        rootSurface = 0.0;
        rootLength = 0.0;
        rootLengthDensity = 0.0;
        volElement = 0.0;
        //CoNumber = 0.0;

        // iterate over all elements
        for (const auto& element : elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

            auto dofGlobalIdx = this->model().dofMapper().subIndex(element, 0, dofCodim);
            GlobalPosition globalPos = element.geometry().center();
            PrimaryVariables Sources;
            // only call the source function if we are not initializing
            if (!(this->timeManager().time() < 0))
            {
                PrimaryVariables priVars;
                initial_(priVars, globalPos);
                //ratioPressure[dofGlobalIdx] = elemVolVars[0].pressure()/priVars[conti0EqIdx];

                ratioMassFraction[dofGlobalIdx] = elemVolVars[0].massFraction(transportCompIdx)/priVars[transportEqIdx];
                diffMassCon[dofGlobalIdx] = (elemVolVars[0].massFraction(transportCompIdx) - priVars[transportEqIdx])*elemVolVars[0].density();

                Scalar initialPc_ = pnRef_-priVars[conti0EqIdx];
                Scalar initialSw_ = MaterialLaw::sw(this->spatialParams().materialLawParams(globalPos),initialPc_);
                //ratioWaterContent[dofGlobalIdx] = elemVolVars[0].saturation(phaseIdx)/initialSw_;
                //diffWaterContent[dofGlobalIdx] = (elemVolVars[0].saturation(phaseIdx) - initialSw_)*elemVolVars[0].porosity();

                this->scvPointSources(Sources, element, fvGeometry, 0, elemVolVars);
                Sources *= fvGeometry.subContVol[0].volume
                                             * this->boxExtrusionFactor(element, fvGeometry, 0);
                sourceP[dofGlobalIdx] += Sources[conti0EqIdx];

                sourceC[dofGlobalIdx] += Sources[transportEqIdx];

                rootSurface[dofGlobalIdx] = this->couplingManager().lowDimSurface(element);
                rootLength[dofGlobalIdx] = this->couplingManager().lowDimLength(element);
                rootLengthDensity[dofGlobalIdx] = this->couplingManager().lowDimLength(element)/fvGeometry.subContVol[0].volume;
                volElement[dofGlobalIdx] = fvGeometry.subContVol[0].volume;
            }
        }

        // attach data to the vtk output
        this->resultWriter().attachDofData(sourceP, "waterSource_[kg/s]", isBox);
        this->resultWriter().attachDofData(sourceC, "soluteSource_[kg/s]", isBox);
        //this->resultWriter().attachDofData(ratioPressure, "pw/pw0_[-]", isBox);
        //this->resultWriter().attachDofData(ratioWaterContent, "theta/theta0_[-]", isBox);
        this->resultWriter().attachDofData(ratioMassFraction, "w/w0_[-]", isBox);
        //this->resultWriter().attachDofData(diffWaterContent, "diffWaterContent_[-]", isBox);
        this->resultWriter().attachDofData(diffMassCon, "diffMassCon_[kg/m3]", isBox);
        this->resultWriter().attachDofData(rootSurface, "rootSurface_[m2]", isBox);
        this->resultWriter().attachDofData(rootLength, "rootLength_[m]", isBox);
        this->resultWriter().attachDofData(rootLengthDensity, "rootLengthDensity_[m/m3]", isBox);
        this->resultWriter().attachDofData(volElement, "voxelVol_[m3]", isBox);
    }

    void preTimeStep()
    {
        //this->couplingManager().resetCouplingSources();
        if ((startEpisode_) and (startEpisodeWithTimeStepSize_ > 0))
        {
            this->timeManager().setTimeStepSize(startEpisodeWithTimeStepSize_);
            std::cout<<" SOIL set time step "<< this->timeManager().timeStepSize()<<"\n";
            startEpisode_ = false;
        }
        //update new born root segments in soil voxels
        //this->couplingManager().resetCouplingSources();

        if ((mimicGrowth_ ) and (this->timeManager().time() < mimicGrowthEndTime_))
        {
        	this->couplingManager().updateLowDimSurfaceInBulkVolume();
        }
        //update initialMassFractionInRhizospheres_ values due to root growth
        initialMassFractionInRhizospheres_.clear();
        for (const auto& element : elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);
            auto soilEIdx = this->elementMapper().index(element);
            Scalar initialTotalMassConcentrationInRhizosphere = 0;
            if (this->couplingManager().newbornRhizoVolumnInBulkElement(soilEIdx) >  0)
            {
                initialTotalMassConcentrationInRhizosphere = (elemVolVars[0].totalMassConcentration()*element.geometry().volume()
                                        - this->couplingManager().rhizoMassInBulkElement(soilEIdx))/
                    (this->couplingManager().newbornRhizoVolumnInBulkElement(soilEIdx));
            }
            initialMassFractionInRhizospheres_[soilEIdx] = elemVolVars[0].dissolveMassFraction(initialTotalMassConcentrationInRhizosphere);

        }
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

        PrimaryVariables mFlowUpBC (0.0);
        PrimaryVariables mFlowBotBC (0.0);

        areaBotBoundary_ = 0.0;
        // iterate over all elements
        bool topSoilSaturated = false;
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

            //GlobalPosition globalPos = intersection.geometry().center();
            //int eIdx = this->elementMapper().index(element);

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (!intersection.boundary())
                    continue;

                GlobalPosition globalPos = intersection.geometry().center();
                if ((globalPos[2] <= this->bBoxMax()[2]- eps_) and (globalPos[2] >= this->bBoxMin()[2]+ eps_))
                    continue;

                unsigned bfIdx = intersection.indexInInside();

                PrimaryVariables BCFlux(0.0);
                this->neumann(BCFlux, element,fvGeometry, intersection, 0, bfIdx);
                BCFlux *= intersection.geometry().volume();

                if (globalPos[2] > this->bBoxMax()[2] - eps_)
                {
                    mFlowUpBC += BCFlux;

                    if (!topSoilSaturated)
                    {
                        ElementVolumeVariables elemVolVars;
                        elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);
                        if (elemVolVars[0].saturation(phaseIdx) > 0.98)
                                            topSoilSaturated = true;
                    }
                }
                else if (globalPos[2] < this->bBoxMin()[2] + eps_)
                {
                    mFlowBotBC += BCFlux;
                    areaBotBoundary_ += intersection.geometry().volume();
                }
            }
        }

        topSoilSaturated_ = topSoilSaturated;
        //if (topSoilSaturated_)
        //    std::cout << "Top soil is saturated. \n";
        std::cout << "Integrated mass source (3D): " << totalSourceP << std::endl;
        std::cout << "Integrated concentration source (3D): " << totalSourceC << std::endl;
        mFluxOut_ = mFlowBotBC;
        mFluxOut_[0] += -totalSourceP;
        mFluxOut_ /=areaBotBoundary_;
        //std::cout<<mFluxOut_[0]<<" "<<areaBotBoundary_<<"\n";
        logFile_.open(nameUptakeRate_, std::ios::app);
        logFile_ << this->timeManager().time() + this->timeManager().timeStepSize()
                << " " << -totalSourceP
                << " " << -totalSourceC
                << " " << mFlowUpBC[0] //up
                << " " << mFlowUpBC[1] //up
                << " " << mFlowBotBC[0] //bot
                << " " << mFlowBotBC[1] //bot
                << std::endl;
        logFile_.close();
        //this->couplingManager().resetCouplingSources();
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
            std::cout<<" SOIL set next time step "<< this->timeManager().timeStepSize()<<"\n";
            startEpisode_ = true;
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
        Scalar sw =  tree.template get<Scalar>("Soil.InitialConditions.InitialSoilWaterContent");
        GlobalPosition soiSurface {0,0,0};
        if (useSoilData_)
        {	//priVars[massOrMoleFracIdx]= soilData_.getValueAtDepth("NutrientConcentration(mg/l)",-globalPos[2])
	        //       /(1 + soilData_.getValueAtDepth("bufferCapacity",-globalPos[2]));

            // mass fraction [convert mg/l to kg/m3 then take ratio with mass water in 1 m3 (1000kg)]
            priVars[massOrMoleFracIdx]= this->spatialParams().initialMassFraction(globalPos);//*1e-6;

        	sw /= this->spatialParams().porosity(soiSurface);
        }
        else
        {	priVars[massOrMoleFracIdx]= tree.template get<Scalar>("Soil.InitialConditions.InitialSoluteMassFracInSoil");
    		sw /= tree.template get<Scalar>("Soil.SpatialParams.Porosity");
        }
        Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParams(soiSurface),sw);

        bool isMultiRoot = tree.template get<bool>("RootSystem.BoundaryConditions.MultiRoot");

        if (isMultiRoot)
            priVars[pressureIdx] = 1e5 - pc;
        else
            priVars[pressureIdx] = 1e5 - pc - 1000*9.81*globalPos[2];

        //priVars[massOrMoleFracIdx]= tree.template get<Scalar>("Soil.InitialConditions.InitialSoluteMassFracInSoil");
	//priVars[massOrMoleFracIdx]= soilData_.getValueAtDepth("NutrientConcentration",-globalPos[2])
        //    /(1 + soilData_.getValueAtDepth("bufferCapacity",-globalPos[2]));
    };

    std::ofstream logFile_;
    const Scalar eps_ = 1e-9;
    Scalar episodeTime, temperature_, pnRef_, pc_, sw_, pw_, uptakeC_num_, uptakeC_ana_, wiltingPointPressure_, startEpisodeWithTimeStepSize_;
    std::string name_, nameUptakeRate_;
    std::shared_ptr<CouplingManager> couplingManager_;
    CSVReader weatherData_;
    bool useWeatherData_, useSoilData_, mimicGrowth_, useRhizoScaleModel_, bottomInFlow_, bottomFreeDrainageFlow_,
                    neglectRootCollarSegment_, reuseCouplingTerms_, periodicVerticalFlow_, startEpisode_, outputEveryTimeStep_, diurnalFluctuation_,
                    rootHydraulicRedistribution_, topSoilSaturated_;

    Scalar previousTimeStepSize_, radialTimeStepRatio_, mimicGrowthEndTime_;
    PrimaryVariables mFluxOut_;
    Scalar areaBotBoundary_;
    //Scalar DiffCoef_;

    mutable std::map<int, Scalar> initialMassFractionInRhizospheres_;
};

} //end namespace

#endif

