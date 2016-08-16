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
#ifndef DUMUX_SOIL_RICHARD2C_TEST_PROBLEM_HH
#define DUMUX_SOIL_RICHARD2C_TEST_PROBLEM_HH

#include <cmath>

// explicitly set the fvgeometry as we currently use another one as in stable
#include <dumux/multidimension/cellcentered/fvelementgeometry.hh>
#include <dumux/multidimension/box/fvelementgeometry.hh>

//#include <dumux/porousmediumflow/richards/implicit/model.hh>

//#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/fluidsystems/liquidphase.hh>

//#include <dumux/porousmediumflow/1p2c/implicit/model.hh>
#include <dumux/porousmediumflow/richards2c/richards2cmodel.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

//#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/h2ono3.hh>

//! get the properties needed for subproblems
#include <dumux/multidimension/subproblemproperties.hh>

//#include "soiltestspatialparams1p2c.hh"
 #include "richardstestspatialparams.hh"

#define NONISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class SoilRichardsTwoCTestProblem;

namespace Properties
{
NEW_TYPE_TAG(SoilRichardsTwoCTestProblem, INHERITS_FROM(RichardsTwoC, RichardsTestSpatialParams)); //RichardsTwoC from properties file
NEW_TYPE_TAG(SoilRichardsTwoCTestBoxProblem, INHERITS_FROM(BoxModel, SoilRichardsTwoCTestProblem));
NEW_TYPE_TAG(SoilRichardsTwoCTestCCProblem, INHERITS_FROM(CCModel, SoilRichardsTwoCTestProblem));

// Set the wetting phase
//SET_PROP(RichardsTestProblem, WettingPhase)
//{
//private:
//    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//public:
//    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
//};
// Set fluid configuration
SET_TYPE_PROP(SoilRichardsTwoCTestProblem,
              FluidSystem,
              Dumux::FluidSystems::H2ONO3<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

// Set the grid type
SET_TYPE_PROP(SoilRichardsTwoCTestProblem, Grid, Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3> >);
//SET_TYPE_PROP(RichardsTestProblem, Grid, Dune::UGGrid<3>);
//SET_TYPE_PROP(RichardsTestProblem, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::conforming>);

// explicitly set the fvgeometry as we currently use another one as in stable
SET_TYPE_PROP(SoilRichardsTwoCTestBoxProblem, FVElementGeometry, Dumux::MultiDimensionBoxFVElementGeometry<TypeTag>);
SET_TYPE_PROP(SoilRichardsTwoCTestCCProblem, FVElementGeometry, Dumux::MultiDimensionCCFVElementGeometry<TypeTag>);

// Set the problem property
SET_TYPE_PROP(SoilRichardsTwoCTestProblem, Problem, Dumux::SoilRichardsTwoCTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(SoilRichardsTwoCTestProblem, SpatialParams, Dumux::RichardsTestSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(SoilRichardsTwoCTestProblem, ProblemEnableGravity, false);

// Enable velocity output
SET_BOOL_PROP(SoilRichardsTwoCTestProblem, VtkAddVelocity, true);

// Set the grid parameter group
SET_STRING_PROP(SoilRichardsTwoCTestProblem, GridParameterGroup, "SoilGrid");

//SET_BOOL_PROP(SoilRichardsTwoCTestProblem, UseHead, false);
//SET_BOOL_PROP(SoilOnePTwoCTestProblem, UseMoles, true);
SET_BOOL_PROP(SoilRichardsTwoCTestProblem, UseMoles, false);

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
class SoilRichardsTwoCTestProblem : public ImplicitPorousMediaProblem<TypeTag>
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
        //contiEqIdx = Indices::contiEqIdx,
        //hIdx = Indices::hIdx,
        pressureIdx = Indices::pwIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx
#endif
    };
     enum {
        // index of the transport equation
        conti0EqIdx = Indices::contiEqIdx,
        transportEqIdx = Indices::transportEqIdx,
        wPhaseIdx = Indices::wPhaseIdx,
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

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag) GlobalProblemTypeTag;
    typedef typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager) CouplingManager;

public:
    SoilRichardsTwoCTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "-soil";
        //stating in the console whether mole or mass fractions are used
        //episodeTime = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
        //                                           Scalar,
        //                                           TimeManager,
        //                                           EpisodeTime);
        //this->timeManager().startNextEpisode(episodeTime);// time episode

        if(useMoles)
        {
            std::cout<<"problem uses mole fractions"<<std::endl;
        }
        else
        {
            std::cout<<"problem uses mass fractions"<<std::endl;
        }
        pnRef_ = 1e5;
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
        // compute source at every integration point
        // needs convertion of  units of 1d pressure if pressure head in richards is used
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[pressureIdx];
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[pressureIdx] ;

        const auto& spatialParams = this->couplingManager().lowDimProblem().spatialParams();
        const unsigned int rootEIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = spatialParams.Kr(rootEIdx);
        const Scalar rootRadius = spatialParams.radius(rootEIdx);

        PrimaryVariables sourceValues;
        sourceValues=0.0;
        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        sourceValues[conti0EqIdx] = 2* M_PI *rootRadius * Kr *(pressure1D - pressure3D)
                                   *elemVolVars[scvIdx].density(wPhaseIdx);
        // needs concentrations in soil and root
        Scalar c1D;
        if(useMoles)
            c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];//*elemVolVars[0].molarDensity();
        else
            c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];//*1000;//elemVolVars[0].density();
        //std::cout << "concentrations c1D " <<c1D<< std::endl;
        Scalar c3D;
        if(useMoles)
            c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];//*elemVolVars[0].molarDensity();
        else
            c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];//*1000;//elemVolVars[0].density();
        //std::cout << "concentrations c3D " <<c3D<< std::endl;
        //Diffussive flux term of transport
        Scalar DiffValue;
        Scalar DiffCoef_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.DiffCoeffRootSurface);
        DiffValue = 2* M_PI *rootRadius *DiffCoef_*(c1D - c3D)*elemVolVars[scvIdx].density(wPhaseIdx);

        //Advective flux term of transport
        Scalar AdvValue;
        //if (sourceValues[conti0EqIdx]>0)
        //    //AdvValue = 0;
        //    AdvValue = sourceValues[conti0EqIdx]*c1D;
        //else
        //    AdvValue = sourceValues[conti0EqIdx]*c3D;
        AdvValue = 0;

        sourceValues[transportEqIdx] = (AdvValue + DiffValue)*source.quadratureWeight()*source.integrationElement();
        //std::cout << "SOIL transportEqIdx " << transportEqIdx <<" " << sourceValues[transportEqIdx]<< std::endl;
        sourceValues[conti0EqIdx] *= source.quadratureWeight()*source.integrationElement();
        //std::cout << "SOIL conti0EqIdx " << conti0EqIdx <<" " << sourceValues[conti0EqIdx]<< std::endl;
        //sourceValues[transportEqIdx] =-1e-9*source.quadratureWeight()*source.integrationElement();
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
        //values.setAllNeumann();

        //values.setAllDirichlet();
        if(globalPos[2] > this->bBoxMax()[2] - eps_)
        {
            values.setAllNeumann();
        }
        else
        {
            values.setAllDirichlet();
            //values.setOutflow(transportEqIdx);
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
    //    initial_(values, globalPos);
        Scalar sw_ = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilSaturation);
        Scalar pc_ = MaterialLaw::pc(this->spatialParams().materialLawParams(globalPos),sw_);

        priVars[pressureIdx] = pnRef_ - pc_;

        priVars[massOrMoleFracIdx] = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilFracNO3);
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
        sourceP = 0.0;
        sourceC = 0.0;

        // iterate over all elements
        for (const auto& element : elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

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
                    sourceC[dofGlobalIdx] += values[transportEqIdx] * fvGeometry.subContVol[scvIdx].volume
                                            * this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                }
            }
        }

        const auto totalSourceP = std::accumulate(sourceP.begin(), sourceP.end(), 0);
        const auto totalSourceC = std::accumulate(sourceC.begin(), sourceC.end(), 0);

        std::cout << "Integrated mass source (3D): " << totalSourceP << std::endl;
        std::cout << "Integrated concentration source (3D): " << totalSourceC << std::endl;

        // attach data to the vtk output
        this->resultWriter().attachDofData(sourceP, "mass source (kg/s)", isBox);
        this->resultWriter().attachDofData(sourceC, "concentration source (kg/s)", isBox);
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
        Scalar sw_ = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilSaturation);
        Scalar pc_ = MaterialLaw::pc(this->spatialParams().materialLawParams(globalPos),sw_);

        priVars[pressureIdx] = pnRef_ - pc_;

        priVars[massOrMoleFracIdx] = GET_RUNTIME_PARAM(TypeTag,
                                   Scalar,
                                   BoundaryConditions.InitialSoilFracNO3);

//        priVars[pressureIdx] = GET_RUNTIME_PARAM(TypeTag,
//                                   Scalar,
//                                   BoundaryConditions.InitialSoilPressure); // initial condition for the pressure
//
//        priVars[massOrMoleFracIdx] = GET_RUNTIME_PARAM(TypeTag,
//                                   Scalar,
//                                   BoundaryConditions.InitialSoilFracN2);  // initial condition for the N2 molefraction
//    //    #if !NONISOTHERMAL
//    //        priVars[temperatureIdx] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Temperature);
//    //    #endif
};

    const Scalar eps_ = 1e-9;
    Scalar episodeTime, temperature_, pnRef_, pc_, sw_;
    std::string name_;
    std::shared_ptr<CouplingManager> couplingManager_;
    //Scalar DiffCoef_;
};

} //end namespace

#endif
