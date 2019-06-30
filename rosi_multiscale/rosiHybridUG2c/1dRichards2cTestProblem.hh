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
#ifndef DUMUX_SOIL_RICHARD2C_BUFFER_RADIALLY_SYMMETRIC_TEST_PROBLEM_HH
#define DUMUX_SOIL_RICHARD2C_BUFFER_RADIALLY_SYMMETRIC_TEST_PROBLEM_HH

#include <cmath>
#include <stdlib.h>
#include <math.h>
//#include "math_expint.hh"
#include <boost/math/special_functions/expint.hpp>

#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/richards2cCylindrical1d/richards2cbuffermodel.hh>

#include <dumux/material/fluidsystems/h2osolute.hh>
#include "1dRichardsbuffertestspatialparams.hh"

#define NONISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class SoilRichardsTwoCBufferRadiallySymmetricTestProblem;

namespace Properties
{
NEW_TYPE_TAG(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, INHERITS_FROM(RichardsTwoCBufferRadiallySymmetric, RichardsBufferRadiallySymmetricTestSpatialParams)); //RichardsTwoC from properties file
NEW_TYPE_TAG(SoilRichardsTwoCBufferRadiallySymmetricTestBoxProblem, INHERITS_FROM(BoxModel, SoilRichardsTwoCBufferRadiallySymmetricTestProblem));
NEW_TYPE_TAG(SoilRichardsTwoCBufferRadiallySymmetricTestCCProblem, INHERITS_FROM(CCModel, SoilRichardsTwoCBufferRadiallySymmetricTestProblem));

SET_TYPE_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OSOLUTE<TypeTag>);
              //Dumux::H2OXFluidSystem<TypeTag>);

// Set the grid type
SET_TYPE_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, Grid, Dune::YaspGrid<1, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 1> >);

// Set the problem property
SET_TYPE_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, Problem, Dumux::SoilRichardsTwoCBufferRadiallySymmetricTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, SpatialParams, Dumux::RichardsBufferRadiallySymmetricTestSpatialParams<TypeTag>);
// Enable gravity
SET_BOOL_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, ProblemEnableGravity, false);

// Enable velocity output
SET_BOOL_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, VtkAddVelocity, true);

// Set the grid parameter group
SET_STRING_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, GridParameterGroup, "Rhizosphere.Grid");

//SET_BOOL_PROP(SoilOnePTwoCTestProblem, UseMoles, true);
SET_BOOL_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, UseMoles, false);

SET_TYPE_PROP(SoilRichardsTwoCBufferRadiallySymmetricTestProblem, LinearSolver, UMFPackBackend<TypeTag>);

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
class SoilRichardsTwoCBufferRadiallySymmetricTestProblem : public ImplicitPorousMediaProblem<TypeTag>
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

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    const Dune::ParameterTree &tree = ParameterTree::tree();

public:
    SoilRichardsTwoCBufferRadiallySymmetricTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();
        name_ = tree.template get<std::string>("radialProblemName");
        restartTime_ = ParameterTree::tree().template get<Scalar>("radialTimeManager.RestartTime");

        Vmax_ = tree.template get<Scalar>("RootSystem.SpatialParams.Vmax");
        Km_ = tree.template get<Scalar>("RootSystem.SpatialParams.Km");
        Cmin_ = tree.template get<Scalar>("RootSystem.SpatialParams.Cmin",0);
        rootRadius_ = ParameterTree::tree().template get<std::vector<Scalar>>("Rhizosphere.Grid.Positions0")[0];
        cylinderRadius_ = ParameterTree::tree().template get<std::vector<Scalar>>("Rhizosphere.Grid.Positions0")[2];
        diffustionRootMembrane_ = tree.template get<Scalar>("RootSystem.SpatialParams.DiffCoeffRootMembrane");
        sigma_ = tree.template get<Scalar>("RootSystem.SpatialParams.PartitionCoeff");
        radialRootHydraulicConductivity_ = tree.template get<Scalar>("RootSystem.SpatialParams.Kr");
        shouldWriteOutput_ = ParameterTree::tree().template get<bool>("radialTimeManager.WriteOutput", false);

        pressureRoot_ = tree.template get<Scalar>("RootSystem.BoundaryConditions.RootPressure");
        pressureSoil_ = tree.template get<Scalar>("Soil.BoundaryConditions.InitialSoilPressure");
        saturationMacro_ = MaterialLaw::sw(this->spatialParams().materialLawParams(this->bBoxMax()), (1e5 - pressureSoil_));
        initialSoluteMassFracInSoil_ = tree.template get<Scalar>("Soil.BoundaryConditions.InitialSoluteMassFracInSoil");
        compFracRoot_ = tree.template get<Scalar>("RootSystem.InitialConditions.InitialSoluteMassFracInRoot");
        outerFluxesWater_ = tree.template get<Scalar>("Rhizosphere.BoundaryConditions.OuterFluxesWater");
        outerFluxesNutrient_ = tree.template get<Scalar>("Rhizosphere.BoundaryConditions.OuterFluxesNutrient");
        //std::cout<<"!!!! radialProblem "<< name_<<" "<<tree.template get<std::vector<Scalar>>("Rhizosphere.Grid.Positions0")[0]<<" "
        //    <<tree.template get<std::vector<Scalar>>("Rhizosphere.Grid.Positions0")[1]<<"\n";
        pnRef_ = 1e5;
        uptakeRate_ = 0.0;
        cumulativeUptake_ = 0.0;
        cumulativeTime_ = 0.0;
        wiltingPointPressure_ = tree.template get<Scalar>("RootSystem.BoundaryConditions.CriticalCollarPressure", -1.4e6);
        rootHydraulicRedistribution_ = tree.template get<bool>("RootSystem.BoundaryConditions.RootHydraulicRedistribution", true);
        couplingInitialTime_ = tree.template get<bool>("radialTimeManager.CouplingInitialTime", true);
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
        //soilTemp_ = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.SoilTemperature);
        soilTemp_ = 10;//tree.template get<Scalar>("BoundaryConditions.SoilTemperature");
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
        //if (globalPos[0] < this->bBoxMin()[0] + eps_)
        //    values.setAllNeumann();
        //else
        //    values.setAllDirichlet();
        values.setAllNeumann();
        if (globalPos[0] > this->bBoxMax()[0] - eps_)
        //if (globalPos[0] < this->bBoxMin()[0] + eps_)
            values.setDirichlet(conti0EqIdx);
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
        priVars[pressureIdx] = pressureSoil_;
        //priVars[massOrMoleFracIdx] = tree.template get<Scalar>("BoundaryConditions.InitialSoluteMassFracInSoil");
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
        auto posOuterRadius = std::max(element.geometry().corner(0)[0],element.geometry().corner(1)[0]);
        auto posInnerRadius = std::min(element.geometry().corner(0)[0],element.geometry().corner(1)[0]);
        if (posOuterRadius > this->bBoxMax()[0] - eps_)
        {
            values[conti0EqIdx] = outerFluxesWater_;
            values[transportEqIdx] = outerFluxesNutrient_;
            values *=posOuterRadius;
        }
        else
        {
            Scalar pressureSoil = pressureSoil_;//elemVolVars[scvIdx].pressure();
            if (((pressureRoot_ < pressureSoil)  and  (pressureSoil > wiltingPointPressure_))
                     or ((pressureRoot_ > pressureSoil) and (rootHydraulicRedistribution_)))
                values[conti0EqIdx] = radialRootHydraulicConductivity_*(pressureSoil - pressureRoot_)*elemVolVars[scvIdx].density();

            Scalar c3D = useMoles ? elemVolVars[scvIdx].moleFraction(transportCompIdx) : elemVolVars[scvIdx].massFraction(transportCompIdx);

            Scalar passive_uptake = (values[conti0EqIdx]>0) ? values[conti0EqIdx]*c3D : 0; //uptake by advection

            passive_uptake += diffustionRootMembrane_*(c3D - compFracRoot_); //uptake by diffussion

            passive_uptake *= useMoles ? elemVolVars[scvIdx].molarDensity() : elemVolVars[scvIdx].density();

            Scalar active_uptake = 0.0;
            if (c3D*elemVolVars[scvIdx].molarDensity() - Cmin_ > 0)
                active_uptake =  useMoles ? Vmax_*(c3D*elemVolVars[scvIdx].molarDensity()-Cmin_)/(Km_+c3D*elemVolVars[scvIdx].molarDensity()-Cmin_)
                                                : Vmax_*(c3D*elemVolVars[scvIdx].density()-Cmin_)/(Km_+c3D*elemVolVars[scvIdx].density()-Cmin_);

            values[transportEqIdx] = (sigma_*passive_uptake + (1-sigma_)*active_uptake);
            values *=posInnerRadius;
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

    /*!
     * \brief If we should write output
     */
    bool shouldWriteRestartFile()
    {
        return (shouldWriteOutput_ and this->timeManager().willBeFinished());
    }

    void addOutputVtkFields()
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        unsigned numDofs = this->model().numDofs();

        // create required scalar fields for the vtk output
        //ScalarField& ratioPressure = *(this->resultWriter().allocateManagedBuffer(numDofs));
        //ScalarField& ratioWaterContent = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ratioMassFraction = *(this->resultWriter().allocateManagedBuffer(numDofs));

        //ratioPressure = 0.0;
        //ratioWaterContent = 0.0;
        ratioMassFraction = 0.0;

        if (!(this->timeManager().time() < 0))
        {
	        // iterate over all elements
	        for (const auto& element : elements(this->gridView()))
	        {
	            FVElementGeometry fvGeometry;
	            fvGeometry.update(this->gridView(), element);

	            ElementVolumeVariables elemVolVars;
	            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

	            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
	            {
		            // only call the source function if we are not initializing
	            	auto dofGlobalIdx = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);

		            GlobalPosition globalPos;
	                if (isBox)
	                        globalPos = element.geometry().corner(scvIdx);
	                else
	                        globalPos = element.geometry().center();

		            PrimaryVariables priVars;
		            initial_(priVars, globalPos);

	                //ratioPressure[dofGlobalIdx] = elemVolVars[scvIdx].pressure()/priVars[conti0EqIdx];
	                ratioMassFraction[dofGlobalIdx] = elemVolVars[scvIdx].massFraction(transportCompIdx)/priVars[transportEqIdx];

	                //Scalar initialPc_ = pnRef_-priVars[conti0EqIdx];
	                //Scalar initialSw_ = MaterialLaw::sw(this->spatialParams().materialLawParams(globalPos),initialPc_);
	                //ratioWaterContent[dofGlobalIdx] = elemVolVars[scvIdx].saturation(phaseIdx)/initialSw_;
	          	}
	         }
        }
        // attach data to the vtk output
        //this->resultWriter().attachDofData(ratioPressure, "pw/pw0", isBox);
        //this->resultWriter().attachDofData(ratioWaterContent, "theta/theta0", isBox);
        this->resultWriter().attachDofData(ratioMassFraction, "w/w0", isBox);
    }

    void preTimeStep()
    {
        if (this->timeManager().time()>=restartTime_)
        //if (this->timeManager().willBeFinished())
        {
            Scalar previousNutrientFraction = 0;
            Scalar previousPosition = 0;
            Scalar MaxNutrientGradient = 0;
            Scalar NutrientGradient = 0;
            Scalar rhizosphereRadius = 0;
            Scalar concentrationAtRootSurface = 0;
            Scalar porousDiffCoeff = 0;
            Scalar waterPressureAtRootSurface = 0;
            Scalar waterContentAtRootSurface = 0;
            Scalar PeclectAtRootSurface = 0;

            int ElementIdx = 0;
            // iterate over all elements
            for (const auto& element : elements(this->gridView()))
            {
                FVElementGeometry fvGeometry;
                fvGeometry.update(this->gridView(), element);

                ElementVolumeVariables elemVolVars;
                elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

                auto posInnerVertex = std::min(element.geometry().corner(0)[0],element.geometry().corner(1)[0]);
                if (posInnerVertex < this->bBoxMin()[0] + eps_)
                {
                    Scalar pressureSoil = pressureSoil_;//elemVolVars[0].pressure();
                    uptakeRate_[conti0EqIdx] = radialRootHydraulicConductivity_*(pressureSoil - pressureRoot_)*elemVolVars[0].density();
                    waterPressureAtRootSurface = elemVolVars[0].pressure();
                    waterContentAtRootSurface = elemVolVars[0].waterContent(0);
                    PeclectAtRootSurface = PelectNumber( element, fvGeometry, elemVolVars);

                    Scalar c3D = useMoles ? elemVolVars[0].moleFraction(transportCompIdx) : elemVolVars[0].massFraction(transportCompIdx);

                    //uptake by advection
                    Scalar passive_uptake = (uptakeRate_[conti0EqIdx]>0) ? uptakeRate_[conti0EqIdx]*c3D : 0;
                    passive_uptake += diffustionRootMembrane_*(c3D - compFracRoot_);

                    ////uptake by diffussion
                    passive_uptake *= useMoles ? elemVolVars[0].molarDensity() : elemVolVars[0].density();
                    Scalar active_uptake = 0;// =  useMoles ? Vmax_*c3D*elemVolVars[0].molarDensity()/(Km_+c3D*elemVolVars[0].molarDensity())
                                                       // : Vmax_*c3D*elemVolVars[0].density()/(Km_+c3D*elemVolVars[0].density());
                    if (c3D*elemVolVars[0].molarDensity() - Cmin_ > 0)
                        active_uptake =  useMoles ? Vmax_*(c3D*elemVolVars[0].molarDensity()-Cmin_)/(Km_+c3D*elemVolVars[0].molarDensity()-Cmin_)
                                                                    : Vmax_*(c3D*elemVolVars[0].density()-Cmin_)/(Km_+c3D*elemVolVars[0].density()-Cmin_);

                    uptakeRate_[transportEqIdx] = (sigma_*passive_uptake + (1-sigma_)*active_uptake);

                    concentrationAtRootSurface = useMoles ? c3D*elemVolVars[0].molarDensity() : c3D*elemVolVars[0].density();

                    cumulativeUptake_[0] += uptakeRate_[0]*this->timeManager().timeStepSize();
                    cumulativeUptake_[1] += uptakeRate_[1]*this->timeManager().timeStepSize();
                    cumulativeTime_ += this->timeManager().timeStepSize();

                    if (this->timeManager().time()==restartTime_)
                    {
                        initialUptakeRate_ = uptakeRate_;
                        //porousDiffCoeff at root soil interface
                        FluxVariables fluxVars;
                        fluxVars.update(*this, element, fvGeometry, 0 /*scvfIdx*/, elemVolVars, false /* oldSol? */);
                        porousDiffCoeff = fluxVars.porousDiffCoeff();
                    }
                }

                //calculate the extension of nutrient gradient (assumed to have monotonic shape)
                Scalar c3D = useMoles ? elemVolVars[0].moleFraction(transportCompIdx) : elemVolVars[0].massFraction(transportCompIdx);

                if (ElementIdx == 1)
                {
                    MaxNutrientGradient = std::abs(c3D - previousNutrientFraction)/(element.geometry().center()[0] - previousPosition);
                    if (MaxNutrientGradient == 0)
                    {
                        rhizosphereRadius = 0;
                        break;
                    }
                }
                else if (ElementIdx != 0)
                {
                    NutrientGradient = std::abs(c3D - previousNutrientFraction)/(element.geometry().center()[0] - previousPosition);
                    if (NutrientGradient <= MaxNutrientGradient*0.01) break;
                }

                previousNutrientFraction = c3D;
                previousPosition = element.geometry().center()[0];
                rhizosphereRadius = std::max(element.geometry().corner(0)[0],element.geometry().corner(1)[0]);
                waterContentMacro_ = saturationMacro_*elemVolVars[0].porosity();
                ElementIdx ++;
            }

            //if (this->timeManager().time() == 0)
            //{   logFile_.open(name_ + "_uptakeRate.log", std::ofstream::out | std::ofstream::trunc);
            //    logFile_ << "Time[s]"
            //                        << " " <<"WaterUptakeFlux[kg/s]"
            //                        << " " <<"NutrientUptakeFlux[kg/s]"
            //                        << " " <<"ConcentrationAtRootSurface[kg/m3]"
            //                        << " " <<"RhizosphereRadius[m]"
            //                        << " " <<"DepletionZone[m]"
            //                        << " " <<"DepletionRatio[-]"
            //                        << std::endl;
            //}
            //else
            if (this->timeManager().time()==restartTime_)
            {
                Scalar DepletionZone = (rhizosphereRadius - rootRadius_) > 0 ? (rhizosphereRadius - rootRadius_) : 0;
                logFile_.open(name_ + "_uptakeRate.log", std::ios::app);
                logFile_ << this->timeManager().time()
                    << " " <<uptakeRate_[0]
                    << " " <<uptakeRate_[1]
                    << " " <<concentrationAtRootSurface
                    << " " <<rhizosphereRadius
                    << " " <<DepletionZone
                    << " " <<DepletionZone/(cylinderRadius_ - rootRadius_)
                    << " " <<cylinderRadius_
                    << " " <<porousDiffCoeff
                    << " " <<waterPressureAtRootSurface
                    << " " <<waterContentAtRootSurface
                    << " " <<pressureSoil_
                    << " " <<waterContentMacro_
                    << " " <<PeclectAtRootSurface
                    << " " <<waterContentAtRootSurface/waterContentMacro_
                    << std::endl;
                logFile_.close();
            }
        }
    }
    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        if (this->timeManager().willBeFinished())
        {
            if (!shouldWriteOutput_)
                for (const auto& solution : this->model().curSol())
                    Solutions_.push_back(solution);

            //Calculate total nutrient mass in the rhizosphere
            Scalar nutrientMassPerLength = 0;
            Scalar massConcentrationAtBoundary = 0;
            Scalar outerBoundary = 0;
            for (const auto& element : elements(this->gridView()))
            {
                FVElementGeometry fvGeometry;
                fvGeometry.update(this->gridView(), element);

                ElementVolumeVariables elemVolVars;
                elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

                auto posInnerVertex = std::min(element.geometry().corner(0)[0],element.geometry().corner(1)[0]);
                auto posOuterVertex = std::max(element.geometry().corner(0)[0],element.geometry().corner(1)[0]);
                auto massConcentration = elemVolVars[0].totalMassConcentration();
                nutrientMassPerLength += M_PI*(posOuterVertex*posOuterVertex - posInnerVertex*posInnerVertex)*massConcentration;
                if (posOuterVertex > this->bBoxMax()[0] - eps_)
                    {
                        outerBoundary = posOuterVertex;
                        massConcentrationAtBoundary = massConcentration;
                    }
            }
            //
            PrimaryVariables rhizoBC;
            rhizoBC[0] = outerBoundary;
            rhizoBC[1] = massConcentrationAtBoundary;
            Solutions_.push_back(rhizoBC);
            //
            PrimaryVariables massInRhizospherePerLength;
            massInRhizospherePerLength[0] = 0;
            massInRhizospherePerLength[1] = nutrientMassPerLength;
            Solutions_.push_back(massInRhizospherePerLength); //
            //
            if (couplingInitialTime_)
                Solutions_.push_back(initialUptakeRate_);
            else
            {
                PrimaryVariables averagedUptakeRate;
                averagedUptakeRate[0] =cumulativeUptake_[0]/cumulativeTime_;
                averagedUptakeRate[1] =cumulativeUptake_[1]/cumulativeTime_;
                Solutions_.push_back(averagedUptakeRate);
            }
            //std::cout<<"RHIZO uptake: "<<Solutions_.end()[0]<<" "<<Solutions_.end()[1]<<" "<<"\n";
        }
    }

    auto uptakeRate() const
    {
        return Solutions_ ;
    }

    Scalar PelectNumber( Element element, FVElementGeometry fvGeometry,
                                                 ElementVolumeVariables elemVolVars)
    {
        FluxVariables fluxVars;
        fluxVars.update(*this, element, fvGeometry, 0/*scvfIdx*/, elemVolVars, false /* oldSol? */);

        Scalar advFlux = 0;
        advFlux = fluxVars.radialVolumeFlux(phaseIdx) * elemVolVars[0].density() * elemVolVars[0].massFraction(1);

        Scalar diffFlux = 0;
        diffFlux = (fluxVars.radialConcentrationGrad(transportCompIdx)*fluxVars.face().normal)*fluxVars.porousDiffCoeff();

        return (advFlux==0)? 0:advFlux/diffFlux;
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
        return (((this->timeManager().willBeFinished()) or (this->timeManager().time()<0)) and
                        shouldWriteOutput_);
    }


private:
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        priVars[pressureIdx] = pressureSoil_;
        priVars[massOrMoleFracIdx] = initialSoluteMassFracInSoil_;
    };

    std::ofstream logFile_;
    const Scalar eps_ = 1e-9;
    Scalar pnRef_, Vmax_, Km_, Cmin_, rootRadius_, cylinderRadius_, diffustionRootMembrane_, compFracRoot_, sigma_, pressureRoot_, pressureSoil_,
                    initialSoluteMassFracInSoil_;
    Scalar radialRootHydraulicConductivity_, restartTime_, outerFluxesWater_, outerFluxesNutrient_, wiltingPointPressure_;
    std::string name_;
    PrimaryVariables uptakeRate_,initialUptakeRate_, cumulativeUptake_;
    Scalar cumulativeTime_;
    std::vector<PrimaryVariables> Solutions_;
    bool shouldWriteOutput_, rootHydraulicRedistribution_, couplingInitialTime_;
    Scalar saturationMacro_, waterContentMacro_;

//    std::shared_ptr<CouplingManager> couplingManager_;
    //Scalar DiffCoef_;
};

} //end namespace

#endif
