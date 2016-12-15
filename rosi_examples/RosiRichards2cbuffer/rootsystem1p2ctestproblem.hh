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

//#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/fluidsystems/liquidphase.hh>
// fluidsystem
//#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/fluidsystems/liquidphase.hh>
//#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/h2ono3.hh>
//! get the properties needed for subproblems
#include <dumux/multidimension/subproblemproperties.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include "rootsystemtestspatialparams.hh"

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
              Dumux::FluidSystems::H2ONO3<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

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

public:
    RootsystemOnePTwoCTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "-root";
        this->spatialParams().setParams();
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
        rootTemp_ = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.RootTemperature);
     return rootTemp_;
    } // 10C
//    Scalar temperature() const
//    { temperatureRoot_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.TemperatureRoot);
//    return temperatureRoot_; } // 10C

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    void boundaryTypesAtPos (BoundaryTypes &values,
                             const GlobalPosition &globalPos ) const
    {

        //////GlobalPosition globalPos = intersection.geometry().center();
        if (globalPos[2] + eps_ >  this->bBoxMax()[2] )
            {
                values.setAllDirichlet();
        ////         outflow condition for the transport equation at upper boundary
        ////        values.setOutflow(transportEqIdx);
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
//                                               BoundaryConditions.CriticalCollarPressure);
//
//          //values[pIdx] *= 1.01;
//    }

    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.CollarPressure);
        values[massOrMoleFracIdx] = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.InitialRootFracNO3);
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
    //    ElementVolumeVariables elemVolVars;
    //    elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);
//
    //    GlobalPosition globalPos = intersection.geometry().center();
    //   if (globalPos[2] + eps_ >  this->bBoxMax()[2] )
    //       {
    //       values[conti0EqIdx] = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.TranspirationRate);
    //       values[transportEqIdx] = values[conti0EqIdx]*elemVolVars[scvIdx].moleFraction(massOrMoleFracIdx);
    //       }
    //   else
    //       {
    //       values[conti0EqIdx] = 0;
    //       values[transportEqIdx] = 0;
    //       }
        values[conti0EqIdx] = 0;
        values[transportEqIdx] = 0;
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
        priVars[pressureIdx] =  GET_RUNTIME_PARAM(TypeTag,
                                           Scalar,
                                           BoundaryConditions.InitialRootPressure);
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
        // compute source at every integration point
        const SpatialParams &spatialParams = this->spatialParams();
        const Scalar Kr = spatialParams.Kr(element, fvGeometry, scvIdx);
        const Scalar rootRadius = spatialParams.rootRadius(element, fvGeometry, scvIdx);

        // convert units of 3d pressure if pressure head is used !!!
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[pressureIdx];
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[pressureIdx];

        PrimaryVariables sourceValues;
        sourceValues=0.0;
        // sink defined as radial flow Jr [m^3 s-1]*density
        sourceValues[conti0EqIdx] = 2 * M_PI *rootRadius * Kr *(pressure3D - pressure1D)
                                   * elemVolVars[scvIdx].density();

        // needs concentrations in soil and root
        Scalar c1D;
        if(useMoles)
            c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];
        else
            c1D = this->couplingManager().lowDimPriVars(source.id())[massOrMoleFracIdx];
        //std::cout << "concentrations c1D " <<c1D<< std::endl;
        Scalar c3D;
        if(useMoles)
            c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];
        else
            c3D = this->couplingManager().bulkPriVars(source.id())[massOrMoleFracIdx];

        //Difussive flux term of transport
        const Scalar DiffValue = 0.0;
        //2* M_PI *rootRadius *DiffCoef_*(c1D - c3D)*elemVolVars[scvIdx].density(/*phaseIdx=*/0);

        //Advective flux term of transport
        Scalar AdvValue;
        if (sourceValues[conti0EqIdx]>0)
            AdvValue = 2 * M_PI *rootRadius * Kr *(pressure3D - pressure1D)
                                   * elemVolVars[scvIdx].density()*c3D;
        else
            AdvValue = 2 * M_PI *rootRadius * Kr *(pressure3D - pressure1D)
                                   * elemVolVars[scvIdx].density()*c1D;


        //Active flux - active uptake based on Michaeles Menten
        Scalar ActiveValue;
        ActiveValue = 0;
        const Scalar Vmax = spatialParams.Vmax();
        const Scalar Km = spatialParams.Km();
        ActiveValue = 2 * M_PI*rootRadius*Vmax*c3D*elemVolVars[scvIdx].density()/(Km+c3D* elemVolVars[scvIdx].density());

        Scalar sigma;
        sigma = spatialParams.PartitionCoeff();

        sourceValues[transportEqIdx] = (sigma*(AdvValue + DiffValue) + (1-sigma)*ActiveValue)*source.quadratureWeight()*source.integrationElement();
        sourceValues[conti0EqIdx] *= source.quadratureWeight()*source.integrationElement();

        //std::cout << "ROOT transportEqIdx " << transportEqIdx <<" "<< sourceValues[transportEqIdx]<< std::endl;
        //std::cout << "ROOT conti0EqIdx " << conti0EqIdx <<" "<< sourceValues[conti0EqIdx]<< std::endl;

        //sourceValues[transportEqIdx] =1e-9*source.quadratureWeight()*source.integrationElement();
        source = sourceValues;
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
                    sourceP[dofGlobalIdx] += values[pressureIdx] * fvGeometry.subContVol[scvIdx].volume
                                             * this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                    sourceC[dofGlobalIdx] += values[massOrMoleFracIdx] * fvGeometry.subContVol[scvIdx].volume
                                             * this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                }
            }
        }
        auto totalSourceP = std::accumulate(sourceP.begin(), sourceP.end(), 0);
        auto totalSourceC = std::accumulate(sourceC.begin(), sourceC.end(), 0);

        std::cout << "Integrated mass source (1D): " << totalSourceP << std::endl;
        std::cout << "Integrated concentration source (1D): " << totalSourceC << std::endl;

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
    std::string name_;

    Scalar sinkSum;
    const Scalar eps_ = 1e-9;
    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif
