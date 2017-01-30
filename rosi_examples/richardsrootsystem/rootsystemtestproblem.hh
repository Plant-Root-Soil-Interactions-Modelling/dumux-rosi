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
#ifndef DUMUX_ROOTSYSTEM_TEST_PROBLEM_HH
#define DUMUX_ROOTSYSTEM_TEST_PROBLEM_HH

#include <cmath>

#include <dumux/porousmediumflow/1d/rootsystem/model.hh>
#include <dumux/porousmediumflow/1d/rootsystem/problem.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

//! get the properties needed for subproblems
#include <dumux/multidimension/subproblemproperties.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include "rootsystemtestspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class RootsystemTestProblem;

namespace Properties
{
NEW_TYPE_TAG(RootsystemTestProblem, INHERITS_FROM(Rootsystem));
NEW_TYPE_TAG(RootsystemTestBoxProblem, INHERITS_FROM(BoxModel, RootsystemTestProblem));
NEW_TYPE_TAG(RootsystemTestCCProblem, INHERITS_FROM(CCModel, RootsystemTestProblem));

SET_PROP(RootsystemTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(RootsystemTestProblem, Grid, Dune::FoamGrid<1, 3>);

// Set the problem property
SET_TYPE_PROP(RootsystemTestProblem, Problem, Dumux::RootsystemTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(RootsystemTestProblem, SpatialParams, Dumux::RootsystemTestSpatialParams<TypeTag>);

// Linear solver settings
#if HAVE_UMFPACK
SET_TYPE_PROP(RootsystemTestProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif

// Enable gravity
SET_BOOL_PROP(RootsystemTestProblem, ProblemEnableGravity, true);

// write newton convergence to vtk
SET_BOOL_PROP(RootsystemTestProblem, NewtonWriteConvergence, false);

// Enable velocity output
SET_BOOL_PROP(RootsystemTestProblem, VtkAddVelocity, true);

//SET_BOOL_PROP(RootsystemTestProblem, GrowingGrid, false);
}

/*!
 * \ingroup OneDRootSystem
 * \ingroup ImplicitTestProblems
 * \brief TODO
 */
template <class TypeTag>
class RootsystemTestProblem : public RootsystemProblem<TypeTag>
{
    typedef RootsystemProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

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
        pIdx = Indices::pIdx
    };
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

    typedef typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag) GlobalProblemTypeTag;
    typedef typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager) CouplingManager;

public:
    RootsystemTestProblem(TimeManager &timeManager, const GridView &gridView)
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
    { return 273.15 + 10; } // 10C


    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    void boundaryTypesAtPos (BoundaryTypes &values,
                             const GlobalPosition &globalPos ) const
    {
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
    void dirichlet(PrimaryVariables &values,
                   const Intersection &intersection) const
    {
        GlobalPosition globalPos = intersection.geometry().center();
        if (globalPos[2] + eps_ >  this->bBoxMax()[2] )
              values[pIdx] = GET_RUNTIME_PARAM(TypeTag,
                                               Scalar,
                                               BoundaryConditions.CriticalCollarPressure);

          //values[pIdx] *= 1.01;
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
        GlobalPosition globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        if (globalPos[2] + eps_ > this->bBoxMax()[2] ){
             values = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.TranspirationRate);
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
        priVars[pIdx] =  GET_RUNTIME_PARAM(TypeTag,
                                           Scalar,
                                           BoundaryConditions.InitialRootPressure);
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
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[pIdx];
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[pIdx];

        // sink defined as radial flow Jr [m^3 s-1]*density
        const Scalar sourceValue = 2 * M_PI *rootRadius *  Kr *(pressure3D - pressure1D)
                                   * elemVolVars[scvIdx].density();
        source = sourceValue*source.quadratureWeight()*source.integrationElement();
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

     void addOutputVtkFields()
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        unsigned numDofs = this->model().numDofs();

        // create required scalar fields for the vtk output
        ScalarField& source = *(this->resultWriter().allocateManagedBuffer(numDofs));
        source = 0.0;

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
                if (!(this->timeManager().time() < 0))
                {
                    PrimaryVariables values;
                    this->scvPointSources(values, element, fvGeometry, scvIdx, elemVolVars);
                    source[dofGlobalIdx] += values[pIdx] * fvGeometry.subContVol[scvIdx].volume
                                            * this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                }
            }
        }


        // attach data to the vtk output
        this->resultWriter().attachDofData(source, "source (kg/s)", isBox);
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
