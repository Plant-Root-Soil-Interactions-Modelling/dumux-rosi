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
//#include <dumux/multidimension/cellcentered/fvelementgeometry.hh>
//#include <dumux/multidimension/box/fvelementgeometry.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
//#include "./richards2cCylindrical1d/richards2cbuffermodel.hh"
#include <dumux/porousmediumflow/richards2cbuffer/richards2cbuffermodel.hh>

#include <dumux/material/fluidsystems/h2ok.hh>
//#include "h2oX.hh"
//! get the properties needed for subproblems
//#include <dumux/multidimension/subproblemproperties.hh>

//#include "soiltestspatialparams1p2c.hh"
#include "richardsbuffertestspatialparams.hh"
#include <dumux/linear/amgbackend.hh>

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
//SET_PROP(RichardsTestProblem, WettingPhase)
//{
//private:
//    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//public:
//    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
//};
// Set fluid configuration
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OK<typename GET_PROP_TYPE(TypeTag, Scalar), false>);
              //Dumux::H2OXFluidSystem<TypeTag>);

// Set the grid type
//SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::YaspGrid<1, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 1> >);
//SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3> >);
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::UGGrid<2>);
//SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::YaspGrid<2>);
//SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Grid, Dune::ALUGrid< 2, 2, Dune::cube, Dune::nonconforming>);
//SET_TYPE_PROP(RichardsTestProblem, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::conforming>);

// explicitly set the fvgeometry as we currently use another one as in stable
//SET_TYPE_PROP(SoilRichardsTwoCBufferTestBoxProblem, FVElementGeometry, Dumux::MultiDimensionBoxFVElementGeometry<TypeTag>);
//SET_TYPE_PROP(SoilRichardsTwoCBufferTestCCProblem, FVElementGeometry, Dumux::MultiDimensionCCFVElementGeometry<TypeTag>);


// Set the problem property
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, Problem, Dumux::SoilRichardsTwoCBufferTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, SpatialParams, Dumux::RichardsBufferTestSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(SoilRichardsTwoCBufferTestProblem, ProblemEnableGravity, false);

// Use the AMG backend to allow parallel computation
SET_TYPE_PROP(SoilRichardsTwoCBufferTestProblem, LinearSolver, AMGBackend<TypeTag>);

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

    //typedef typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag) GlobalProblemTypeTag;
    //typedef typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager) CouplingManager;

public:
    SoilRichardsTwoCBufferTestProblem(TimeManager &timeManager, const GridView &gridView)
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
        Scalar rootRadius = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.rootRadius);
        if (sqrt((globalPos[0]*globalPos[0])+(globalPos[1]*globalPos[1])) < rootRadius + eps_)
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
    //    initial_(values, globalPos);
        //Scalar sw_ = GET_RUNTIME_PARAM(TypeTag,
        //                                Scalar,
        //                                BoundaryConditions.InitialSoilSaturation);
        //Scalar pc_ = MaterialLaw::pc(this->spatialParams().materialLawParams(globalPos),sw_);
        Scalar pw_ = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilPressure);
        //priVars[pressureIdx] = pnRef_ - pc_;
        priVars[pressureIdx] = pw_;

        priVars[massOrMoleFracIdx] = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilFracK);
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
        Scalar Vmax = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        SpatialParams.Vmax);
        Scalar Km = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        SpatialParams.Km);
        Scalar rootRadius = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        SpatialParams.rootRadius);
        Scalar c3D;
        if(useMoles)
            c3D = elemVolVars[scvIdx].moleFraction(1);
        else
            c3D = elemVolVars[scvIdx].massFraction(1);

        Scalar active_uptake = 0;
        active_uptake = 0;
        active_uptake =  Vmax*c3D*elemVolVars[scvIdx].density()
                                /(Km+c3D*elemVolVars[scvIdx].density());

        //std::cout<<"        "<<this->timeManager().time()<<" Active_uptake "<< active_uptake <<" "<<c3D
        //        <<" "<<std::endl;
        values[conti0EqIdx] = 0.0;
        values[transportEqIdx] = active_uptake;
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

        Scalar Vmax = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        SpatialParams.Vmax);
        Scalar Km = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        SpatialParams.Km);

        // create required scalar fields for the vtk output
        ScalarField& sourceP = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& sourceC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ratioCNum = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ratioCAna = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& relErrC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& relErrS = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& Neumann = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& Dirichlet = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& DiffCoeff = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& EffDiffCoeff = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& AnalyticalC = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& AnalyticalS = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& Distance = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ut_num = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& ut_ana = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& diffFlux = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& diffFluxMass = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& advFlux = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& storage = *(this->resultWriter().allocateManagedBuffer(numDofs));
        ScalarField& storageAna = *(this->resultWriter().allocateManagedBuffer(numDofs));

        sourceP = 0.0;
        sourceC = 0.0;
        ratioCNum = 0.0;
        ratioCAna = 0.0;
        relErrC = 0.0;
        relErrS = 0.0;
        Neumann = 0.0;
        Dirichlet = 0.0;
        DiffCoeff = 0.0;
        EffDiffCoeff = 0.0;
        AnalyticalC = 0.0;
        AnalyticalS = 0.0;
        Distance = 0.0;
        ut_ana = 0.0;
        ut_num = 0.0;
        diffFlux = 0.0;
        diffFluxMass = 0.0;
        advFlux = 0.0;
        storage = 0.0;
        storageAna = 0.0;

        Scalar Km_=GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Km);
        Scalar Vmax_=GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Vmax);
        Scalar rootRadius_=GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.rootRadius);
        Scalar Frac0=GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.InitialSoilFracK);

        Scalar Uptake_num, Uptake_ana;
        Uptake_num = 0.0;
        Uptake_ana = 0.0;

        // iterate over all elements
        for (const auto& element : elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            ElementBoundaryTypes bcTypes;
            bcTypes.update(*this, element, fvGeometry);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this, element, fvGeometry, false /* oldSol? */);

            FluxVariables fluxVars;
            fluxVars.update(*this, element, fvGeometry, 0, elemVolVars, false /* oldSol? */);

            //FluxVariables fluxVars;
            //fluxVars.update(*this, element, fvGeometry);

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
                    sourceC[dofGlobalIdx] += values[transportEqIdx] * fvGeometry.subContVol[scvIdx].volume
                                            * this->boxExtrusionFactor(element, fvGeometry, scvIdx);
                    ratioCNum[dofGlobalIdx] += elemVolVars[scvIdx].massFraction(1)/Frac0;
                    //const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
                    if (bcTypes.hasNeumann())
                        Neumann[dofGlobalIdx] = 1;
                    if (bcTypes.hasDirichlet())
                        Dirichlet[dofGlobalIdx] = 1;
                    DiffCoeff[dofGlobalIdx] += elemVolVars[scvIdx].diffCoeff();
                    EffDiffCoeff[dofGlobalIdx] += elemVolVars[scvIdx].effDiffCoeff();

                    Scalar x = pow((pow(element.geometry().center()[0],2) + pow(element.geometry().center()[1],2)),0.5);// 2D
                    //Scalar x = pow((pow(element.geometry().corner(0)[0],2) + pow(element.geometry().corner(0)[1],2)),0.5);
                    //Scalar x = pow(pow(element.geometry().center()[0],2),0.5);
                    //Scalar x = pow(pow(element.geometry().corner(0)[0],2),0.5);
                    //std::cout << element.geometry().center() << std::endl;
                    if (x<rootRadius_)
                        x=rootRadius_;
                    Distance[dofGlobalIdx] = x;

                    Scalar lambda = Vmax_*rootRadius_/(elemVolVars[scvIdx].effDiffCoeff()*Km_);
                            lambda /=(elemVolVars[scvIdx].saturation(phaseIdx)* elemVolVars[scvIdx].porosity()+elemVolVars[scvIdx].buffer());

                    Scalar L = lambda/2*log(4*exp(-0.5772)*elemVolVars[scvIdx].effDiffCoeff()*pow(rootRadius_,(-2))*(this->timeManager().time()+this->timeManager().timeStepSize())+1);

                    AnalyticalC[dofGlobalIdx] += (Cinf_-Cinf_*lambda/(1+Cinf_dl+L+sqrt(4*Cinf_dl+pow((1-Cinf_dl+L),2)))*
                                                boost::math::expint(1,pow(x,2)/(4*elemVolVars[scvIdx].effDiffCoeff()*(this->timeManager().time()+this->timeManager().timeStepSize()))))
                                                /elemVolVars[scvIdx].density();

                    ratioCAna[dofGlobalIdx] += AnalyticalC[dofGlobalIdx]/Frac0;
                    if (AnalyticalC[dofGlobalIdx] != 0)
                        relErrC[dofGlobalIdx]=std::abs((elemVolVars[scvIdx].massFraction(transportCompIdx)-AnalyticalC[dofGlobalIdx])/AnalyticalC[dofGlobalIdx]);

                    //if (sourceC[dofGlobalIdx] != 0)
                    AnalyticalS[dofGlobalIdx] += 2*Vmax_*Cinf_dl/(1+Cinf_dl+L+sqrt(4*Cinf_dl+pow((1-Cinf_dl+L),2)));
                    //AnalyticalS[dofGlobalIdx] +=-2*Vmax_*Cinf_dl/
                    //                                        (1+Cinf_dl+(lambda/2*log(4*exp(-0.5772)*elemVolVars[scvIdx].effDiffCoeff()*pow(rootRadius_,(-2))*this->timeManager().time()+1))
                    //                                    +sqrt(4*Cinf_dl+pow((1-Cinf_dl+(lambda/2*log(4*exp(-0.5772)*elemVolVars[scvIdx].effDiffCoeff()*pow(rootRadius_,(-2))*this->timeManager().time()+1))),2)));
                    if (AnalyticalS[dofGlobalIdx] != 0)
                    relErrS=std::abs((sourceC[dofGlobalIdx])/AnalyticalS[dofGlobalIdx]);

                    ut_num[dofGlobalIdx] = Vmax*elemVolVars[scvIdx].massFraction(1)*elemVolVars[scvIdx].density()
                                /(Km+elemVolVars[scvIdx].massFraction(1)*elemVolVars[scvIdx].density());

                    ut_ana[dofGlobalIdx] = Vmax*AnalyticalC[dofGlobalIdx]*elemVolVars[scvIdx].density()
                                /(Km+AnalyticalC[dofGlobalIdx]*elemVolVars[scvIdx].density());

                    if (Uptake_num < ut_num[dofGlobalIdx])
                        Uptake_num = ut_num[dofGlobalIdx];
                    if (Uptake_ana < ut_ana[dofGlobalIdx])
                        Uptake_ana = ut_ana[dofGlobalIdx];

                    diffFlux[dofGlobalIdx] = -(fluxVars.moleFractionGrad(transportCompIdx)*fluxVars.face().normal)
                                            *fluxVars.porousDiffCoeff() * fluxVars.molarDensity();
                    diffFlux[dofGlobalIdx] =std::abs(diffFlux[dofGlobalIdx]);


                    diffFluxMass[dofGlobalIdx] = -(fluxVars.massFractionGrad(transportCompIdx)*fluxVars.face().normal)
                                            *fluxVars.porousDiffCoeff()*elemVolVars[scvIdx].density();
                    diffFluxMass[dofGlobalIdx] =std::abs(diffFluxMass[dofGlobalIdx]);

                    advFlux[dofGlobalIdx] = fluxVars.volumeFlux(phaseIdx) * fluxVars.molarDensity() *
                                            elemVolVars[scvIdx].moleFraction(transportCompIdx) / elemVolVars[scvIdx].viscosity();
                    storage[dofGlobalIdx] = elemVolVars[scvIdx].density() * elemVolVars[scvIdx].massFraction(transportCompIdx) *
                                            (elemVolVars[scvIdx].saturation(phaseIdx)*elemVolVars[scvIdx].porosity()+elemVolVars[scvIdx].buffer());
                    storageAna[dofGlobalIdx] = elemVolVars[scvIdx].density() * AnalyticalC[dofGlobalIdx] *
                                            (elemVolVars[scvIdx].saturation(phaseIdx)*elemVolVars[scvIdx].porosity()+elemVolVars[scvIdx].buffer());

                //    std::cout <<" TIME !!! "<< this->timeManager().time() <<" "<< elemVolVars[scvIdx].massFraction(1)
                //            <<" "<< AnalyticalC[dofGlobalIdx]<<" "<<ut_num[dofGlobalIdx]
                //            <<" "<<ut_ana[dofGlobalIdx]<<" "<<AnalyticalS[dofGlobalIdx]<<" "<< x<<  std::endl;
                //    std::cout<<"      lambda "<<lambda<<" L "<< L <<" "<<
                //            (1+Cinf_dl+L+sqrt(4*Cinf_dl+pow((1-Cinf_dl+L),2)))<<" "<<
                //            boost::math::expint(1,pow(x,2)/(4*elemVolVars[scvIdx].effDiffCoeff()*this->timeManager().time()))<<" "<<std::endl;
                //    std::cout <<"               "<< elemVolVars[scvIdx].density()<<" "<<elemVolVars[scvIdx].saturation(phaseIdx)<<" "
                //                <<elemVolVars[scvIdx].porosity()<<" "<<elemVolVars[scvIdx].buffer()<<" "<<fluxVars.porousDiffCoeff()<<std::endl;
                }
            }

    //        //print flux.vols
    //        int fIdxInner = 0;
    //        for (const auto& intersection : intersections(this->gridView(), element))
    //            {
    //                int fIdx = intersection.indexInInside();
	//				if (intersection.boundary())
    //                {
    //
    //                    FluxVariables fluxVars;
    //                    fluxVars.update(*this,
    //                                    element,
    //                                    fvGeometry,
    //                                    fIdx,
    //                                    elemVolVars,true);
    //                    std::cout<< fIdxInner<< " Results: moleFractionGrad_ bound"<<fluxVars.moleFractionGrad(transportCompIdx) << std::endl;
    //                    //scvfFluxes[fIdx] = fluxVars.volumeFlux(phaseIdx);
    //                }
    //            }
    //    std::cout<< " " << std::endl;
        }

        const auto totalSourceP = std::accumulate(sourceP.begin(), sourceP.end(), 0);
        const auto totalSourceC = std::accumulate(sourceC.begin(), sourceC.end(), 0);

        std::cout << "Integrated mass source (3D): " << totalSourceP << std::endl;
        std::cout << "Integrated concentration source (3D): " << totalSourceC << std::endl;

        // attach data to the vtk output
        this->resultWriter().attachDofData(sourceP, "water_uptake(kg/s)", isBox);
        this->resultWriter().attachDofData(sourceC, "uptake_num(kg/s)", isBox);
        this->resultWriter().attachDofData(ratioCNum, "ratioCNum", isBox);
        this->resultWriter().attachDofData(ratioCAna, "ratioCAna", isBox);
        this->resultWriter().attachDofData(Neumann, "BC_Neumann", isBox);
        this->resultWriter().attachDofData(Dirichlet, "BC_Dirichlet", isBox);
        this->resultWriter().attachDofData(DiffCoeff, "DiffCoeff", isBox);
        this->resultWriter().attachDofData(EffDiffCoeff, "EffDiffCoeff", isBox);
        this->resultWriter().attachDofData(AnalyticalC, "AnalyticalC", isBox);
        this->resultWriter().attachDofData(relErrC, "relErrC", isBox);
        this->resultWriter().attachDofData(relErrS, "relErrS", isBox);
        this->resultWriter().attachDofData(AnalyticalS, "uptake_analytical(kg/s)", isBox);
        this->resultWriter().attachDofData(Distance, "Distance", isBox);
        this->resultWriter().attachDofData(ut_num, "ut_num", isBox);
        this->resultWriter().attachDofData(ut_ana, "ut_ana", isBox);
        //this->resultWriter().attachDofData(advFlux, "advFlux", isBox);
        //this->resultWriter().attachDofData(diffFlux, "diffFlux", isBox);
        this->resultWriter().attachDofData(diffFluxMass, "diffFluxMass", isBox);
        this->resultWriter().attachDofData(storage, "storage", isBox);
        this->resultWriter().attachDofData(storageAna, "storageAna", isBox);

        logFile_.open(this->name() + ".log", std::ios::app);
        logFile_ << "time = " << this->timeManager().time()+this->timeManager().timeStepSize() << " uptake_num = " << ut_num[0]
                << " uptake_ana = " << ut_ana[0] << std::endl;
        logFile_.close();
    }

    //! Set the coupling manager
//    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
//    { couplingManager_ = cm; }

    //! Get the coupling manager
//    const CouplingManager& couplingManager() const
//    { return *couplingManager_; }

private:
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        Scalar pw_ = GET_RUNTIME_PARAM(TypeTag,
                                        Scalar,
                                        BoundaryConditions.InitialSoilPressure);
        //priVars[pressureIdx] = pnRef_ - pc_;
        priVars[pressureIdx] = pw_;

        //priVars[massOrMoleFracIdx] = 0;
        //if((globalPos[0] < this->bBoxMin()[0] + eps_))
        priVars[massOrMoleFracIdx]= GET_RUNTIME_PARAM(TypeTag,
                                       Scalar,
                                       BoundaryConditions.InitialSoilFracK);
        //Scalar sw_ = GET_RUNTIME_PARAM(TypeTag,
        //                                Scalar,
        //                                BoundaryConditions.InitialWaterContent)/
        //             GET_RUNTIME_PARAM(TypeTag,
        //                                Scalar,
        //                                SpatialParams.Porosity);
        //Scalar pc_ = MaterialLaw::pc(this->spatialParams().materialLawParams(globalPos),sw_);
        //std::cout << "Pressure  " << pnRef_-pc_ << std::endl;
        //std::cin.ignore(100000, '\n');
};
    std::ofstream logFile_;
    const Scalar eps_ = 1e-9;
    Scalar episodeTime, temperature_, pnRef_, pc_, sw_, pw_;
    std::string name_;
//    std::shared_ptr<CouplingManager> couplingManager_;
    //Scalar DiffCoef_;
};

} //end namespace

#endif
