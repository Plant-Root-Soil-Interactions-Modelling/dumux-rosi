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
 * \brief Base class for all fully implicit Rootsystem problems
 */
#ifndef DUMUX_ROOTSYSTEM_PROBLEM_1P2C_HH
#define DUMUX_ROOTSYSTEM_PROBLEM_1P2C_HH

#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/growth/gridgrowth.hh>
//#include <dumux/dumux/porousmediumflow/1p2c/implicit/properties.hh>
#include "properties1p2c.hh"

namespace Dumux
{
/*!
 * \ingroup RootsystemModel
 * \ingroup ImplicitBaseProblems
 * \brief Base class for all fully implicit Rootsystem problems
 *
 * For a description of the Rootsystem model, see Dumux::RootsystemModel
 */
template<class TypeTag>
class RootsystemOnePTwoCProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    //add
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables

        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
//#if NONISOTHERMAL
//        temperatureIdx = Indices::temperatureIdx
//#endif
    };
    enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
//#if NONISOTHERMAL
//        energyEqIdx = Indices::energyEqIdx
//#endif
    };
    //enum { growingGrid = GET_PROP_VALUE(TypeTag, GrowingGrid) };

    //typedef ImplicitGridGrowth<TypeTag, growingGrid> GridGrowthModel;

    typedef typename GridView::template Codim<0>::Entity Element;
//    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    enum { growingGrid = GET_PROP_VALUE(TypeTag, GrowingGrid) };
    typedef ImplicitGridGrowth<TypeTag, growingGrid> GridGrowthModel;

    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief The constructor.
     *
     * The overloaded class must allocate all data structures
     * required, but _must not_ do any calls to the model, the
     * jacobian assembler, etc inside the constructor.
     *
     * If the problem requires information from these, the
     * ImplicitProblem::init() method be overloaded.
     *
     * \param timeManager The TimeManager which keeps track of time
     * \param gridView The GridView used by the problem.
     */
    RootsystemOnePTwoCProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //Add
        //FluidSystem::init();

        waterStress_ = false;
        switchBC_ = false;

        // if we are calculating on an growing grid get the grid growth model
        if (growingGrid)
            gridGrowth_ = Dune::make_shared<GridGrowthModel>(*static_cast<Implementation*>(this));

        //Add
    //    //stating in the console whether mole or mass fractions are used
    //    if(useMoles)
    //    {
    //        std::cout<<"problem uses mole fractions"<<std::endl;
    //    }
    //    else
    //    {
    //        std::cout<<"problem uses mass fractions"<<std::endl;
    //    }

    }

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        ParentType::init();

        if (growingGrid)
            gridGrowth().init();
    }

    /*!
     * \brief Called by the time manager after the time integration.
     */
    void preTimeStep()
    {
        ParentType::preTimeStep();
        // If gridGrowth is used, this method grows the grid.
        // Remeber to call the parent class function if this is overwritten
        // on a lower problem level when using an adaptive grid
        if (growingGrid && this->timeManager().timeStepIndex() > 0){
            this->gridGrowth().growGrid();
            this->resultWriter().gridChanged();
        }
        preSol_ = this->model().curSol();
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
            values.setOutflow(transportEqIdx);
            Scalar criticalCollarPressure = GET_RUNTIME_PARAM(TypeTag,
                                         Scalar,
                                         BoundaryConditions.CriticalCollarPressure);
            //get element index Eid of root segment at root colar
            int Eid=-1;
            for (const auto& element : elements(this->gridView()))
            {
                Eid ++;
                auto posZ = std::max(element.geometry().corner(0)[2],element.geometry().corner(1)[2]);
                if (posZ + eps_ > this->bBoxMax()[2])
                    break;
            }
            if (this->timeManager().time()>=0)
            {
                if ((preSol_[Eid][conti0EqIdx] < criticalCollarPressure ))
                {
                    std::cout<<"Collar pressure: "<<preSol_[Eid][conti0EqIdx]<<" < " <<criticalCollarPressure<<"\n";
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
        else
            values.setAllNeumann();
    }

    void setWaterStress()
    {
//       const SolutionVector& curSol = this->model().curSol();
//       FVElementGeometry fvGeometry;
//       bool waterStress_old = waterStress_;
//       double Hcrit = GET_RUNTIME_PARAM(TypeTag,
//                                        Scalar,
//                                        BoundaryConditions.CriticalCollarPressure);

//       ElementIterator eIt = this->gridView().template begin<0>();
//       ElementIterator eEndIt = this->gridView().template end<0>();
//       for (; eIt != eEndIt; ++eIt){

//           fvGeometry.update(this->gridView(), *eIt);

//           for (unsigned int i=0; i<eIt->geometry().corners();i++){
//               GlobalPosition globalPos = eIt->geometry().corner(i);

//               if ( globalPos[2] < eps_ ) {
//                   if( curSol[0] <=  Hcrit){
//                       waterStress_ = true;
//                   }
//                   else
//                       waterStress_ = false;
//               }
//           }
//       }
//       if (waterStress_ !=  waterStress_old)
//           switchBC_= true;
    }

    bool getWaterStress() const
    {
        return waterStress_;
    }

    void resetSwitchBC()
    {
       switchBC_= false;
    }

    bool getSwitchBC()
    {
       return switchBC_;
    }

     /*!
      * \brief Returns growth model used for the problem.
      */
     GridGrowthModel& gridGrowth()
     {
         return *gridGrowth_;
     }

     /*!
      * \brief Returns growth model used for the problem.
      */

     const GridGrowthModel& gridGrowth() const
     {
         return *gridGrowth_;
     }

    /*!
     * \brief Capability to introduce problem-specific routines at the
     * beginning of the grid growth
     *
     * Function is called at the beginning of the standard grid
     * modification routine, GridGrowth::growGrid(.
     */
    void preGrowth()
    {}

    /*!
     * \brief Capability to introduce problem-specific routines after grid growth
     *
     * Function is called at the end of the standard grid
     * modification routine, GridGrowth::growGrid() , to allow
     * for problem-specific output etc.
     */
    void postGrowth()
    {}

private:
    bool waterStress_ = false;
    bool switchBC_ = false;
    const Scalar eps_ = 1e-6;
    Dune::shared_ptr<GridGrowthModel> gridGrowth_;
    SolutionVector preSol_;
};

}

#include <dumux/implicit/growth/gridgrowthpropertydefaults.hh>

#endif
