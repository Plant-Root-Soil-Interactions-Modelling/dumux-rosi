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
 * \brief A test problem for the one-phase root model:
 * Sap is flowing through a 1d network root xylem.
 */
#ifndef DUMUX_ROOT_PROBLEM_HH
#define DUMUX_ROOT_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/io/gnuplotinterface.hh>

#include "rootspatialparams_.hh"

namespace Dumux {

/*!
 * \brief a root problem file
 */
template <class TypeTag>
class RootProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    // copy some indices for convenience
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    RootProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<CouplingManager> couplingManager,
                const GlobalPosition& domainSize)
    : ParentType(fvGridGeometry)
    , couplingManager_(couplingManager)
    , domainSize_(domainSize)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_1d";
        criticalCollarPressure_ = getParam<Scalar>("BoundaryConditions.CriticalCollarPressure", -1.4e6); // in Pa
        useCyclicTranspiration_ = getParam<bool>("BoundaryConditions.UseCyclicTranspiration", true);
        dailyTranspirationRate_ = getParam<Scalar>("BoundaryConditions.DailyTranspirationRate", 5.33); // mm/day
        dailyTranspirationRate_ *= domainSize_[0]*domainSize_[1]/86400; // kg/s

        // push first plot point
        tPotPlot_.push_back(0);
        tActPlot_.push_back(0);
        tCumulPlot_.push_back(0);
        timePlot_.push_back(0);
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here makes extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto radius = this->spatialParams().rootParams(element).radius;
        return M_PI*radius*radius;
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
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Return the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return 273.15 + 10.0; }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The global position
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }


    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    template<class ElementVolumeVariables>
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        // set transpiration rate at the root collar
        const auto globalPos = scvf.center();
        if (globalPos[2] + eps_ > this->fvGridGeometry().bBoxMax()[2])
        {
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];

            // if the plant has water stress reduce the transpiration rate (imposing Dirichlet boundary condition weakly)
            const auto p = volVars.pressure(0);
            if (p <= criticalCollarPressure_*0.5)
            {
                const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
                const Scalar Kx = this->spatialParams().Kx(eIdx);
                const auto dist = (globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm();
                const Scalar criticalTranspiration = volVars.density(0)*Kx*(volVars.pressure(0) - criticalCollarPressure_)/dist;
                values[Indices::conti0EqIdx] = std::min(potentialTranspirationRate(), criticalTranspiration);
            }
            else
                values[Indices::conti0EqIdx] = potentialTranspirationRate();

            values /= volVars.extrusionFactor() * scvf.area(); // convert from kg/s to kg/(s*m^2)
        }

        return values;
    }

    //! Compute potential transpiration rate in kg/s
    Scalar potentialTranspirationRate() const
    {
        // possibly make the transpiration rate dependent on the current root length for growth
        const auto reductionFactor = 1.0;// rootLength_ < 60 ? rootLength_/60 : 1.0;
        if (useCyclicTranspiration_)
            return reductionFactor*(dailyTranspirationRate_*std::sin(time_*2*M_PI / 86400 - M_PI/2.0) + dailyTranspirationRate_);
        else
            return reductionFactor*dailyTranspirationRate_;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

     /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
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
     * \param pointSource A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub-control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // compute source at every integration point
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];

        const auto lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = this->spatialParams().Kr(lowDimElementIdx);
        const Scalar rootRadius = this->spatialParams().radius(lowDimElementIdx);

        // relative soil permeability
        const auto krel = 1.0;//this->couplingManager().relPermSoil(pressure3D);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        const auto density = 1000;
        const Scalar sourceValue = 2* M_PI *krel*rootRadius * Kr *(pressure3D - pressure1D)*density;
        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(1e5); }

    // \}

    //! Called after every time step
    //! Output the total global exchange term
    void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars) const
    {
        NumEqVector source(0.0);
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                source += pointSources;
            }
        }

        std::cout << "Global integrated source (root): " << source << " (kg/s) / "
                  <<                           source*3600*24*1000 << " (g/day)" << '\n';
    }

    void updateRootVolume()
    {
        rootVolume_ = 0.0, rootLength_ = 0.0;
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto& params = this->spatialParams().rootParams(element);
            if (params.rootId >= 0) // exclude shoot
            {
                const auto radius = params.radius;
                const auto length = element.geometry().volume();
                rootVolume_ += length*M_PI*radius*radius;
                rootLength_ += length;
            }
        }
    }

    //! compute the actual transpiration rate
    Scalar computeActualTranspirationRate(const SolutionVector& sol, const GridVariables& gridVars) const
    {
        NumEqVector transpirationRate(0.0);
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (const auto& scvf : scvfs(fvGeometry))
                if (scvf.boundary())
                    transpirationRate += this->neumann(element, fvGeometry, elemVolVars, scvf)*scvf.area()*elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        }

        std::cout << "Actual transpiration rate:       " << transpirationRate << " (kg/s) / "
                  << transpirationRate[0]*86400*1000 << " (g/day) / "
                  << transpirationRate[0]/domainSize_[0]/domainSize_[1]*86400 << " (mm/day)\n"
                  << "Potential transpiration rate:    " << potentialTranspirationRate() << " (kg/s) / "
                  << potentialTranspirationRate()*86400*1000 << " (g/day) / "
                  << potentialTranspirationRate()/domainSize_[0]/domainSize_[1]*86400 << " (mm/day)\n";

        return transpirationRate[0];
    }

    //! plot the transpiration rate
    void plotTranspirationRate(double dt, const SolutionVector& sol, const GridVariables& gridVars)
    {
        transpirationPlot_.resetPlot();

        transpirationPlot_.setXlabel("time [days]");
        transpirationPlot_.setYlabel("transpiration rate [mm/days]");
        transpirationPlot_.setOption("set y2label \"cumulative transpiration [mm]\"");

        const auto tAct = computeActualTranspirationRate(sol, gridVars)/domainSize_[0]/domainSize_[1]*86400; // convert to mm/day
        const auto tPot = potentialTranspirationRate()/domainSize_[0]/domainSize_[1]*86400; // convert to mm/day

        tPotPlot_.push_back(tPot);
        tActPlot_.push_back(tAct);
        tCumulPlot_.push_back(tCumulPlot_.back() + tAct*dt/86400);
        timePlot_.push_back(time_/86400);

        transpirationPlot_.addDataSetToPlot(timePlot_, tPotPlot_, "potentialtranspiration", "with lines axes x1y1 lw 2 lc rgb 'black'");
        transpirationPlot_.addDataSetToPlot(timePlot_, tActPlot_, "actualtranspiration", "with lines axes x1y1 lw 3");
        transpirationPlot_.addDataSetToPlot(timePlot_, tCumulPlot_, "cumulativetranspiration", "with lines axes x1y2 lw 3");
        transpirationPlot_.setOption("set ytics nomirror");
        transpirationPlot_.setOption("set y2tics");

        transpirationPlot_.setYRange(0, 12);
        transpirationPlot_.setOption("set autoscale x");
        transpirationPlot_.setOption("set y2range [0 : 30]");

        transpirationPlot_.setOption("set title \"Plant transpiration\"");
        transpirationPlot_.plot("transpiration");
    }

    //! set the current time for evaluation of time-dependent boundary conditions
    void setTime(Scalar t)
    { time_= t; }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar criticalCollarPressure_, dailyTranspirationRate_;
    bool useCyclicTranspiration_;
    Scalar time_;

    GnuplotInterface<Scalar> transpirationPlot_;
    std::vector<Scalar> tPotPlot_, tActPlot_, tCumulPlot_, timePlot_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
    const GlobalPosition domainSize_; // soil domain size to compute the correct transpiration rates

    Scalar rootVolume_ = 0.0;
    Scalar rootLength_ = 0.0;
};

} //end namespace Dumux

#endif
