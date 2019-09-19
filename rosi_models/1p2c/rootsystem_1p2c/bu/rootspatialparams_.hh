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
 * \ingroup OneTests
 * \brief The spatial parameters class blood flow problem
 */
#ifndef DUMUX_ROOT_SPATIALPARAMS_HH
#define DUMUX_ROOT_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/growth/rootparameters.hh>

#include <RootSystem.h>

namespace Dumux {

/*!
 * \brief Root spatial params
 */
template<class FVGridGeometry, class Scalar>
class RootSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar, RootSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = RootSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using RootParams = GrowthModule::RootParams<Scalar>;
    using DGFParamIndices = GrowthModule::DGFParamIndices;
public:
    // export permeability type
    using PermeabilityType = Scalar;

    RootSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        porosity_ = getParam<Scalar>("Root.SpatialParams.Porosity", 0.4);
        higherRadialConductivities_ = getParam<bool>("RootSystem.EnableHigherRadialConductivities", false);

        // output root conductivities
        if (getParam<bool>("Output.PlotRootConductivities", false))
        {
            // plot for primary and secondary roots
            plotRootConductivities_(0);
            plotRootConductivities_(1);
        }
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \param ipGlobal The integration point
     * \note Kx has units [m^4/(Pa*s)] so we have to divide by the cross-section area
     *       and multiply with a characteristic viscosity
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    { return permeability(element); }

    //! simpler interface
    PermeabilityType permeability(const Element& element) const
    {
        const Scalar r = rootParams(element).radius;
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        return this->Kx(eIdx) / (M_PI*r*r) * Components::SimpleH2O<Scalar>::liquidViscosity(285.15, 1e5);
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    //! Return the root parameters
    const RootParams& rootParams(const Element &element) const
    {
        const auto& elementMapper = this->fvGridGeometry().elementMapper();
        const auto eIdx = elementMapper.index(element);
        int rootId = rootId_[eIdx];
        if (rootId < 0)
            return shootParams_;
        else
            return rootParams_[rootId];
    }

    Scalar previousLength(std::size_t eIdx) const
    {
        return previousLength_[eIdx];
    }

    Scalar age(std::size_t eIdx) const
    {
        return ages_[eIdx];
    }

    Scalar radius(std::size_t eIdx) const
    {
        int rootId = rootId_[eIdx];
        if (rootId < 0)
            return shootParams_.radius;
        else
            return rootParams_[rootId].radius;
    }

    Scalar Kr(std::size_t eIdx) const
    {
        return Kr_[eIdx];
    }

    Scalar Kx (std::size_t eIdx) const
    {
        return Kx_[eIdx];
    }

    //! Access to the radii vector for output
    const std::vector<Scalar>& radii() const
    { return radii_; }

    //! Access to the orders vector for output
    const std::vector<Scalar>& orders() const
    { return orders_; }

    //! Access to the ages vector for output
    const std::vector<Scalar>& ages() const
    { return ages_; }

    //! Access to the axial conductivities vector for output
    const std::vector<Scalar>& axialConductivities() const
    { return Kx_; }

    //! Access to the radial conductivities vector for output
    const std::vector<Scalar>& radialConductivities() const
    { return Kr_; }

    //! Read parameters from the a crootbox rootsystem
    void initParameters(const CRootBox::RootSystem& rs)
    {
        // be careful: the total number of roots and the grown number of roots is not the same
        rootParams_.resize(rs.getNumberOfRoots(/*allRoots=*/true));
        // this gets the root that already have at least one segment
        const auto roots = rs.getRoots();
        for (auto&& r : roots)
        {
            auto& rp = rootParams_[r->id];

            // compute the root order
            rp.order = 0;
            auto* r_ = r;
            while (r_->parent != nullptr)
            {
                ++rp.order;
                r_=r_->parent;
            }

            rp.radius = r->param.a * 0.01; // conversion from cm to m
            rp.rootId = r->id;
            rp.plantId = 0;
        }

        // set the shootParams
        shootParams_.order = 0;
        shootParams_.radius = 0.01;
        shootParams_.rootId = -1;
        shootParams_.plantId = 0;

        const auto shootSegments = rs.getShootSegments();
        const auto segmentRoots = rs.getSegmentsOrigin();
        const auto segmentAges = rs.getNETimes();
        const auto& gridView = this->fvGridGeometry().gridView();
        rootId_.resize(gridView.size(0));
        radii_.resize(gridView.size(0));
        orders_.resize(gridView.size(0));
        previousLength_.resize(gridView.size(0));
        ages_.resize(gridView.size(0));
        Kx_.resize(gridView.size(0));
        Kr_.resize(gridView.size(0));
        for (const auto& element : elements(gridView))
        {
            const auto& elementMapper = this->fvGridGeometry().elementMapper();
            const auto eIdx = elementMapper.index(element);

            // initialize the previous length to the current length
            previousLength_[eIdx] = element.geometry().volume();

            // subtract the number of shoot segments to get the normal segment index
            const int segmentIdx = eIdx - shootSegments.size();
            // shoot segments
            if (segmentIdx < 0)
            {
                rootId_[eIdx] = -1;
                ages_[eIdx] = rs.getSimTime();
                const auto& rp = rootParams(element);
                radii_[eIdx] = rp.radius;
                orders_[eIdx] = rp.order;
                Kx_[eIdx] = 1e-10; // very high value to not contribute to resistance
                Kr_[eIdx] = 0.0; // no exchange with the soil
            }
            // regular segements
            else
            {
                // we hope that the crootbox indices and the dune indices on initialization are the same!
                rootId_[eIdx] = segmentRoots[segmentIdx]->id;
                ages_[eIdx] = std::max(0.0, rs.getSimTime()-segmentAges[segmentIdx]);
                const auto& rp = rootParams_[rootId_[eIdx]];
                radii_[eIdx] = rp.radius;
                orders_[eIdx] = rp.order;

                // compute age dependent conductivities
                Kx_[eIdx] = computeAxialConductivity_(orders_[eIdx], ages_[eIdx]);
                Kr_[eIdx] = computeRadialConductivity_(orders_[eIdx], ages_[eIdx]);
            }
        }
    }

    //! Output an analysis of the root system
    void analyseRootSystem(const CRootBox::RootSystem& rs) const
    {
        // compute some statistics of the initial root system
        Scalar totalLength = 0.0, totalLengthTop = 0.0, totalLengthBottom = 0.0;
        Scalar totalLengthPrimary = 0.0, totalLengthSecondary = 0.0;
        Scalar rootVolume = 0.0;
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            if (rootId_[eIdx] >= 0) // exclude shoot
            {
                const auto geo = element.geometry();
                const auto length = geo.volume();
                const auto r = radii_[eIdx];
                totalLength += length;
                rootVolume += length*M_PI*r*r;

                if (geo.center()[2] > -0.42)
                    totalLengthTop += length;
                else
                    totalLengthBottom += length;

                if (orders_[eIdx] == 0)
                    totalLengthPrimary += length;
                else
                    totalLengthSecondary += length;
            }
        }

        std::cout << ".........................................................\n"
                  << "-- Root system age:            " << rs.getSimTime() << " days\n"
                  << "-- Total length:               " << totalLength << " m\n"
                  << "-- Total length (top 42 cm):   " << totalLengthTop << " m\n"
                  << "-- Total length (below 42 cm): " << totalLengthBottom << " m\n"
                  << "-- Total length (primary):     " << totalLengthPrimary << " m\n"
                  << "-- Total length (secondary):   " << totalLengthSecondary << " m\n"
                  << "-- Total volume:               " << rootVolume << " m³\n"
                  << ".........................................................\n";
    }

    //! update the root parameters that aren't set by growth directly
    void updateParameters(const CRootBox::RootSystem& rs)
    {
        // update radii for output
        const auto& gridView = this->fvGridGeometry().gridView();
        radii_.resize(gridView.size(0));
        orders_.resize(gridView.size(0));
        Kx_.resize(gridView.size(0));
        Kr_.resize(gridView.size(0));
        for (const auto& element : elements(gridView))
        {
            const auto& elementMapper = this->fvGridGeometry().elementMapper();
            const auto eIdx = elementMapper.index(element);
            const auto& rp = rootParams(element);
            radii_[eIdx] = rp.radius;
            orders_[eIdx] = rp.order;

            // compute age dependent conductivities
            Kx_[eIdx] = computeAxialConductivity_(orders_[eIdx], ages_[eIdx]);
            Kr_[eIdx] = computeRadialConductivity_(orders_[eIdx], ages_[eIdx]);
        }
    }

    //! resize the segment parameters
    void resizeSegmentParams(std::size_t size)
    {
        rootId_.resize(size);
        previousLength_.resize(size);
        ages_.resize(size);
    }

    //! set root index (needed by growth if eIdx of an element changed after growth)
    void setAge(unsigned int eIdx, double age)
    { ages_[eIdx] = age; }

    //! set root index (needed by growth if eIdx of an element changed after growth)
    void setRootId(unsigned int eIdx, int rootId)
    { rootId_[eIdx] = rootId; }

    //! set the previous length
    void setPreviousLength(unsigned int eIdx, Scalar length)
    { previousLength_[eIdx] = length; }

    //! Resize root params to add more roots
    void rootParamsResize(std::size_t size)
    { rootParams_.resize(size); }

    //! Return a reference to the root parameters of a root id to change them (called from growth module)
    RootParams& getRootParams(int rootId)
    { return rootParams_[rootId]; }

private:

    void plotRootConductivities_(int order)
    {
        const std::string type = order == 0 ? "primary root" : "secondary root";
        GnuplotInterface<double> gnuplot;
        gnuplot.resetPlot();
        gnuplot.setXlabel("age [days]");
        gnuplot.setYlabel("segment radial cond. [cm/hPa/days]");
        gnuplot.setOption("set y2label \"segment axial cond. [cm^4/hPa/days]\"");

        constexpr int size = 300;
        std::vector<double> age(size);
        std::vector<double> kr(size);
        std::vector<double> highkr(size);
        std::vector<double> kx(size);

        for (int i = 0; i <= size-1; ++i)
        {
            age[i] = double(i)/double(size-1)*70.0;
            kr[i] = computeRadialConductivity_(order, age[i])/1e-4*86400; // from m/Pa/s to cm/hPa/day
            highkr[i] = kr[i]*10;
            kx[i] = computeAxialConductivity_(order, age[i])/1e-10*86400; // from m^4/Pa/s to cm^4/hPa/day
        }
        gnuplot.setOption("set y2tics");
        gnuplot.setOption("set ytics nomirror");
        gnuplot.setXRange(0, 70);

        if (order == 0)
        {
            gnuplot.setYRange(0, 2e-4);
            gnuplot.setOption("set y2range [0 : 5]");
        }
        else if (order == 1)
        {
            gnuplot.setYRange(0, 2e-3);
            gnuplot.setOption("set y2range [0 : 2e-3]");
        }
        else
        {
            gnuplot.setOption("set autoscale y");
            gnuplot.setOption("set autoscale y2");
        }

        gnuplot.setOption("set format y \"%.1e\"");
        gnuplot.setOption("set format y2 \"%.1e\"");

        gnuplot.addDataSetToPlot(age, kr, "kradial" + std::to_string(order), "with lines axes x1y1 lw 3");
        gnuplot.addDataSetToPlot(age, kx, "kaxial" + std::to_string(order), "with lines axes x1y2 lw 3");
        if (order == 1)
            gnuplot.addDataSetToPlot(age, highkr, "kradialhigh" + std::to_string(order), "with lines axes x1y1 lw 3 dt '-'");

        gnuplot.setOption("set title \"Root segment conductivities " + type + "\"");
        gnuplot.plot("rootconductivities" + std::to_string(order));
    }

    //! compute the radial conductivity (m/Pa/s) given the segment age in days
    double computeRadialConductivity_(int order, double age)
    {
        constexpr double conversion = 1e-4/86400; // from cm/hPa/day to m/Pa/s
        if (order == 0)
        {
            if (age < 5)
                return 1.8e-4*conversion;
            else if (age >= 5 && age < 10)
                return (-(age-5)*1.2e-4/5.0 + 1.8e-4)*conversion;
            else if (age >= 10 && age < 15)
                return 0.6e-4*conversion;
            else if (age >= 15 && age < 20)
                return (-(age-15)*0.42e-4/5.0 + 0.6e-4)*conversion;
            else
                return 0.18e-4*conversion;
        }
        else
        {
            const double factor = higherRadialConductivities_ ? 10.0 : 1.0;
            if (age < 10)
                return 1.8e-4*conversion*factor;
            else if (age >= 10 && age < 15)
                return (-(age-10)*1.62e-4/5.0 + 1.8e-4)*conversion*factor;
            else
                return 0.18e-4*conversion*factor;
        }
    }

    //! compute the axial conductivity in (m^4/Pa/s) given the segment age in days
    double computeAxialConductivity_(int order, double age)
    {
        constexpr double conversion = 1e-10/86400; // from cm^4/hPa/day to m^4/Pa/s
        if (order == 0)
        {
            if (age < 1)
                return 0.01*conversion;
            else if (age >= 1 && age < 3)
                return ((age-1)*0.29/2.0 + 0.01)*conversion;
            else if (age >= 3 && age < 4)
                return 0.3*conversion;
            else if (age >= 4 && age < 5)
                return ((age-4)*4 + 0.3)*conversion;
            else
                return 4.3*conversion;
        }
        else
        {
            if (age < 5)
                return 0.01e-3*conversion;
            else if (age >= 5 && age < 10)
                return ((age-5)*0.09e-3/5.0 + 0.01e-3)*conversion;
            else if (age >= 10 && age < 12)
                return ((age-10)*0.5e-3/2.0 + 0.1e-3)*conversion;
            else if (age >= 12 && age < 20)
                return 0.6e-3*conversion;
            else if (age >= 20 && age < 22)
                return ((age-20)*1.1e-3/2.0 + 0.6e-3)*conversion;
            else
                return 1.7e-3*conversion;
        }
    }

    //! Parameters for every root
    std::vector<RootParams> rootParams_;
    RootParams shootParams_;

    //! Segment paramters
    std::vector<int> rootId_; //! the root id for each segment
    std::vector<Scalar> previousLength_; //! the length of the element at the last time step (new elements return 0 length)
    std::vector<Scalar> ages_; //! the current age of the element
    std::vector<Scalar> Kx_;  //! axial conductivities
    std::vector<Scalar> Kr_;  //! radial conductivities

    std::vector<Scalar> radii_; //! radii for output
    std::vector<Scalar> orders_;  //! root orders for output

    Scalar porosity_;

    bool higherRadialConductivities_;
};

} // end namespace Dumux

#endif
