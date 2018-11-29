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

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv1p.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/growth/rootparameters.hh>

#include <RootSystem.h>

namespace Dumux {

/*!
 * \brief Root spatial params
 */
template<class FVGridGeometry, class Scalar>
class RootSpatialParams: public FVSpatialParamsOneP<FVGridGeometry, Scalar, RootSpatialParams<FVGridGeometry, Scalar>>
{
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, RootSpatialParams<FVGridGeometry, Scalar>>;

    using ElementMapper = typename FVGridGeometry::ElementMapper;
    using Water = Components::SimpleH2O<Scalar>;
	using DGFParamIndices = GrowthModule::DGFParamIndices;
	using RootParams = GrowthModule::RootParams<Scalar>;

public:

    // export permeability type
    using PermeabilityType = Scalar;

    RootSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry) :
        ParentType(fvGridGeometry)
    {
        kr_ = getParam<std::vector<Scalar>>("RootSystem.Conductivity.Kr");
		krTable_ = kr_.size()>1;
		if (krTable_) {
            krT_ = getParam<std::vector<Scalar>>("RootSystem.Conductivity.KrAge");
		}

        kx_ = getParam<std::vector<Scalar>>("RootSystem.Conductivity.Kx");
		kxTable_ = kx_.size()>1;
		if (kxTable_) {
            kxT_ = getParam<std::vector<Scalar>>("RootSystem.Conductivity.KxAge");
		}
    }

	/**
	 * Root radius (m)
	 */
    Scalar radius(const Element &element) const {
    	auto eIdx = this->fvGridGeometry().elementMapper().index(element);
    	int rootId = rootId_[eIdx];
    	if (rootId < 0)
    		return shootParams_.radius;
    	else
    		return rootParams_[rootId].radius;
    }

	/**
	 * Root radial conductivity (m^2 s / kg)
	 */
	Scalar radialConductivity(const Element &element) const {
		if (krTable_) { // table
			auto eIdx = this->fvGridGeometry().elementMapper().index(element);
			auto t = age_[eIdx];
            auto i = locate(t, krT_);
			return kr_.at(i);
		} else { // const
			return kr_[0];
		}
	}

	/**
	 * Root axial conductivity (m^5 s / kg)
	 */
	Scalar axialConductivity(const Element &element) const {
		if (krTable_) { // table
			auto eIdx = this->fvGridGeometry().elementMapper().index(element);
			auto t = age_[eIdx];
            auto i = locate(t, kxT_);
			return kx_.at(i);
		} else { // const
			return kx_[0];
		}
	}

	/*!
	 * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
	 */
	template<class ElementSolution>
	Scalar permeability(const Element& element,
			const SubControlVolume& scv,
			const ElementSolution& elemSol) const {
		Scalar mu = Water::liquidViscosity(285.15, 1e5); // temperature, pressure
		Scalar a = this->radius(element);
		Scalar kz = this->axialConductivity(element);
		return kz*mu/(a*a*M_PI); 		// a^2 * k / mu = kz  --> k = kz/a^2*mu
	}

	/*!
	 * \brief Return the root parameters
	 */
    RootParams& rootParams(const Element &element) const
    {
        const auto& elementMapper = this->fvGridGeometry().elementMapper();
        const auto eIdx = elementMapper.index(element);
        int rootId = rootId_[eIdx];
        if (rootId < 0)
            return shootParams_;
        else
            return rootParams_[rootId];
    }

    //! Read parameters from the a crootbox rootsystem
    void initParameters(const CRootBox::RootSystem& rs)
    {
        rootParams_.resize(rs.getNumberOfRoots(/*allRoots=*/true)); // be careful: the total number of roots and the grown number of roots is not the same

        const auto roots = rs.getRoots(); // this gets the root that already have at least one segment

        // copy parameters given per root into rootParams_
        for (auto&& r : roots) {
            auto& rp = rootParams_[r->id];
            // compute the root order
            rp.order = 0;
            auto* r_ = r;
            while (r_->parent != nullptr) {
                ++rp.order;
                r_=r_->parent;
            }
            rp.radius = r->param.a * 0.01; // conversion from cm to m
            rp.rootId = r->id;
            // rp.axialPerm = constantKx_;
            // rp.radialPerm = constantKr_;
            rp.plantId = 0;
        }

        // define shootParams_
        shootParams_.order = 0;
        shootParams_.radius = 0.01;
        shootParams_.rootId = -1;
        // shootParams_.axialPerm = constantKx_*1e20;
        // shootParams_.radialPerm = 0.0;
        shootParams_.plantId = 0;

        const auto shootSegments = rs.getShootSegments();
        const auto segmentRoots = rs.getSegmentsOrigin();
        const auto& gridView = this->fvGridGeometry().gridView();
        rootId_.resize(gridView.size(0));
        radii_.resize(gridView.size(0));
        orders_.resize(gridView.size(0));
        previousLength_.resize(gridView.size(0));
        for (const auto& element : elements(gridView)) {
            const auto& elementMapper = this->fvGridGeometry().elementMapper();
            const auto eIdx = elementMapper.index(element);

            // initialize the previous length to the current length
            previousLength_[eIdx] = element.geometry().volume();

            // subtract the number of shoot segments to get the normal segment index
            const int segmentIdx = eIdx - shootSegments.size();
            // shoot segments
            if (segmentIdx < 0) {
                rootId_[eIdx] = -1;
                const auto& rp = rootParams(element);
                radii_[eIdx] = rp.radius;
                orders_[eIdx] = rp.order;
            } else { // regular segments
                // we hope that the crootbox indices and the dune indices on initialization are the same!
                rootId_[eIdx] = segmentRoots[segmentIdx]->id;
                const auto& rp = rootParams_[rootId_[eIdx]];
                radii_[eIdx] = rp.radius;
                orders_[eIdx] = rp.order;
            }
        }
    }

//    //! update the root parameters that aren't set by growth directly
//    void updateParameters(const ::RootSystem& rs)
//    {
//        const auto roots = rs.getRoots();
//        for (auto&& r : roots)
//        {
//            auto& rp = rootParams_[r->id];
//            if (rp.isUpdated)
//                continue;
//
////            // TODO: compute the permeabilities dependent on age/radius --> On the segment level
////            rp.axialPerm = constantKx_;
////            rp.radialPerm = constantKr_;
//
//            rp.isUpdated = true;
//        }
//
//        // update radii for output
//        const auto& gridView = this->fvGridGeometry().gridView();
//        radii_.resize(gridView.size(0));
//        orders_.resize(gridView.size(0));
//        for (const auto& element : elements(gridView))
//        {
//            const auto& elementMapper = this->fvGridGeometry().elementMapper();
//            const auto eIdx = elementMapper.index(element);
//            const auto& rp = rootParams(element);
//            radii_[eIdx] = rp.radius;
//            orders_[eIdx] = rp.order;
//        }
//    }

    //! resize the segment parameters
    void resizeSegmentParams(std::size_t size)
    {
        rootId_.resize(size);
        previousLength_.resize(size);
    }

    //! set root index (needed by growth if eIdx of an element changed after growth)
    void setRootId(unsigned int eIdx, int rootId)
    {
    	rootId_[eIdx] = rootId;
    }

    //! Returns segment length of previous time step
    Scalar previousLength(std::size_t eIdx) const
    {
        return previousLength_[eIdx];
    }

    //! set the previous length
    void setPreviousLength(unsigned int eIdx, Scalar length)
    {
    	previousLength_[eIdx] = length;
    }

    //! Resize root params to add more roots
    void rootParamsResize(std::size_t size)
    {
    	rootParams_.resize(size);
    }

    //! Return a reference to the root parameters of a root id to change them (called from growth module)
    RootParams& getRootParams(int rootId)
    {
    	return rootParams_[rootId];
    }

    template<class ElementSol>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSol& elemSol) const
    {
    	return 0.4;
    }

    //! 1d table look up: xx is ascending, returns the index i , so that x>=xx[i] and x<xx[i+1]
    static size_t locate(Scalar x, const std::vector<Scalar>& xx) {
        unsigned int jr, jm, jl;
        jl = 0;
        jr = xx.size() - 1;
        while (jr - jl > 1) {
            jm = (jr + jl) >> 1; // thats a divided by two
            if (x >= xx[jm])
                jl = jm;
            else
                jr = jm;
        }
        // std::cout << "locate(" << x << ", " << xx[0] << ", " << xx[1] << ") = " << jl << "\n";
        return jl; // left index
    }

    //! returns linearly interpolated values of a 1-D function at specific query point x. Vector xx contains the sample points, and vv contains the corresponding values
    static Scalar interp1(Scalar x, const std::vector<Scalar>& vv, const std::vector<Scalar>& xx) {
        assert(xx.size() == vv.size());
        if (xx.size() == 1) {
            return vv[0];
        }
        size_t i = locate(x, xx);
        i = std::min(i, xx.size() - 1);
        Scalar t = (x - xx[i]) / (xx[i + 1] - xx[i]);
        t = std::min(std::max(t, 0.), 1.);
        Scalar v = vv[i] * (1. - t) + vv[i + 1] * t;
        // std::cout << "interp1(" << x << ", " << vv[0] << ", " << vv[1] << ", " << xx[0] << ", " << xx[1] << ")=" << v << ", " << t << "\n";
        return v;
    }

private:

    //! Parameters per root
    std::vector<RootParams> rootParams_;
    RootParams shootParams_;

    //! Parameters per segment
    std::vector<int> rootId_; //! the root id for each segment
    std::vector<Scalar> previousLength_; //! the length of the element at the last time step (new elements return 0 length)
    std::vector<Scalar> radii_; //! radii for output
    std::vector<Scalar> orders_;  //! root orders for output
    std::vector<Scalar> age_; // TODO

    // Kr, Kt could be age dependent
	std::vector<Scalar> kr_;
	std::vector<Scalar> krT_ = std::vector<Scalar>(0);
	bool krTable_;

	std::vector<Scalar> kx_;
	std::vector<Scalar> kxT_ = std::vector<Scalar>(0);
	bool kxTable_;

};

} // end namespace Dumux

#endif
