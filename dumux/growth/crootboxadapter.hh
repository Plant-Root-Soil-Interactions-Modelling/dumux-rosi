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
 * \brief This file contains a generic grid growth class
 */
#ifndef DUMUX_CROOTBOX_ADAPTER_HH
#define DUMUX_CROOTBOX_ADAPTER_HH

#include <Organism.h>
#include <Root.h>
#include <mymath.h>

#include "growthinterface.hh"

namespace Dumux {

namespace GrowthModule {

using namespace CRootBox;

/**
 * Adapter for CRootBox: Converts types, units, and naming conventions
 */
template<class GlobalPosition>
class CRootBoxAdapter :public GrowthInterface<GlobalPosition> {

public:

    CRootBoxAdapter(CRootBox::RootSystem& rs) :rootsystem_(rs) {
        this->root2dune = std::vector<size_t>(rs.getNumberOfNodes());
        std::iota(this->root2dune.begin(), this->root2dune.end(), 0);
    };

    virtual ~CRootBoxAdapter() { }; // nothing to do, rootsystem_ is someone elses problem

    void simulate(double dt) override {
        rootsystem_.simulate(dt/3600/24, false);
    }

    std::vector<size_t> updatedNodeIndices() const override {
        auto ni = rootsystem_.getUpdatedNodeIndices();
        return std::vector<size_t>(ni.begin(), ni.end());
    }

    std::vector<GlobalPosition> updatedNodes() const override {
        std::vector<Vector3d> n = rootsystem_.getUpdatedNodes();
        auto p = std::vector<GlobalPosition>(n.size());
        for (size_t i=0; i<n.size(); i++) {
            p[i][0] = n[i].x/100.; // convert to m
            p[i][1] = n[i].y/100.;
            p[i][2] = n[i].z/100.;
        }
        return p;
    }

    std::vector<size_t> newNodeIndices() const override {
        auto v = std::vector<size_t>(rootsystem_.getNumberOfNewNodes());
        std::iota(v.begin(), v.end(), rootsystem_.getNumberOfNodes()-rootsystem_.getNumberOfNewNodes());
        return v;
    }

    std::vector<GlobalPosition> newNodes() const override {
        std::vector<Vector3d> n = rootsystem_.getNewNodes();
        auto p = std::vector<GlobalPosition>(n.size());
        for (size_t i=0; i<n.size(); i++) {
            p[i][0] = n[i].x/100.; // convert to m
            p[i][1] = n[i].y/100.;
            p[i][2] = n[i].z/100.;
        }
        return p;
    }

    std::vector<std::array<size_t, 2>> newSegments() const override {
        std::vector<Vector2i> s = rootsystem_.getNewSegments();
        auto seg = std::vector<std::array<size_t, 2>>(s.size());
        for (size_t i=0; i<s.size(); i++) {
            seg[i] = std::array<size_t, 2>{size_t(s[i].x), size_t(s[i].y)};
        }
        return seg;
    }

    std::vector<double> segmentCreationTimes() const override {
        auto sct = rootsystem_.getNewSegmentCTs();
        std::transform(sct.begin(), sct.end(), sct.begin(), std::bind1st(std::multiplies<double>(), 24.*3600.)); // convert to s
        return sct;
    }

    std::vector<int> segmentOrders() const override {
        std::vector<Organ*> roots = rootsystem_.getNewSegmentOrigins();
        auto orders = std::vector<int>(roots.size());
        for (size_t i=0; i<roots.size(); i++) {
            Organ* r = roots[i];
            int o = 0;
            while (r->getParent() != nullptr) {
                o++;
                r = (Organ*)r->getParent();
            }
            orders[i] = o;
        }
        return orders;
    }

    /**
     * Currently only for roots (e.g. more general to roots[i]->getParamter("radius") )
     */
    std::vector<double> segmentRadii() const override {
        std::vector<Organ*> roots = rootsystem_.getNewSegmentOrigins();
        auto radii = std::vector<double>(roots.size());
        for (size_t i=0; i<roots.size(); i++) {
            radii[i] = ((Root*)roots[i])->param()->a / 100.; // convert to m
        }
        return radii;
    }

private:

    CRootBox::Organism& rootsystem_;

};

} // end namespace GridGrowth

} // end namespace Dumux

#endif
