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
 * \brief An adpater for growth problems
 */
#ifndef DUMUX_GROWTH_PROBLEM_ADAPTER_HH
#define DUMUX_GROWTH_PROBLEM_ADAPTER_HH

namespace Dumux {
namespace GrowthModule {

/*!
 * \brief a root problem file
 */
template <class Base>
class GrowthProblemAdapter : public Base
{
public:
    using Base::Base;

    //! return the pre growth volume of an scv
    //! needed for the residual for elements changing size by growth
    template<class Element, class SubControlVolume, class VolumeVariables>
    auto preGrowthVolume(const Element& element,
                         const SubControlVolume& scv,
                         const VolumeVariables& volVars)
    {
        if (isNew(element))
            return 0.0;
        else
        {
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            return volVars.extrusionFactor() * this->spatialParams().previousLength(eIdx);
        }
    }

    //! return the post growth volume of an scv
    //! needed for the residual for elements changing size by growth
    template<class Element, class SubControlVolume, class VolumeVariables>
    auto postGrowthVolume(const Element& element,
                          const SubControlVolume& scv,
                          const VolumeVariables& volVars)
    {
        return volVars.extrusionFactor()*scv.volume();
    }

    void clearIsNewMarkers()
    {
        isNew_.assign(this->fvGridGeometry().gridView().size(0), false);
    }

    template<class Element>
    bool isNew(const Element &element) const
    {
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        return isNew_[eIdx];
    }

    template<class Element>
    void setIsNewMarker(const Element& element, bool isNew)
    {
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        isNew_[eIdx] = isNew;
    }

private:
    std::vector<bool> isNew_; //! Marker which element are currently new
};

} // end namespace GrowthModule
} //end namespace Dumux

#endif
