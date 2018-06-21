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
 * \ingroup EmbeddedCoupling
 * \ingroup CCModel
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain. Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */

#ifndef DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_HH
#define DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_HH

#include <dumux/common/properties.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Coupled root and soil domain
 */
template<class MDTraits>
class SoilRootCouplingManager
: public EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::average>
{
    using ParentType = EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::average>;
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<i>;
    template<std::size_t i> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<i>, Problem);

    static constexpr auto bulkIdx = typename MDTraits::template DomainIdx<0>();

public:
    using ParentType::ParentType;

    Scalar relPermSoil(Scalar p) const
    {
        using MaterialLaw = typename Problem<bulkIdx>::SpatialParams::MaterialLaw;
        const Scalar pc = std::max(0.0, this->problem(bulkIdx).nonWettingReferencePressure() - p);
        const auto& materialParams = this->problem(bulkIdx).spatialParams().materialLawParams();
        const Scalar sw = MaterialLaw::sw(materialParams, pc);
        return MaterialLaw::krw(materialParams, sw);
    }
};

} // end namespace Dumux

#endif
