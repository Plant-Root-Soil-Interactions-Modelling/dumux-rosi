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
 * \ingroup Properties
 * \ingroup CCProperties
 * \ingroup CCModel
 * \file
 *
 * \brief Default properties for cell centered models
 */
#ifndef DUMUX_MPFA_PROPERTY_DEFAULTS_HH
#define DUMUX_MPFA_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/cellcentered/elementboundarytypes.hh>
#include "mpfaelementvolumevariables.hh"
#include "mpfafvelementgeometry.hh"
#include "mpfaassembler.hh"
#include "mpfahelper.hh"
#include "mpfalocalresidual.hh"
#include "mpfaproperties.hh"
#include "mpfaimplicitdarcyfluxvariables.hh"
#include "mpfaimplicitlocaljacobian.hh"



namespace Dumux {

struct BoundaryLayers
{
public:
    // actual node on boundary
    static const unsigned int boundary = 0;
    // node on second half edge
    static const unsigned int intermediate = 1;
    // interior node
    static const unsigned int interior = 2;
};

struct CouplingStrategies
{
public:
    // actual node on boundary
    static const unsigned int simpleDirichlet = 0;
    // node on second half edge
    static const unsigned int complexDirichlet = 1;
    // interior node
    static const unsigned int fluxCoupling = 2;
};

// forward declarations
template<class TypeTag> class MpfaLocalResidual;
template<class TypeTag> class CCElementBoundaryTypes;
template<class TypeTag> class MpfaElementVolumeVariables;
template<class TypeTag> class MpfaFVElementGeometry;
template<class TypeTag> class MpfaAssembler;
template<class TypeTag, int method, int dimension> class MPFAHelper;

namespace Properties {
//! Set the default for the FVElementGeometry
SET_TYPE_PROP(MpfaModel, FVElementGeometry, Dumux::MpfaFVElementGeometry<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(MpfaModel, ElementBoundaryTypes, Dumux::CCElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(MpfaModel, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(MpfaModel, BaseLocalResidual, Dumux::MpfaLocalResidual<TypeTag>);

//! An array of secondary variable containers
SET_TYPE_PROP(MpfaModel, ElementVolumeVariables, Dumux::MpfaElementVolumeVariables<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(MpfaModel, JacobianAssembler, Dumux::MpfaAssembler<TypeTag>);

//! The local jacobian operator
SET_TYPE_PROP(MpfaModel, LocalJacobian, Dumux::MpfaImplicitLocalJacobian<TypeTag>);

//! indicate that this is no box discretization
SET_BOOL_PROP(MpfaModel, ImplicitIsBox, false);

//! by default only calculate the flux over the half faces with the o method
SET_BOOL_PROP(MpfaModel, FullMpfaOBoundary, false);

//! Definition of indices for different layers at boundary
SET_TYPE_PROP(MpfaModel, BoundaryLayers, Dumux::BoundaryLayers);

//! set a bool indicating that we use the mpfa
SET_BOOL_PROP(MpfaModel, ImplicitIsMpfa, true);

//! set a bool indicating if we want to couple a lower dimensional model on the element facets
SET_BOOL_PROP(MpfaModel, FacetCoupling, false);

//! set a bool indicating if we want to use flux coupling in the coupled model
SET_TYPE_PROP(MpfaModel, CouplingStrategies, Dumux::CouplingStrategies);

//! default for the coupling strategy atually used in the coupled model
SET_INT_PROP(MpfaModel, CouplingStrategy, Dumux::CouplingStrategies::simpleDirichlet);

//! set the factor for the coupling conditions on internal flux faces when facet coupling is active
SET_SCALAR_PROP(MpfaModel, XiFactor, 1.0);

//! a type to extract the indices for the different mpfa methods
SET_PROP(MpfaModel, MpfaMethods)
{
     // MPFA O-method
    static const int oMethod = 0;
    // MPFA L-method
    static const int lMethod = 1;
};

//! set the mpfa-o method by default
SET_INT_PROP(MpfaModel, MpfaMethod, GET_PROP(TypeTag, MpfaMethods)::oMethod);

//! extract interactionvolume type from the helper
SET_PROP(MpfaModel, MpfaInteractionVolume)
{
private:
    typedef typename Dumux::MpfaHelper<TypeTag, GET_PROP_TYPE(TypeTag, GridView)::dimension> MpfaHelper;
public:
    typedef typename MpfaHelper::InteractionVolume type;
};

//! set boundary interactionVolumeType
// by default, the O-method is applied on the boundaries
SET_TYPE_PROP(MpfaModel, MpfaBoundaryInteractionVolume, Dumux::MpfaO2DInteractionVolume<TypeTag>);

//! extract interactionvolume manager type from the helper
SET_PROP(MpfaModel, MpfaInteractionVolumeManager)
{
private:
    typedef typename Dumux::MpfaHelper<TypeTag, GET_PROP_TYPE(TypeTag, GridView)::dimension> MpfaHelper;
public:
    typedef typename MpfaHelper::InteractionVolumeManager type;
};

//! set the interaction volume manager type for the boundary volume
// by default, the mpfa-o method is applied on the boundaries
SET_TYPE_PROP(MpfaModel, MpfaBoundaryInteractionVolumeManager, Dumux::MpfaO2DManager<TypeTag>);

//! extract interactionvolume container type from the helper
SET_PROP(MpfaModel, MpfaInteractionVolumeContainer)
{
private:
    typedef typename Dumux::MpfaHelper<TypeTag, GET_PROP_TYPE(TypeTag, GridView)::dimension> MpfaHelper;
public:
    typedef typename MpfaHelper::InteractionVolumeContainer type;
};

//! set the class for the face types
SET_PROP(MpfaModel, MpfaFaceTypes)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MpfaBoundaryInteractionVolume) InteractionVolume;
public:
    typedef typename InteractionVolume::FaceTypes type;
};

//! extract flux calculator type from the helper
SET_PROP(MpfaModel, MpfaFluxCalculator)
{
private:
    typedef typename Dumux::MpfaHelper<TypeTag, GET_PROP_TYPE(TypeTag, GridView)::dimension> MpfaHelper;
public:
    typedef typename MpfaHelper::FluxCalculator type;
};

//! extract flux calculator type for the boundary volumes
// by default, the mpfa-o method is used on the boundaries
SET_TYPE_PROP(MpfaModel, MpfaBoundaryFluxCalculator, Dumux::MpfaO2DFluxCalculator<TypeTag>);

//! by default, continuity in the interactionvolumes is enforced on the intersection centers
SET_SCALAR_PROP(MpfaModel, MpfaContinuityPoint, 0.);



} // namespace Properties
} // namespace Dumux

#endif
