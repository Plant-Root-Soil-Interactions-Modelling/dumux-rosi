// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS_PROPERTIES_CYL_HH
#define RICHARDS_PROPERTIES_CYL_HH

#include <dumux/common/properties.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richardsCylindrical1d/model.hh>
#include "../soil_richards/richardsproblem.hh" // the problem class
#include <dumux/discretization/extrusion.hh>
//where do i add teh local residual file?
//do i add the spatial params?

namespace Dumux {
namespace Properties {

namespace TTag { // Create new type tags
struct RichardsTT { using InheritsFrom = std::tuple<Richards>; };
struct RichardsBox { using InheritsFrom = std::tuple<RichardsTT, BoxModel>; };
struct RichardsCC { using InheritsFrom = std::tuple<RichardsTT, CCTpfaModel>; };
}

template<class TypeTag> // Set grid type
struct Grid<TypeTag, TTag::RichardsCylFoamTT> { using type = Dune::FoamGrid<1,1>; }; //  Dune::SPGrid<GetPropType<TypeTag, Properties::Scalar>, 1>



// As mentioned before, we modify the areas and volumes used
// for integrals over the control volume and the control volume faces by changing the extrusion type.
// Here, we pass these traits to the grid geometry of the box scheme (the scheme
// that we use here) and specialize the `GridGeometry` property accordingly.
// [[codeblock]]
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::OnePRotSym>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    // We take the default traits as basis and exchange the extrusion type
    // The first axis (x-axis) is the radial axis, hence the zero. That means we rotate about the second axis (y-axis).
    struct GGTraits : public BoxDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

public:
    // Pass the above traits to the box grid geometry such that it uses the
    // rotation-symmetric sub-control volumes and faces.
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsTT> { using type = RichardsProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsTT> {
    using type = RichardsParams<GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>>;
};

} // end namespace properties
} // end namespace DUMUX

#endif
