// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_ROOTS_PROPERTIES_NOCOUPLING_HH
#define DUMUX_ROOTS_PROPERTIES_NOCOUPLING_HH

//#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>
#include <dumux/multidomain/embedded/integrationpointsource.hh>

namespace Dumux {
namespace Properties {

/**
 * no functionality (but Dumux wants its bindings)
 * ugly, but I found no other option...
 */
// The point source type (not used)
template<class TypeTag>
struct PointSource<TypeTag, TTag::Roots> {
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using type = IntegrationPointSource<Dune::FieldVector<double, 3>, NumEqVector>;
};
/// Dummy types
class DummyPointSourceDataR {
public:
    double lowDimElementIdx() { throw 1; };
};
class DummySpatialR {
public:
    double kr(int i) const { throw 1; };
    double radius(int i) const { throw 1; };
};
class DummyProblemR {
public:
    DummySpatialR spatialParams() { throw 1; };
};
class DummyCouplingManagerR {
public:
    std::vector<double> bulkPriVars(int i) { throw 1; };
    std::vector<double>  lowDimPriVars(int i) { throw 1; };
    DummyProblemR& problem(int i) { throw 1; };
    DummyPointSourceDataR& pointSourceData(int i) { throw 1; };
    std::vector<double>& lowDimPointSources() { throw 1; };
    std::vector<double>& bulkPointSources() { throw 1; };
};
// For a dummy manager
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Roots> {
    using type = DummyCouplingManagerR;
};


} // namespace Properties
} // namespace Dumux

#endif
