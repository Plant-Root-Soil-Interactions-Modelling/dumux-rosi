#ifndef DUMUX_SOIL_PROPERTIES_NC_NOCOUPLING_HH
#define DUMUX_SOIL_PROPERTIES_NC_NOCOUPLING_HH

namespace Dumux {
namespace Properties {

/**
 * no functionality (but Dumux wants its static bindings)
 * ugly, but I found no other option...
 */
// The point source type (not used)
template<class TypeTag>
struct PointSource<TypeTag, TTag::RichardsTT> {
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using type = IntegrationPointSource<Dune::FieldVector<double, 3>, NumEqVector>;
};
/// Dummy types
class DummyPointSourceData {
public:
    double lowDimElementIdx() { throw 1; };
};
class DummySpatial {
public:
    double kr(int i) const { throw 1; };
    double radius(int i) const { throw 1; };
};
class DummyProblem {
public:
    DummySpatial spatialParams() { throw 1; };
};
class DummyCouplingManager {
public:
    std::vector<double> bulkPriVars(int i) { throw 1; };
    std::vector<double>  lowDimPriVars(int i) { throw 1; };
    DummyProblem& problem(int i) { throw 1; };
    DummyPointSourceData& pointSourceData(int i) { throw 1; };
    std::vector<double>& lowDimPointSources() { throw 1; };
    std::vector<double>& bulkPointSources() { throw 1; };
};
// For a dummy manager
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RichardsTT> {
    using type = DummyCouplingManager;
};

} // end namespace properties
} // end namespace DUMUX

#endif
