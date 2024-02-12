// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:

#ifndef RICHARDS_PARAMETERS_HH
#define RICHARDS_PARAMETERS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabsdefaultpolicy.hh>
#include <dumux/material/fluidmatrixinteractions/2p/materiallaw.hh>

#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/io/inputfilefunction.hh>


namespace Dumux {

/*!
 * The SpatialParams class of RichardsProblem
 *
 * supports multiple soil layers (in z-direction),
 * with different VG parameters sets
 */
template<class GridGeometry, class Scalar>
class RichardsParams : public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, RichardsParams<GridGeometry, Scalar>>
{
public:

    using ThisType = RichardsParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Water = Components::SimpleH2O<Scalar>;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;
    using BasicParams = typename PcKrSwCurve::BasicParams;
    using EffToAbsParams = typename PcKrSwCurve::EffToAbsParams;
    using RegularizationParams = typename PcKrSwCurve::RegularizationParams;
    // using MaterialLaw = FluidMatrix::TwoPMaterialLaw<Scalar, FluidMatrix::VanGenuchten, FluidMatrix::VanGenuchtenRegularization<Scalar>, FluidMatrix::TwoPEffToAbsDefaultPolicy>;

    enum { dimWorld = GridView::dimensionworld };

    using PermeabilityType = Scalar; // export permeability type

    RichardsParams(std::shared_ptr<const GridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {

        /* SimpleH2O is constant in regard to temperature and reference pressure */
        Scalar mu = Water::liquidViscosity(0.,0.); // Dynamic viscosity: 1e-3 [Pa s]
        Scalar rho = Water::liquidDensity(0.,0.);  // Density: 1000 [kg/m³]

        /* Get Van Genuchten parameters from the input file */
        std::vector<Scalar> qr = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Qr"); // [1]
        std::vector<Scalar> qs = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Qs"); // [1]
        std::vector<Scalar> alpha = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Alpha");  // [1/cm]
        std::vector<Scalar> n = getParam<std::vector<Scalar>>("Soil.VanGenuchten.N"); // [1]
        kc_ = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Ks"); // hydraulic conductivity [cm/day]
        std::transform(kc_.begin (), kc_.end (), kc_.begin (), std::bind1st(std::multiplies<Scalar>(), 1./100./24./3600.)); // convert from [cm/day] to [m/s]
        homogeneous_ = qr.size()==1; // more than one set of VG parameters?
        three_ = false;

        phi_.resize(qr.size());
        // Qr, Qs, alpha, and n goes to the PcKrSwCurve VanGenuchten
        for (int i=0; i<qr.size(); i++) {
            phi_[i] =  qs.at(i); // Richards equation is independent of phi [1]
            basicParams_.push_back(BasicParams(0.,0.));
            effToAbsParams_.push_back(EffToAbsParams());
            regularizationParams_.push_back(RegularizationParams());

            effToAbsParams_.at(i).setSwr(qr.at(i)/phi_[i]); // Qr
            effToAbsParams_.at(i).setSnr(1.-qs.at(i)/phi_[i]); // Qs

            Scalar a = alpha.at(i) * 100.; // from [1/cm] to [1/m]
            basicParams_.at(i).setAlpha(a/(rho*g_)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
            basicParams_.at(i).setN(n.at(i)); // N
            k_.push_back(kc_.at(i)*mu/(rho*g_)); // Convert to intrinsic permeability

            // Regularisation parameters
            double eps = 1.e-4; // with 1.e-9 benchmark 3 does not work anymore (and everything becomes slower)
            regularizationParams_.at(i).setPcLowSwe(eps);
            regularizationParams_.at(i).setPcHighSwe(1. - eps);
            regularizationParams_.at(i).setKrnLowSwe(eps);
            regularizationParams_.at(i).setKrwHighSwe(1 - eps);

            materialLaw_.push_back(PcKrSwCurve(basicParams_.at(i), effToAbsParams_.at(i), regularizationParams_.at(i)));
        }
        layerIdx_ = getParam<int>("Soil.Grid.layerIdx", 1);
        layer_ = InputFileFunction("Soil.Layer", "Number", "Z", layerIdx_, 0); // [1]([m])

        // std::cout << "RichardsParams created: homogeneous " << homogeneous_ << " " << "\n" << std::endl;
    }

    /*!
     * \brief \copydoc GridGeometry::porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const {
        return phi_.at(index_(element));
    }

    /*!
     * \brief \copydoc GridGeometry::porosity
     * simper interface
     */
    Scalar porosity(const Element& element) const {
        return phi_.at(index_(element));
    }

    /*!
     * \brief \copydoc FVSpatialParamsOneP::permeability
     * [m^2]\
     */
    template<class ElementSolution>
    Scalar permeability(const Element& element,
        const SubControlVolume& scv, const ElementSolution& elemSol) const {
        return permeability(element);
    }

    //! simpler interface
    Scalar permeability(const Element& element) const {
        return k_.at(index_(element));
    }

    /*
     * \brief Hydraulic conductivities [m/s], called by the problem for conversions
     */
    const Scalar hydraulicConductivity(const Element& element) const {
        return kc_.at(index_(element));
    }

    //! set of VG parameters for the element
    const BasicParams& basicParams(const Element& element) const {
        return basicParams_.at(index_(element));
    }

    //! set of VG parameters for the element
    const EffToAbsParams& effToAbsParams(const Element& element) const {
        return effToAbsParams_.at(index_(element));
    }

    //! set of VG parameters for the element
    const RegularizationParams& regularizationParams(const Element& element) const {
        return regularizationParams_.at(index_(element));
    }

    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const {
        //return makeFluidMatrixInteraction(PcKrSwCurve(basicParams_.at(index_(element)), effToAbsParams_.at(index_(element)), regularizationParams_.at(index_(element))));
        return materialLaw_.at(index_(element));
    }

    //! pointer to the soils layer input file function
    InputFileFunction* layerIFF() {
        return &layer_;
    }

    /**
     * Call to change default setting (of 1.e-6 for both)
     *
     * pcEps    capillary pressure regularisation
     * krEps 	relative permeabiltiy regularisation
     */
    void setRegularisation(double pcEps, double krEps) {
    	for (int i =0; i<regularizationParams_.size(); i++) {
    	    regularizationParams_.at(i).setPcLowSwe(pcEps);
    	    regularizationParams_.at(i).setPcHighSwe(1. - pcEps);
    	    regularizationParams_.at(i).setKrnLowSwe(krEps);
    	    regularizationParams_.at(i).setKrwHighSwe(1 - krEps);
    	    materialLaw_.at(i) = PcKrSwCurve(basicParams_.at(i), effToAbsParams_.at(i), regularizationParams_.at(i)); // update material Law
    	}
    }

    void addVanGenuchtenDomain(double minx, double miny, double minz, double maxx, double maxy, double maxz, int layerIndex) {
        homogeneous_ = false;
        three_ = true;
        std::vector<double> minmax{ minx, miny, minz, maxx, maxy, maxz };
        boxes.push_back(minmax);
        layerIndices.push_back(layerIndex);
    }

    void changeVanGenuchtenSet(int vgIndex, double qr, double qs, double alpha, double n, double ks) {
        Scalar mu = Water::liquidViscosity(0.,0.); // Dynamic viscosity: 1e-3 [Pa s]
        Scalar rho = Water::liquidDensity(0.,0.);  // Density: 1000 [kg/m³]
        int i = vgIndex; // rename
        phi_.at(i) =  qs; // Richards equation is independent of phi [1]
        effToAbsParams_.at(i).setSwr(qr/phi_[i]); // Qr
        effToAbsParams_.at(i).setSnr(1.-qs/phi_[i]); // Qs
        Scalar a = alpha * 100.; // from [1/cm] to [1/m]
        basicParams_.at(i).setAlpha(a/(rho*g_)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
        basicParams_.at(i).setN(n); // N
        k_.push_back(ks*mu/(rho*g_)); // Convert to intrinsic permeability
        materialLaw_.at(i) = PcKrSwCurve(basicParams_.at(i), effToAbsParams_.at(i), regularizationParams_.at(i)); // update material Law
    }

private:

    //! returns the index of the soil layer
    size_t index_(const Element& element) const {
        if (homogeneous_) {
            return 0;
        } else if (three_){ // use 3d min max boxes
            auto mid = element.geometry().center();
            int c = 0;
            for (auto& v : boxes) {
                if ((mid[0]>v[0]) and (mid[1]>v[1]) and(mid[2]>v[2]) and (mid[0]<v[3])and (mid[1]<v[4]) and (mid[2]<v[5])) {
                    return layerIndices[c];
                }
                c++;
            }
            return 0;
        } else { // use input file function
            auto eIdx = this->gridGeometry().elementMapper().index(element);
            Scalar z = element.geometry().center()[dimWorld - 1];
            //std::cout << z << "\n";
            return size_t(layer_.f(z, eIdx)-1); // layer number starts with 1 in the input file
        } // add 3D things
    }

    std::vector<Scalar> phi_; // porosity
    std::vector<std::vector<double>> boxes;
    std::vector<int> layerIndices;

    bool homogeneous_; // soil is homogeneous
    bool three_; // 3d
    InputFileFunction layer_;
    int layerIdx_; // index of layer data within the grid
    std::vector<Scalar> k_; // permeability [m²]
    std::vector<Scalar> kc_; // hydraulic conductivity [m/s]

    std::vector<BasicParams> basicParams_;
    std::vector<EffToAbsParams> effToAbsParams_;
    std::vector<RegularizationParams> regularizationParams_;
    std::vector<PcKrSwCurve> materialLaw_;

    static constexpr Scalar g_ = 9.81; // cm / s^2

};

} // end namespace Dumux

#endif
