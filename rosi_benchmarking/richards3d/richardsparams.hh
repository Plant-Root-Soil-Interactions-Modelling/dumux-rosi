/*!
 * \file
 *
 * \brief spatial parameters for the RichardsLensProblem
 */
#ifndef RICHARDS_PARAMETERS_HH
#define RICHARDS_PARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>

#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/richards/implicit/model.hh>
#include <dumux/porousmediumflow/richards/implicit/properties.hh>

#include <vector>



namespace Dumux
{

template<class TypeTag> class RichardsParams; // forward declaration



/**
 *  RegularizedVanGenuchten is not working! VanGenuchten is not working, but this does...
 */
template <class ScalarT, class ParamsT = VanGenuchtenParams<ScalarT> >
class MyVanGenuchten : public VanGenuchten<ScalarT,ParamsT>
{
public:
    static double pc(const ParamsT &params, double swe)
    {
        if (swe<=0) {
            return HUGE_VAL; // something reasonable
        }
        if (swe>=1) {
            return 0;
        }
        return pow(pow(swe, -1.0/params.vgm()) - 1, 1.0/params.vgn())/params.vgAlpha();
    }
};



namespace Properties
{

// The spatial parameters TypeTag
NEW_TYPE_TAG(RichardsParams);

// Set the spatial parameters defined as TTAG of type RichardsParams and attached to TTAG RichardsParam
SET_TYPE_PROP(RichardsParams, SpatialParams, Dumux::RichardsParams<TypeTag>);

// Set the material law
SET_PROP(RichardsParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef MyVanGenuchten<Scalar> EffectiveLaw; //RegularizedVanGenuchten
public:
    typedef EffToAbsLaw<EffectiveLaw> type; // define the material law parameterized by absolute saturations
};

}



/**
 * Solves the Richard equation in 1D with Van Genuchten model
 *
 * Van Genuchten parameters are passed using the .input file:
 * [VanGenuchten] # silt
 * Qr = 0.034
 * Qs = 0.46
 * alpha = 0.016 # [1/cm]
 * n = 1.37
 * Ks = 6.94444e-7 # [m/s]
 *
 * Initial values are passed by the .dgf file. The file can be produced by the Matlab script create1Ddgf
 *
 * Boundary values can be chosen in the .input file:
 * [BC_Top]
 * type = 2 # 1 constant pressure head, 2 constant flux, more will follow
 * value = 0
 *
 * [BC_Bot]
 * type = 2 # 1 constant pressure head, 2 constant flux, more will follow
 * value = 0
 *
 * Output times can be chosen in the .input file:
 * [TimeManager]
 * Episodes = 60 120 2000 # s
 */
template<class TypeTag> class RichardsParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid; // defined in myrichards problem
    typedef typename Grid::ctype CoordScalar;  //typedef double CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar; // defined in the depths of dumux //typedef double Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView; // defined at the problem constructor, in start.hh
    typedef typename GridView::template Codim<0>::Entity Element;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry; // CC-, or BoxFVElementGeometry
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WP;
    static const bool useHead = GET_PROP_VALUE(TypeTag, UseHead); // head is still not working...
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator; //typedef Dumux::GridCreator<TypeTag> GridCreator;  // not too sure about that

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw; // defined above
    typedef typename MaterialLaw::Params MaterialLawParams; // type of MaterialLaw

public:

    /*!
     * \brief Constructor
     *
     * \param gridView The DUNE GridView representing the spatial domain of the problem
     *
     */
    RichardsParams(const GridView& gridView) : ParentType(gridView)
    {
        phi_ = 1; // Richards equation is independent of phi
        const double g = 9.81; // TODO fetch values from problem class

        const double temp =  273.15 + 10;  // -> 10°C
        const double pnRef = 1e5;
        const double mu = WP::viscosity(temp,pnRef); // h2o: Dynamic viscosity of water 0.001 Pa·s (independent of temp and p)
        const double rho = WP::density(temp,pnRef);  // h2o: 1000 kg/m³ (independent of temp and p)

        // get van genuchten parameters from the input file
        std::vector<double> Qr = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<double>, VanGenuchten, Qr);
        std::vector<double> Qs = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<double>, VanGenuchten, Qs);
        std::vector<double> alpha = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<double>, VanGenuchten, alpha);
        std::vector<double> n = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<double>, VanGenuchten, n);
        Kc_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<double> , VanGenuchten, Ks); // hydraulic conductivity

        more_ = Qr.size()>1; // more than one set of VG parameters?

        // Qr, Qs, alpha, and n goes to the MaterialLaw VanGenuchten
        for (int i=0; i<Qr.size(); i++) {
            materialParams_.push_back(MaterialLawParams());
            // QR
            materialParams_.at(i).setSwr(Qr.at(i)/phi_);
            // QS
            materialParams_.at(i).setSnr(1.-Qs.at(i)/phi_);
            // ALPHA
            double a = alpha.at(i) * 100.; // from [1/cm] to [1/m]
            materialParams_.at(i).setVgAlpha(a/(rho*g)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
            // N
            materialParams_.at(i).setVgn(n.at(i));
            // Kc
            K_.push_back(Kc_.at(i)*mu/(rho*g)); // convert to intrinsic permeability (from the hydraulic conductivity)
            // Debug
            std::cout << "\nVan Genuchten Parameters are " << Qr.at(i) << ", " << Qs.at(i) <<
                    ", "<< alpha.at(i) << ", "<< n.at(i) << ", "<< Kc_.at(i) << "\n";
        }
    }

    /*!
     * \brief Function for defining the intrinsic (absolute) permeability. [m²]
     *
     * \return intrinsic (absolute) permeability
     * \param globalPos The position of the center of the element
     */
    const DimWorldMatrix intrinsicPermeability(const Element &element, const FVElementGeometry &fvGeometry, int scvIdx) const
    {
        return K_.at(getDI(element));
    }

    /*
     * \brief Function for defining the (absolute) hydraulic conductivity. [m/s]
     */
    const Scalar hydraulicConductivity(const Element &element, const FVElementGeometry &fvGeometry, int scvIdx) const
    {
        return Kc_.at(getDI(element));
    }

    /*!
     * \brief Function for defining the porosity.
     *
     * \return porosity
     * \param globalPos The position of the center of the element
     */
    Scalar porosity(const Element &element, const FVElementGeometry &fvGeometry, int scvIdx) const
    {
        return phi_; // porosity*saturation = watercontent
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \return the material parameters object
     * \param globalPos The position of the center of the element
     */
    const MaterialLawParams& materialLawParams(const Element &element, const FVElementGeometry &fvGeometry, int scvIdx) const
    {
        return materialParams_.at(getDI(element));
    }

private:

    /**
     * returns the domain index
     */
    int getDI(const Element &element) const
    {
        if (more_) {
            int i = GridCreator::parameters(element).at(1)-1; // starting from 1
            return i;
        } else {
            return 0;
        }
    }

    bool more_;
    double phi_;
    std::vector<double> K_; // permeability [m²]
    std::vector<double> Kc_; // hydraulic conductivity [m/s]
    std::vector<MaterialLawParams> materialParams_;

};

} // end namespace

#endif
