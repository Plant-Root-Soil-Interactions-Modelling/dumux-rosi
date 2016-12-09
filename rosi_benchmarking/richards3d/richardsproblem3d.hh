#ifndef MYRICHARDSPROBLEM_HH
#define MYRICHARDSPROBLEM_HH



#include <dune/grid/io/file/dgfparser.hh>

#include <dumux/porousmediumflow/richards/implicit/model.hh>
// richards model -> implicit model ->  properties
// richards model -> richards problem -> properties
// a lot happens in in dumux/implict/model.hh

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "richardsparams.hh"



namespace Dumux
{

template <class TypeTag> class RichardsProblem3d;

/*
 * the property definitions are distributed over various header files.
 * most relevant properties are defined below, and in richardsparams.hh
 */
namespace Properties
{

NEW_TYPE_TAG(RichardsProblem3d, INHERITS_FROM(Richards, RichardsParams));
NEW_TYPE_TAG(RichardsBoxProblem3d, INHERITS_FROM(BoxModel, RichardsProblem3d));   // box is the discretisation method
NEW_TYPE_TAG(RichardsCCProblem3d, INHERITS_FROM(CCModel, RichardsProblem3d));     // cc = cell centered is the discretisation method

SET_TYPE_PROP(RichardsProblem3d, Grid, Dune::YaspGrid<3>);

// attach the type class Dumux::MyRichardsProblem to the tag Properties::MyRichardsProblem
SET_TYPE_PROP(RichardsProblem3d, Problem, Dumux::RichardsProblem3d<TypeTag>);

// Set the wetting phase
SET_PROP(RichardsProblem3d, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(RichardsProblem3d, ProblemEnableGravity, true);

SET_BOOL_PROP(RichardsProblem3d, UseHead, false); // head is not working, i dont know why

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
template <class TypeTag>
class RichardsProblem3d : public RichardsProblem<TypeTag>
{
    typedef RichardsProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry; // BoxFVElementGeometry
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables; // RichardsVolumeVariables
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables; // BoxElementVolumeVariables
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;  // defined in richardsparam (spatialParams)
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar; // no idea where this is defined (always double anyway)
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager; // in propertydefaults.hh (no idea, where this is included) // typedef Dumux::TimeManager<TypeTag> TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // copy some indices for convenience
        pwIdx = Indices::pwIdx,
        hIdx = Indices::hIdx,
        contiEqIdx = Indices::contiEqIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    static const bool useHead = GET_PROP_VALUE(TypeTag, UseHead); // false, because true is not working
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator; //typedef Dumux::GridCreator<TypeTag> GridCreator;  // not too sure about that
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase)  WP;

public:

    /*!
     * \brief Constructor
     *
     * \param timeManager The Dumux TimeManager for simulation management.
     * \param gridView The grid view on the spatial domain of the problem
     */
    RichardsProblem3d(TimeManager &timeManager, const GridView &gridView) : ParentType(timeManager, gridView)
{
        pnRef_ = 1e5; // reference pressure if Pascal are used

        /*
         * read relevant run time parameter
         */
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        bcTop_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, BC_Top, Type);
        bcBot_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, BC_Bot, Type);

        bcTopValue_ = 0; // default value
        if ((bcTop_==1) || (bcTop_==2)) { // read constant pressure if dirchlet head (in [cm]) or constant flux (in [ kg/(m² s)])
            bcTopValue_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BC_Top, Value);
        }

        bcBotValue_ = 0; // default value
        if ((bcBot_==1) || (bcBot_==2)) { // read constant pressure head (in [cm]) or constant flux (in [ kg/(m² s)])
            bcBotValue_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BC_Bot, Value);
        }

        if (bcTop_==4) {
            precTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<double>, Climate, Times);
            precData_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<double>, Climate, Precipitation); // in [cm / s]
        }

        /*
         * episode story
         */
        episodeNumber_=1;  // episode counter
        try { // Episodes defined
            episodeTimes_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<double>, TimeManager, Episodes);
            auto it = episodeTimes_.begin();
            episodeTimes_.insert(it,0.); // add initial
            double tend = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, TimeManager, TEnd);
            episodeTimes_.push_back(tend); // add final

        } catch(std::exception& e) {
            try { // Number of episodes defined
                int noe =  GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, TimeManager, NubmerOfEpisodes);
                double tend = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, TimeManager, TEnd);
                double dt;
                if (noe>1) {
                    dt = tend/(noe-1);
                } else {
                    dt = 0;
                }
                for (int i = 0; i<noe; i++) {
                    episodeTimes_.push_back(i*dt);
                }
            } catch(std::exception& e) { // undefined
                std::cout << "\nOptionally, specific output times in the input file can be defined in [TimeManager] Episodes = 200 1000 # [s] \n";
                episodeTimes_  = std::vector<double>(); // create empty vector
                auto it = episodeTimes_.begin();
                episodeTimes_.insert(it,0.); // add initial
                double tend = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, TimeManager, TEnd);
                episodeTimes_.push_back(tend); // add final
            }
        }
        double dt = episodeTimes_.at(episodeNumber_)-episodeTimes_.at(episodeNumber_-1);
        episodeNumber_++;
        this->timeManager().startNextEpisode(dt);
}

    /**
     * deactivate the use of restart file
     */
    bool shouldWriteRestartFile()
    {
        return false;
    }

    /**
     * write only initial, end, or episode times
     */
    bool shouldWriteOutput() const
    {
        if (this->timeManager().time() == 0 ||
                this->timeManager().willBeFinished() ||
                this->timeManager().episodeWillBeOver())
        {
            return true;
        } else {
            return false;
        }
    }

    /**
     * nothing really ends
     */
    void episodeEnd()
    {
        if (!this->timeManager().willBeFinished()) {
            double dt = episodeTimes_.at(episodeNumber_)-episodeTimes_.at(episodeNumber_-1);
            episodeNumber_++;
            this->timeManager().startNextEpisode(dt);
        }
    }

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    {
        return name_;
    }

    /*!
     * \brief Returns the temperature [K] within a finite volume
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 273.15 + 10; // -> 10°C // this function is called, and I don't know what temperature does in the model (probably nothing)
    };

    /*!
     * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar referencePressure(const Element &element, const FVElementGeometry &fvGeometry, const int scvIdx) const
    {
        return pnRef_;; // reference pressure when using Pascal
    };

    /*!
     * \brief Return the sources within the domain.
     *
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0;
    }

    /**
     * before each time step
     */
    void preTimeStep()
    {
       ParentType::preTimeStep();
       std::cout << "\npreTimeStep()\n\n";
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the boundary type is set
     */
    void boundaryTypesAtPos(BoundaryTypes &values, const GlobalPosition &globalPos) const
    {
//        cout << "\n\nBoundariesAtPos\n\n";

        if (globalPos[2]==topZ) { // top bc
            switch (bcTop_) {
            case 1: // constant pressure head
                //std::cout << "top Dirichlet \n";
                values.setAllDirichlet();
                break;
            case 2: // constant flux
                //std::cout << "top constant flux \n";
                values.setAllNeumann();
                break;
            case 4: // atmospheric boundary condition (with surface run-off)
                //std::cout << "top atmospheric \n";
                values.setAllNeumann();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Top boundary type not implemented");
            }
        } else if (globalPos[2]==botZ) { // bot bc
            switch (bcBot_) {
            case 1: // constant pressure head
                //std::cout << "bot Dirichlet \n";
                values.setAllDirichlet();
                break;
            case 2: // constant flux
                //std::cout << "bot constant flux \n";
                values.setAllNeumann();
                break;
            case 5: // free drainage
                //std::cout << "bot free drainage \n";
                values.setAllNeumann();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type not implemented");
            }
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        if (globalPos[2]==topZ) { // top bc
            switch (bcTop_) {
            case 1: // constant pressure
                values[hIdx] = pnRef_ - toPa_(bcTopValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Top boundary type Dirichlet: unknown error");
            }
        } else if (globalPos[2]==botZ) { // bot bc
            switch (bcBot_) {
            case 1: // constant pressure
                values[hIdx] = pnRef_ - toPa_(bcBotValue_);
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Dirichlet: unknown error");
            }
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local subcontrolvolume index
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    void solDependentNeumann(PrimaryVariables &values, const Element &element, const FVElementGeometry &fvGeometry, const Intersection &intersection,
            const int scvIdx, const int boundaryFaceIdx, const ElementVolumeVariables &elemVolVars) const
    {
        const double temp =  273.15 + 10;  // -> 10°C
        const double pnRef = 1e5; // Pa
        const double rho = WP::density(temp,pnRef); // kg/m³
        const double g = abs(this->gravity()[0]); // 1D
        double const atm = 1e5/(rho*g); // atmospheric pressure [Pa]

        GlobalPosition pos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        if ((pos[2]!=botZ) && (bcBot_==5)) { // free drainage at bottom boundary

            double Kc = this->spatialParams().hydraulicConductivity(element,fvGeometry,scvIdx);
            VolumeVariables  v0 = elemVolVars[0];
            VolumeVariables  v1 = elemVolVars[1];
            double swe = 0.5*(v0.saturation(wPhaseIdx) + v1.saturation(wPhaseIdx)); // TODO i take the mean because i don't know better
            double krw = MaterialLaw::krw(this->spatialParams().materialLawParams(element,fvGeometry,scvIdx), swe);
            values[contiEqIdx] = krw*Kc*rho; // * 1 [m]

        } else if((pos[2]==topZ) && (bcTop_==4)) { // atmospheric boundary condition (with surface run-off) at top

            double Kc = this->spatialParams().hydraulicConductivity(element,fvGeometry,scvIdx);
            VolumeVariables  v0 = elemVolVars[0];
            VolumeVariables  v1 = elemVolVars[1];
            double swe;
            //swe = 0.5*(v0.saturation(wPhaseIdx) + v1.saturation(wPhaseIdx)); // TODO i take the mean because i don't know better
            swe = v0.saturation(wPhaseIdx); // 0 is better....

            double krw = MaterialLaw::krw(this->spatialParams().materialLawParams(element,fvGeometry,scvIdx), swe);
            double h = MaterialLaw::pc(this->spatialParams().materialLawParams(element,fvGeometry,scvIdx), swe);
            h = - h/(rho*g); // from Pa -> m pressure head
            double dz = fvGeometry.elementVolume; // 1D

            double prec = getPrec_(this->timeManager().time()); // precipitation or evaporation
            if (prec<0) { // precipitation
                double imax = rho*Kc*((h-0.)/dz -1.); // maximal infiltration
                double v = (std::max(prec,imax));
                values[contiEqIdx] = v;
                std::cout << "\n its precipitation: "<< prec << ", max inf " << imax << " Swe "<< swe << " Pressurehead "<< h << " values " << v << " at time " << this->timeManager().time() <<"\n";
            } else { // evaporation
                double emax = rho*krw*Kc*((h-atm)/dz -1.); // maximal evaporation
                double v  = std::min(prec,emax);
                values[contiEqIdx] = v;
                std::cout << "\n its evaporation: "<< prec << ", max eva " << emax << " Swe "<< swe << " Pressurehead "<< h <<" values " << v << " at time " << this->timeManager().time() << "\n";
            }

        } else { // forward it to the interface without the volume variables
            this->neumann(values,element,fvGeometry,intersection,scvIdx,boundaryFaceIdx);
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * \param values The neumann values for the conservation equations
     * \param globalPos The position for which the Neumann value is set
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {

        values[contiEqIdx] = 0; return;

        if (globalPos[2]==0) { // top bc
            switch (bcTop_) {
            case 2: // constant flux
                std::cout << " top flux " << bcTopValue_ << " ";
                values[contiEqIdx] = -10*bcTopValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] * rho
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Top boundary type Neumann: unknown error");
            }
        } else { // bot bc
            switch (bcBot_) {
            case 2: // constant flux
                std::cout << " bot flux " << bcBotValue_<< " ";
                values[contiEqIdx] = -10*bcBotValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] *rho
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Neumann: unknown error");
            }
        }
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * initial usually calls initialAtPos, that must be overwritten in that case
     */
    void initial(PrimaryVariables &values, const Element &element, const FVElementGeometry &fvGeometry, const int scvIdx) const
    {

        values[hIdx] = 0.5; return;

        double iv = GridCreator::parameters(element).at(0);
        if (useHead) {
            values[hIdx] = iv;
        } else {
            values[hIdx] = pnRef_ - toPa_(iv);
        }
        //std::cout << values[hIdx] << "\n";
    }



private:

    // pressure head to pascal
    double toPa_(double ph) const
    {
        return -ph*10.*abs(this->gravity()[0]); // 1D
    }

    /*
     * returns the precipitation of the following data point.
     * e.g. (day_n, prec_n), means $t \in [day_{n-1}..day_n]$ will return prec_n
     * this makes sense, since normally data are given as mean precipitation per day
     */
    double getPrec_(double t) const
    {
        return getPrec_(t,0,precData_.size()-1);
    }

    // table look up (binary search O(log(n)))
    double getPrec_(double t, int l, int r) const
    {
        if ((t<=precTime_.at(l))||(l==r)) {
            return precData_.at(l);
        } else {
            int i = ceil((double(l)+double(r))/2.);
            if (t<=precTime_.at(i)) {
                return getPrec_(t,l+1,i);
            } else {
                return getPrec_(t,i,r);
            }
        }
    }

    double topZ = 1;
    double botZ = 0;

    int bcTop_;
    int bcBot_;
    double bcTopValue_;
    double bcBotValue_;
    std::vector<double> precTime_;
    std::vector<double> precData_;

    double pnRef_; // reference pressure
    std::string name_; // problem name
    std::vector<double> episodeTimes_;
    int episodeNumber_;

};

} //end namespace

#endif
