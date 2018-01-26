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
 * \ingroup RichardsTests
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards box model.
 */
#ifndef RICHARDS_PROBLEM1D_HH
#define RICHARDS_PROBLEM1D_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <rosi_benchmarking/richards1d/richardsparams.hh>
#include <dune/foamgrid/foamgrid.hh>


namespace Dumux
{
template <class TypeTag>
class RichardsProblem1d;


// Specify the properties for the lens problem
namespace Properties
{
NEW_TYPE_TAG(RichardsProblem1d, INHERITS_FROM(Richards, RichardsParams));
NEW_TYPE_TAG(RichardsBoxProblem1d, INHERITS_FROM(BoxModel, RichardsProblem1d));
NEW_TYPE_TAG(RichardsCCProblem1d, INHERITS_FROM(CCTpfaModel, RichardsProblem1d));

// Set the grid type
SET_TYPE_PROP(RichardsProblem1d, Grid, Dune::FoamGrid<1,1>);
// apparently Dune::OneDGrid does not read the simplex parameters from the dgf (without warning),
// no idea how to pass initial conditions otherwise. FoamGrid works fine!
// SGrid is no longer available, or I don't know how to use it

// Set the physical problem to be solved
SET_TYPE_PROP(RichardsProblem1d, Problem, RichardsProblem1d<TypeTag>);
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
 * [TimeLoop]
 * Episodes = 60 120 2000 # s
 */
template <class TypeTag>
class RichardsProblem1d : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Element = typename GridView::template Codim<0>::Entity;
    enum {
        // copy some indices for convenience
        pressureIdx = Indices::pressureIdx,
        conti0EqIdx = Indices::conti0EqIdx,
		wPhaseIdx = Indices::wPhaseIdx,
        bothPhases = Indices::bothPhases
    };
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Water = SimpleH2O<Scalar>;
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables); // from fvproblem.hh
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace); // from fvproblem.hh
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The Dumux TimeManager for simulation management.
     * \param gridView The grid view on the spatial domain of the problem
     */
    RichardsProblem1d(std::shared_ptr<const FVGridGeometry> fvGridGeometry): ParentType(fvGridGeometry) {
		pnRef_ = 1e5; // reference pressure if Pascal are used

		/*
		 * read relevant run time parameter
		 */
		name_ = getParam<std::string>("Problem.Name");
		bcTop_ = getParam<int>("BC_Top.Type");
		bcBot_ = getParam<int>("BC_Bot.Type");

		bcTopValue_ = 0; // default value
		if ((bcTop_==1) || (bcTop_==2)) { // read constant pressure if dirchlet head (in [cm]) or constant flux (in [ kg/(m² s)])
			bcTopValue_ = getParam<double>("BC_Top.Value");
		}

		bcBotValue_ = 0; // default value
		if ((bcBot_==1) || (bcBot_==2)) { // read constant pressure head (in [cm]) or constant flux (in [ kg/(m² s)])
			bcBotValue_ = getParam<double>("BC_Bot.Value");
		}

		if (bcTop_==4) {
			precTime_ = getParam<std::vector<double>>("Climate.Times");
			precData_ = getParam<std::vector<double>>("Climate.Precipitation"); // in [cm / s]
		}
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const {
    	return name_;
    }

    /*!
     * \brief Returns the temperature [K] within a finite volume
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const {
    	return 273.15 + 10; // -> 10°C
    };

    /*!
     * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonWettingReferencePressure() const {
    	return 1.0e5; // [Pa]
    };

    // \}

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

	/*!
	 * \copydoc FVProblem::sourceAtPos
	 */
    ResidualVector sourceAtPos(const GlobalPosition &globalPos) const {
        //! As a default, i.e. if the user's problem does not overload any source method
        //! return 0.0 (no source terms)
        return ResidualVector(0.0);
	}

	/*!
	 * \copydoc FVProblem::boundaryTypesAtPos
	 */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const	{
		//        cout << "\n\nBoundariesAtPos\n\n";
        BoundaryTypes bcTypes;
    	if (globalPos[0]==0) { // top bc
			switch (bcTop_) {
			case 1: // constant pressure head
				bcTypes.setAllDirichlet();
				break;
			case 2: // constant flux
				bcTypes.setAllNeumann();
				break;
			case 4: // atmospheric boundary condition (with surface run-off)
				bcTypes.setAllNeumann();
				break;
			default:
				DUNE_THROW(Dune::InvalidStateException,"Top boundary type not implemented");
			}
		} else { // bot bc
			switch (bcBot_) {
			case 1: // constant pressure head
				bcTypes.setAllDirichlet();
				break;
			case 2: // constant flux
				bcTypes.setAllNeumann();
				break;
			case 5: // free drainage
				bcTypes.setAllNeumann();
				break;
			default:
				DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type not implemented");
			}
		}
    	return bcTypes;
	}

	/*!
	 * \copydoc FVProblem::dirichletAtPos
	 */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
	{
    	PrimaryVariables values;
		if (globalPos[0]==0) { // top bc
			switch (bcTop_) {
			case 1: // constant pressure
				values[pressureIdx] = pnRef_ - toPa_(bcTopValue_);
				break;
			default:
				DUNE_THROW(Dune::InvalidStateException,"Top boundary type Dirichlet: unknown error");
			}
		} else { // bot bc
			switch (bcBot_) {
			case 1: // constant pressure
				values[pressureIdx] = pnRef_ - toPa_(bcBotValue_);
				break;
			default:
				DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Dirichlet: unknown error");
			}
		}
        values.setState(bothPhases); // wtf?
		return values;
	}

	/*!
	 * see FVProblem::neumann
	 */
    ResidualVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf) const
	{
		const double rho = Water::liquidDensity(this->temperature(),pnRef_); // h2o: 1000 kg/m³ Density of water(independent of temp and p)
//		const double g = abs(this->gravity()[0]); // 1D
//		double const atm = 1e5/(rho*g); // atmospheric pressure [Pa]

		GlobalPosition pos = scvf.ipGlobal();

		ResidualVector values;
		if (pos[0]!=0) { // bottom boundary
			switch (bcBot_) {
			case 2: { // constant flux
				//std::cout << " bot flux " << bcBotValue_<< " ";
				values[conti0EqIdx] = -10*bcBotValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] *rho
				break;
			}
			case 5: {// free drainage
				double Kc = this->spatialParams().hydraulicConductivity(element);
				VolumeVariables  v0 = elemVolVars[0];
				VolumeVariables  v1 = elemVolVars[1];
				double swe = 0.5*(v0.saturation(wPhaseIdx) + v1.saturation(wPhaseIdx)); // TODO i take the mean because i don't know better
				ElementSolutionVector esv;
				SubControlVolume scv;
				MaterialLawParams params = this->spatialParams().materialLawParams(element, scv, esv);
				double krw = MaterialLaw::krw(params, swe);
				values[conti0EqIdx] = krw*Kc*rho; // * 1 [m]
				break;
			}
			default:
				DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Neumann: unknown error");
			}
		} else if(pos[0]==0) { // top boundary

			switch (bcTop_) {
			case 2: { // constant flux
				//std::cout << " top flux " << bcTopValue_ << " ";
				values[conti0EqIdx] = -10*bcTopValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] * rho
				break;
			}
			case 4: { // atmospheric boundary condition (with surface run-off)
//				double Kc = this->spatialParams().hydraulicConductivity(element);
//				ElementVolumeVariables  v0 = elemVolVars[0];
//				ElementVolumeVariables  v1 = elemVolVars[1];
//				double swe;
//				//swe = 0.5*(v0.saturation(wPhaseIdx) + v1.saturation(wPhaseIdx)); // TODO i take the mean because i don't know better
//				swe = v0.saturation(0); // todo is it wPhaseIdx?
//
//				double krw = MaterialLaw::krw(this->spatialParams().materialLawParams(element), swe);
//				double h = MaterialLaw::pc(this->spatialParams().materialLawParams(element), swe);
//				h = - h/(rho*g); // from Pa -> m pressure head
//				double dz = 0.01; // m // todo no idea how this works ... fvGeometry.elementVolume; // 1D
//
//				double prec = getPrec_(0); // precipitation or evaporation  todo how do I get the time?
//				if (prec<0) { // precipitation
//					double imax = rho*Kc*((h-0.)/dz -1.); // maximal infiltration
//					double v = (std::max(prec,imax));
//					values[conti0EqIdx] = v;
//					std::cout << "\n its precipitation: "<< prec << ", max inf " << imax << " Swe "<< swe << " Pressurehead "<< h << " values " << v
//							<< " at time " << 0 <<"\n"; // todo time
//				} else { // evaporation
//					double emax = rho*krw*Kc*((h-atm)/dz -1.); // maximal evaporation
//					double v  = std::min(prec,emax);
//					values[conti0EqIdx] = v;
//					std::cout << "\n its evaporation: "<< prec << ", max eva " << emax << " Swe "<< swe << " Pressurehead "<< h <<" values " << v
//							<< " at time " << 0 << "\n"; // todo time
//				}
				break;
			}
			default:
				DUNE_THROW(Dune::InvalidStateException,"Top boundary type Neumann: unknown error");
			}
		}
		return values;
	}

	/*!
	 * see FVProblem::initial
	 */
    template<class Entity>
    PrimaryVariables initial(const Entity& entity) const
	{
    	std::cout << "initial()\n";
		PrimaryVariables values;
		double iv = GridCreator::parameters(entity).at(0); // todo thats wrong, is it? wahts entity and how i get my element from it
		values[pressureIdx] = pnRef_ - toPa_(iv);
		std::cout << values[pressureIdx] << "\n";
        values.setState(bothPhases); // wtf?
		return values;
	}


private:

	/**
	 *  pressure head to pascal
	 *
	 *  @param ph           pressure head [cm]
	 */
	double toPa_(double ph) const {
		return -ph*10.*abs(this->gravity()[0]);
	}

	/**
	 * returns the precipitation of the following data point.
	 * e.g. (day_n, prec_n), means $t \in [day_{n-1}..day_n]$ will return prec_n
	 * this makes sense, since normally data are given as mean precipitation per day
	 */
	double getPrec_(double t) const {
		return getPrec_(t,0,precData_.size()-1);
	}

	/**
	 * table look up (binary search O(log(n)))
	 */
	double getPrec_(double t, int l, int r) const {
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

	int bcTop_;
	int bcBot_;
	double bcTopValue_;
	double bcBotValue_;
	std::vector<double> precTime_;
	std::vector<double> precData_;

	double pnRef_; // reference pressure

	std::string name_; // problem name
};

} //end namespace Dumux

#endif
