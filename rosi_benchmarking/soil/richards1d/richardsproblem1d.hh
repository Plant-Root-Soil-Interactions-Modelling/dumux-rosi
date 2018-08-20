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

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "richardsparams.hh"

#include <dumux/common/timeloop.hh>
#include <iostream>
#include <fstream>


namespace Dumux {

template <class TypeTag>
class RichardsProblem1d;

// Specify the properties for the lens problem
namespace Properties {

NEW_TYPE_TAG(RichardsProblem1d, INHERITS_FROM(Richards));
NEW_TYPE_TAG(RichardsBoxProblem1d, INHERITS_FROM(BoxModel, RichardsProblem1d));
NEW_TYPE_TAG(RichardsCCProblem1d, INHERITS_FROM(CCTpfaModel, RichardsProblem1d));

// Set the grid type
SET_TYPE_PROP(RichardsProblem1d, Grid, Dune::FoamGrid<1,1>);

// Set the physical problem to be solved
SET_TYPE_PROP(RichardsProblem1d, Problem, RichardsProblem1d<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(RichardsProblem1d, SpatialParams, RichardsParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                                               typename GET_PROP_TYPE(TypeTag, Scalar)>);
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
 * Episodes = 60 120 2000
 */
template <class TypeTag>
class RichardsProblem1d : public PorousMediumFlowProblem<TypeTag>
{
	using ParentType = PorousMediumFlowProblem<TypeTag>;
	using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
	using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
	using MaterialLaw = typename GET_PROP_TYPE(TypeTag, SpatialParams)::MaterialLaw;
	using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
	using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
	using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
	using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
	using Element = typename GridView::template Codim<0>::Entity;
	enum {
		// copy some indices for convenience
		pressureIdx = Indices::pressureIdx,
		conti0EqIdx = Indices::conti0EqIdx,
		bothPhases = Indices::bothPhases
	};
	using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
	using SubControlVolume = typename FVElementGeometry::SubControlVolume;
	using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
	using SourceValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
	static constexpr int dimWorld = GridView::dimensionworld;
	using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
	using Water = Components::SimpleH2O<Scalar>;
	using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
	using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
	using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
	using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace; // from fvproblem.hh
	using MaterialLawParams = typename MaterialLaw::Params;
	using Grid = typename GET_PROP_TYPE(TypeTag, Grid);

public:
	/*!
	 * \brief Constructor
	 *
	 * \param timeManager The Dumux TimeManager for simulation management.
	 * \param gridView The grid view on the spatial domain of the problem
	 */
	RichardsProblem1d(std::shared_ptr<const FVGridGeometry> fvGridGeometry, GridManager<Grid>* gm): ParentType(fvGridGeometry) {

	    gridmanager_ = gm;

		bcTop_ = getParam<int>("BC_Top.Type");
		bcBot_ = getParam<int>("BC_Bot.Type");

		bcTopValue_ = 0; // default value
		if ((bcTop_==1) || (bcTop_==2)) { // read constant pressure head (in [cm]) or constant flux (in [ kg/(m² s)])
			bcTopValue_ = getParam<Scalar>("BC_Top.Value");
		}

		bcBotValue_ = 0; // default value
		if ((bcBot_==1) || (bcBot_==2)) { // read constant pressure head (in [cm]) or constant flux (in [ kg/(m² s)])
			bcBotValue_ = getParam<Scalar>("BC_Bot.Value");
		}

		if (bcTop_==4) {
			precTime_ = getParam<std::vector<Scalar>>("Climate.Times");
			precData_ = getParam<std::vector<Scalar>>("Climate.Precipitation"); // in [cm / s]
		}

		 myfile_.open(this->name()+".txt");
	}

	~RichardsProblem1d() {
	    myfile_.close();
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
		return pRef_; // [Pa]
	};

	/*!
	 * \copydoc FVProblem::boundaryTypesAtPos
	 */
	BoundaryTypes boundaryTypesAtPos(const GlobalPosition &pos) const	{
		BoundaryTypes bcTypes;
		if (onUpperBoundary_(pos)) { // top bc
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
	PrimaryVariables dirichletAtPos(const GlobalPosition &pos) const
	{
		PrimaryVariables values;
		if (onUpperBoundary_(pos)) { // top bc
			switch (bcTop_) {
			case 1: // constant pressure
				values[pressureIdx] =toPa_(bcTopValue_);
				break;
			default:
				DUNE_THROW(Dune::InvalidStateException,"Top boundary type Dirichlet: unknown boundary type");
			}
		} else { // bot bc
			switch (bcBot_) {
			case 1: // constant pressure
				values[pressureIdx] = toPa_(bcBotValue_);
				break;
			default:
				DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Dirichlet: unknown boundary type");
			}
		}
		values.setState(bothPhases);
		return values;
	}

	/*!
	 * \copydoc FVProblem::neumann // [kg/(m²*s)]
	 */
	ResidualVector neumann(const Element& element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const SubControlVolumeFace& scvf) const {

		ResidualVector values;
		GlobalPosition pos = scvf.center();

		if (onUpperBoundary_(pos)) { // top boundary
			switch (bcTop_) {
			case 2: { // constant flux
				values[conti0EqIdx] = -10*bcTopValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] * rho
				break;
			}
			case 4: { // atmospheric boundary condition (with surface run-off)
				Scalar Kc = this->spatialParams().hydraulicConductivity(element);
				Scalar mS = 0;
				auto numScv = fvGeometry.numScv();
				for (auto i = 0; i<numScv; i++) {
					mS += (elemVolVars[i].saturation()/numScv);
				}
				MaterialLawParams params = this->spatialParams().materialLawParams(element);
				Scalar krw = MaterialLaw::krw(params, mS);
				Scalar p = MaterialLaw::pc(params, mS)+pRef_;
				Scalar h = - toHead_(p)/100.; // from Pa -> m pressure head

				GlobalPosition ePos = element.geometry().center();
				Scalar dz = 2 * std::abs(ePos[2]- pos[2]); // 0.01; // m // fvGeometry.geometry().volume()?;
				Scalar prec = getPrec_(time_); // precipitation or evaporation
				gridmanager_
				if (prec<0) { // precipitation
					Scalar imax = rho_*Kc*((h-0.)/dz -1.); // maximal infiltration
					Scalar v = std::max(prec,imax);
					values[conti0EqIdx] = v;
					std::cout << "\nprecipitation: "<< prec << ", max inf " << imax << " S "<< mS << " Pressurehead "<< h << " values " << v << " at time " << time_ << "dz" << dz;
				} else { // evaporation
					Scalar emax = rho_*krw*Kc*((h-(-100))/dz -1.); // maximal evaporation (-100 m = -10.000 cm)
					Scalar v  = std::min(prec,emax);
					values[conti0EqIdx] = v;
					std::cout << "\nevaporation: "<< prec << ", max eva " << emax << " S "<< mS << " Pressurehead "<< h <<" values " << v << " at time " << time_<< "dz" << dz;
				}
		    	// hack for benchmark 4
		    	myfile_ << time_ << ", "; //
		    	myfile_ << values[conti0EqIdx] << "\n";
				break;
			}
			default:
				DUNE_THROW(Dune::InvalidStateException,"Top boundary type Neumann: unknown error");
			}
		} else { // bottom boundary
			switch (bcBot_) {
			case 2: { // constant flux
				values[conti0EqIdx] = -10*bcBotValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] *rho
				break;
			}
			case 5: {// free drainage
				Scalar Kc = this->spatialParams().hydraulicConductivity(element);
				Scalar mS = 0; // mean saturation
				auto numScv = fvGeometry.numScv();
				for (auto i = 0; i<numScv; i++) {
					mS += (elemVolVars[i].saturation()/numScv);
				}
				MaterialLawParams params = this->spatialParams().materialLawParams(element);
				Scalar krw = MaterialLaw::krw(params, mS);
				values[conti0EqIdx] = krw*Kc*rho_; // * 1 [m]
				break;
			}
			default:
				DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Neumann: unknown error");
			}
		}
		return values;
	}

	/*!
	 * see FVProblem::initial
	 */
	template<class Entity>
	PrimaryVariables initial(const Entity& entity) const {
		PrimaryVariables values;
		Scalar iv = gridmanager_->getGridData()->parameters(entity).at(0); // TODO
		values[pressureIdx] = toPa_(iv);
		values.setState(bothPhases);
		return values;
	}

	/**
	 * Set current simulation time (within the simulation loop)
	 * for atmospheric look up
	 */
	void setTime(Scalar t) {
		time_ = t;
	}

private:

	// cm -> Pa
	Scalar toPa_(Scalar ph) const {
		return pRef_ +  ph/100.*rho_*g_;
	}

	Scalar toHead_(Scalar p) const {
		return (p-pRef_) * 100./rho_/g_;
	}


	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_;
	}

	/*
	 * returns the precipitation of the following data point.
	 * e.g. (day_n, prec_n), means $t \in [day_{n-1}..day_n]$ will return prec_n
	 * this makes sense, since normally data are given as gridmanager_mean precipitation per day
	 */
	Scalar getPrec_(Scalar t) const {
		return getPrec_(t,0,precData_.size()-1);
	}

	// table look up (binary search O(log(n)))
	Scalar getPrec_(Scalar t, int l, int r) const {
		if ((t<=precTime_.at(l))||(l==r)) {
			return precData_.at(l);
		} else {
			int i = ceil((Scalar(l)+Scalar(r))/2.);
			if (t<=precTime_.at(i)) {
				return getPrec_(t,l+1,i);
			} else {
				return getPrec_(t,i,r);
			}
		}
	}

	static constexpr Scalar g_ = 9.81; // cm / s^2
	static constexpr Scalar rho_ = 1.e3; // kg / m^3
	static constexpr Scalar pRef_ = 1.e5; // Pa

	static constexpr Scalar eps_ = 1e-8;

	int bcTop_;
	int bcBot_;
	Scalar bcTopValue_;
	Scalar bcBotValue_;
	std::vector<Scalar> precTime_;
	std::vector<Scalar> precData_;

	Scalar time_ = 0.;

	mutable std::ofstream myfile_; // file for flux over time

	GridManager<Grid>* gridmanager_;

};

} //end namespace Dumux

#endif
