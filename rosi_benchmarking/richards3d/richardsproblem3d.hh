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
#ifndef RICHARDS_PROBLEM3D_HH
#define RICHARDS_PROBLEM3D_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <rosi_benchmarking/richards3d/richardsparams.hh>

namespace Dumux
{
/*!
 * \ingroup RichardsTests
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards box model.
 */
template <class TypeTag>
class RichardsProblem3d;


// Specify the properties for the lens problem
namespace Properties
{
NEW_TYPE_TAG(RichardsProblem3d, INHERITS_FROM(Richards, RichardsParams));
NEW_TYPE_TAG(RichardsBoxProblem3d, INHERITS_FROM(BoxModel, RichardsProblem3d));
NEW_TYPE_TAG(RichardsCCProblem3d, INHERITS_FROM(CCTpfaModel, RichardsProblem3d));

// Use 2d YaspGrid
SET_TYPE_PROP(RichardsProblem3d, Grid, Dune::YaspGrid<3>);

// Set the physical problem to be solved
SET_TYPE_PROP(RichardsProblem3d, Problem, RichardsProblem3d<TypeTag>);
} // end namespace Dumux

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 *
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards model.
 *
 * The domain is box shaped. Left and right boundaries are Dirichlet
 * boundaries with fixed water pressure (fixed Saturation \f$S_w = 0\f$),
 * bottom boundary is closed (Neumann 0 boundary), the top boundary
 * (Neumann 0 boundary) is also closed except for infiltration
 * section, where water is infiltrating into an initially unsaturated
 * porous medium. This problem is very similar the the LensProblem
 * which uses the TwoPBoxModel, with the main difference being that
 * the domain is initally fully saturated by gas instead of water and
 * water instead of a %DNAPL infiltrates from the top.
 *
 * This problem uses the \ref RichardsModel
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxrichards -parameterFile test_boxrichards.input -TimeManager.TEnd 10000000</tt>
 * <tt>./test_ccrichards -parameterFile test_ccrichards.input -TimeManager.TEnd 10000000</tt>
 *
 * where the initial time step is 100 seconds, and the end of the
 * simulation time is 10,000,000 seconds (115.7 days)
 */
template <class TypeTag>
class RichardsProblem3d : public PorousMediumFlowProblem<TypeTag>
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
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Water = SimpleH2O<Scalar>;
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables); // from fvproblem.hh
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace; // from fvproblem.hh
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
	RichardsProblem3d(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry) {
		name_ = getParam<std::string>("Problem.Name");
		initial_ = getParam<Scalar>("Grid.Initial");
		bcTop_ = getParam<int>("BC_Top.Type");
		bcBot_ = getParam<int>("BC_Bot.Type");

		bcTopValue_ = 0; // default value
		if ((bcTop_==1) || (bcTop_==2)) { // read constant pressure if dirchlet head (in [cm]) or constant flux (in [ kg/(m² s)])
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
	}

	/*!
	 * \brief Returns the reference pressure [Pa] of the non-wetting
	 *        fluid phase within a finite volume
	 *
	 * This problem assumes a constant reference pressure of 1 bar.
	 */
	Scalar nonWettingReferencePressure() const {
		return 1.0e5;
	}

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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &pos) const	{
		//        cout << "\n\nBoundariesAtPos\n\n";
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann(); // default
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
		} else if (onLowerBoundary_(pos)) { // bot bc
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
				values[pressureIdx] = nonWettingReferencePressure() - toPa_(bcTopValue_);
				break;
			default:
				DUNE_THROW(Dune::InvalidStateException,"Top boundary type Dirichlet: unknown error");
			}
		} else if (onLowerBoundary_(pos)) { // bot bc
			switch (bcBot_) {
			case 1: // constant pressure
				values[pressureIdx] = nonWettingReferencePressure() - toPa_(bcBotValue_);
				break;
			default:
				DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Dirichlet: unknown error");
			}
		}
        values.setState(bothPhases);
		return values;
	}

	/*!
	 * see FVProblem::neumann
	 */
    ResidualVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf) const
    /*
     * dumux/discretization/box: fvelementgeometry.h (fe local basis), elementvolumevariables.h, subcontrolvolume.h, subcontrolvolumeface.h
     * dumux/porousmediumflow/richards/: volumevariables.h
     * Element (e.g. dune/grid/yaspgrid) Geometry (e.g. dune/grid/common)
     */
	{
		const Scalar rho = Water::liquidDensity(this->temperature(),nonWettingReferencePressure()); // h2o: 1000 kg/m³ Density of water(independent of temp and p)
		const Scalar g = 9.81; // abs(this->gravity()[dimWorld-1]);
		Scalar const atm = nonWettingReferencePressure()/(rho*g); // atmospheric pressure [Pa]
		GlobalPosition pos = scvf.center();
		ResidualVector values;
		values[conti0EqIdx] = 0; // default

		if (onUpperBoundary_(pos)) { // top boundary

			switch (bcTop_) {
			case 2: { // constant flux
				//std::cout << " top flux " << bcTopValue_ << " ";
				values[conti0EqIdx] = -10*bcTopValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] * rho
				break;
			}
			case 4: { // atmospheric boundary condition (with surface run-off)
				Scalar Kc = this->spatialParams().hydraulicConductivity(pos);
				Scalar mS = 0;
				auto numScv = fvGeometry.numScv();
				for (auto i = 0; i<numScv; i++) {
					mS += (elemVolVars[i].saturation(wPhaseIdx)/numScv);
				}
				MaterialLawParams params = this->spatialParams().materialLawParamsAtPos(pos);
				Scalar krw = MaterialLaw::krw(params, mS);
				Scalar p = MaterialLaw::pc(params, mS);
				Scalar h = - p/(rho*g); // from Pa -> m pressure head

				GlobalPosition ePos = element.geometry().center();
				Scalar dz = 2 * std::abs(ePos[2]- pos[2]); // 0.01; // m // fvGeometry.geometry().volume()?;
				Scalar prec = getPrec_(time_); // precipitation or evaporation
				if (prec<0) { // precipitation
					Scalar imax = rho*Kc*((h-0.)/dz -1.); // maximal infiltration
					Scalar v = std::max(prec,imax);
					values[conti0EqIdx] = v;
					std::cout << "\nprecipitation: "<< prec << ", max inf " << imax << " S "<< mS << " Pressurehead "<< h << " values " << v << " at time " << time_ ;
				} else { // evaporation
					Scalar emax = rho*krw*Kc*((h-atm)/dz -1.); // maximal evaporation
					Scalar v  = std::min(prec,emax);
					values[conti0EqIdx] = v;
					std::cout << "\nevaporation: "<< prec << ", max eva " << emax << " S "<< mS << " Pressurehead "<< h <<" values " << v << " at time " << time_;
				}
				break;
			}
			default:
				DUNE_THROW(Dune::InvalidStateException,"Top boundary type Neumann: unknown error");
			}

		} else if (onLowerBoundary_(pos)) { // bot boundary

//			auto eIndex = fvGeometry.fvGridGeometry().elementMapper().index(element);
//			auto scvIndex = scvf.insideScvIdx();
//			GlobalPosition ePos = element.geometry().center();
//			std::cout << "Top boundary condition: Center of subcontrolvolume face with local index " << scvIndex << " is (" << pos << "), at element with index " <<
//					eIndex << " and vertex position (" << ePos << "), Number of subdcontrolvolumes:  " << fvGeometry.numScv() << "\n";

//			SubControlVolume scv = fvGeometry.scv(scvIndex);
//			GlobalPosition dofPos = scv.dofPosition();

			switch (bcBot_) {
			case 2: { // constant flux
				//std::cout << " bot flux " << bcBotValue_<< " ";
				values[conti0EqIdx] = -10*bcBotValue_/(24.*60.*60.); // [kg/(m²*s)] = 1/10 [cm/s] *rho
				break;
			}
			case 5: {// free drainage
				Scalar Kc = this->spatialParams().hydraulicConductivity(pos);
				Scalar mS = 0;
				auto numScv = fvGeometry.numScv();
				for (auto i = 0; i<numScv; i++) {
					mS += (elemVolVars[i].saturation(wPhaseIdx)/numScv);
				}
				MaterialLawParams params = this->spatialParams().materialLawParamsAtPos(pos);
				Scalar krw = MaterialLaw::krw(params, mS);
				values[conti0EqIdx] = krw*Kc*rho; // * 1 [m]
				break;
			}
			default:
				DUNE_THROW(Dune::InvalidStateException,"Bottom boundary type Neumann: unknown error");
			}

		}
		return values;
	}

	/*!
	 * \name Volume terms
	 */
	// \{

	/*!
	 * see FVProblem::initial
	 */
	template<class Entity>
	PrimaryVariables initial(const Entity& entity) const
	{
		PrimaryVariables values;
		if (initial_<0) {
			Scalar iv = initial_;
			values[pressureIdx] = nonWettingReferencePressure() - toPa_(iv);
			// std::cout << values[pressureIdx] << "\n";
			values.setState(bothPhases);
		} else { // hard coded for example jan2
			GlobalPosition pos = entity.geometry().center();
			Scalar iv = -54 + pos[2];
			values[pressureIdx] = nonWettingReferencePressure() - toPa_(iv);
			values.setState(bothPhases);
		}
		return values;
	}

	void setTime(Scalar t) {
		time_ = t;
	}

	// \}

private:

	/**
	 *  pressure head to pascal
	 *
	 *  @param ph           pressure head [cm]
	 */
	Scalar toPa_(Scalar ph) const {
		const Scalar g = 9.81; // TODO stupid gravity bug
		return -ph*10.*g;
	}

	bool onLowerBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld-1] < this->fvGridGeometry().bBoxMin()[dimWorld-1] + eps_;
	}

	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_;
	}

	/**
	 * returns the precipitation of the following data point.
	 * e.g. (day_n, prec_n), means $t \in [day_{n-1}..day_n]$ will return prec_n
	 * this makes sense, since normally data are given as mean precipitation per day
	 */
	Scalar getPrec_(Scalar t) const {
		return getPrec_(t,0,precData_.size()-1);
	}

	/**
	 * table look up (binary search O(log(n)))
	 */
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

	int bcTop_;
	int bcBot_;
	Scalar bcTopValue_;
	Scalar bcBotValue_;
	std::vector<Scalar> precTime_;
	std::vector<Scalar> precData_;

	std::string name_; // problem name

	Scalar time_ = 0.;

	Scalar initial_ = -200; // cm pressure head
	Scalar depth_ = -1; // m // todo useg gridview. bBoxMin Max instead
	Scalar eps_ = 1.e-5; // m

};

} //end namespace Dumux

#endif
