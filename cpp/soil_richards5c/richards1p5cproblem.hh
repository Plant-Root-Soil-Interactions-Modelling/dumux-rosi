// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef RICHARDS1P5C_PROBLEM_HH
#define RICHARDS1P5C_PROBLEM_HH
#include <algorithm>
#include <vector>
#include <dumux/porousmediumflow/problem.hh> // base class
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>

#include "../soil_richards/richardsparams.hh"

#include <dune/common/exceptions.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>


// getDofIndices, getPointIndices, getCellIndices
#include <dune/grid/utility/globalindexset.hh>

namespace Dumux {

/*!
 * RichardsProblem:
 * Uses Dumux as an easy to use Richards equation solver,
 * where most parameters can be set dynamically
 */
template <class TypeTag>
class Richards1P5CProblem : public PorousMediumFlowProblem<TypeTag>
{
public:
    using ParentType = PorousMediumFlowProblem<TypeTag>;
	
	// exports, used by the binding
	using Grid = GetPropType<TypeTag, Properties::Grid>;
	using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
	using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
	using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
	using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
	using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;

	// other
	using GridView = typename FVGridGeometry::GridView;
	using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
	
	using NumEqVector = typename Dumux::NumEqVector<PrimaryVariables>;
	using FVElementGeometry = typename FVGridGeometry::LocalView;
	using SubControlVolume = typename FVGridGeometry::SubControlVolume;
	using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
	using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
	//using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
	using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

	using Scalar = GetPropType<TypeTag, Properties::Scalar>;
	using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
	using Element = typename GridView::template Codim<0>::Entity;
	using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

																													   
	using PcKrSwCurve = Dumux::FluidMatrix::VanGenuchtenDefault<Scalar>;
	// VanGenuchtenNoReg<double>  // both of type TwoPMaterialLaw // using MaterialLaw = typename GetPropType<TypeTag, Properties::SpatialParams>::MaterialLaw;
	using BasicParams = typename PcKrSwCurve::BasicParams;
    using EffToAbsParams = typename PcKrSwCurve::EffToAbsParams;
    using RegularizationParams = typename PcKrSwCurve::RegularizationParams;

	using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;

	using PointSource = GetPropType<TypeTag, Properties::PointSource>;
	using CouplingManager= GetPropType<TypeTag, Properties::CouplingManager>;
	using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
	using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
																   
	using EffectiveDiffusivityModel = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
	static constexpr bool isBox = FVGridGeometry::discMethod == DiscretizationMethods::box;

	
	
    static constexpr bool useMoles = true;
	
    static constexpr int numFluidComps = FluidSystem::numComponents;
    static constexpr int numSolidComps = SolidSystem::numComponents;
    static constexpr int numInertSolidComps =  SolidSystem::numInertComponents;

	enum {
		// W elements
		pressureIdx = 0, // index of primary variables
		h2OIdx = FluidSystem::liquidPhaseIdx, // fluid index
		
		//in W phase
		soluteIdx = 1, // 1st solute index
		//mucilIdx = 3,//not used for now. no advection for eqIdx == mucilIdx
		//in solid phase
		CSS2Idx = 2,//adsorbed C
		decayIdx = 3,//decayed C
				
		//don t think EqIndx != pvdIdx when we have just 1 phase
		conti0EqIdx = pressureIdx, // indices of the equations

		dimWorld = GridView::dimensionworld,
		
		
		numFluidSolutes = numFluidComps - 1,//all minus water
		numSolidSolutes = numSolidComps -  1,// all minus sil
		numSolutes = numFluidSolutes + numSolidSolutes,


		//!!!!!!!!!!!!
		numComponents_ = numFluidComps + numSolidComps - numInertSolidComps//ignore the soil as component as it is inert
		//!!!!!!!!!!!!
	};

	enum dzScalingType {
		dx_2 = 1,
		dx = 2
	};
	enum BCTypes {
		constantPressure = 1,
		constantConcentration = 1,
		constantFlux = 2,
		constantFluxCyl = 3,
		atmospheric = 4,
		freeDrainage = 5,
		outflow = 6,
		linear = 7,
		michaelisMenten = 8,
		managed = 9,
		outflowCyl = 10
	};

	enum GridParameterIndex {
		materialLayerNumber = 0
	};
	
	
	
		
	/*!
	 * \brief Constructor: constructed in the main file
	 */
	Richards1P5CProblem(std::shared_ptr<const FVGridGeometry> gridGeometry)
	: PorousMediumFlowProblem<TypeTag>(gridGeometry) {
		
		gravityOn_ = Dumux::getParam<bool>("Problem.EnableGravity", (dimWorld > 1));
		temperatureK = Dumux::getParam<double>("Problem.temperatureK", 273.15+10);
		verbose_local_residual = Dumux::getParam<bool>("Problem.verbose_local_residual", verbose_local_residual);

		source_.resize(numComponents_); // numComponents_ equations (currently hard coded, where can I get the value?)
		
		 verbose =  getParam<int>("Problem.verbose", 0);
		 toFile =  getParam<bool>("Problem.toFile", false);
		doSoluteFlow =   getParam<bool>("Problem.doSoluteFlow", doSoluteFlow);
		reactionExclusive=   getParam<bool>("Problem.reactionExclusive",(dimWorld > 1));
				dzScaling = getParam<int>("Soil.BC.dzScaling", 2); 
		dobioChemicalReaction_ = getParam<bool>("Problem.dobioChemicalReaction",true);
		doDecay_ = getParam<bool>("Problem.doDecay",true);
		
		// Uptake params
		vMax =  getParam<Scalar>("RootSystem.Uptake.Vmax", 6.2e-11/1e4*(24.*3600.))*1e4/(24.*3600.); //  [mol cm-2 day-1] -> [mol m-2 s-1]
		km = getParam<Scalar>("RootSystem.Uptake.Km", 3.1e-9 /1e6 )*1e6;  // [mol cm-3] -> [mol m-3]
		
		for(int i = 0; i < numComponents_; i++)//all components except h2o
		{
			//std::cout<<"PorousMediumFlowProblem ID"<<i<<" ";
			
			source_[i] = nullptr;
			if(i ==h2OIdx)
			{
				
				// BC			
				bcTopType_ = getParam<int>("Soil.BC.Top.Type", outflow); 
				bcBotType_ = getParam<int>("Soil.BC.Bot.Type", outflow);
				
				//const auto& myParams = Parameters::paramTree();
				//myParams.report();
				bcTopValues_.at(i) =  getParam<double>("Soil.BC.Top.Value", 0.);
				bcBotValues_.at(i) =  getParam<double>("Soil.BC.Bot.Value", 0.);
				//std::cout<<"start h2o "<<bcTopValues_.at(i)<<" "<<bcBotValues_.at(i)<<" ";
				
				//IC
				//std::cout<<"initialSoil_"<<std::endl;
				initialSoil_.at(i) = InputFileFunction("Soil.IC", "P", "Z", 0., this->spatialParams().layerIFF()); // [cm]([m]) pressure head, conversions hard coded
				// Precipitation & Evaporation, and solute input
				if (bcTopType_==atmospheric) {
					componentInput_.at(i) = InputFileFunction("Climate", "Precipitation", "Time", 0.); // cm/day (day)
				}
			}else{
				
				// BC			
				bcSTopType_.at(i - soluteIdx) = getParam<int>("Soil.BC.Top.C"+std::to_string(i)+"Type", outflow); 
				bcSBotType_.at(i - soluteIdx) = getParam<int>("Soil.BC.Bot.C"+std::to_string(i)+"Type", outflow);
				bcSTopValue_.at(i - soluteIdx) = getParam<Scalar>("Soil.BC.Top.C"+std::to_string(i)+"Value", 0.);
				bcSBotValue_.at(i - soluteIdx) = getParam<Scalar>("Soil.BC.Bot.C"+std::to_string(i)+"Value", 0.);
				
				//IC
				initialSoil_.at(i) = InputFileFunction("Soil.IC", "C"+std::to_string(i), "C"+std::to_string(i)+"Z", 
													0., this->spatialParams().layerIFF()); // kg/kg or mol/mol soil
				if (bcSTopType_.at(i - soluteIdx)==managed) {
					componentInput_.at(i) = InputFileFunction(std::to_string(i)+".Component.Managed", "Input", "Time", 0.); // cm/day (day)
					
				}
				
			}
			
			
			componentInput_.at(i).setVariableScale(1./(24.*60.*60.)); //day-> s  
			Scalar g2kg = 1/1000 ;
			Scalar m2_2_cm2 = 10000;
			Scalar unitConversion = useMoles ? m2_2_cm2 : m2_2_cm2 * g2kg; //something else needed? 
			componentInput_.at(i).setFunctionScale(unitConversion/(24.*60.*60.)); // g/(cm2 day) or mol/(cm2 day)  -> kg/(m²*s) or mol/(m²*s) 
			//std::cout<<std::endl;
		}

		criticalPressure_ = getParam<double>("Soil.CriticalPressure", -1.e4); // cm
		criticalPressure_ = getParam<double>("Climate.CriticalPressure", criticalPressure_); // cm
		sourceSlope = getParam<double>("Soil.SourceSlope", -1.); // cm, negative value disables regularisation
		
		double m3_2_cm3 = 1e6;// cm3/m3
		f_sorp = getParam<double>("Soil.f_sorp", f_sorp);//[-]
		k_sorp = getParam<double>("Soil.k_sorp", k_sorp) ;// mol/cm3 water or mol
		css1Function = getParam<double>("Soil.css1Function", 0);// 0: no css1 , 1: linear css1
		CSSmax = getParam<double>("Soil.CSSmax", CSSmax)*m3_2_cm3;//mol/cm3 => [mol/m3]	
		//alpha = getParam<double>("Soil.alpha", alpha)/(24.*60.*60.);//[d-1] => [s-1]
		//kads = getParam<double>("Soil.kads", kads)  /m3_2_cm3;//[cm3/mol] => [m3/mol] 
		kdes = getParam<double>("Soil.kdes", kdes)/(24.*60.*60.);//[d-1] => [s-1]         
		vmax_decay = getParam<double>("Soil.vmax_decay", vmax_decay);//#[mol C / m^3 scv / s]
		km_decay = getParam<double>("Soil.km_decay", km_decay);//[mol C / m^3 scv] 
		
		// Buffer power
		b_ = getParam<Scalar>("Component.BufferPower", 0.);
		freundlichN_ = getParam<Scalar>("Component.FreundlichN", 0.);
		freundlichK_ = getParam<Scalar>("Component.FreundlichK", 0.);
		
		// todo: check if need to re-implement the lines below
		// if (css1Function != 8){
			// CSSmax = CSSmax * m3_2_cm3; //mol/cm3 zone 1 to //mol/m3 zone 1
			// k_sorp = k_sorp * m3_2_cm3;//mol/cm3 water to mol/m3 water 
		// }
		// if(css1Function == 9)
		kads = getParam<double>("Soil.kads", kads)  /(24.*60.*60.);//[cm3/mol/d] => [cm3/mol/s] or [1/d] => [1/s] 
         if(css1Function == 9)
         {
             kads *= 1/m3_2_cm3;//[cm3/mol/s] / [cm3/m3] => [m3/mol/s] 
         }      


		///
		computedCellVolumesCyl = false;
		if(dimWorld == 1)
		{
			segLength  = Dumux::getParam<double>("Problem.segLength")/100;
			this->setVolumesCyl(computeCellVolumesCyl_());
			computedCellVolumesCyl = true;
		}else{
			segLength  = -1.;
		}

	}

	double getCellVolumesCyl(int dofIndex) const
	{
		
		if(!computedCellVolumesCyl) // if we want to change the segLEngth without recreating the object?
		{
			DUNE_THROW(Dune::InvalidStateException, "getCellVolumesCyl: cellVolumesCyl not set");
		}
		return cellVolumesCyl.at(dofIndex);
	}

    void setVolumesCyl(std::vector<double> cellvoles) const
	{
		const_cast<std::vector<double>&>(cellVolumesCyl) = cellvoles;
    }

    /**
     * The volume [m3] of each element (vtk cell)
     *
     * This is done for a single process, gathering and mapping is done in Python.
     */
    std::vector<double> computeCellVolumesCyl_() const
	{
        std::vector<double> vols;
		auto points = this->getPoints_();//get the vertices == faces of 1D domain
        for (int i = 0; i < (points.size()-1); i++) {
			double rIn = points.at(i).at(0);
			double rOut = points.at(i + 1).at(0);
            vols.push_back((rOut*rOut - rIn*rIn)* M_PI * segLength);
        }
        return vols;
    }

    /**
     * Returns the Dune vertices (vtk points) of the grid for a single mpi process.
     * Gathering and mapping is done in Python.
     */
    std::vector<std::array<double,1>> getPoints_() const
	{
		//auto eIdx = this->gridGeometry().elementMapper().index(entity);
		//Scalar z = entity.geometry().center()[dimWorld - 1];
		
        std::vector<std::array<double,1>> points;
        points.reserve(this->gridGeometry().gridView().size(dimWorld));
        for (const auto& v : vertices(this->gridGeometry().gridView())) {
            auto p = v.geometry().center();
            std::array<double,1> vp;
            for (int i=0; i<dimWorld; i++) { // found no better way
                vp[i] = p[i];
            }
            points.push_back(vp);
        }
        return points;
    }
	/**
	 * \brief Eventually, closes output file
	 */
	~Richards1P5CProblem() {
		//std::cout << "closing file \n";
		//myfile_.close();
	}

	/*!
	 * \brief Temperature [K] within a finite volume. This problem assumes a temperature of 10 degrees Celsius.
	 *
	 * called EnergyVolumeVariablesImplementation::updateTemperature(...) in porousmediumflow/nonisothermal/volumevariables.hh,
	 * included by porousmediumflow/volumevariables.hh,
	 *
	 * todo this makes very little sense for isothermal!
	 *
	 * overwrites PorousMediumFlowProblem::temperature (compiles without, throws exception of base class)
	 */
	Scalar temperature() const {
		return temperatureK; //273.15 + 10; // -> 10°C
	}

	/*!
	 * \brief Reference pressure [Pa] of the non-wetting. This problem assumes a constant reference pressure of 1 bar.
	 *
	 * called by porousmediumflow/richards/volumevariables.hh
	 */
	Scalar nonWettingReferencePressure() const {return pRef_;}
	Scalar nonwettingReferencePressure() const {
		return nonWettingReferencePressure();
	}
	
	Scalar massOrMoleDensity(const auto& volVars, const int compIdx, const bool isFluid) const
	{
		return isFluid ? (useMoles ? volVars.molarDensity(compIdx) : volVars.density(compIdx) ):
				(useMoles ? volVars.solidComponentMolarDensity(compIdx) : volVars.solidComponentDensity(compIdx) ); 
	}

	Scalar massOrMoleFraction(const auto& volVars, const int phaseIdx, const int compIdx, const bool isFluid) const
	{
		return isFluid ?( useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx) ): 
				(useMoles ? volVars.solidMoleFraction(compIdx) : volVars.solidMassFraction(compIdx)); 
	}

	/**
	 * The buffer power for a scv for a volVar (linear in this implementation), equals $\rho_b K_d$ in Eqn (4) in phosphate draft
	 *
	 * used by my the modified localresidual.hh (see dumux-rosi/dumux/porousmediumflow/compositional)
	 */
	Scalar bufferPower(const VolumeVariables& volVars, int compIdx, int dofIndex) const {
		
		switch(compIdx)
		{
			case h2OIdx:{
				DUNE_THROW(Dune::InvalidStateException, "bufferPower used for water");
				break;
			}
			case soluteIdx:{
				
					double C_S_W = massOrMoleDensity(volVars, h2OIdx, true) * std::max(massOrMoleFraction(volVars,h2OIdx, compIdx, true), 0.);//mol C/m3 soil water

					return this->computeB(volVars, C_S_W, dofIndex);
				// if (getB_(volVars, dofIndex)>0.) {//b_
					// return b_;
				// } else {
					// if (getFreundlichK_(volVars, dofIndex)==0.) {//freundlichK_
						// return 0.;
					// }
					
					
				// }
			}
			default:
			{
				return 0.;
			}
		}
		
	}
	
	double getFreundlichK_(const VolumeVariables& volVars, int dofIndex) const 
	{
		double Css = massOrMoleFraction(volVars,0, 0, false);
		if ((Css > 0.)||((source_[2] != nullptr)&&(source_[2]->at(dofIndex) != 0.))) { return 0.;}
		else{return freundlichK_;}
	}
	
	double getB_(const VolumeVariables& volVars, int dofIndex) const 
	{
		double Css = massOrMoleFraction(volVars,0, 0, false);
		if ((Css > 0.)||((source_[2] != nullptr)&&(source_[2]->at(dofIndex) != 0.))) { return 0.;}
		else{return b_;}
	}

	double computeB(const VolumeVariables& volVars, double C_S_W, int dofIndex) const
	{
		//double CSS1 = this->computeCSS1(C_S_W, theta); // [ mol / m^3 zone 1]
		//return 1+(CSS1*f_sorp)/(C_S_W*theta);// = 1+f * (Cssmax*csw)/(csw+ksrop)^2
		
		//mol/m3
			double solidDensity = massOrMoleDensity(volVars, SolidSystem::mainCompIdx , false);
			double solVolFr = (1. - volVars.porosity());
			double bulkSoilDensity = solidDensity * solVolFr;
			double C_S_S2 = std::max(massOrMoleFraction(volVars,0, CSS2Idx - numFluidComps, false) //mol C/mol solid soil
							, 0.);
			if((C_S_W <= 0)||(C_S_S2>0.)){return 0.;}
			 
		return computeCSS1(bulkSoilDensity, C_S_W, dofIndex)/C_S_W;//CSSmax*C_S_W/(k_sorp+C_S_W);//bulkDensity_*freundlichK_*std::pow(C_S_W, freundlichN_)/C_S_W;
	}
	
	double computeCSS1(double bulkSoilDensity, double C_S_W, int dofIndex) const
	{// mol		
	
		//double Css = massOrMoleFraction(volVars,0, 0, false);
		//if ((Css > 0.)||((source_[2] != nullptr)&&(source_[2]->at(dofIndex) != 0.)) { return 0.;}
		if (((source_[2] != nullptr)&&(source_[2]->at(dofIndex) != 0.))||(C_S_W <= 0)) { return 0.;}
		switch(css1Function) {
		  case 0:
			return 0.;//none
		  case 1:
			return CSSmax*C_S_W/(k_sorp+C_S_W);//langumuir
		  case 2:
			// [mol/m3] = [mol/m3] * [cm^{3*n}/mol^n] * [mol^n/cm^{3*n}]
			// see DeBauw et al. 2020 supplementary, Eq. 2
			// to do: re-check
			return bulkSoilDensity * freundlichK_*std::pow(C_S_W, freundlichN_);  // freundlich
		  case 4:
			return b_*C_S_W;
		  case 5:
			return (kads * C_S_W * CSSmax)/(kads * C_S_W + kdes);
		  case 3:
			return 0.;//only pde
		  case 9:
			return 0.;//only pde
		  default:
			DUNE_THROW(Dune::InvalidStateException, "css1Function not recognised (0 or 1) "+ std::to_string(css1Function));
		}
		return 1.;
	}

    std::function<double(double,double)> computeDtCSS2 =
        std::bind(&Richards1P5CProblem::computeDtCSS2_, this, std::placeholders::_1, std::placeholders::_2); 
        
    double computeInitCSS2_(double CSS1, double CSW) // mol/m3
    {
        if(css1Function == 1){return CSS1;}
        if(css1Function == 5){return CSS1;}
        if(css1Function == 3)
        {
            if(verbose>1)
            {
                std::cout<<"computeInitCSS2_ "<<CSSmax<<" "<<CSW<<" "<<k_sorp<<" "<<( CSSmax * CSW/(CSW+k_sorp))<<std::endl;
            }
            return CSSmax * CSW/(CSW+k_sorp);
        }
        if(css1Function == 9){return (kads * CSW * CSSmax)/(kads * CSW + kdes);}
        if(css1Function == 5){return (kads * CSW * CSSmax)/(kads * CSW + kdes);}
        DUNE_THROW(Dune::InvalidStateException, "css1Function not recognised (0, 1, or 2)"+ std::to_string(css1Function));
        return 0.;
    }
    double computeDtCSS2_(double CSW, double CSS2) // mol/m3/s
    {
        if(css1Function == 9){return kads * CSW * (CSSmax - CSS2) - kdes * CSS2;}
        if(css1Function == 5){return 0.;}
        //if(css1Function == 10){return kads * (CSSmax - CSS2) - kdes * CSS2;}
        
        DUNE_THROW(Dune::InvalidStateException, "css1Function not recognised (0, 1, or 2)"+ std::to_string(css1Function));
        return 0.;
    }
    
	/*!
	 * \copydoc FVProblem::initial
	 *
	 * called by FVProblem::applyInitialSolution(...)
	 */
	template<class Entity>
	PrimaryVariables initial(const Entity& entity) const {
		auto eIdx = this->gridGeometry().elementMapper().index(entity);
		Scalar z = entity.geometry().center()[dimWorld - 1];
		PrimaryVariables v(0.0);
		//Scalar h_init = initialSoil_[h2OIdx].f(z,eIdx)
		//Scalar h_mucil = std::max(h_init, 0.);
		//if()
		//{
			
		//}else{
			v[pressureIdx] = toPa_(initialSoil_[h2OIdx].f(z,eIdx));
		//}
		if(verbose>1)
		{
			 std::cout<<"PrimaryVariables initial(1p5Cproblem) "<<z<<" "<<v[pressureIdx];
		}
		for(int i = soluteIdx; i<numComponents_;i++)//solutes
		{
			v[i] = initialSoil_.at(i).f(z,eIdx);
			if(verbose>1)
			{
				 std::cout<<" result : "<<v[i];
			}
		}
		if(verbose>1)
		{
			 std::cout<<std::endl;
		}
		
		return v;
	}

	/*!
	 * @copydoc FVProblem::boundaryTypesAtPos
	 *
	 * discretization dependent, e.g. called by BoxElementBoundaryTypes::boundaryTypes(...)
	 * when?
	 */
	BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const {
		BoundaryTypes bcTypes; ///pressureIdx, bcTopType_, bcSTopType_, ...
		bcTypes.setAllNeumann();
        if (onUpperBoundary_(globalPos)) { // top, bot, or outer bc
            if (bcTopType_ == constantPressure) {
                bcTypes.setDirichlet(pressureIdx);
			 
														  
												  
            }
			for(int i = soluteIdx; i<numComponents_;i++)//solutes
			{
				if (bcSTopType_.at(i - soluteIdx) == constantConcentration) 
				{
					bcTypes.setDirichlet(i);
				}
			}
            
        } else if (onLowerBoundary_(globalPos)) { // top, bot, or outer bc
            if (bcBotType_ == constantPressure) {
                bcTypes.setDirichlet(pressureIdx);//,conti0EqIdx
            }
			for(int i = soluteIdx; i<numComponents_;i++)//solutes
			{
				if (bcSBotType_.at(i - soluteIdx) == constantConcentration) 
				{
					bcTypes.setDirichlet(i);
				}
			}
            
        }
		return bcTypes;
	}

    /*!
     * \copydoc FVProblem::dirichletAtPos
     *
     * dirchlet(...) is called by the local assembler, e.g. BoxLocalAssembler::evalDirichletBoundaries
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const {
		DUNE_THROW(Dune::InvalidStateException, "Do not use Dirichlet BC please");
        PrimaryVariables values;
        if (onUpperBoundary_(globalPos)) 
		{ // top bc
            switch (bcTopType_) 
			{
				case constantPressure: values[pressureIdx] = toPa_(bcTopValues_[pressureIdx]); break;
				default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Dirichlet: unknown boundary type");
            }
			for(int i = soluteIdx; i<numComponents_;i++)//solutes
			{
				switch (bcSTopType_.at(i - soluteIdx)) 
				{
					case constantConcentration: values[i] = bcSTopValue_.at(i - soluteIdx)*rho_; break;
					default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Dirichlet: unknown boundary type");
				}
            }
        } else if (onLowerBoundary_(globalPos)) 
		{ // bot bc
            switch (bcBotType_) 
			{
				case constantPressure: values[pressureIdx] = toPa_(bcBotValues_[pressureIdx]); break;
				default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Dirichlet: unknown boundary type");
            }
			for(int i = soluteIdx; i<numComponents_;i++)//solutes
			{
				switch (bcSBotType_.at(i - soluteIdx)) 
				{
					case constantConcentration: values[i] = bcSBotValue_.at(i - soluteIdx)*rho_; break;
					default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Dirichlet: unknown boundary type");
				}
            }
        }
        // values.setState(Indices::bothPhases); /// <--?????
        return values;
    }

	
   /*!
	 * \copydoc FVProblem::neumann // [kg/(m²*s)] or [mol/(m²*s)]
	 *
	 * called by BoxLocalResidual::evalFlux,
	 * negative = influx, mass flux in \f$ [ (kg or mol) / (m^2 \cdot s)] \f$// effective diffusion coefficient !!!!!!!!!!!!
	 */
	NumEqVector neumann(const Element& element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const ElementFluxVariablesCache&  fluxCache,
			const SubControlVolumeFace& scvf) const {
		
		NumEqVector flux;
		GlobalPosition pos = scvf.center();
		auto& volVars = elemVolVars[scvf.insideScvIdx()];
		

		/*
		 *  WATER
		 */
		double f = 0.; // return value [kg m-2 s-1)] or [mol m-2 s-1]
		double pos0 = 1;		double scvf_area = scvf.area();
		if(dimWorld == 1){
			pos0 =pos[0]; 
			scvf_area = 2 * M_PI * pos0 * segLength;//m2
		}
		if ( onUpperBoundary_(pos) || onLowerBoundary_(pos) ) {

			//Scalar s = volVars.saturation(h2OIdx);
			// m2 * [kg/m^3] * [m/s^2]/[Pa*s] = m/s
			Scalar kc = volVars.permeability() * volVars.density(h2OIdx) * g_/volVars.viscosity(h2OIdx);//this->spatialParams().hydraulicConductivity(element); //  [m/s], 
			Scalar h = volVars.pressureHead(0);//-toHead_(p); // in cm todo why minus -pc?
			GlobalPosition ePos = element.geometry().center();
			// multiply distance by 2 to better limit the flow (because we can have sharp decrease of wat. pot. 
            //between cell center and cell face) which cannot be taken into account because of the discretisation.
			Scalar dz = 100 *  std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]); //2 *// m->cm
            if (dzScaling == dx){dz = dz * 2;}
            
			Scalar krw = volVars.relativePermeability();//materialLaw_ .krw(s);//	The relative permeability for the wetting phase [between 0 and 1]
			
			//useMole fraction or mass fraction? 
			//[kg/m3] or [mol/m3]
			Scalar rhoW = useMoles ? volVars.molarDensity(h2OIdx) : volVars.density(h2OIdx) ;
			double cm2m = 1./100. ; //[m/cm]
			double unitConversion = cm2m; //something else needed? 
			if (onUpperBoundary_(pos)) { // top bc
				switch (bcTopType_) {
                case constantPressure: {
                    f = rhoW * kc * ((h - bcTopValues_[pressureIdx]) / dz - gravityOn_)*pos0 *unitConversion; // maximal inflow
                    //std::cout << "!";
                    break;
                }
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcTopValues_[pressureIdx]*rhoW/(24.*60.*60.) * unitConversion; // cm/day -> kg/(m²*s) or 
					if (f < 0) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_); // maximal inflow
						f = std::max(f, imax)*pos0;
					} else { // outflow
                    // mol/m2/s = mol/m3 * m/s * - * (cm -cm )/cm
						Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
						f = std::min(f, omax)*pos0;
					}
					break;
				}
				case constantFluxCyl: { // with switch for maximum in- or outflow
					//useMoles:  [cm/d] * [mol/m^3] * [d/s] * [m/cm] = [m/s] *  [mol/m^3] = [mol /(m^2 * s)]
					//useMass:  [cm/d] * [kg/m^3] * [d/s] * [m/cm] = [m/s] *  [kg/m^3] = [kg /(m^2 * s)]
					f = -bcTopValues_[pressureIdx]*rhoW/(24.*60.*60.) * unitConversion * pos[0];
					if (f < 0) { // inflow
					//[mol/m^3] * [m/s] * [cm/cm] = mol/(m2 * s)
					
					//kg/m3 * m2 
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_)*pos[0]; // maximal inflow
						f = std::max(f, imax);
                        if ((f!= 0)&&verbose)
						{
							std::cout<<"onupperBoundary_constantFluxCyl, f: "<<bcTopValues_[pressureIdx]<<" "<<
                            f<<", imax: "<<imax<<", std::max(f, imax): "<<(std::max(f, imax))
							<<", krw: "<<krw<<", kc: "<<kc<<", h: "<<h<<" rho "<<rhoW<<" pos[0] "<<pos[0]<<" dz "<<dz<<std::endl;
						}
					} else { // outflow
                    //mol/m3 * [m/s] *[-] *[cm/cm] = mol/m2/s * pos0 for axysimmetric scaling
						Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
                        
						if ((f!= 0)&&verbose)
						{
							std::cout<<"onUpperBoundary_constantFluxCyl, finput: "<<bcTopValues_[pressureIdx]
                            <<", f: "<<f<<", omax: "<<omax<<", std::min(f, omax): "<<(std::min(f, omax))
							<<", krw: "<<krw<<", kc: "<<kc<<", h: "<<h<<" pos[0] "<<pos[0]<<" dz "<<dz<<std::endl;
						}
						f = std::min(f, omax);
					}
					break;
				}
				case atmospheric: { // atmospheric boundary condition (with surface run-off) // TODO needs testing & improvement
					Scalar prec = -componentInput_[h2OIdx].f(time_);//-precipitation_.f(time_);
					if (prec < 0) { // precipitation
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_); // maximal infiltration
						f = std::max(prec, imax)*pos0;
					} else { // evaporation
						Scalar emax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal evaporation
						f = std::min(prec, emax)*pos0;
					}
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann (water) unknown type: "+std::to_string(bcTopType_));
				}
			} else if (onLowerBoundary_(pos)) { // bot bc
				switch (bcBotType_) {
                case constantPressure: {
                    f = rhoW * kc * ((h - bcBotValues_[pressureIdx]) / dz - gravityOn_)* pos0; // maximal inflow
//                    Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
//                    f = std::min(f, omax); rho_
                    break;
                }
				case constantFlux: { // with switch for maximum in- or outflow
					f = -bcBotValues_[pressureIdx]*rhoW/(24.*60.*60.) * unitConversion; // cm/day -> kg/(m²*s) or mol/(m2*s)
					if (f < 0) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_); // maximal inflow
						f = std::max(f, imax)*pos0;
					} else { // outflow
						Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_); // maximal outflow (evaporation)
						f = std::min(f, omax)*pos0;
					}
					break;
				}
				case constantFluxCyl: { // with switch for maximum in- or outflow
					f = -bcBotValues_[pressureIdx]*rhoW/(24.*60.*60.) * unitConversion * pos[0];
					
					//std::cout<<" water flow "<<f<<" "<<bcBotValues_[pressureIdx]<<" rhoW "<<rhoW<<" useMoles "<<useMoles<<" molarDensity "<<volVars.molarDensity(h2OIdx)
					//<<" density "<<volVars.density(h2OIdx)<<" unit conversion "<<unitConversion<<" "<<pos[0]<<std::endl;
					if (f < 0) { // inflow
						Scalar imax = rhoW * kc * ((h - 0.) / dz - gravityOn_)* pos[0]; // maximal inflow
						if ((f!= 0)&&verbose)
						{
							std::cout<<"onLowerBoundary_constantFluxCyl, f: "<<bcBotValues_[pressureIdx]<<" "<<
                            f<<", imax: "<<imax<<", std::max(f, imax): "<<(std::max(f, imax))
							<<", krw: "<<krw<<", kc: "<<kc<<", h: "<<h<<" rho"<<rhoW<<" pos[0] "<<pos[0]<<" dz "<<dz<<std::endl;
						}
						f = std::max(f, imax);
					} else { // outflow
						Scalar omax = rhoW * krw * kc * ((h - criticalPressure_) / dz - gravityOn_)* pos[0]; // maximal outflow (evaporation)
						// std::cout << " f " << f*1.e9  << ", omax "<< omax << ", value " << bcBotValue_.at(0) << ", crit "  << criticalPressure_ << ", " << pos[0] << "\n";
						if ((f!= 0)&&verbose)
						{
							std::cout<<"onLowerBoundary_constantFluxCyl, f: "<<bcBotValues_[pressureIdx]<<" "<<
                            f<<", omax: "<<omax<<", std::min(f, omax): "<<(std::min(f, omax))
							<<", krw: "<<krw<<", kc: "<<kc<<", h: "<<h<<" rho "<<rhoW<<" pos[0] "<<pos[0]<<" dz "<<dz
							<<" criticalPressure_ "<<criticalPressure_<<" gravityOn_ "<<gravityOn_<<std::endl;
						}
						f = std::min(f, omax);
					}
					break;
				}
				case freeDrainage: { // holds when useMoles?
					f = krw * kc * rhoW *pos0; // * 1 [m]
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann (water) unknown: "+std::to_string(bcBotType_));
				}
			}
			
		}
		
		flux[h2OIdx] = f;// [mol /(m^2 * s)] * pos0

		/*
		 * SOLUTES
		 */
		
			//[kg/m3] or [mol/m3]
		Scalar rhoW = useMoles ? volVars.molarDensity(h2OIdx) : volVars.density(h2OIdx) ;
		double g2kg = 1./1000. ;
		double m2_2_cm2 = 10000;
		double unitConversion = useMoles ? m2_2_cm2 : m2_2_cm2 * g2kg; //something else needed? 
		for(int i_s = 0; i_s < numFluidSolutes; i_s++) //for(int i = soluteIdx;i<numComponents_;i++)
		{
        
			int i = i_s + 1;//int i_s = i - soluteIdx;//for vectors which do not have a value for the H2O primary variable
            
		if(verbose>1)
		{
			std::cout << "neumann_solutes() ";
			std::cout<<i_s <<" "<<i<<" "<<bcSBotType_.size()<<" "<<(i_s < bcSBotType_.size());
            std::cout<<" bc? "<<onUpperBoundary_(pos)<<" "<<onLowerBoundary_(pos)<<" ";
            std::cout<<" type "<<bcSTopType_.at(i_s)<<" "<<bcSTopValue_.at(i_s)<<std::endl;
		}
			Scalar massOrMolFraction = useMoles? volVars.moleFraction(0, i) : volVars.massFraction(0, i);
			if (onUpperBoundary_(pos)) { // top bc Solute
				//std::cout<<"neumann solute, upper BC "<<bcSTopType_.at(i_s)<<" ";
				switch (bcSTopType_.at(i_s)) {
				case constantConcentration: {
					GlobalPosition ePos = element.geometry().center();
					Scalar dz =  std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]);
                    if (dzScaling == dx){dz = dz * 2;}
					//!!! att component param set
					Scalar de = volVars.effectiveDiffusionCoefficient(conti0EqIdx,soluteIdx,0);
					flux[i] = (de * (massOrMolFraction*rhoW-bcSTopValue_.at(i_s)*rhoW) / dz + f * massOrMolFraction)*pos0;
					break;

				}
				case constantFlux: {
					//flux[i] = -bcSTopValue_.at(i_s)*rhoW/(24.*60.*60.)/100*pos0; // cm/day -> kg/(m²*s)
					//usemoles:
					// mol/(cm2 * d)  * [d/s]  * cm2/m2 =  mol/(m2 * s)
					// use mass:
					// g/(cm2/d) *[d/s] * [kg/g] * cm2/m2 = kg/(m2 * s)
					flux[i] = -bcSTopValue_.at(i_s)/(24.*60.*60.)* unitConversion; // g/cm2/day || mol/cm2/day  -> kg/(m²*s) || mol/(m²*s)
					break;
				}
				case constantFluxCyl: {// massOrMolFraction
                
                    double soluteFlow =  -bcSTopValue_.at(i_s)/(24.*60.*60.)*unitConversion; // [mols/(m2 * s)]
                    if (soluteFlow > 0) { // outflow
                        GlobalPosition ePos = element.geometry().center();
                        Scalar dz =  std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]);// m
                        if (dzScaling == dx){dz = dz * 2;}
                        //!!! att component param set
                        Scalar de = volVars.effectiveDiffusionCoefficient(conti0EqIdx,soluteIdx,0);
						Scalar minCvalue = 0.;// mol/mol
                        // [m^2/s] * ([mols / molw]* [molw/m3] + [mols / molw]* [molw/m3]) /m + [molw /(m^2 * s)] * [mols / molw]
                        //  [mols/m2/s]
                        Scalar omax = (de * (massOrMolFraction*rhoW - minCvalue*rhoW) / dz + flux[h2OIdx] * massOrMolFraction);
						if (verbose)
						{
							std::cout<<"onUpperBoundary_constantFluxCyl, Fs: "<<bcSTopValue_.at(i_s)<<" soluteFlow "<<
                            soluteFlow<<", omax: "<<omax<<" min "<<std::min(soluteFlow, omax)<<", massOrMolFraction: "<<massOrMolFraction
							<<" unitConversion "<<unitConversion<<" pos[0] "<<pos[0]<<std::endl;
						}
						soluteFlow = std::min(soluteFlow, omax);
                        
                    }
                    
					flux[i] = soluteFlow*pos[0]; // g/cm2/day || mol/cm2/day -> kg/(m²*s) || mol/(m²*s)
					//std::cout<<"constantfluxCyl "<<flux[i];
					break;
				}
				case outflow: {
					// std::cout << "f " << f << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n";
					flux[i] = f * massOrMolFraction;//*pos0;
					break;
				}
				case outflowCyl: {
					// std::cout << "f " << f << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n";
					flux[i] = f * massOrMolFraction *pos0;
					break;
				}
				case managed: {
					Scalar input = componentInput_.at(i).f(time_);
					flux[i] = input*pos0;
					break;
				}
				case michaelisMenten: {	
					// [mol m-2 s-1] * [mols / molw] * [molw/m3] / ([mol m-3] + [mols / molw] * [molw/m3])
					flux[i] = vMax * std::max(massOrMolFraction,0.)*rhoW/(km + std::max(massOrMolFraction,0.)*rhoW)*pos0;
						if (verbose)
						{
							std::cout<<"onUpperBoundary_michaelisMenten, vMax: "<<vMax<<" massOrMolFraction "<<
                            massOrMolFraction<<", rhoW: "<<rhoW<<" km "<<km<<" pos0 "<<pos0<<std::endl;
						}
					break;
				}
				default:
					DUNE_THROW(Dune::InvalidStateException, "Top boundary type Neumann (solute) unknown: "+std::to_string(bcSTopType_.at(i_s)));
				}
			} else if (onLowerBoundary_(pos)) { // bot bc Solute
				switch (bcSBotType_.at(i_s)) {
				case constantConcentration: {
					GlobalPosition ePos = element.geometry().center();
					Scalar dz = std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]);
                    if (dzScaling == dx){dz = dz * 2;}
					Scalar de = volVars.effectiveDiffusionCoefficient(conti0EqIdx,soluteIdx,0);
					//diffusion + advection (according to wat flux computed above)
					flux[i] =(de * (massOrMolFraction*rhoW-bcSBotValue_.at(i_s)*rhoW) / dz + f * massOrMolFraction)*pos0;
					//[kg_solute/(m²*s)] = [m2 / s] * ([kg_solute/kg_tot] * [kg_tot / m^3_tot] - kg_solute/ m^3_tot)/m + [kg_tot/(m²*s)] * [kg_solute/kg_tot]
					// std::cout << d*1.e9 << ", "<< de*1.e9 << ", " << volVars.massFraction(0, i) << ", " << bcSBotValue_ << ", " << flux[i]*1.e9  << "\n";
					break;
				}
				case constantFlux: {
					//flux[i] = -bcSBotValue_.at(i_s)*rhoW/(24.*60.*60.)/100*pos0; // cm/day -> kg/(m²*s)
					flux[i] = -bcSBotValue_.at(i_s)/(24.*60.*60.)*unitConversion; // g/cm2/day -> kg/(m²*s)
					//[kg_solute/(m²*s)] = [cm_solute/day] * (kg_tot / m^3_tot)* s/day * m/cm ??
					break;
				}
				case constantFluxCyl: {
                        
                    double soluteFlow =  -bcSBotValue_.at(i_s)/(24.*60.*60.)*unitConversion; // [mols/(m2 * s)]
                    if (soluteFlow > 0) { // outflow
                        GlobalPosition ePos = element.geometry().center();
                        Scalar dz =  std::fabs(ePos[dimWorld - 1] - pos[dimWorld - 1]);// m
                        if (dzScaling == dx){dz = dz * 2;}
                        //!!! att component param set
                        Scalar de = volVars.effectiveDiffusionCoefficient(conti0EqIdx,soluteIdx,0);
						Scalar minCvalue = 0.;// mol/mol
                        // [m^2/s] * ([mols / molw]* [molw/m3] + [mols / molw]* [molw/m3]) /m + [molw /(m^2 * s)] * [mols / molw]
                        //  [mols/m2/s]
                        Scalar omax = (de * (massOrMolFraction*rhoW - minCvalue*rhoW) / dz + flux[h2OIdx] * massOrMolFraction);
						if (verbose)
						{
							std::cout<<"onLowerBoundary_constantFluxCyl, solute, Fs: "<<bcSBotValue_.at(i_s)<<" soluteFlow "<<
                            soluteFlow<<", omax: "<<omax<<" min "<<std::min(soluteFlow, omax)<<", massOrMolFraction: "<<massOrMolFraction
							<<" unitConversion "<<unitConversion<<" pos[0] "<<pos[0]<<std::endl;
						}
						soluteFlow = std::min(soluteFlow, omax);
                        
                    }
                    
					flux[i] = soluteFlow*pos[0]; // g/cm2/day || mol/cm2/day -> kg/(m²*s) || mol/(m²*s)
                        
					break;
				}
				case outflow: {//?? free outflow??
					// std::cout << "f " << f*1.e6 << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n";
					flux[i] = f * massOrMolFraction;
					break;
				}
				case outflowCyl: {
					// std::cout << "f " << f << ", "  << volVars.massFraction(0, i) << "=" << f*volVars.massFraction(0, i) << "\n"; rho_
					flux[i] = f * massOrMolFraction*pos0;
					break;
				}
				case michaelisMenten: {	
					// [mol m-2 s-1] * [mols / molw] * [molw/m3] / ([mol m-3] + [mols / molw] * [molw/m3])
					flux[i] = vMax * std::max(massOrMolFraction,0.)*rhoW/(km + std::max(massOrMolFraction,0.)*rhoW)*pos0;
						if (verbose)
						{
							std::cout<<"onLowerBoundary_michaelisMenten, vMax: "<<vMax<<" massOrMolFraction "<<
                            massOrMolFraction<<", rhoW: "<<rhoW<<" km "<<km<<" pos0 "<<pos0<<std::endl;
						}
					break;
				}
				default: DUNE_THROW(Dune::InvalidStateException, "Bottom boundary type Neumann (solute): unknown error");
				}
			} else {
				flux[i] = 0.;
			}
		  
							 
		}

		return flux;
	}

	/*!
	 * \copydoc FVProblem::source
	 *
	 * called by FVLocalResidual:computeSource(...)
	 *
     * E.g. for the mass balance that would be a mass rate in \f$ [ mol / (m^3 \cdot s)] \f
     */
	NumEqVector source(const Element &element, const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars,
			const SubControlVolume &scv) const {
		NumEqVector source;
		GlobalPosition pos = scv.center();
		bool dobioChemicalReaction = dobioChemicalReaction_; //by default, do biochemical reactions
		double pos0; 
        double svc_volume;
        int dofIndex = scv.dofIndex();
		if (dimWorld == 1)//1daxissymmetric model
		{
			pos0 = pos[0];
            svc_volume = getCellVolumesCyl(dofIndex);//with 1d model, need to evaluate manually the volume of the cell.
                            // for simplicity, we give directly source as [ mol / (m^3 \cdot s)] for now
		}else{ // dimWorld == 3
			pos0 = 1.;
            svc_volume = scv.volume();
        }
  
															  
		auto eIdx = this->spatialParams().gridGeometry().elementMapper().index(element);
        const auto& volVars = elemVolVars[scv];
		for(int i = 0;i < numComponents_;i++)
		{			
			if (source_[i] != nullptr) {
				source[i] = source_[i]->at(eIdx)/svc_volume * pos0;// [ mol / (m^3 \cdot s)]
                if(i == h2OIdx)
                {
                
                    if (sourceSlope>=0.) {
                        Scalar h = volVars.pressureHead(0);// cm
                        if (h<criticalPressure_) {
                            source[i] = 0.;
                        } else if (h<=criticalPressure_+sourceSlope) { //  h in [crit, crit+slope]
                            double ratio = (h - criticalPressure_)/sourceSlope;
                            // std::cout << "source(): " << h << ", "<< theta << "\n" << std::flush;
                            source[i] = ratio * source_[i]->at(eIdx)/svc_volume * pos0;
                        } 
                    }
                }
				if(dobioChemicalReaction&&(i> h2OIdx)&&(source[i]> 0)&&reactionExclusive)
				{//if we have a reaction in the cell via source and source-biochemreaction mutually exclusive, 
					//biochem is disabled
					dobioChemicalReaction = false;
				}
			}else{source[i] = 0.;}												 
		}		
		if(((dimWorld == 1)&&(reactionExclusive))||((dimWorld > 1)&&(!reactionExclusive)))
		{
			std::cout<<"(((dimWorld == 1)&&(reactionExclusive))||((dimWorld > 1)&&(!reactionExclusive))) "
			<<dimWorld<<" "<<reactionExclusive<<std::endl;
			DUNE_THROW(Dune::InvalidStateException, "(((dimWorld == 1)&&(reactionExclusive))||((dimWorld > 1)&&(!reactionExclusive))) ");
			
		}
		if(dobioChemicalReaction)
		{
			bioChemicalReaction(element, source, volVars, pos0, scv);
		}
	
		if(verbose>1)
		{
			std::cout << "source() ";
			for(int i = 0;i < numComponents_;i++)
			{
				std::cout<<source[i]<<" ";
			}			
			std::cout << std::endl;
		}
		
		return source;
	}
	
	
		/*!
	 *
     * E.g. for the mol balance that would be a mass rate in \f$ [ mol / (m^3 \cdot s)] \f
     */
	void bioChemicalReaction(const Element &element, 
								NumEqVector &q, const VolumeVariables &volVars, double pos0, 
								const SubControlVolume &scv) const
	{
		
		// (mol Soil / m3 soil)  
		double solidDensity = massOrMoleDensity(volVars, SolidSystem::mainCompIdx , false);//soilIdx -  numFluidComps
		// m3 soil/m3 scv
		double solVolFr = (1. - volVars.porosity());
		// (mol soil / m3 scv) = (mol Soil / m3 soil)  * ([m3 space - m3 pores]/m3 scv)
		double bulkSoilDensity = solidDensity * solVolFr;
		bool doDecay = doDecay_; //by default, do decay
        double svc_volume;
        int dofIndex = scv.dofIndex();
		
		if (dimWorld == 1)//1daxissymmetric model
		{
            svc_volume = getCellVolumesCyl(dofIndex);//with 1d model, need to evaluate manually the volume of the cell.
                            // for simplicity, we give directly source as [ mol / (m^3 \cdot s)] for now
		}else{ // dimWorld == 3
            svc_volume = scv.volume();
        }
        
		std::vector<double> WorCorN(numComponents_);
		
        WorCorN[h2OIdx] = volVars.saturation(h2OIdx) * volVars.porosity(); //m3 water / m3 scv		
		
		for(int Idx = 1; Idx <= numFluidSolutes; Idx++)
		{
			WorCorN[Idx] = massOrMoleDensity(volVars, h2OIdx, true) * 
					std::max(massOrMoleFraction(volVars,0, Idx, true),					//mol C or N/mol soil water
					0.);																//mol C or N/m3 soil water
		}
		
		for(int Idx = numFluidSolutes +1 ; Idx <= numSolutes; Idx++)
		{
			//[mol C / m3 scv] = [mol scv/m3 scv] * [mol C/mol scv]
			WorCorN[Idx]  = bulkSoilDensity * 
						std::max(massOrMoleFraction(volVars,0, Idx - numFluidComps, false) //mol C/mol solid soil
						, 0.); //mol C/m3 scv						
		}
		
        double CSS2 = bulkSoilDensity * std::max(massOrMoleFraction(volVars,0, CSS2Idx - numFluidComps, false), 0.) ; // mol C / m3 scv zone 2
		double C_SfrW = std::max(massOrMoleFraction(volVars,0, soluteIdx, true), 0.);					//mol C/mol soil water
		double C_S_W = massOrMoleDensity(volVars, h2OIdx, true) * C_SfrW;								//mol C/m3 soil water
        double dtCSS2 =  computeDtCSS2( C_S_W, CSS2);
		//[mol solute / m3 scv/s] 
		double decay = 0.;
		if (doDecay)
		{
			decay = vmax_decay * WorCorN[soluteIdx]/(km_decay + WorCorN[soluteIdx]);
			//std::cout<<WorCorN[soluteIdx]<<"WorCorN";
			
		}//else{
			//std::cout<<"noDecay";
		//	q[soluteIdx] += -dtCSS2* pos0;
		//}
		q[soluteIdx] += - ( decay + dtCSS2) * pos0; //0* pos0 ; 
		q[decayIdx]  += decay * pos0;
		q[CSS2Idx] +=  dtCSS2 * pos0 ;
			
	}

	/*!
	 * \brief Applies a vector of point sources. The point sources
	 *        are possibly solution dependent.
	 *
	 * \param pointSources A vector of Dumux::PointSource s that contain
              source values for all phases and space positions.
	 *
	 * For this method, the \a values method of the point source
	 * has to return the absolute mass rate in kg/s. Positive values mean
	 * that mass is created, negative ones mean that it vanishes.
	 */
	template<class PointSource>
	void addPointSources(std::vector<PointSource>& pointSources) const {}

	/*!
	 * \brief Evaluate the point sources (added by addPointSources) TODO
	 *        for all phases within a given sub-control-volume.
	 *
	 * This is the method for the case where the point source is
	 * solution dependent and requires some quantities that
	 * are specific to the fully-implicit method.
	 *
	 * For this method, the \a values() method of the point sources returns
	 * the absolute rate mass generated or annihilate in kg/s. Positive values mean
	 * that mass is created, negative ones mean that it vanishes.
	 */
	 //not used??
	template<class ElementVolumeVariables>
	void pointSource(PointSource& source,
			const Element &element,
			const FVElementGeometry& fvGeometry,
			const ElementVolumeVariables& elemVolVars,
			const SubControlVolume &scv) const {

		PrimaryVariables sourceValue(0.);
		std::cout<<"template<class ElementVolumeVariables>"<<std::endl;
		DUNE_THROW(Dune::InvalidStateException, "template<class ElementVolumeVariables>");

	}

	/*!
	 * Sets the current simulation time (within the simulation loop) for atmospheric look up [s]
	 *
	 * eventually, called in the main file (example specific, richards.cc)
	 */
	void setTime(Scalar t, Scalar dt) {
		time_ = t;
		dt_ = dt; // currently unused
	}

	/*!
	 * Source per element index \f$ [ kg / s)] \f$
	 *
	 * eventually, called in the main file (example specific, richards.cc)
	 */
	void setSource(std::shared_ptr<std::vector<double>> s, int eqIdx = 0) {
		source_[eqIdx] = s;
	}

	//! Set the coupling manager
	void setCouplingManager(CouplingManager* cm) {
		couplingManager_ = cm;
	}

	/*!
	 * sets the critical pressure for evaporation [cm] (default = -10000 cm)
	 *
	 *  eventually, called in the main file (example specific, richards.cc)
	 */
	void criticalPressure(Scalar p) {
		criticalPressure_ = p;
	}

	/**
	 * Sets boundary fluxes according to the last solution TODO
	 */
	void postTimeStep(const SolutionVector& sol, const GridVariables& gridVars) {
		double bc_flux_upper = 0.;
		double bc_flux_lower = 0.;
		int uc = 0;
		int lc = 0;
		for (const auto& e :elements(this->gridGeometry().gridView())) {
			auto fvGeometry = localView(this->gridGeometry());
			fvGeometry.bindElement(e);
			auto elemVolVars = localView(gridVars.curGridVolVars());
			elemVolVars.bindElement(e, fvGeometry, sol);
			for (const auto& scvf :scvfs(fvGeometry)) { // evaluate root collar sub control faces
				auto p = scvf.center();
				if (onUpperBoundary_(p)) { // top
					bc_flux_upper += neumann(e, fvGeometry, elemVolVars, scvf);
					uc++;
				} else if (onLowerBoundary_(p)) { // bottom
					bc_flux_lower += neumann(e, fvGeometry, elemVolVars, scvf);
					lc++;
				}
			}
		}
		bc_flux_upper /= uc;
		bc_flux_lower /= lc;
	}

	/*!
	 * Writes the actual boundary fluxes (top and bottom) into a text file. Call postTimeStep before using it.
	 */
	void writeBoundaryFluxes() {}

	/**
	 * debug info TODO make meaningful for 2c
	 */
	void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars) const {
		NumEqVector source(0.0);
		for (const auto& element : elements(this->gridGeometry().gridView())) {
			auto fvGeometry = localView(this->gridGeometry());
			fvGeometry.bindElement(element);
			auto elemVolVars = localView(gridVars.curGridVolVars());
			elemVolVars.bindElement(element, fvGeometry, sol);
			for (auto&& scv : scvs(fvGeometry)) {
				auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
				pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
				source += pointSources;
			}
		}
		std::cout << "Global integrated source (soil): " << source[h2OIdx] << " (kg/s) / " << source[h2OIdx]*3600*24*1000 << " (g/day)" << '\n';
	}

	/**
	 * Forwards to spatialParams
	 *
     * Call to change default setting (of 1.e-6 for both)
     *
     * pcEps    capillary pressure regularisation
     * krEps 	relative permeabiltiy regularisation
     */
    void setRegularisation(double pcEps, double krEps) {
    	this->spatialParams().setRegularisation(pcEps,krEps);
    }

    /**
     * Forwards to spatialParams
     */
    void addVanGenuchtenDomain(double minx, double miny, double minz, double maxx, double maxy, double maxz, int layerIndex) {
        this->spatialParams().addVanGenuchtenDomain(minx, miny, minz, maxx, maxy, maxz, layerIndex);
    }

    /**
     * Forwards to spatialParams
     */
    void changeVanGenuchtenSet(int vgIndex, double qr, double qs, double alpha, double n, double ks){
        this->spatialParams().changeVanGenuchtenSet(vgIndex, qr, qs, alpha, n, ks);
    }
	
	void setFaceGlobalIndexSet(std::map<int,int>  faceIdx_)
	{
		faceIdx = faceIdx_;											   
	}

	// BC, direct access for Python binding (setTopBC, setBotBC, in richards.hh)
	int bcTopType_;
	int bcBotType_;
	std::vector<double> bcTopValues_ = std::vector<double>(1);
	std::vector<double> bcBotValues_ = std::vector<double>(1);

	std::vector<int> bcSTopType_ = std::vector<int>(numSolutes); // one solute
	std::vector<int> bcSBotType_ = std::vector<int>(numSolutes);
	std::vector<double> bcSTopValue_ = std::vector<double>(numSolutes);
	std::vector<double> bcSBotValue_= std::vector<double>(numSolutes);
	
	bool doSoluteFlow = true;
    
	int verbose;
    int dzScaling;
    
    bool RFmethod2 = false;
	bool verbose_local_residual = false;
	double segLength;//m
	std::vector<double> cellVolumesCyl;
	int nCells_all;
	int nFaces;
	bool computedCellVolumesCyl = false;
	std::map<int,int>  faceIdx;
    
    
private:



	//! cm pressure head -> Pascal
	Scalar toPa_(Scalar ph) const {
		return pRef_ + ph / 100. * rho_ * g_;
	}

	//! Pascal -> cm pressure head
	Scalar toHead_(Scalar p) const {
		return (p - pRef_) * 100. / rho_ / g_;
	}

	//! true if on the point lies on the left boundary
	bool onLeftBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[0] > this->gridGeometry().bBoxMin()[0] - eps_;
	}
	//! true if on the point lies on the right boundary
	bool onRightBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
	}
	
	//! true if on the point lies on the upper boundary
	bool onUpperBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] > this->gridGeometry().bBoxMax()[dimWorld - 1] - eps_;
	}

	//! true if on the point lies on the upper boundary
	bool onLowerBoundary_(const GlobalPosition &globalPos) const {
		return globalPos[dimWorld - 1] < this->gridGeometry().bBoxMin()[dimWorld - 1] + eps_;
	}

	// Initial
	std::vector<InputFileFunction> initialSoil_ = std::vector<InputFileFunction>(numComponents_);
								 

	bool gravityOn_;
	bool dobioChemicalReaction_;
	bool doDecay_;
	// Source
	std::vector<std::shared_ptr<std::vector<double>>> source_; // [kg/s]
	CouplingManager* couplingManager_ = nullptr;

	//InputFileFunction precipitation_;
	std::vector<InputFileFunction> componentInput_ = std::vector<InputFileFunction>(numComponents_);
	Scalar criticalPressure_; // cm
	Scalar time_ = 0.;
	Scalar dt_ = 0.;

	std::ofstream myfile_;
	bool toFile;
	bool reactionExclusive;
	bool doNeffect = true;
	
	
	static constexpr Scalar eps_ = 1.e-7;
	double temperatureK;
	static constexpr Scalar g_ = 9.81; // m / s^2 (for type conversions)
	static constexpr Scalar rho_ = 1.e3; // kg / m^3 (for type conversions)
	static constexpr Scalar pRef_ = 1.e5; // Pa

	
	// default value in CPB units
	double  v_maxL = 5e5; //Maximum reaction rate of enzymes targeting large polymers [d-1]
	double  K_L = 10e-3 ; //Half-saturation coefficients of enzymes targeting large polymers [mol cm-3 soil]
	std::vector<double>  m_max{0.01,0.001}; //[d-1] Maximum maintenance rate coefficient for the corresponding microbial groups
	std::vector<double>  micro_max{1,10}; //[d-1] Maximum growth rate coefficient for the corresponding microbial group
	std::vector<double>  beta{0.1,0.3}; //[-] Reduction factor of maintenance requirements in dormant state for the corresponding microbial group
	std::vector<double>  k_S{ 50000, 30000}; // [mol soil / mol C soil / d] Specific substrate affinity to small polymers for the corresponding microbial group
	std::vector<double>  C_S_W_thres{ 0.001/(12*1000), 0.01/(12*1000)}; //[mol C/cm3 soil water] Threshold concentration for reactivation and deactivation for the corresponding microbial groups
	double  k_phi= 0.1; //[-] Sharpness parameter for the switch function from active to dormancy
	std::vector<double>  k_D{ 1, 5}; //[d-1] Deactivation rate coefficient for the corresponding microbial group
	std::vector<double>  k_R{ 1, 5}; //[d-1] Reactivation rate coefficient for the corresponding microbial group
	double  k_decay= 0.75; //[-] Maintenance yield
	double  k_decay2= 0.5;//[-] Proportion of large polymers formed from dead microbial biomass due to maintenance
	std::vector<double>  k_growth{ 0.8, 0.6};//[-] Growth yield on small polymers for the corresponding microbial groups
	double k_sorp = 0.2;//mol/cm3 or mol
	double CSSmax = 1.75;// [mol/cm3 soil scv zone 1] or mol, max sorption capacity
	double alpha = 0.1; //[d-1]
    double kads = 23265;//[cm3/mol]
    double kdes=4;//[-]
	double vmax_decay=7.32e-5;//[mol C / m^3 scv / s]
	double km_decay=10.5;//[mol C / m^3 scv]
	Scalar b_; // buffer power
	double freundlichK_;// [mol^(1-n)*cm^{3*(n-1)}]
	double freundlichN_;// [-]
	Scalar sourceSlope = -1.; // slope for regularization, negative values disables regularisation
	
	// active root solute uptake
	Scalar vMax; // Michaelis Menten Parameter [mol m-2 s-1]
	Scalar km;  // Michaelis Menten Parameter  [mol m-3]
	
	std::vector<double>  C_aLim ; //[mol C/m3] Minimum concentraiton of activated organisme (to be able to restart from a solute concentration of 0)
	
    double psikPa_opt = -3.;//[kPa]
    double psikPa_th = -15800.;//[kPa]
    double alpha_A = 1.47;//[s-1]
    double psiMPa_A2D = -0.46;//[MPa] psiMPa_A2D
    double tau_DA = 0.39;//[-]
    double w_DA = 3.38;//[-]
    int css1Function;
	double f_sorp = 0;
	

};

    
} //end namespace Dumux

#endif
