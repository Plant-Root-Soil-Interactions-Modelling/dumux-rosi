// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:

#ifndef RICHARDS_PARAMETERS10C_HH
#define RICHARDS_PARAMETERS10C_HH

#include "../soil_richards/richardsparams.hh"


namespace Dumux {

/*!
 * The SpatialParams class of RichardsProblem
 *
 * supports multiple soil layers (in z-direction),
 * with different VG parameters sets
 */
template<class GridGeometry, class Scalar>
class RichardsParams10C : public RichardsParams<GridGeometry, Scalar, RichardsParams<GridGeometry, Scalar>>
{
public:

	RichardsParams10C(std::shared_ptr<const GridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
		
		mucilCAffectsW = getParam<bool>("Soil.mucilCAffectsW", mucilCAffectsW);
		mucilCMolarMass = getParam<Scalar>("Soil.mucilCMolarMass", mucilCMolarMass); // only used if interactions hydraulic conductivity <=> [mucilage] implemented
	}
	
	
    /*
     * decreases the soil hydraulic conductivity caused by mucilage
	 * from Landl (2021, doi:10.3389/fagro.2021.622367)
     */
    template<class ElementSolution>
	Scalar relativePermeabilityMucilage(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol, int idx_mucil) const {
		Scalar d = 1.4
		Scalar Nu = 566
		Scalar cw = gravimetricMucilageConcentration(scv, elemSol, idx_mucil); // g/g
		return 1/(1+ Nu * cw^d);
	}
	
    /*
     * gravimetric mucilage concentration in water [g/g]
     */
    template<class ElementSolution>
	Scalar gravimetricMucilageConcentration(
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol, int idx_mucil) const {
        const auto& priVars = elemSol[scv.localDofIndex()];
		Scalar Cmolar = priVars[idx_mucil]); // molar concentration [mol mucil-C/mol water]
		Scalar Cgravimetric = Cmolar / (Water::molarMass() * 1000) * mucilCMolarMass; // gravimetric concentration [mol mucil-C/mol water] * [mol water/g water] * [g mucil / mol mucil-C]
		return Cgravimetric;
	}
	
    /*
     * gravimetric mucilage concentration in water [g/g]
     */
    template<class ElementSolution>
	Scalar gravimetricMucilageConcentration(
                    const SubControlVolumeFace& scvf,
                    const ElementVolumeVariables& elemVolVars, int idx_mucil) const {
		auto& volVars = elemVolVars[scvf.insideScvIdx()];
        Scalar Cmolar = volVars.moleFraction(0, idx_mucil); // molar concentration [mol mucil-C/mol water]
		Scalar Cgravimetric = Cmolar / (Water::molarMass() * 1000) * mucilCMolarMass; // gravimetric concentration [mol mucil-C/mol water] * [mol water/g water] * [g mucil / mol mucil-C]
		return Cgravimetric;
	}
	
private:

	Scalar mucilCMolarMass = 30; // g mucil/ mol mucilage-C
	bool mucilCAffectsW = false; // does the mucilage-C concentration affect the water flow?
};

} // end namespace Dumux

#endif
