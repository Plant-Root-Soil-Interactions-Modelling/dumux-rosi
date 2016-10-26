#ifndef DISPERSIONCLASS
#define DISPERSIONCLASS
#include "definitions.hpp"

class ParticleClass;
/** Class for storing dispersion data.
 * This class is used for storing dispersion parameters
 * between calls to the subroutines ParticleClass::Dispersion
 * and ElementsClass::MoveParticle.
 */
class DispersionClass {
public:
  vector3 v;     ///< velocity vector.
  vector3 v1;    ///< vector perpendicular to v.
  vector3 v2;    ///< vector perpendicular to v and v1.
  vector3 ds;    ///< vector containing the dispersion step.
  vector3 rand;  ///< vector containing random numbers.
  double vNorm; ///< norm of the velocity.
  double wc;    ///< water content for the element.
	double dmc;
  double* parameter; ///< pointer to the dispersion parameters.
	
	//<<MB
	vector3 v_from;
	vector3 v_to;
	vector3 v1_from;
	vector3 v1_to;
	vector3 v2_from;
	vector3 v2_to;  
	vector3 ds_from;
	vector3 ds_to;
	vector3 ds_to_correlated;
	vector3 ds_from_correlated;
	vector3 ds_to_uncorrelated;
	vector3 dref_from;		///< vector containing the expression needed for the reflection barrier in x,y,z 
	vector3 dref_to;
	vector3 d_from;
	vector3 d_to;
	//vector3 rand_to;
	//vector3 rand_from;
	//vector3 rand_XYZ;
	double dmc_from;   	///< double containing the diffusion coefficient in x,y,z 
	double dmc_to;
	vector3 convection_from;
	vector3 convection_to;
  double vNorm_from;
	double vNorm_to;
	double wc_from;
	double wc_to;
	//>>MB
  ParticleClass* particle; ///< Pointer to the particle class.
  DispersionClass() {}
  DispersionClass(ParticleClass* p): particle(p) {}
  virtual ~DispersionClass() {}
};

#endif
