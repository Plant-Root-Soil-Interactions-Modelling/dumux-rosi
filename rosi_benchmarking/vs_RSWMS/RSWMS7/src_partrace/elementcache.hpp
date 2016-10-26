#ifndef DISPERSIONCLASS
#define DISPERSIONCLASS
#include "definitions.hpp"
#include "hash.hpp"

class ParticleClass;

/** Class for storing dispersion data.
 * This class is used for storing dispersion parameters
 * between calls to the subroutines ParticleClass::Dispersion
 * and ElementsClass::MoveParticle.
 */
class ElementCacheClass: public HashClassData {
public:
  long e, ex,ey,ez;
  vector3 velo;        ///< velocity vector.
  vector3 ds_conv;     ///< vector containing the convection step.
  vector3 ds_disp;     ///< vector containing the dispersion step.
  vector3 rand;        ///< vector containing random numbers.
  vector3 reflect;     ///< value used for the refection of the dispersion.
  vector3 x1;          ///< start coordinates of the element.
  double retar;        ///< retardation.
  double wc;           ///< water content for the element.
  double diffusion;    ///< value for the diffusion.
  double* parameter;   ///< pointer to the dispersion parameters.
  double interpol[24]; ///< interpolation coefficients.
  long nodes[8];       ///< nodes of the element.
  int finitevolindex;  ///< index for finite volume defined velocity.

  ParticleClass* particle;
  ElementCacheClass(long ele);
  virtual ~ElementCacheClass();
private:
  ElementCacheClass() {}
};

#endif
