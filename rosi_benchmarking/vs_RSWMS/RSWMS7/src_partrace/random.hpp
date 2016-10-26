#ifndef RANDOMCLASS
#define RANDOMCLASS
#include <cmath>

/**
  * Class for generating random numbers.
  * Methods for equally spaced and gauss distributed random
  * numbers are included.
  */
class RandomClass {

protected:
  const double m;
  const double a;
  double seed;
  int fast_iset;
  double fast_v1;
  double fast_v2;
  double fast_rsq;
  double fast_fac;
  // 
  int N, M;
  unsigned long ptgfsr[624];
  unsigned long MATRIX_A;
  unsigned long UPPER_MASK;
  unsigned long LOWER_MASK;
  unsigned long MAX32_MASK;
  unsigned long TEMPERING_MASK_B;
  unsigned long TEMPERING_MASK_C;
  unsigned long mag01[2];
  unsigned long y;
  int k, kk;
  double max32;
  double infinity;
  void sgenrand(unsigned long seed); /* seed should not be 0 */

public:
  RandomClass();
  RandomClass(double s);
  ~RandomClass() {  }

  inline void set_seed(double s);

  inline double draw() { 
    // Der multiplikative lineare Kongruenzgenerator (MLKG)
    // Erzeugung einer gleichverteilten Zufallsvariablen 0..m-1
    // m Primzahl, a primitive Wurzel von m ==> Periode m
    // Verfahren: s2=(a*s1) mod m
    seed=fmod(a*seed,m);
    return seed; 
  }

  inline double draw(double lower, double upper) { 
    // s2=(a*s1) mod m
    return lower+draw()/(m-1.0)*(upper-lower);
  }

  inline double gauss(double mean=0, double std=1) { 
    // Gauss-Verteilung
    // Erzeuge 2 gleichverteilte 0..1 Zufallszahlen z1,z2,
    // Erzeuge Radius r=sqrt(-2*ln(z1))
    // Erzeuge Winkel w=2*pi*z2
    // r,w sind Polarkoordinaten eines zufaelligen Punktes
    // karth. Koordinaten x=r*sin(w) und y=r*cos(w) sind 0/1 normalverteilt
    return mean + std * sqrt(-2e0*log(draw(0e0,1e0))) * sin(2e0*M_PI*draw(0e0,1e0));
  }

  inline unsigned long fastdraw_max() const { return MAX32_MASK; }
  double fastdraw();
  double fastgauss();
};

#endif



