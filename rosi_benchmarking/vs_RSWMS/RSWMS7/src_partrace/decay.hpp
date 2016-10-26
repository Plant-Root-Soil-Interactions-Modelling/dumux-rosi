#ifndef DECAYCLASS
#define DECAYCLASS

#include <string>

class ParticleClass;
class PartraceClass;
class ConstData;

/** The base class for calculating mikrobial decay.
 * This base class introduces the interface and does not calculate
 * any decay.
 */
class DecayClass {

protected:
  ParticleClass* pclass;   ///< pointer to particle data.
  PartraceClass* partrace; ///< pointer to global partrace data.
  std::string decaytype;   ///< decay name.
  ConstData *decayfield;

public:
  DecayClass(PartraceClass* pt, ParticleClass* part); ///< constructor.
  virtual ~DecayClass() {  } ///< destructor.

  virtual void CalculateMicrobialDecay() {}
  virtual double DecayProbability(double decay, double conc, double dt) { return 0; }
  std::string DecayType() { return decaytype; }
  inline ConstData *data() { return decayfield; }

};

/** Class for calculating first order mikrobial decay.
 * Model: p = (decay*dt+conc)/conc
 *
 * \li p = probability to decay
 * \li decay = decay rate
 * \li dt = time step
 * \li conc = concentration
 *
 */
class DecayClass1: public DecayClass {
public:
  DecayClass1(PartraceClass* pt, ParticleClass* part, double rate,
	      std::string& DecayField, int d);
  virtual ~DecayClass1();
  virtual double DecayProbability(double decay, double conc, double dt, double wc);
  void CalculateMicrobialDecay();
};


/** Class for calculating second order mikrobial decay.
 * Model: p = exp(decay*dt)
 *
 * \li p = probability to decay
 * \li decay = decay rate
 * \li dt = time step
 */
class DecayClass2: public DecayClass1 {
public:
  DecayClass2(PartraceClass* elem, ParticleClass* part, double rate,
	      std::string& DecayField, int d);
  virtual ~DecayClass2() { }
  virtual double DecayProbability(double decay, double conc, double dt, double wc);
};


// /** Class for calculating Michaelis-Menten decay.
//  * Model: p = (V_max/(K_m + conc)+f) * decay
//  */
// class DecayClass3: public DecayClass1 {
// public:
//   DecayClass3(PartraceClass* elem, ParticleClass* part, double rate,
// 	      std::string& DecayField, int d);
//   virtual ~DecayClass3() { }
//   virtual double DecayProbability(double decay, double conc, double dt, double wc);
// };



// /** Class for calculating decay equal to "sink terms".
//  * Model: p =  decay
//  */
// class DecayClass4: public DecayClass1 {
// public:
//   DecayClass4(PartraceClass* elem, ParticleClass* part, double rate,
// 	      std::string& DecayField, int d);
//   virtual ~DecayClass4() { }
//   virtual double DecayProbability(double decay, double conc, double dt, double wc);
// };

#endif
