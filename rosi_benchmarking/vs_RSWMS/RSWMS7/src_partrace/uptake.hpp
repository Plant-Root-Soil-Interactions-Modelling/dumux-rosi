#ifndef UPTAKECLASS
#define UPTAKECLASS

#include <string>

class ParticleClass;
class PartraceClass;
class ConstData;

/** The base class for calculating solute uptake.
 * This base class introduces the interface and does not calculate
 * any uptake.
 */
class UptakeClass {

protected:
  ParticleClass* pclass;   ///< pointer to particle data.
  PartraceClass* partrace; ///< pointer to global partrace data.
  std::string uptaketype;   ///< uptake name.
  ConstData *uptakefield;

public:
  UptakeClass(PartraceClass* pt, ParticleClass* part); ///< constructor.
  virtual ~UptakeClass() {  } ///< destructor.

  virtual void CalculateUptake() {}
  virtual double UptakeProbability(double uptake, double conc, double dt) { return 0; }
  std::string UptakeType() { return uptaketype; }
  inline ConstData *data() { return uptakefield; }

};


/** Class for calculating passiv uptake.
 * Model: p = uptake
 *
 * \li p = probability to uptake
 * \li uptake = uptake rate
 * \li dt = time step
 * \li conc = concentration
 *
 */
class UptakeClass1: public UptakeClass {
public:
  UptakeClass1(PartraceClass* pt, ParticleClass* part, double rate,
	       std::string& UptakeField,int d);
  virtual ~UptakeClass1();
  virtual double UptakeProbability(double uptake, double conc, double dt, double wc);
  void CalculateUptake();
};


/** Class for calculating Michaelis-Menten uptake (active)
 * Model: p = (V_max/(K_m + conc)+f) * uptake
 */
class UptakeClass2: public UptakeClass1 {
public:
  UptakeClass2(PartraceClass* elem, ParticleClass* part, double rate,
	       std::string& UptakeField, int d);
  virtual ~UptakeClass2() { }
  virtual double UptakeProbability(double uptake, double conc, double dt, double wc);
};


/** Class for calculating concentration dependent uptake
 * Model: p = uptake
 */
class UptakeClass3: public UptakeClass1 {
public:
  UptakeClass3(PartraceClass* elem, ParticleClass* part, double rate,
	       std::string& UptakeField,int d);
  virtual ~UptakeClass3() { }
  virtual double UptakeProbability(double uptake, double conc, double dt, double wc);
};


#endif
