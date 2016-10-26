#ifndef PARTICLECLASS_H
#define PARTICLECLASS_H
#include "particle.hpp"
#include "definitions.hpp"
#include "elements.hpp"
#include "dispersion.hpp"
#include <string>
#include <fstream>

class Event;
class InteractClass;
class ParticleListClass;
class RunTimeClass;
class ConstData;
class DecayClass;
class UptakeClass;


/** Class for the injection of particles.
 * The inserting of particles into the geometry is described by injections.
 */
class InjectionClass {
public:
  // Start distribution of the particles 
  int StartModel;
  vector3 StartPos;
  vector3 StartRadius;
  // injection time information
  double Start;
  double End;
  Event *startEvent; ///< event for the begin of the injection.
  Event *endEvent;   ///< event for the end of the injection.
  // injection mass
  ulong ParticlesToInject; ///< number of particles to inject.
  double InitConcentration;///< initial concentration for the whole volume.
  int ProcLast;            ///< for distributing particles on PEs.
  InjectionClass();
  ~InjectionClass();
  /// returns true if all particles are injected.
  inline bool ready() { return ParticlesToInject<=0; }
  void init(RunTimeClass& RunTime);
};


/** The base particle type.
 * This class contains all properties of a konservative particle.
 * All none konservative particles have to inherit from this class.
 */
class ParticleClass {

  friend class InteractClass;
  friend class ParticleListClass;

protected:
  std::string Name, Abbrev; // Name of particle
  ulong MaxNo;              // max number of particles
  PartraceClass* partrace;  // pointer to global Partrace data
  ParticleClass* next;      // pointer to next ParticleClass
  long NoElxyz;
  // injection information
  int injectionNo;                           // number of injections
  InjectionClass *injection;                 // array of injections
  // Concentration
  double MassPerPart; 
  double C0;      // norm value for concentration
  ConstData *concentration; // contains: conc in water for TP, whole conc for Retardation
  ParticleList* ParticlesInElement; // array of lists of particles
  ParticleList* ParticlesInElementSorbed; // array of lists of sorbed particles
  ParticleSimpleList* MovedParticles;  // temporary array for moved particles

  // Output filenames
  std::string fnBTC, fnBTCSeq, fnConc, fnMom, fnParticles;
  std::ofstream BTCout;
  std::ofstream *BTCSeqout;
  // decay information
  bool decay;
  // uptake information
  bool uptake;
 
  // processor information
  int ProcMax; // no. of processors
  int ProcNo;  // processor number 0..ProcMax-1
  // adsorption parameters
  ConstData *adsorption;

public:
  bool AlwaysCalculateConcentration;
  Event *SaveConcTime;
  Event *SaveMomentsTime;
  Event *SaveParticlesTime;
  // Input filenames
  std::string RestartFromConcFile; ///< Initial concentration from conc file.
  DecayClass *Decay;
  UptakeClass *Uptake;

  // Constructor
  ParticleClass(PartraceClass* pt, long no, double m, 
		std::string &name, std::string &abbrev,
		int injno, InjectionClass *inj,	double c0);
  // Destructor
  virtual ~ParticleClass();

  Particle* new_particle();
  Particle* new_particle(vector3 &pos);
  bool InjectParticlePart();
  void Distribute(int& pno, int& ProcLast);
  void MoveParticles();
  void StartDirac(ulong pno, InjectionClass &inj);  
  void StartCube(ulong pno, InjectionClass &inj);   
  void StartEllipse(ulong pno, InjectionClass &inj);
  void InitDecay(std::string& DecayField, int order, int useField, double rate);
  void InitUptake(std::string& UptakeField, int order, int useField, double rate);
  void SaveConcentration(double simtime);
  void SaveParticlePositions(double simtime);
  void OpenBTC(double simtime);
  void WriteBTC(double simtime);
  void SaveParticles(std::fstream& res);
  void SaveParticlesList(std::fstream& res, ParticleList* plist);
  void LoadParticles(std::fstream& res);
  void LoadParticlesList(std::fstream& res, ParticleList* plist);
  void setInitialConcentration(int &pno, InjectionClass &inj);
  inline double massPerPart() const { return MassPerPart; }
  inline double Czero() const { return C0; }
  inline bool has_decay() const { return decay; }
  inline bool has_uptake() const { return uptake; }
  inline ParticleList* firstParticleList() const { return ParticlesInElement; }
  inline ParticleList* firstSorbedParticleList() const { return ParticlesInElementSorbed; }
  virtual double Concentration(long ElementNo) const;
  virtual double ConcentrationSorbed(long ElementNo) const { return 0.0; };
  virtual double Retardation(long ElementNo) const { return 1.0; }
  virtual void CalculateRetardation() {  }
  virtual void InitSorption( std::string &fn,vector3& variance, int predefined) {  }
  virtual void CalculateSorption() {  }
  virtual void CalculateConcentration(double simtime);
  void CalculateConcentration2(ConstData* conc, ParticleList* plist);
  virtual void CalculateConcentrationInit();
  virtual void WriteMoments(double simtime);
  virtual std::string SorptionType();
  std::string DecayType();
  std::string UptakeType();
  const std::string & ParticleName() const { return Name; }
};


/** The base class for sorption by transition probability.
 * This class contains all properties of a particle with a sorption
 * by transition probability. Every sorption method that wants to
 * implement sorption by transition probability has to inherit from this class.
 */
class ParticleClassTP: public ParticleClass {
protected:
  // Sorption
  int sorparamno;
  ConstData *concentrationSorbed;
public:
  ParticleClassTP(PartraceClass* pt, long no, double m,
		  std::string &name, std::string &abbrev, 
		  int injno, InjectionClass *inj, double c0);
  virtual ~ParticleClassTP();
  virtual void CalculateConcentrationInit() {}
  virtual void CalculateConcentration(double simtime);
  void CalculateSorption();
  virtual void InitSorption(std::string &fn, vector3& sorpparam, int predefined);
  virtual std::string SorptionType();
  virtual double ConcentrationSorbed(long ElementNo) const;
  virtual void CalcSorpProb(const double* soilprop, const double* sorparm,
			    double c, double wc, double& p1, double& p2);
  virtual void WriteMoments(double simtime);
};


/** Linear Sorption by transition probability.
 * Implements the linear sorption model: Cs = kd * Cw
 *
 * \li Cs = sorbed concentration [M/M]
 * \li Cw = concentration in water [M/L**3]
 * \li kd = distribution coefficient [L**3/M]
 */
class ParticleClassTPLinear: public ParticleClassTP {
public:
  ParticleClassTPLinear(PartraceClass* pt, long no, double m,
			std::string &name, std::string &abbrev,
			int injno, InjectionClass *inj,	double c0):
    ParticleClassTP(pt, no, m, name, abbrev, injno, inj, c0) { sorparamno=1; }
  virtual ~ParticleClassTPLinear() { }
  virtual std::string SorptionType();
  virtual void CalcSorpProb(const double* soilprop, const double* sorparm,
			    double c, double wc, double& p1, double& p2);
};


/** Linear Sorption by transition probability (non-equilibrium).
 * see: ParticleClassTPLinear
 */
class ParticleClassTPLinearNE: public ParticleClassTP {
public:
  ParticleClassTPLinearNE(PartraceClass* pt, long no, double m,
			  std::string &name, std::string &abbrev, 
			  int injno, InjectionClass *inj, double c0):
    ParticleClassTP(pt, no, m, name, abbrev, injno, inj, c0) { sorparamno=2; }
  virtual ~ParticleClassTPLinearNE() {}
  virtual std::string SorptionType();
  virtual void CalcSorpProb(const double* soilprop, const double* sorparm,
			    double c, double wc, double& p1, double& p2);
};


/** Freundlich Sorption by transition probability.
 * Implements the freundlich sorption model: Cs = kd * Cw**n
 *
 * \li Cs = mass of solute sorbed per dry unit weight of solid [M/M]
 * \li Cw = concentration of solute in solution [M/L**3]
 * \li kd = Freundlich sorption coefficient [L**3/M]
 * \li  n = constant
 */
class ParticleClassTPFreundlich: public ParticleClassTP {
public:
  ParticleClassTPFreundlich(PartraceClass* pt, long no, double m,
			    std::string &name, std::string &abbrev, 
			    int injno, InjectionClass *inj, double c0):
    ParticleClassTP(pt, no, m, name, abbrev, injno, inj, c0) { sorparamno=2; }
  virtual ~ParticleClassTPFreundlich() {}
  virtual std::string SorptionType();
  virtual void CalcSorpProb(const double* soilprop, const double* sorparm,
			    double c, double wc, double& p1, double& p2);
};


/** Freundlich Sorption by transition probability (non-equilibrium).
 * see: ParticleClassTPFreundlich
 */
class ParticleClassTPFreundlichNE: public ParticleClassTP {
public:
  ParticleClassTPFreundlichNE(PartraceClass* pt, long no, double m,
			      std::string &name, std::string &abbrev, 
			      int injno, InjectionClass *inj, double c0):
    ParticleClassTP(pt, no, m, name, abbrev, injno, inj, c0) { sorparamno=3; }
  virtual ~ParticleClassTPFreundlichNE() {}
  virtual std::string SorptionType();
  virtual void CalcSorpProb(const double* soilprop, const double* sorparm,
			    double c, double wc, double& p1, double& p2);
};


/** Langmuir Sorption by transition probability.
 * Implements the langmuir sorption model: Cs = (K*Kmax*Cw) / (1+k*Cw)
 *
 * \li Cs = mass of solute sorbed per dry unit weight of solid [M/M]
 * \li Cw = concentration of solute in solution [M/L**3]
 * \li K  = Langmuir constant [L**3/M]
 * \li Kmax = concentration of sorption sites or the maximum sorption capacity
 */
class ParticleClassTPLangmuir: public ParticleClassTP {
public:
  ParticleClassTPLangmuir(PartraceClass* pt, long no, double m,
			  std::string &name, std::string &abbrev, 
			  int injno, InjectionClass *inj, double c0):
    ParticleClassTP(pt, no, m, name, abbrev, injno, inj, c0) { sorparamno=2; }
  virtual ~ParticleClassTPLangmuir() {}
  virtual std::string SorptionType();
  virtual void CalcSorpProb(const double* soilprop, const double* sorparm,
			    double c, double wc, double& p1, double& p2);
};


// Sorption by transition probability (TP) langmuir non-equilibrium
/** Langmuir Sorption by transition probability (non-equilibrium).
 * see: ParticleClassTPLangmuir
 */
class ParticleClassTPLangmuirNE: public ParticleClassTP {
public:
  ParticleClassTPLangmuirNE(PartraceClass* pt, long no, double m,
			    std::string &name, std::string &abbrev, 
			    int injno, InjectionClass *inj, double c0):
    ParticleClassTP(pt, no, m, name, abbrev, injno, inj, c0) { sorparamno=3; }
  virtual ~ParticleClassTPLangmuirNE() {}
  virtual std::string SorptionType();
  virtual void CalcSorpProb(const double* soilprop, const double* sorparm,
			    double c, double wc, double& p1, double& p2);
};


/** Freundlich Sorption by retardation.
 * Instead of using transition probabilities the sorption is described
 * by a retardation factor. This factor is calculated in every element
 * by solving the equation: c = cw + cs = cw + k*cw**n for the unknown
 * variable cw with a Newton iteration scheme.
 * 
 * see: ParticleClassTPFreundlich
 */
class ParticleClassR: public ParticleClass {
protected:
  // Sorption
  int sorparamno;
  ConstData *retardation;
public:
  ParticleClassR(PartraceClass* pt, long no, double m,
		 std::string &name, std::string &abbrev, 
		 int injno, InjectionClass *inj, double c0);
  virtual ~ParticleClassR();
  virtual void InitSorption(std::string &fn, vector3& sorpparam, int predefined);
  virtual void CalculateConcentrationInit() {}
  virtual double Retardation(long ElementNo) const;
  virtual double Concentration(long ElementNo) const;
  virtual double ConcentrationSorbed(long ElementNo) const;
  virtual void CalculateConcentration(double simtime);
  virtual void CalculateRetardation();
  virtual void WriteMoments(double simtime);
  virtual std::string SorptionType();
};


/** Freundlich Sorption with interactions.
 * Before calculating the freundlich sorption of the particles Pa, Pb and Pab 
 * the interaction between them is calcutated. At the moment only the following
 * model is available:
 *
 * Cab = k * Ca**na * Cb**nb
 * 
 * \li Ca  = concentration in solution for particle Pa
 * \li Cb  = concentration in solution for particle Pb
 * \li Cab = concentration in solution for particle Pab 
 * \li k = interaction parameter
 * \li na = interaction parameter
 * \li nb = interaction parameter
 * 
 * see: ParticleClassTPFreundlich and InteractClass
 */
class ParticleClassInterF: public ParticleClassTP {
public:
  ParticleClassInterF(PartraceClass* pt, long no, double m,
		      std::string &name, std::string &abbrev, 
		      int injno, InjectionClass *inj, double c0);
  virtual ~ParticleClassInterF() { }
  virtual std::string SorptionType();
  virtual void CalculateSorption() {  }
};

#endif
