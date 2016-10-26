#ifndef INTERACTCLASS
#define INTERACTCLASS

#include <string>

class ParticleClass;
class PartraceClass;

/** Class for interactions.
 *  This class defines the interactions between particles.
 */
class InteractClass {

protected:
  std::string InteractName;
  ParticleClass *TargetParticle; // pointer to target particle data
  ParticleClass *Particle1, *Particle2; // pointer to source particle data
  // internal list of interactions
  static InteractClass *top, *bottom;
  InteractClass* next;
  // parameters for interaction
  enum { paramNo=10 };
  double param[paramNo];
  const int resultsize;
  PartraceClass *partrace; ///< pointer to global partrace data.

public:
  InteractClass(std::string &name, ParticleClass* target, ParticleClass* part1,
		ParticleClass* part2, double* p, PartraceClass *pt);
  virtual ~InteractClass() { }

  static InteractClass* TopOfList() { return top; }
  inline InteractClass* NextInList() const { return next; }

  static void CalculateAllInteractions();
  void CalculateInteraction();
  inline virtual std::string InteractionName() const { return InteractName; }
  void CreateMolecules(long ele, int noab);
  void CrackMolecules(long ele, int noab);
  void Sorption(long ele, ParticleClass* particles, double x);

};

#endif


