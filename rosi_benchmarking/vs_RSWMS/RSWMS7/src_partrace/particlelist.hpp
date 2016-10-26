#ifndef PARTICLELISTCLASS
#define PARTICLELISTCLASS

#include "particleclass.hpp"
#include "definitions.hpp"
#include "runtime.hpp"

class PartraceClass;

/** Class for handling all particle types.
 * This class contains a list of all particle types and routines
 * to iterate over this list.
 */
class ParticleListClass {

protected:
  ParticleClass *top;     // pointer to first ParticleClass
  ParticleClass *bottom;  // pointer to last ParticleClass
  void DeleteList();
  const int RestartCode;
  bool FirstRestart;
  std::string MoveCommand;
  std::string DeleteCommand;

public:
  std::string RestartFile;
  std::string RestartFileOut;
  bool LoadRestart;
  bool WriteRestart;
  PartraceClass *partrace;
  double SaveTime;
  double SaveTimeStep;

  // Constructor
  ParticleListClass();
  // Destructor
  virtual ~ParticleListClass();

  // methods
  inline ParticleClass* TopOfList() const { return top; }
  inline ParticleClass* NextInList(ParticleClass* p) const { return p->next; }
  ParticleClass* FindParticleType(std::string &name);
  void AddInList(ParticleClass* p);
  void CalculateAllMicrobialDecay();
  void CalculateAllUptake();
  bool InjectAllParticles();
  void MoveAllParticles();
  void CalculateAllConcentrations(double simtime);
  void SaveAllConcentrations(double simtime);
  void OpenAllBTCs(double simtime);
  void WriteAllBTCs(double simtime);
  void SetRestart(std::string &restartfilein, std::string &restartfileout,
		  double firsttime, double sd, int me);
  void WriteRestartFile(bool ready);
  void LoadRestartFile();
};

#endif
