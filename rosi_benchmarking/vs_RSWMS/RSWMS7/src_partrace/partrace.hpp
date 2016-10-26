#ifndef PARTRACE_H
#define PARTRACE_H
#include <string>
#include <fstream>
#include "runtime.hpp"
#include "particlelist.hpp"

class Parallel;
class ElementsClass;
class RandomClass;
class DecayClass;
class UptakeClass;

class PartraceClass {
protected:
  bool debugging;
  PartraceClass(): debugging(false) {}
  
public:
  static std::string NewHandlerText;
  int PartraceCaller;   ///< 0=Partrace  1=ParSWMS  2=muphi
  int distributedData;  ///< distribute data over the processors.
  bool PrintMass; ///< Print also masses in the concentration output.
  Parallel *parallel;
  ElementsClass *elements;
  RandomClass *Random;
  RunTimeClass RunTime;
  ParticleListClass particles;
  DecayClass *Decay;
  UptakeClass *Uptake;
  std::ofstream debug;
  // member functions
  PartraceClass(int caller, Parallel *para);
  virtual ~PartraceClass();
  inline bool Debugging() const { return debugging; }
  int run();
  void input(std::string& name);
  void error(const std::string &text);
  void SetNewHandlerText(const std::string &text);
  void Report(const std::string &text);
  void Debug(const std::string &text);
  void start_debug(bool deb=false);
  void SetFNSequence(std::string& fn, int no, unsigned int w);
  void Strip(std::string &s);
  void appendSlash(std::string &path);
};

void PartraceNewHandler();

#endif
