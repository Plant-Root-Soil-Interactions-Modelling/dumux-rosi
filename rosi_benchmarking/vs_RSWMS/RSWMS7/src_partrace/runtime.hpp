#ifndef RUNTIMECLASS
#define RUNTIMECLASS

#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
class Event;
class PartraceClass;

/** Class for handling the simulation time.
 * Contains routines for setting up the simulation time, testing if
 * the simulation time reaches some specified time, incrementing
 * the simulation time and so on.
 * Contains also the routine Seconds() for getting the CPU-time
 * in seconds.
 */
class RunTimeClass {
protected:
  typedef std::list<Event *> eventlist;
  eventlist EventList;
  PartraceClass *partrace;
public:
  unsigned long TimeStep;
  double dt;
  double dt_max;
  double dt_external;      ///< Externally defined time step size.
  double SystemTime;
  double Time;
  double oldTime;
  double SimulationTime;
  double nextEventTime;
  double time_eps;
  int Progress;
  double ProgressStep;
  double ProgressTime;

  RunTimeClass();
  virtual ~RunTimeClass();

  inline void SetTime(const double t) { Time=t; }
  inline double GiveTime() const { return Time; }
  inline double GiveRunTime() const { return SimulationTime; }
  inline void SetRunTime(double t, double dt_ext) { SimulationTime=t; dt_external=dt_ext; }
  inline bool NotReady() const { return Time<SimulationTime-time_eps; }
  inline bool Reached(const double t) const { return Time>=t-time_eps; }
  inline bool NotReached(const double t) const { return Time<t-time_eps; }
  inline void setMaxTimeStep(double dtmax) { dt_max=dtmax; }
  inline double remainingRunTime() const { return SimulationTime-Time; }
  inline bool firstTimeStep() const { return oldTime==0.0; }
  double Seconds() const;
  void Init(PartraceClass *pt, double simtime);
  void Increment();
  bool progress(int &pro);
  void insert(Event *e);
  void resetEvents();
};

#endif
