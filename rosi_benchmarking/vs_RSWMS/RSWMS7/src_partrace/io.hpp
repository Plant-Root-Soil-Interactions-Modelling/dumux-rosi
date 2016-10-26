#ifndef IOCLASS
#define IOCLASS

#include "definitions.hpp"

class interdata;
class boundarydata;
class InjectionClass;

/// Temporary class used for sending the interaction data via MPI.
class interdata {
public:
  string200 name, target, part1, part2;
  enum { paramNo=10 };
  double param[paramNo];
};

/// Temporary class used for sending the main input data via MPI.
class indata {
public:
  string200 InpPath, OutPath, DispField, VelField;
  string200 fn_conc, fn_btc, fn_mom, fn_wc, CoordFile;
  string200 RestartFile, RestartFileOut, fnParticles;
  string200 fnReflectInside;
  double simtime;
  double SaveDistanceConc, SaveDistanceMom, SaveDistanceRestart;
  double FirstTimeRestart, SaveParticlesTimeStep;
  triple NumberOfElements;
  int GeometryType, DiffusionType, DistributeData, ConvectionMethod, DispersionMethod, CoefAtBoundary, SecDispersiveNew;
  int ReflectionMethod, NoSurfaceReflection, EleVxVy0, TimeSplitting, Recalculate, sub, RandomNumber;
  int UseDispField,UseVelField, PorosityHom, UseWC, PrintMass;
  int ReflectInside;
  vector3 DispMean, DispWeight, VelMean, SideLength; 
  vector3 VelVariance;
  double VelTimeStepFactor;
  double MaxTime;
  double DispDiffusion, Porosity, BulkDensity;
  int CompPartPerPro,CompRand,CompSeed, OutputType, RestartType;
  int PartTypes, VelType, WCType;
  int BTCNo;
  vector3* BTC_List;
  int InterNo;
  interdata* Inter_List;
  int BoundariesNo;
  boundarydata* Boundaries;
  int ElementSequenceNo;
  // constructor, destructor
  indata();
  ~indata();
};

/// Temporary class used for sending the particle input data via MPI.
class partdata {
public:
  string200 PartName,PartAbbrev,fn_sorption,DecayField,UptakeField,RestartFromConcFile;
  int PartSorpType,PartUseSorpField;
  int UseDecayField,DecayOrder;
  int UseUptakeField,UptakeOrder;
  vector3 PartMean;
  long PartNumber;
  double PartMass,PartC0;
  vector3 PartKdNB;
  double DecayRate,UptakeRate;
  int injectionNo;
  InjectionClass *injection; // pointer to injection(s)
};

#endif
