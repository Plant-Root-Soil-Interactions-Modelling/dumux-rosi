#ifndef PARTRACE_INTERFACE
#define PARTRACE_INTERFACE

class PartraceClass;
class ConstData;
class Parallel;

extern "C" {
  int *parswms2partrace;         ///< convert the index from parswms to partrace.
  /// \a Partrace object for handling all partrace types.
  PartraceClass *partrace;
  /// Parallel class used for \a Partrace.
  Parallel *partrace_parallel;
  /// contains the water content on the nodes.
  ConstData *watercontentOnNodes;
  void initPartrace();  ///< constructor.
  void closePartrace(); ///< destructor.
  void runPartrace(double &t);
  void setPartraceFields(int &inodes, double v1[], double v2[], double v3[], double theta[], double sinkE[]);
  /// not used in the current version of the interface.
  void setConcentration(int &iElm, double concPar[]);
};
#endif
