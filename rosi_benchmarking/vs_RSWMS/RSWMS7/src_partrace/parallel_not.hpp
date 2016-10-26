#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include <string>

class MPIWindow
{
protected:
  void *windata;     ///< Pointer to the data.
  int dim;           ///< Number of elements in windata.
  int typesize;      ///< length of one element in windata.
  MPIWindow() {}
public:
  MPIWindow(int nel, int len);
  virtual ~MPIWindow();
  void get(void *buffer, int from, int size, int pe) {}
  void put(void *buffer, int from, int size, int pe) {}
  void sum(int  *buffer, int from, int size, int pe) {}
  void fence(int assert=0) {}
  inline void* data() const { return windata; }
};

class Parallel
{
protected:
  int numpes;
  int mype;
  char *processor_name;
  double startwtime, endwtime;
  int initMPI;
public:
  Parallel();
  Parallel(int &argc, char** &argv);
  virtual ~Parallel();
  void output(const std::string &s);
  inline int cpus() const { return numpes; }
  inline int mycpu() const { return mype; }
  inline void sync() {}
  static void abort(const std::string &s, int returncode=1);
  double getTime();
  inline void send(int to_pe, int *data, int size) {}
  inline void receive(int from_pe, int *data, int size) {}
  inline void send(int to_pe, double *data, int size) {}
  inline void receive(int from_pe, double *data, int size) {}
  void sum(double *in, double *out, int len, int pe);
  void sum(double *in, double *out, int len);
  void sum(  long *in,   long *out, int len);
  inline void sum(double &inout) {}
  inline void sum(int    &inout) {}
  inline void sum(double *inout, int len) {}
  inline void sum(  long *inout, int len) {}
  inline void broadcast(double *inout, int len, int root=0) {}
  inline void broadcast(  long *inout, int len, int root=0) {}
  inline void broadcast(   int *inout, int len, int root=0) {}
  inline void broadcast(  void *inout, int len, int root=0) {}
  inline void broadcast(std::string &s, int root=0) {}
  void gather(int*    inbuffer, int*    outbuffer, int len);
  void gather(double* inbuffer, double* outbuffer, int len);
};

#endif
