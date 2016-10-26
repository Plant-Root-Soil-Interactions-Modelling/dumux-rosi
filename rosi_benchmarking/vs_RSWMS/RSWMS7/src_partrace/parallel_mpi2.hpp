#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include <mpi.h>
#include <string>

class MPIWindow
{
protected:
  MPI::Win win;      ///< MPI window to the data on the other PEs.
  void *windata;     ///< Pointer to the data.
  int dim;           ///< Number of elements in windata.
  int typesize;      ///< length of one element in windata.
  bool mpi_memory;   ///< allocate mpi memory in constructor true/false.
public:
  MPIWindow(int nel, int len);
  MPIWindow(void *wdata, int nel, int len);
  virtual ~MPIWindow();
  void get(void *buffer, int from, int size, int pe);
  void put(void *buffer, int from, int size, int pe);
  void sum(int  *buffer, int from, int size, int pe);
  void fence(int assert=0);
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
  void sync();
  static void abort(const std::string &s, int returncode=1);
  double getTime();
  void send(int to_pe, int *data, int size);
  void receive(int from_pe, int *data, int size);
  void send(int to_pe, double *data, int size);
  void receive(int from_pe, double *data, int size);
  void sum(double *in, double *out, int len, int pe);
  void sum(double *in, double *out, int len);
  void sum(  long *in,   long *out, int len);
  void sum(double &inout);
  void sum(int    &inout);
  void sum(double *inout, int len);
  void sum(  long *inout, int len);
  void broadcast(double *inout, int len, int root=0);
  void broadcast(  long *inout, int len, int root=0);
  void broadcast(   int *inout, int len, int root=0);
  void broadcast(  void *inout, int len, int root=0);
  void broadcast(std::string &s, int root=0);
  void gather(int* inbuffer, int* outbuffer, int len);
  void gather(double* inbuffer, double* outbuffer, int len);
};

#endif
