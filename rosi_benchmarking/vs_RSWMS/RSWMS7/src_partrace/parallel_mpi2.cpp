#include "parallel.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
using namespace std;


Parallel::Parallel():
  numpes(1), mype(0), processor_name(0), initMPI(0)
{
  numpes = MPI::COMM_WORLD.Get_size();
  mype   = MPI::COMM_WORLD.Get_rank();
  startwtime = MPI::Wtime();
}


Parallel::Parallel(int &argc, char** &argv):
  numpes(1), mype(0), processor_name(0), initMPI(1)
{
  MPI::Init(argc,argv);
  numpes = MPI::COMM_WORLD.Get_size();
  mype   = MPI::COMM_WORLD.Get_rank();
  processor_name=new char[MPI::MAX_PROCESSOR_NAME];
  if(mype==0) cout<<"Starting "<<argv[0]<<" on "<<numpes<<" PEs."<<endl<<flush;
  int  namelen;
  MPI::Get_processor_name(processor_name, namelen);
  startwtime = MPI::Wtime();
  ostringstream os;
  os << "Process " << mype << " of " << numpes
     << " is on " <<  processor_name;
  output(os.str());
}


Parallel::~Parallel()
{
  if(initMPI) {
    delete [] processor_name;
    MPI::Finalize();
  }
}


void Parallel::output(const std::string &s)
{
  sync();
  for(int p=0;p<numpes;p++) {
    if(p==mype) cout<<"PE"<<mype<<": "<<s<<endl<<flush;
    sync();
  }
}


void Parallel::abort(const std::string &s, int returncode)
{
  cerr<<s<<endl;
  MPI::COMM_WORLD.Abort(returncode);
}


void Parallel::sync()
{
  MPI::COMM_WORLD.Barrier();
}


double Parallel::getTime()
{
  //Time in sec.
  return MPI::Wtime()-startwtime;
}


void Parallel::send(int to_pe, int *data, int size)
{
  const int tag=1;
  // send data to processor to_pe.
  MPI::COMM_WORLD.Send(data, size, MPI::INT, to_pe, tag);
}


void Parallel::receive(int from_pe, int *data, int size)
{
  const int tag=1;
  // receive data from processor from_pe.
  MPI::COMM_WORLD.Recv(data, size, MPI::INT, from_pe, tag);
}


void Parallel::send(int to_pe, double *data, int size)
{
  const int tag=2;
  // send data to processor to_pe.
  MPI::COMM_WORLD.Send(data, size, MPI::DOUBLE, to_pe, tag);
}


void Parallel::receive(int from_pe, double *data, int size)
{
  const int tag=2;
  // receive data from processor from_pe.
  MPI::COMM_WORLD.Recv(data, size, MPI::DOUBLE, from_pe, tag);
}


void Parallel::sum(double *in, double *out, int len, int pe)
{ 
  // Summation to local node pe
  MPI::COMM_WORLD.Reduce((void *) in, (void *) out, len, MPI::DOUBLE, MPI::SUM, pe); 
}


void Parallel::sum(double *in, double *out, int len)
{ 
  // global summation
  MPI::COMM_WORLD.Allreduce((void *) in, (void *) out, len, MPI::DOUBLE, MPI::SUM); 
}


void Parallel::sum(double &in)
{ 
  // global summation
  double out;
  MPI::COMM_WORLD.Allreduce((void *) &in, (void *) &out, 1, MPI::DOUBLE, MPI::SUM); 
  in=out;
}


void Parallel::sum(long *in, long *out, int len)
{ 
  // global summation
  MPI::COMM_WORLD.Allreduce((void *) in, (void *) out, len, MPI::LONG, MPI::SUM); 
}


void Parallel::sum(long *in, int len)
{ 
  // global summation, automatic allocation of temporary buffer
  long *buffer=new long[len];
  MPI::COMM_WORLD.Allreduce((void *) in, (void *) buffer, len, MPI::LONG, MPI::SUM); 
  for(int i=0;i<len;i++) in[i]=buffer[i];
  delete [] buffer;
}


void Parallel::sum(int &in)
{ 
  // global summation.
  int buffer;
  MPI::COMM_WORLD.Allreduce((void *) &in, (void *) &buffer, 1, MPI::INT, MPI::SUM); 
  in=buffer;
}


void Parallel::sum(double *in, int len)
{ 
  // global summation
  double *buffer=new double[len];
  MPI::COMM_WORLD.Allreduce((void *) in, (void *) buffer, len, MPI::DOUBLE, MPI::SUM); 
  for(int i=0;i<len;i++) in[i]=buffer[i];
  delete [] buffer;
}


void Parallel::broadcast(int* buffer, int count, int root)
{ 
  MPI::COMM_WORLD.Bcast(buffer, count, MPI::INT, root);
}


void Parallel::broadcast(long* buffer, int count, int root)
{ 
  MPI::COMM_WORLD.Bcast(buffer, count, MPI::LONG, root);
}


void Parallel::broadcast(double* buffer, int count, int root)
{ 
  MPI::COMM_WORLD.Bcast(buffer, count, MPI::DOUBLE, root);
}


void Parallel::broadcast(void* buffer, int count, int root)
{ 
  MPI::COMM_WORLD.Bcast(buffer, count, MPI::BYTE, root);
}


void Parallel::broadcast(std::string &s, int root)
{ 
  const char *buffer;
  int len=s.size()+1;
  MPI::COMM_WORLD.Bcast(&len, 1, MPI::INT, root);
  if(mype==root) buffer=s.c_str();
  else buffer=new char[len];
  MPI::COMM_WORLD.Bcast((void*) buffer, len, MPI::CHAR, root);
  if(mype!=root) {
    s=buffer;
    delete [] buffer;
  }
}


void Parallel::gather(int* inbuffer, int* outbuffer, int len)
{
  // put inbuffer from cpu i into array outbuffer at position i
  MPI::COMM_WORLD.Allgather((void *) inbuffer,  len, MPI::INT,
		            (void *) outbuffer, len, MPI::INT);
}


void Parallel::gather(double* inbuffer, double* outbuffer, int len)
{
  // put inbuffer from cpu i into array outbuffer at position i
  MPI::COMM_WORLD.Allgather((void *) inbuffer,  len, MPI::DOUBLE,
		            (void *) outbuffer, len, MPI::DOUBLE);
}


MPIWindow::MPIWindow(int nel, int len):
  dim(nel), typesize(len), mpi_memory(true)
{
  int winsize_bytes=dim*typesize;
  windata=MPI::Alloc_mem(winsize_bytes, MPI::INFO_NULL);
  win=MPI::Win::Create(windata, winsize_bytes, typesize,
		       MPI::INFO_NULL, MPI::COMM_WORLD);
}


MPIWindow::MPIWindow(void *wdata, int nel, int len):
  dim(nel), typesize(len), mpi_memory(false)
{
  int winsize_bytes=dim*typesize;
  win=MPI::Win::Create(wdata, winsize_bytes, typesize,
		       MPI::INFO_NULL, MPI::COMM_WORLD);
}


MPIWindow::~MPIWindow()
{
  win.Free();
  if(mpi_memory) MPI::Free_mem(windata);
}


void MPIWindow::get(void *buffer, int from, int size, int pe)
{
  const int assert=0;
  win.Lock(MPI::LOCK_SHARED, pe, assert);
  size*=typesize;
  win.Get(buffer, size, MPI::BYTE, pe, from, size, MPI::BYTE);
  win.Unlock(pe);
}


void MPIWindow::put(void *buffer, int from, int size, int pe)
{
  const int assert=0;
  win.Lock(MPI::LOCK_SHARED, pe, assert);
  size*=typesize;
  win.Put(buffer, size, MPI::BYTE, pe, from, size, MPI::BYTE);
  win.Unlock(pe);
}


void MPIWindow::sum(int *buffer, int from, int size, int pe)
{
  const int assert=0;
  win.Lock(MPI::LOCK_SHARED, pe, assert);
  win.Accumulate(buffer, size, MPI::INT, pe, from, size, MPI::INT, MPI::SUM);
  win.Unlock(pe);
}


void MPIWindow::fence(int assert)
{
  win.Fence(assert);
}
