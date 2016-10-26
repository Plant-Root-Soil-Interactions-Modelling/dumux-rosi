#include "parallel.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include <ctime>
#include <algorithm>
using namespace std;


Parallel::Parallel():
  numpes(1), mype(0), processor_name(0), initMPI(0)
{
  startwtime = (double)time(0);
}


Parallel::Parallel(int &argc, char** &argv):
  numpes(1), mype(0), processor_name(0), initMPI(0)
{
  const int max_processor_name=256;
  processor_name=new char[max_processor_name];
  processor_name[0]='\n';
  gethostname(processor_name, max_processor_name);
  cout<<"Starting "<<argv[0]<<" on "<<numpes<<" PE."<<endl<<flush;
  startwtime = (double)time(0);
  ostringstream os;
  os << "Process " << mype << " of " << numpes
     << " is on " <<  processor_name;
  output(os.str());
}


Parallel::~Parallel()
{
  delete [] processor_name;
}


void Parallel::output(const std::string &s)
{
  cout<<s<<endl<<flush;
}


void Parallel::abort(const std::string &s, int returncode)
{
  cerr<<s<<endl;
  exit(returncode);
}


double Parallel::getTime()
{
  //Time in sec.
  return (double)time(0);
}


void Parallel::sum(double *in, double *out, int len, int pe)
{ 
  // non parallel global summation is only copying
  if(in!=out) copy(in, in+len, out);
}


void Parallel::sum(double *in, double *out, int len)
{ 
  // non parallel global summation is only copying
  if(in!=out) copy(in, in+len, out);
}


void Parallel::sum(long *in, long *out, int len)
{ 
  // non parallel global summation is only copying
  if(in!=out) copy(in, in+len, out);
}


void Parallel::gather(int* inbuffer, int* outbuffer, int len)
{
   copy(inbuffer, inbuffer+len, outbuffer);
}


void Parallel::gather(double* inbuffer, double* outbuffer, int len)
{
   copy(inbuffer, inbuffer+len, outbuffer);
}


MPIWindow::MPIWindow(int nel, int len)
{
  dim=nel;
  typesize=len;
  windata=new char[dim*typesize];
}


MPIWindow::~MPIWindow()
{
  delete [] static_cast<char*>(windata);
}

