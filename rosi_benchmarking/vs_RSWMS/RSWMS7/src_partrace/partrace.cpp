#include <iostream>
#include <sstream>
#include <iomanip>
#include "partrace.hpp"
#include "parallel.hpp"
#include "elements.hpp"
#include "random.hpp"
#include "interact.hpp"
#include "distdata.hpp"

using namespace std;

string PartraceClass::NewHandlerText;

PartraceClass::PartraceClass(int caller, Parallel *para):
  PartraceCaller(caller), parallel(para), elements(0), Random(0)
{
  particles.partrace=this;
}


PartraceClass::~PartraceClass()
{
  delete elements;
  delete Random;
}


int PartraceClass::run()
{ 
  int progress;
  int me=parallel->mycpu();
  bool inject, boundary_inject, concentration_inject;
	//ElementsClass *elements=this->elements;

  cout<<"run partrace at "<<RunTime.Time<<" until "<<RunTime.SimulationTime<<endl;
  //Timeloop 
  using namespace std;

  parallel->sync();
  if(Debugging()) {debug<<"start of time loop "<<RunTime.Seconds()<<"secs"<<endl;  }
  if(PartraceCaller && elements->DispersionMethod==3) elements->FiniteVolToNodes(); // interpolate FiniteVolVelocity to Nodes
  
  while(RunTime.NotReady()) {
    if(!PartraceCaller) { elements->ReadFlowFields();
      if(this->elements->VelType==1) {elements->ElementsToNodes(); }
    }
    RunTime.Increment();
    parallel->sync();
    if(Debugging()) debug<<"time="<<RunTime.Time<<" dt="<<RunTime.dt<<" SimTime="<<RunTime.SimulationTime<<endl;
    inject=boundary_inject=concentration_inject=false;
    concentration_inject=elements->setAllBoundaryConcentration(RunTime, 1);
    if(Debugging()) debug<<"setAllBoundaryConcentration"<<RunTime.Seconds()<<"secs"<<endl;
    boundary_inject=elements->doAllBoundaryInjections(RunTime);
    if(Debugging()) debug<<"doAllBoundaryInjections"<<RunTime.Seconds()<<"secs"<<endl;
    inject=particles.InjectAllParticles();
    if(inject || boundary_inject || concentration_inject) {
      if(Debugging()) debug<<"InjectAllParticles "<<RunTime.Seconds()<<"secs"<<endl;
      particles.CalculateAllConcentrations(RunTime.oldTime);
      if(Debugging()) debug<<"CalculateAllConcentrations "<<RunTime.Seconds()<<"secs"<<endl;
      particles.SaveAllConcentrations(RunTime.oldTime);
      if(Debugging()) debug<<"SaveAllConcentrations "<<RunTime.Seconds()<<"secs"<<endl;
      InteractClass::CalculateAllInteractions();
      if(me==0 && RunTime.oldTime==0) particles.WriteAllBTCs(RunTime.oldTime);
    }
    particles.CalculateAllMicrobialDecay();
    parallel->sync();
    if(Debugging()) debug<<"CalculateAllMicrobialDecay "<<RunTime.Seconds()<<"secs"<<endl;
    particles.CalculateAllUptake();
    parallel->sync();
    if(Debugging()) debug<<"CalculateAllUptake "<<RunTime.Seconds()<<"secs"<<endl;
    particles.MoveAllParticles();
    parallel->sync();
    elements->setAllBoundaryConcentration(RunTime, 2);
    elements->setAllGlobalConcentration(RunTime);
    if(Debugging()) debug<<"MoveAllParticles "<<RunTime.Seconds()<<"secs"<<endl;
    particles.CalculateAllConcentrations(RunTime.Time);
    parallel->sync();
    if(Debugging()) debug<<"CalculateAllConcentrations "<<RunTime.Seconds()<<"secs"<<endl;
    particles.SaveAllConcentrations(RunTime.Time);
    particles.WriteRestartFile(false);
    parallel->sync();
    if(Debugging()) debug<<"SaveAllConcentrations "<<RunTime.Seconds()<<"secs"<<endl;
    InteractClass::CalculateAllInteractions();
    if(me==0) particles.WriteAllBTCs(RunTime.Time);
    parallel->sync();
    if(Debugging()) debug<<"WriteAllBTCs "<<RunTime.Seconds()<<"secs"<<endl;
    if(me==0 && !PartraceCaller && RunTime.progress(progress))
      std::cerr<<progress<<"%, execution time: "<<RunTime.Seconds()<<"secs"<<std::endl;
  } //Timeloop
  
  return 0;
}


void PartraceClass::error(const string &text)
{
  cerr<<"Uncorrectable error while running Partrace\n";
  cerr<<text<<endl;
  Parallel::abort("Program terminated.");
}


void PartraceClass::SetNewHandlerText(const std::string &text)
{
  NewHandlerText=text;
}


void PartraceClass::Report(const string &text)
{
  if(parallel->mycpu()==0) cout<<text<<endl;
}


void PartraceClass::Debug(const string &text)
{
  if(debugging) debug<<text<<endl;
}


void PartraceClass::start_debug(bool deb)
{
  debugging=deb;
  if(debugging) {
    ostringstream os;
    os<<setw(5)<<setfill('0')<<parallel->mycpu();
    string fn="partrace.debug." + os.str();
    debug.open(fn.c_str());
    if(debug.fail()) {
      cerr<<"Can not open debug file "<<fn<<endl;
      parallel->abort("Program terminated.", 100);
    }
  }
}


void PartraceClass::SetFNSequence(string& fn, int no, unsigned int w)
{
  ostringstream ext;
  ext.width(w);
  ext.fill('0');
  ext<<no;
  unsigned int len=fn.size();
  unsigned int p=fn.find_last_not_of("0123456789");
  if(p==string::npos) p=len-w-1;
  if(len<w || len-p<=w) 
    fn+=ext.str();
  else
    fn.replace(p+1,w,ext.str());
}


void PartraceClass::Strip(string &s)
{
  string::size_type p=s.find_last_not_of(' ');
  if(p==string::npos) s="";
  else s.erase(p+1);
}


void PartraceClass::appendSlash(string &path)
{
  if(path.empty()) return;
  if(path[path.size()-1] != '/') path+='/';
}


void PartraceNewHandler()
{
  cerr<<"Uncorrectable error while running Partrace\n";
  cerr<<"No more memory for: "<<PartraceClass::NewHandlerText<<endl;;
  Parallel::abort("Program terminated.");
}
