#include <cstdlib>
#include <iostream>
#include <fstream>
#include "particlelist.hpp"
#include "partrace.hpp"
#include "parallel.hpp"
#include "decay.hpp"
#include "uptake.hpp"

using namespace std;


// Constructor
ParticleListClass::ParticleListClass(): RestartCode(20060516)
{ 
  top=bottom=0;
  LoadRestart=WriteRestart=false;
  FirstRestart=true;
  SaveTime=0.0;
}


// Destructor
ParticleListClass::~ParticleListClass()
{
  DeleteList();
}


ParticleClass* ParticleListClass::FindParticleType(string &name)
{
  // Calculate Decay
  ParticleClass* p=TopOfList();
  while(p) {
    if(p->Name==name || p->Abbrev==name) return p;
    p=NextInList(p);
  }
  partrace->error(string("Particle '")+name+string("' not found"));
  return 0;
}


void ParticleListClass::AddInList(ParticleClass* p) 
{ 
  if(bottom) {
    bottom->next=p;
    bottom=p;
  }
  else top=bottom=p;
}  


void ParticleListClass::DeleteList()
{
  ParticleClass *p=top, *q;
  while(p) {
    q=p;
    p=p->next;
    delete q;
  }
}


void ParticleListClass::CalculateAllMicrobialDecay()
{
  // Calculate Decay
  ParticleClass* p=TopOfList();
  while(p) {
    p->Decay->CalculateMicrobialDecay();
    p=NextInList(p);
  }
}


void ParticleListClass::CalculateAllUptake()
{
  // Calculate Uptake
  ParticleClass* p=TopOfList();
  while(p) {
    p->Uptake->CalculateUptake(); 
    p=NextInList(p);
  }
}


void ParticleListClass::MoveAllParticles()
{
  // move particle
  ParticleClass* p=TopOfList();
  while(p) {
    p->MoveParticles();
    p=NextInList(p);
  }
}


void ParticleListClass::CalculateAllConcentrations(double simtime)
{
  ParticleClass* p=TopOfList();
  while(p) {
    p->CalculateConcentration(simtime);
    p=NextInList(p);
  }
}


void ParticleListClass::SaveAllConcentrations(double simtime)
{
  ParticleClass* p=TopOfList();
  while(p) {
    p->SaveConcentration(simtime);
    p->WriteMoments(simtime);
    p->SaveParticlePositions(simtime);
    p=NextInList(p);
  }
}


void ParticleListClass::OpenAllBTCs(double simtime)
{
  ParticleClass* p=TopOfList();
  while(p) {
    p->OpenBTC(simtime);
    p=NextInList(p);
  }
}


void ParticleListClass::WriteAllBTCs(double simtime)
{
  ParticleClass* p=TopOfList();
  while(p) {
    p->WriteBTC(simtime);
    p=NextInList(p);
  }
}


bool ParticleListClass::InjectAllParticles()
{
  int injected=0;
  ParticleClass* p=TopOfList();
  while(p) {
    if(p->InjectParticlePart()) injected=1;
    p=NextInList(p);
  }
  partrace->parallel->sum(injected);
  return injected;
}


void ParticleListClass::SetRestart(string &restartfilein,string &restartfileout,
				   double firsttime, double sd, int me)
{
  RestartFile=restartfilein;
  RestartFileOut=restartfileout;
  if(RestartFile!=string("")) {
    LoadRestart=true;
    partrace->SetFNSequence(RestartFile, me, 5);
    if(me==0) cout<<"Using Restart-File "<<RestartFile<<endl;
  }
  if(RestartFileOut!=string("")) {
    WriteRestart=true;
    partrace->SetFNSequence(RestartFileOut, me, 5);
  }
  SaveTime=firsttime*60.0; // savetime in seconds
  if(sd<=0.0) SaveTimeStep=3600.0; // save every hour;
  else SaveTimeStep=sd*60.0;
}


void ParticleListClass::WriteRestartFile(bool ready)
{
  if(WriteRestart==false) return;
  int ior=0;
  double elapsedTime=partrace->RunTime.Seconds();
  partrace->parallel->broadcast(&elapsedTime, 1);
  if(elapsedTime<SaveTime && !ready) return;
  if(FirstRestart) {
     MoveCommand=string("mv ")+RestartFileOut+string(" ")+RestartFileOut+string(".old");
     DeleteCommand=string("rm -f ")+RestartFileOut+string(".old");
  }
  else ior=system(MoveCommand.c_str());
  fstream res(RestartFileOut.c_str(),ios::out);
  if(res.fail()) 
    partrace->error(string("Open error for the Partrace output restart-file ")+RestartFileOut);
  res.write((char *) &RestartCode, sizeof(int));
  if(partrace->elements->RestartType != 2) {
    res.write((char *) &(partrace->RunTime.Time), sizeof(partrace->RunTime.Time));
    res.write((char *) partrace->elements->ParticlesOutside,
	      6*sizeof(*(partrace->elements->ParticlesOutside)));
    res.write((char *) &(partrace->elements->VFieldRestartPos), sizeof(int));
  }

  ParticleClass* p=TopOfList();
  while(p) {
    p->SaveParticles(res);
    p=NextInList(p);
  }
  res.close();
  partrace->parallel->sync();
  if(FirstRestart) FirstRestart=false;
  else ior=system(DeleteCommand.c_str());
  SaveTime+=SaveTimeStep;
}


void ParticleListClass::LoadRestartFile()
{
  string line;
  double time;
  int code;

  fstream res(RestartFile.c_str(),ios::in);
  if(res.fail()) partrace->error(string("Open error for the Partrace restart-file ")+RestartFile);
  res.read((char *) &code, sizeof(int));
  if(code!=RestartCode) partrace->error("Wrong file format for the partrace restart-file");
  if(partrace->elements->RestartType != 2) {
    res.read((char *) &time, sizeof(time));
    partrace->RunTime.SetTime(time);
    res.read((char *) partrace->elements->ParticlesOutside,
	     6*sizeof(*(partrace->elements->ParticlesOutside)));
    res.read((char *) &(partrace->elements->VFieldPos), sizeof(int));
		partrace->RunTime.resetEvents();
    partrace->particles.CalculateAllConcentrations(partrace->RunTime.Time);
    partrace->particles.SaveAllConcentrations(partrace->RunTime.Time);
  }

  ParticleClass* p=TopOfList();
  while(p) {
    p->LoadParticles(res);
    p=NextInList(p);
  }
  res.close();
}
