#include <cmath>
#include "decay.hpp"
#include "runtime.hpp"
#include "elements.hpp"
#include "partrace.hpp"
#include "particleclass.hpp"
#include "random.hpp"
#include "distdata.hpp"

using namespace std;


DecayClass::DecayClass(PartraceClass* pt, ParticleClass* part):
      pclass(part), partrace(pt)
{
  decaytype="No decay";
  decayfield=0;
}


DecayClass1::DecayClass1(PartraceClass* pt, ParticleClass* part, double rate,
			 std::string& DecayField, int d): 
  DecayClass(pt, part)
{
  partrace->SetNewHandlerText("decayfield");
  if(d==1) {
    decaytype="Homogeneous zero order decay";
    decayfield=new ConstData(partrace->parallel, "decayfield", partrace->elements->TotalNumber(), 1);
    decayfield->init(rate);
  }
  else {
    decaytype="Heterogeneous zero order decay";
    int ibuffer[4];
    partrace->elements->ReadField(DecayField, "Decay-rate file", "#decayfile",
				  ibuffer, decayfield);
  }
}

DecayClass1::~DecayClass1()
{
  //delete decayfield;
}


double DecayClass1::DecayProbability(double decay, double conc, double dt, double wc)
{ 
  register double fact;
  fact=(decay*dt+conc) / conc;
  return fact<0.0 ? 0.0 : fact;
}

double DecayClass2::DecayProbability(double decay, double conc, double dt, double wc)
{ 
  return exp(decay*dt);
}

// double DecayClass3::DecayProbability(double decay, double conc, double dt, double wc)
// { 
//   register double fact; 
//   double V_max = 5.;
//   double K_m = 0.5;                                                                                                                                              
//   double f = 0. ;                                                   
                                                                        
//   fact= 1 - ((V_max / (K_m+conc) + f) * decay * dt/wc);                                                                                                                          
//   //if (fact != 1) {cout << fact << endl;}                                                                                                                                       
//   return fact;    

// }

// double DecayClass4::DecayProbability(double decay, double conc, double dt, double wc)
// { 
//   return (1- decay*dt/wc);
// }

void DecayClass1::CalculateMicrobialDecay()
{
  double conc, decay, prob, wc, mass_p, mass_n, mass;
  double dt=partrace->RunTime.dt;
  Particle *particle;
  ParticleList* plist=pclass->firstParticleList();
  ParticleList* plist2=pclass->firstSorbedParticleList();
  int e_to=decayfield->lines();
  ConstData *watercontent=partrace->elements->watercontent;
  mass_p=0.;
  mass_n=0.;
  for(int e=0; e<e_to; e++,plist++) {
    // decay for particles in solution
    if(plist->count()>0) {
      conc=pclass->Concentration(e);
      decay=decayfield->get(e);
      wc = watercontent->get(e);
      particle=plist->start();
      prob=DecayProbability(decay, conc, dt, wc);
      while(particle) {
	mass=particle->get_mass();
	mass_p +=mass;
	particle->reduce_mass(prob);
	mass=particle->get_mass();
	mass_n += mass;
	particle=plist->next(particle);

      }
    }
    // decay for sorbed particles
    if(plist2) {
      if(plist2->count()>0) {
	conc=pclass->ConcentrationSorbed(e);
	decay=decayfield->get(e);
	particle=plist2->start();
	prob=DecayProbability(decay, conc, dt, wc);
	while(particle) {
	  particle->reduce_mass(prob);
	  particle=plist2->next(particle);
	}
      }
      plist2++;
    }
  }
}


DecayClass2::DecayClass2(PartraceClass* pt, ParticleClass* part, 
			 double rate, std::string& DecayField, int d): 
  DecayClass1(pt, part, rate, DecayField, d)
{
  if(d==1) decaytype="Homogeneous first order decay";
  else decaytype="Heterogeneous first order decay";
}

// DecayClass3::DecayClass3(PartraceClass* pt, ParticleClass* part, 
// 			 double rate, std::string& DecayField, int d): 
//   DecayClass1(pt, part, rate, DecayField, d)
// {
//   if(d==1) decaytype="Homogeneous second order decay";
//   else decaytype="Heterogeneous second order decay";
// }

// DecayClass4::DecayClass4(PartraceClass* pt, ParticleClass* part, 
// 			 double rate, std::string& DecayField, int d): 
//   DecayClass1(pt, part, rate, DecayField, d)
// {
//   if(d==1) decaytype="Homogeneous third order decay";
//   else decaytype="Heterogeneous third order decay";
// }
