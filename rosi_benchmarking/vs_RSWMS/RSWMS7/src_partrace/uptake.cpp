#include <cmath>
#include "uptake.hpp"
#include "runtime.hpp"
#include "elements.hpp"
#include "partrace.hpp"
#include "particleclass.hpp"
#include "random.hpp"
#include "distdata.hpp"

using namespace std;


UptakeClass::UptakeClass(PartraceClass* pt, ParticleClass* part):
      pclass(part), partrace(pt)
{
  uptaketype="No solute uptake / Exclusion";
  uptakefield=0; 
}


UptakeClass1::UptakeClass1(PartraceClass* pt, ParticleClass* part, double rate,
			   std::string& UptakeField, int d): 
  UptakeClass(pt, part)
{
  uptaketype="Passive solute uptake";
  partrace->SetNewHandlerText("uptakefield");
  uptakefield=new LocalData(partrace->parallel, "uptakefield", partrace->elements->TotalNumber(), 1);
  uptakefield->init(rate);
}

UptakeClass1::~UptakeClass1()
{
  //delete uptakefield;
}


double UptakeClass1::UptakeProbability(double uptake, double conc, double dt, double wc)
{ 
  return (1- uptake*dt/wc); //passiv uptake
}


double UptakeClass2::UptakeProbability(double uptake, double conc, double dt, double wc)
{ 
  register double fact; 
  const double V_max = 5.;
  const double K_m = 0.5;                                                                                                       
  const double f = 0. ;                                                   
                                        
  fact= 1 - ((V_max / (K_m+conc) + f) * uptake * dt/wc); // Michaelis Menten uptake                          
  //if (fact != 1) {cout << fact << endl;}                                                                                                                                       
  return fact;    

}


double UptakeClass3::UptakeProbability(double uptake, double conc, double dt, double wc)
{ 
  // concentration dependent uptake - calculated in rswms
  //std::cout <<"uptake " << uptake <<   std::endl;      
  return (uptake);
}


void UptakeClass1::CalculateUptake()
{
  double conc, uptake, prob, wc, mass_p, mass_n, mass;
  double dt=partrace->RunTime.dt;
  Particle *particle;
  ParticleList* plist=pclass->firstParticleList();
  int e_to=uptakefield->lines();
  ConstData *watercontent=partrace->elements->watercontent;
  mass_p=0.;
  mass_n=0.;
  for(int e=0; e<e_to; e++,plist++) {
    // uptake ONLY from particles in solution
    if(plist->count()>0) {
      conc=pclass->Concentration(e);
      uptake=uptakefield->get(e);
      wc = watercontent->get(e);
      particle=plist->start();
      prob=UptakeProbability(uptake, conc, dt, wc); 
      //std::cout << " e "<< e << " prob " << prob <<std::endl;   
      while(particle) {
	mass=particle->get_mass();
	mass_p +=mass;
	particle->reduce_mass(prob);
	mass=particle->get_mass();
	mass_n += mass;
	particle=plist->next(particle);
      }
    }
  }
}


UptakeClass2::UptakeClass2(PartraceClass* pt, ParticleClass* part, 
			   double rate, std::string& UptakeField, int d): 
  UptakeClass1(pt, part, rate, UptakeField, d)
{                             
  uptaketype="Active solute uptake / Michaelis Menten"; 
}

UptakeClass3::UptakeClass3(PartraceClass* pt, ParticleClass* part, 
			   double rate, std::string& UptakeField, int d): 
  UptakeClass1(pt, part, rate, UptakeField, d)
{
  uptaketype="Diffusive solute uptake";
}
