#include <iostream>
#include <fstream>
#include <iomanip>
#include "particleclass.hpp"
#include "elements.hpp"
#include "random.hpp"
#include "event.hpp"
#include "partrace.hpp"
#include "parallel.hpp"
#include "distdata.hpp"

using namespace std;

// Sorption models
// 1 = konservative
// 2 = linear adsorbing           cs=ka*cw
// 3 = freundlich isotherm        cs=k*cw**n
// 4 = langmuir isotherm          cs=ka*kb*cw/(1+kb*cw)
// 5 = linear non equilibrium
// 6 = freundlich non equilibrium
// 7 = langmuir non equilibrium
// 8 = freundlich isotherm with retardation
// 9 = interaction with freundlich isotherm

// Sorption with transition-probability
ParticleClassTP::ParticleClassTP(PartraceClass *pt, long no, double m,
	       string &name, string &abbrev, int injno, InjectionClass *inj, double c0):
    ParticleClass(pt, no, m, name, abbrev, injno, inj, c0), sorparamno(1)
{
  partrace->SetNewHandlerText("concentrationSorbed");
  if(partrace->distributedData)
    concentrationSorbed=new DistData( partrace->parallel, "concentrationSorbed", NoElxyz, 1);
  else
    concentrationSorbed=new LocalData( partrace->parallel, "concentrationSorbed", NoElxyz, 1);
  partrace->SetNewHandlerText("ParticlesInElementSorbed");
  ParticlesInElementSorbed=new ParticleList[NoElxyz];
}


ParticleClassTP::~ParticleClassTP() { 
  delete concentrationSorbed; 
  delete [] ParticlesInElementSorbed;
}


string ParticleClassTP::SorptionType()
{
  return string("Sorption by transition probability");
}


double ParticleClassTP::ConcentrationSorbed(long ElementNo) const
{
  return concentrationSorbed->get(ElementNo);
}


void ParticleClassTP::InitSorption(string &fn, vector3& sorpparam, int predefined)
{
  partrace->SetNewHandlerText("adsorption");
  if(predefined==2) {
    // heterogeneous adsorption
    int ibuffer[4];
    partrace->elements->ReadField(fn,"Adsorption file",
			     "#partrace adsorption file",ibuffer,adsorption);
  }
  else {
    // homogeneous adsorption
    adsorption=new ConstData(partrace->parallel, "adsorption", NoElxyz, sorparamno);
    adsorption->set_dataline(0, sorpparam.data());
  }
}


void ParticleClassTP::CalcSorpProb(const double* soilprop, const double* sorparm,
				   double c, double wc, double& p1, double& p2) 
{
  p1=1;
  p2=0;
}


void ParticleClassTP::CalculateConcentration(double simtime)
{
  CalculateConcentration2(concentration, ParticlesInElement);
  CalculateConcentration2(concentrationSorbed, ParticlesInElementSorbed);
  CalculateSorption();
}


void ParticleClassTP::CalculateSorption()
{
  Particle *particle, *prev, *particle2, *ptemp;
  ParticleList* plist=ParticlesInElement;
  ParticleList* plist2=ParticlesInElementSorbed;
  double sorpProb1, sorpProb2;
  double* soilprop=partrace->elements->soilprop->new_buffer();
  double* sorpprop=adsorption->new_buffer();
  double c, wc;
  RandomClass *rand=partrace->Random;

  for(long e=0; e<NoElxyz; e++,plist++,plist2++) {
    if(plist->count()==0 && plist2->count()==0)
      continue; // no particles in element
    partrace->elements->soilprop->get_line(e, soilprop);
    wc=partrace->elements->watercontent->get(e);
    adsorption->get_line(e, sorpprop);
    c=concentration->get(e);
    CalcSorpProb(soilprop, sorpprop, c, wc, sorpProb1, sorpProb2);
    // particles in solution
    particle=plist->start(prev);
    plist2->start(particle2);
    while(particle) {
      if( rand->fastdraw() <= sorpProb1 ) {
	// move particle into list of sorbed particles
	ptemp=particle;
	particle=plist->remove_from_list(particle,prev);
	particle2=plist2->insert_after(ptemp,particle2);
      }
      else particle=plist->next(particle,prev);
    } // while

    // sorbed particles
    prev=particle2;
    particle=prev->next;
    while(particle) {
      if( rand->fastdraw() <= sorpProb2 ) {
	// move particle into list of particles in solution
	ptemp=particle;
	particle=plist2->remove_from_list(particle,prev);
	plist->insert(ptemp);
      }
      else particle=plist2->next(particle,prev);
    } // while
  }
  // delete local buffers
  partrace->elements->soilprop->delete_buffer(soilprop);
  adsorption->delete_buffer(sorpprop);
}


void ParticleClassTPLinear::CalcSorpProb(const double* soilprop,
     const double* sorparm, double c, double wc, double& p1, double& p2)
{
  // probability for not sorbed particle to be sorbed
  // p1 = kd/(kd+1)
  double kd=sorparm[0]*soilprop[4]/wc;
  p1=kd/(kd+1.0);
  // probability for sorbed particle not to be sorbed
  // p2 = 1-p1 = 1/(kd+1)
  p2=1.0-p1;
}


string ParticleClassTPLinear::SorptionType()
{
  return string("Linear-Sorption by transition probability");
}


void ParticleClassTPLinearNE::CalcSorpProb(const double* soilprop,
     const double* sorparm, double c, double wc, double& p1, double& p2)
{
  // probability for not sorbed particle to be sorbed
  // p1 = b*kd*dt
  p1=sorparm[1]*sorparm[0]*soilprop[4]/wc*partrace->RunTime.dt;
  // probability for sorbed particle not to be sorbed
  // p2 = b*dt        
  p2=sorparm[1]*partrace->RunTime.dt;
}


string ParticleClassTPLinearNE::SorptionType()
{
  return string("Linear non-equilibrium Sorption by transition probability");
}


// Freundlich sorption:   Csorp = Kd * Cwat**n
// with:  Csorp = mass of solute sorbed per dry unit weight of solid (mg/kg)
//        Cwat  = concentration of solute in solution in equilibrium with the
//                mass of solute sorbed onto the solid (mg/L)
//           Kd = Freundlich sorption coefficient (L/kg)
//            n = constant
void ParticleClassTPFreundlich::CalcSorpProb(const double* soilprop,
     const double* sorparm, double c, double wc, double& p1, double& p2)
{
  double cs;
  // probability for not sorbed particle to be sorbed
  // p = kd*c**n / (c+kd*c**n)
  if(c>0.0) {
    // Cs = kd * BulkDensity / Porosity * C**n
    cs=sorparm[0]*soilprop[4]/wc * pow(c,sorparm[1]);
    p1=cs/(cs+c);
    // probability for sorbed particle not to be sorbed
    // p2 = 1 - p1 = c/(cs+c)
    p2=1.0-p1;
  }
  else {
    p1=0.0;
    p2=1.0;
  }
}


string ParticleClassTPFreundlich::SorptionType()
{
  return std::string("Freundlich-Sorption by transition probability");
}


void ParticleClassTPFreundlichNE::CalcSorpProb(const double* soilprop,
     const double* sorparm, double c, double wc, double& p1, double& p2)
{
  // probability for not sorbed particle to be sorbed
  // p1 = b*dt*kd*c**(n-1)
  if(c>0.0) p1=sorparm[2]*partrace->RunTime.dt*
	      sorparm[0]*soilprop[4]/wc * pow(c,sorparm[1]-1.0);
  else p1=0.0;
  // probability for sorbed particle not to be sorbed
  // p2 = b*dt
  p2=sorparm[2]*partrace->RunTime.dt;
}


string ParticleClassTPFreundlichNE::SorptionType()
{
  return string("Freundlich non-equilibrium Sorption by transition probability");
}


void ParticleClassTPLangmuir::CalcSorpProb(const double* soilprop,
     const double* sorparm, double c, double wc, double& p1, double& p2)
{
  // probability for not sorbed particle to be sorbed
  // p1 = ka*kb / (1+ka*kb+kb*c)
  double kakb=sorparm[0]*sorparm[1]*soilprop[4]/wc;
  p1=kakb/(1.0+kakb+sorparm[1]*c);
  // probability for sorbed particle not to be sorbed
  // p2 = 1 - p1 = (1-kb*c) / (1+ka*kb+kb*c)
  p2=1.0-p1;
}


string ParticleClassTPLangmuir::SorptionType()
{
  return string("Langmuir-Sorption by transition probability");
}


void ParticleClassTPLangmuirNE::CalcSorpProb(const double* soilprop,
     const double* sorparm, double c, double wc, double& p1, double& p2)
{
  // probability for not sorbed particle to be sorbed
  // p1 = (b*ka*kb*dt) / (1+kb*c)
  p1=sorparm[2]*sorparm[1]*sorparm[0]*soilprop[4]/wc*partrace->RunTime.dt /
    (1.0+sorparm[2]*partrace->RunTime.dt);
  // probability for sorbed particle not to be sorbed
  // p = b*dt
  p2=sorparm[2]*partrace->RunTime.dt;
}


string ParticleClassTPLangmuirNE::SorptionType()
{
  return std::string("Langmuir non-equilibrium Sorption by transition probability");
}


//
// Sorption with retardation
//
ParticleClassR::ParticleClassR(PartraceClass* pt, long no, double m,
			       std::string &name, std::string &abbrev, 
			       int injno, InjectionClass *inj, double c0):
  ParticleClass(pt, no, m,name, abbrev, injno, inj, c0), sorparamno(2)
{
  partrace->SetNewHandlerText("retardation");
  if(partrace->distributedData)
    retardation=new DistData(partrace->parallel, "retardation", NoElxyz, 1);
  else
    retardation=new LocalData(partrace->parallel, "retardation", NoElxyz, 1);
  retardation->init(1.0);
}


ParticleClassR::~ParticleClassR()
{
  delete retardation;
}


string ParticleClassR::SorptionType()
{
  return std::string("Freundlich-Sorption by Retardation");
}


void ParticleClassR::InitSorption(string &fn, vector3& sorpparam, int predefined)
{ 
  partrace->SetNewHandlerText("adsorption");
  if(predefined==2) { // heterogeneous adsorption
    int ibuffer[4];
    partrace->elements->ReadField(fn,"Adsorption file ",
			   "#partrace adsorption file",ibuffer,adsorption);
  }
  else { // homogeneous adsorption
    adsorption=new ConstData(partrace->parallel, "adsorption", NoElxyz, sorparamno);
    adsorption->set_dataline(0, sorpparam.data());
  }
}


double ParticleClassR::Retardation(long ElementNo) const
{
  return retardation->get(ElementNo);
}


double ParticleClassR::Concentration(long ElementNo) const
{
  return concentration->get(ElementNo)*retardation->get(ElementNo);
}


double ParticleClassR::ConcentrationSorbed(long ElementNo) const
{
  return concentration->get(ElementNo)*(1.0-retardation->get(ElementNo));
}


void ParticleClassR::CalculateConcentration(double simtime)
{  
  CalculateConcentration2(concentration, ParticlesInElement);
  CalculateRetardation();
}


void ParticleClassR::CalculateRetardation()
{
  double c, cw, cn, r, k, n;
  const double epsilon=1.0e-5;
  const double null=1.0e-20;
  int j;
  double rdata;
  const double *sorparm;
  const double *soilprop;
  long to=concentration->end();
  long from=concentration->begin();
  double wc;

  // calculate retardation on local elements
  // retardation and concentration are equal distributed
  // c = cw + cs = cw + k*cw**n
  for(long e=from;e<to;e++) {
    c=concentration->get(e);
    rdata=retardation->get(e);
    if(c<null) rdata=1.0;
    else {
      if(rdata>null) cw=rdata*c;
      else cw=0.5*c;
      sorparm=adsorption->get_line(e);
      soilprop=partrace->elements->soilprop->get_line(e);
      wc=partrace->elements->watercontent->get(e);
      k=sorparm[0]*soilprop[4]/wc;
      n=sorparm[1];
      // Newton: f(cw)=cw+k*cw**n-c,  f'(cw)=1+n*k*cw**(n-1)
      for(j=0;j<50;j++) { // max. 50 Iterations
        cn=pow(cw,n);
        r=cw-(cw+k*cn-c)/(n*k*cn/cw+1.0);
        if(r<null) cw/=2.0;
        else {
          if(fabs((cw-r)/cw)<epsilon) {
            cw=r;
            break;
          }
          else cw=r;
        }
      }
      if(j==50) rdata=1.0;
      else rdata=cw/c;
      //      cout<<"k n c cw r "<<k<<' '<<n<<' '<<c<<' '<<cw<<' '<<cw/c<<endl;
    }
    retardation->set_data(e, 0, rdata);
  } // for
}


void ParticleClassTP::WriteMoments(double simtime)
{
  if(!SaveMomentsTime->mustHandle(simtime)) return;

  int i;
  long e;
  double mom[30];
  for(i=0;i<30;i++) mom[i]=0.0;
  // calculate moments
  // 0-9 in solution    10-19 sorbed   20-29 both
  // 0,10,20   = 0. moment = mass
  // 1-3,11-13 = 1. moment
  // 4-6,14,16 = 2. moment
  // 7-9,17-19 = 3. moment

  // calculating mass
  long number_w=0, number_s=0;
  double *p;
  Particle *particle1, *particle2;
  ParticleList* plist1=ParticlesInElement;
  ParticleList* plist2=ParticlesInElementSorbed;
  double mass=MassPerPart;
  for(e=0; e<NoElxyz; e++,plist1++,plist2++) {
    particle1=plist1->start();
    number_w+=plist1->count();
    while(particle1) {
      p=particle1->position.data();
      if(decay) {
	mass=particle1->get_mass();
	mom[0]+=mass;
      }
      mom[1]+=p[0]*mass;
      mom[2]+=p[1]*mass;
      mom[3]+=p[2]*mass;
      mom[4]+=p[0]*p[0]*mass;
      mom[5]+=p[1]*p[1]*mass;
      mom[6]+=p[2]*p[2]*mass;
      mom[7]+=p[0]*p[0]*p[0]*mass;
      mom[8]+=p[1]*p[1]*p[1]*mass;
      mom[9]+=p[2]*p[2]*p[2]*mass;
      particle1=plist1->next(particle1);
    }
    particle2=plist2->start();
    number_s+=plist2->count();
    while(particle2) {
      p=particle2->position.data();
      if(decay) {
	mass=particle2->get_mass();
	mom[10]+=mass;
      }
      mom[11]+=p[0]*mass;
      mom[12]+=p[1]*mass;
      mom[13]+=p[2]*mass;
      mom[14]+=p[0]*p[0]*mass;
      mom[15]+=p[1]*p[1]*mass;
      mom[16]+=p[2]*p[2]*mass;
      mom[17]+=p[0]*p[0]*p[0]*mass;
      mom[18]+=p[1]*p[1]*p[1]*mass;
      mom[19]+=p[2]*p[2]*p[2]*mass;
      particle2=plist2->next(particle2);
    }
  }
  if(!decay) {
    mom[ 0]=number_w*MassPerPart;
    mom[10]=number_s*MassPerPart;
  }
  partrace->parallel->sum( mom, 20);
  partrace->parallel->sum(&number_w, 1);
  partrace->parallel->sum(&number_s, 1);
  mom[20]=mom[0]+mom[10];
  if(ProcNo==0) {
    if(mom[20]>0.0)
      for(i=21;i<30;i++) mom[i]=(mom[i-20]+mom[i-10])/mom[20];
    if(number_w>0)
      for(i=1;i<10;i++) mom[i]/=mom[0];
    if(number_s>0)
      for(i=11;i<20;i++) mom[i]/=mom[10];
    double mom2z[3], mom3z[3];
    double mom2z_w[3], mom3z_w[3];
    double mom2z_s[3], mom3z_s[3];
    // mom2z(x) = var(x) = E( (x-E(x))**2 ) = E(x**2) - E(x)**2
    // mom3z(x) = E( (x-E(x))**3 ) = E(x**3) - 3*E(x**2)*E(x) + 2*E(x)**3
    for(i=0;i<3;i++) {
      mom2z[i]  =mom[i+24]-mom[i+21]*mom[i+21];
      mom3z[i]  =mom[i+27]-3.0*mom[i+24]*mom[i+21]+2.0*mom[i+21]*mom[i+21]*mom[i+21];  
      mom2z_w[i]=mom[i+ 4]-mom[i+ 1]*mom[i+ 1];
      mom3z_w[i]=mom[i+ 7]-3.0*mom[i+ 4]*mom[i+ 1]+2.0*mom[i+ 1]*mom[i+ 1]*mom[i+ 1];  
      mom2z_s[i]=mom[i+14]-mom[i+11]*mom[i+11];
      mom3z_s[i]=mom[i+17]-3.0*mom[i+14]*mom[i+11]+2.0*mom[i+11]*mom[i+11]*mom[i+11];  
    }

    ofstream out;
    if(SaveMomentsTime->get_counter()==0) {
      out.open(fnMom.c_str(), ios::trunc);
      out.close();
    }
    out.open(fnMom.c_str(), ios::app);
    if(out.fail()) partrace->error(string("Can not open moments output file: ") + fnMom);
    out<<"***************************************\n";
    out<<"Moments for "<<Name<<" after "<<simtime<<" days\n";
    out<<"Number of Particles : "<<setw(13)<<number_w+number_s<<'\n';
    out<<"Total Mass          : "<<setw(13)<<mom[20]<<'\n';
    out<<"Center of Mass      : "; for(i=21;i<24;i++) out<<setw(13)<<mom[i]; out<<'\n';
    out<<"Variance            : "; for(i=0;i<3;i++) out<<setw(13)<<mom2z[i]; out<<'\n';
    out<<"Third Moment        : "; for(i=0;i<3;i++) out<<setw(13)<<mom3z[i]; out<<'\n';
    out<<"In solution -------------------------------------------------\n";
    out<<"Number of Particles : "<<setw(13)<<number_w<<'\n';
    out<<"Total Mass          : "<<setw(13)<<mom[0]<<'\n';
    out<<"Center of Mass      : "; for(i=1;i<4;i++) out<<setw(13)<<mom[i];     out<<'\n';
    out<<"Variance            : "; for(i=0;i<3;i++) out<<setw(13)<<mom2z_w[i]; out<<'\n';
    out<<"Third Moment        : "; for(i=0;i<3;i++) out<<setw(13)<<mom3z_w[i]; out<<'\n';
    out<<"sorbed ------------------------------------------------------\n";
    out<<"Number of Particles : "<<setw(13)<<number_s<<'\n';
    out<<"Total Mass          : "<<setw(13)<<mom[10]<<'\n';
    out<<"Center of Mass      : "; for(i=11;i<14;i++) out<<setw(13)<<mom[i];   out<<'\n';
    out<<"Variance            : "; for(i=0;i<3;i++) out<<setw(13)<<mom2z_s[i]; out<<'\n';
    out<<"Third Moment        : "; for(i=0;i<3;i++) out<<setw(13)<<mom3z_s[i]; out<<'\n';
    out<<"Time used: "<<partrace->RunTime.Seconds()<<"secs"<<endl;
    out.close();
  }
  SaveMomentsTime->set_next_time(simtime);
}


void ParticleClassR::WriteMoments(double simtime)
{
  if(!SaveMomentsTime->mustHandle(simtime)) return;

  long i;
  double m, m_w, m_s;
  // calculating mass
  double mass;
  double mass_w;
  double mass_s;
  double vol;
  // mom1   0-2,   mom2  3-5,   mom3  6-8, 
  // mom1w  9-11, mom2w 12-14, mom3w 15-17
  // mom1s 18-20, mom2s 21-23, mom3s 24-26
  // mass_w 27,   mass_s 28
  double mom[29];
  double mom2z[3], mom2z_w[3], mom2z_s[3];
  double mom3z[3], mom3z_w[3], mom3z_s[3];
  for(i=0;i<29;i++) mom[i]=0.0;
  vector3 p;
  int number_w=0;    // no. of particles in water
  int number_s=0;    // no. of particles sorbed
  for(i=0;i<NoElxyz;i++) {
    vol=partrace->elements->volume->get(i) * partrace->elements->watercontent->get(i);
    if(decay) number_w+=ParticlesInElement[i].count();
    m_w=Concentration(i)*vol;
    m_s=ConcentrationSorbed(i)*vol;
    m=m_w+m_s;
    partrace->elements->GetElementsPosition(i, p);
    mom[0]+=p[0]*m;
    mom[1]+=p[1]*m;
    mom[2]+=p[2]*m;
    mom[3]+=p[0]*p[0]*m;
    mom[4]+=p[1]*p[1]*m;
    mom[5]+=p[2]*p[2]*m;
    mom[6]+=p[0]*p[0]*p[0]*m;
    mom[7]+=p[1]*p[1]*p[1]*m;
    mom[8]+=p[2]*p[2]*p[2]*m;
    mom[9]+=p[0]*m_w;
    mom[10]+=p[1]*m_w;
    mom[11]+=p[2]*m_w;
    mom[12]+=p[0]*p[0]*m_w;
    mom[13]+=p[1]*p[1]*m_w;
    mom[14]+=p[2]*p[2]*m_w;
    mom[15]+=p[0]*p[0]*p[0]*m_w;
    mom[16]+=p[1]*p[1]*p[1]*m_w;
    mom[17]+=p[2]*p[2]*p[2]*m_w;
    mom[18]+=p[0]*m_s;
    mom[19]+=p[1]*m_s;
    mom[20]+=p[2]*m_s;
    mom[21]+=p[0]*p[0]*m_s;
    mom[22]+=p[1]*p[1]*m_s;
    mom[23]+=p[2]*p[2]*m_s;
    mom[24]+=p[0]*p[0]*p[0]*m_s;
    mom[25]+=p[1]*p[1]*p[1]*m_s;
    mom[26]+=p[2]*p[2]*p[2]*m_s;
    mom[27]+=m_w;
    mom[28]+=m_s;
  }
  mass_w=mom[27];
  mass_s=mom[28];
  mass=mass_w+mass_s;
  number_w=(long) (mass_w/MassPerPart+0.5);    // no. of particles in water
  number_s=(long) (mass_s/MassPerPart+0.5);    // no. of particles sorbed
  int number=number_w+number_s;                // no. of particles
  if(mass!=0.0) {
    for(i=0;i<9;i++) mom[i]/=mass;
    if(mass_w>0.0) for(i=9;i<18;i++) mom[i]/=mass_w;
    if(mass_s>0.0) for(i=18;i<27;i++) mom[i]/=mass_s;
  }
  // mom2z(x) = var(x) = E( (x-E(x))**2 ) = E(x**2) - E(x)**2
  for(i=0;i<3;i++) {
    mom2z[i]  =mom[i+3 ]-mom[i   ]*mom[i];
    mom2z_w[i]=mom[i+12]-mom[i+ 9]*mom[i+9];
    mom2z_s[i]=mom[i+21]-mom[i+18]*mom[i+18];
  }
  // mom3z(x) = E( (x-E(x))**3 ) = E(x**3) - 3*E(x**2)*E(x) + 2*E(x)**3
  for(i=0;i<3;i++) {
    mom3z[i]  =mom[i+6] -3.0*mom[i+ 3]*mom[i   ]+2.0*mom[i   ]*mom[i   ]*mom[i];
    mom3z_w[i]=mom[i+15]-3.0*mom[i+12]*mom[i+9 ]+2.0*mom[i+ 9]*mom[i+ 9]*mom[i+9];
    mom3z_s[i]=mom[i+24]-3.0*mom[i+21]*mom[i+18]+2.0*mom[i+18]*mom[i+18]*mom[i+18];
  }
  if(ProcNo==0) {
     ofstream out;
     if(SaveMomentsTime->get_counter()==0) out.open(fnMom.c_str(), ios::trunc);
     else out.open(fnMom.c_str(), ios::app);
     if(out.fail()) partrace->error(string("Can not open moments output file: ") + fnMom);
     out<<"*************************************************************\n";
     out<<"Moments for "<<Name<<" after "<<simtime<<" days\n";
     out<<"Number of Particles : "<<setw(13)<<number<<'\n';
     out<<"Total Mass          : "<<setw(13)<<mass<<'\n';
     out<<"Center of Mass      : "; for(i=0;i<3;i++) out<<setw(13)<<mom[i];  out<<'\n';
     out<<"Variance            : "; for(i=0;i<3;i++) out<<setw(13)<<mom2z[i]; out<<'\n';
     out<<"Third Moment        : "; for(i=0;i<3;i++) out<<setw(13)<<mom3z[i]; out<<'\n';
     out<<"In solution -------------------------------------------------\n";
     out<<"Number of Particles : "<<setw(13)<<number_w<<'\n';
     out<<"Total Mass          : "<<setw(13)<<mass_w<<'\n';
     out<<"Center of Mass      : "; for(i=9;i<12;i++) out<<setw(13)<<mom[i];  out<<'\n';
     out<<"Variance            : "; for(i=0;i<3;i++) out<<setw(13)<<mom2z_w[i]; out<<'\n';
     out<<"Third Moment        : "; for(i=0;i<3;i++) out<<setw(13)<<mom3z_w[i]; out<<'\n';
     out<<"sorbed ------------------------------------------------------\n";
     out<<"Number of Particles : "<<setw(13)<<number_s<<'\n';
     out<<"Total Mass          : "<<setw(13)<<mass_s<<'\n';
     out<<"Center of Mass      : "; for(i=18;i<21;i++) out<<setw(13)<<mom[i];  out<<'\n';
     out<<"Variance            : "; for(i=0;i<3;i++) out<<setw(13)<<mom2z_s[i]; out<<'\n';
     out<<"Third Moment        : "; for(i=0;i<3;i++) out<<setw(13)<<mom3z_s[i]; out<<'\n';
     out<<"Time used: "<<partrace->RunTime.Seconds()<<"secs"<<endl;
     out.close();
  }
  SaveMomentsTime->set_next_time(simtime);
}


//
// freundlich sorption with interactions
//
ParticleClassInterF::ParticleClassInterF(PartraceClass* pt, long no, double m,
					 std::string &name, std::string &abbrev, 
					 int injno, InjectionClass *inj, double c0):
  ParticleClassTP(pt, no, m, name, abbrev, injno, inj, c0)
{
  sorparamno=2;
}


string ParticleClassInterF::SorptionType()
{
  return std::string("Freundlich-Sorption with interaction");
}
