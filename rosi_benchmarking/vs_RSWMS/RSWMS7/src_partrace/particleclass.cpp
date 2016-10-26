#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include "particleclass.hpp"
#include "elements.hpp"
#include "random.hpp"
#include "event.hpp"
#include "partrace.hpp"
#include "parallel.hpp"
#include "runtime.hpp"
#include "distdata.hpp"
#include "decay.hpp"
#include "uptake.hpp"

using namespace std;


InjectionClass::InjectionClass(): startEvent(0), endEvent(0)
{
}


InjectionClass::~InjectionClass()
{
  delete startEvent;
  delete endEvent;
}


void InjectionClass::init(RunTimeClass& RunTime)
{
  startEvent=new Event("injection start", Start);
  endEvent=new Event("injection end", End);
  RunTime.insert(startEvent);
  RunTime.insert(endEvent);
  ProcLast=0;
}


// Constructor
ParticleClass::ParticleClass(PartraceClass* pt, long no, double m, string &name,
			     string &abbrev, int injno, InjectionClass *inj, double c0): 
               MaxNo(no), partrace(pt), next(0),
	       injectionNo(injno), injection(inj),
	       MassPerPart(m), C0(c0), BTCSeqout(0),
               Decay(0), Uptake(0), adsorption(0), SaveConcTime(0), 
	       SaveMomentsTime(0), SaveParticlesTime(0)
{
  ProcMax=partrace->parallel->cpus();
  ProcNo=partrace->parallel->mycpu();
  AlwaysCalculateConcentration=true;
  if(C0==0.0) C0=1.0;
  for(int i=0;i<injectionNo;i++) injection[i].init(partrace->RunTime);
  NoElxyz=partrace->elements->TotalNumber();
  partrace->SetNewHandlerText("concentration");
  if(partrace->distributedData)
    concentration=new DistData(partrace->parallel, "concentration", NoElxyz, 1);
  else
    concentration=new LocalData(partrace->parallel, "concentration", NoElxyz, 1);
  Name=name; Abbrev=abbrev;
  fnBTC=partrace->elements->fnBTC+string(".")+Abbrev;
  fnBTCSeq=partrace->elements->fnBTC+string(".seq");
  fnConc=partrace->elements->fnConc+string(".")+Abbrev;
  fnMom=partrace->elements->fnMom+string(".")+Abbrev;
  if(partrace->elements->fnParticles.length()>0) {
    fnParticles=partrace->elements->fnParticles+string(".")+Abbrev;
    partrace->SetFNSequence(fnParticles, ProcNo, 5);
  }
  // allocate and initialize particle lists
  partrace->SetNewHandlerText("ParticlesInElement");
  ParticlesInElement=new ParticleList[NoElxyz];
  partrace->SetNewHandlerText("MovedParticles");
  MovedParticles=new ParticleSimpleList[NoElxyz];
  ParticlesInElementSorbed=0;
}

// Destructor
ParticleClass::~ParticleClass()
{
  delete [] ParticlesInElement;
  delete concentration;
  delete adsorption;
  delete SaveConcTime;
  delete SaveMomentsTime;
  delete SaveParticlesTime;
  delete Decay;
  delete Uptake;
  delete [] BTCSeqout;
  delete [] injection;
  if(ProcNo==0) BTCout.close();
}


Particle* ParticleClass::new_particle()
{
  if(decay || uptake) return new ParticleWithMass(MassPerPart);
  else return new Particle;
}


Particle* ParticleClass::new_particle(vector3 &pos)
{
  if(decay || uptake)  return new ParticleWithMass(pos, MassPerPart);
  else return new Particle(pos);
}


std::string ParticleClass::SorptionType()
{
  return std::string("konservative");
}


std::string ParticleClass::DecayType()
{
  return Decay->DecayType();
}

std::string ParticleClass::UptakeType()
{
  return Uptake->UptakeType();
}

double ParticleClass::Concentration(long ElementNo) const
{
  return concentration->get(ElementNo);
}


void ParticleClass::Distribute(int& pno, int& ProcLast)
{
  if(ProcMax>1) {
    int ino=pno;
    pno=ino/ProcMax;
    int rest=ino%ProcMax;
    if((ProcNo>=ProcLast && ProcNo<ProcLast+rest) || 
       (ProcNo<ProcLast+rest-ProcMax)) pno++;
    ProcLast+=rest;
    if(ProcLast>=ProcMax) ProcLast-=ProcMax;
  }
}


bool ParticleClass::InjectParticlePart()
{
  bool injecting=false;
  int pno;
  for(int i=0;i<injectionNo;i++) {
    if(injection[i].ready()) continue;
    if(injection[i].startEvent->mustHandle(partrace->RunTime.oldTime)) {
      if(partrace->RunTime.Time>=injection[i].End || injection[i].Start==injection[i].End) {
	pno=injection[i].ParticlesToInject;
	injection[i].startEvent->deactivate();
	injection[i].endEvent->deactivate();
      }
      else
	pno=int(partrace->RunTime.dt/(injection[i].End-partrace->RunTime.oldTime)
		*injection[i].ParticlesToInject+0.5);
      injection[i].ParticlesToInject-=pno;
      // particles per injection and processor
      partrace->SetNewHandlerText("ParticlesInElement");
      Distribute(pno, injection[i].ProcLast);
      switch(injection[i].StartModel) {
      case 1: StartDirac(pno, injection[i]); break;
      case 2: StartCube(pno, injection[i]); break;
      case 3: injection[i].StartRadius[2]=injection[i].StartRadius[1]=injection[i].StartRadius[0];
	StartEllipse(pno, injection[i]);
	break;
      case 4: StartEllipse(pno, injection[i]); break;
      case 5: 
      case 6: setInitialConcentration(pno, injection[i]); break;
      default: StartDirac(pno, injection[i]); break;
      }
      cout<<"Processor "<<ProcNo<<" Injecting "<<pno<<" particles of type: "
	  <<SorptionType()<<" at time "<<partrace->RunTime.Time<<endl;
      //     cout<<"StartPos   : "<<StartPos[0]<<' '<<StartPos[1]<<' '<<StartPos[2]<<' '<<endl;
      //     cout<<"StartRadius: "<<StartRadius[0]<<' '<<StartRadius[1]<<' '<<StartRadius[2]<<' '<<endl;
      injecting=true;
    }
  }
  return injecting;
}


void ParticleClass::InitDecay(string& DecayField, int order, int useField, double rate)
{
  decay=1;
  switch(order) {
  case 1: Decay=new DecayClass(partrace, this);
    decay=false;
    break;
  case 2: Decay=new DecayClass1(partrace, this, rate, DecayField, useField);
    decay=true;
    break;
  case 3: Decay=new DecayClass2(partrace, this, rate, DecayField, useField);
    decay=true;
    break;
  default: partrace->error("Invalid microbial decay type");
    break;
  }
}

void ParticleClass::InitUptake(string& UptakeField, int order, int useField, double rate)
{
  uptake=true;
  switch(order) {
  case 1: Uptake=new UptakeClass(partrace, this);
    uptake=false;
    break;
  case 2: Uptake=new UptakeClass1(partrace, this, rate, UptakeField, useField);
    uptake=true;
    break;
  case 3: Uptake=new UptakeClass2(partrace, this, rate, UptakeField, useField);
    uptake=true;
    break;
  case 4: Uptake=new UptakeClass3(partrace, this, rate, UptakeField, useField);
    uptake=true;
    break;
  default: partrace->error("Invalid uptake type");
    break;
  }
}


void ParticleClass::MoveParticles()
{
  long e, ex,ey,ez, ele;
  Particle* particle;
  Particle *p, *prev;
  DispersionClass dispersion(this);


  // move all particles
  ParticleList* plist=ParticlesInElement;
  ElementsClass *elements=partrace->elements;
  for(e=0; e<NoElxyz; e++,plist++) {
    particle=plist->start(prev);
    if(particle==0) continue; // no particles in element
    elements->GetElementXYZ(e,ex,ey,ez);
    //out <<  'e: ' << (int)e << endl;//<< 'NoElx: ' << NoElx << endl;
    //sleep(2);
    while(particle) { // all particles in element e
      ele=elements->MoveParticle(e,ex,ey,ez,particle->position,dispersion);
      if( ele == e ) // particle did not leave element
	particle=plist->next(particle, prev);
      else if(ele<0) // particle has left the simulation area
 	particle=plist->remove(particle,prev); // particle outside
      else { // particle has changed element
	p=particle;
	particle=plist->remove_from_list(particle,prev);
	MovedParticles[ele].insert(p);
      }
    }
  }
  // delete local buffer
  //elements->soilprop->delete_buffer(dispersion.parameter);
  // merge particles in MovedParticles into ParticlesInElement
  for(e=0; e<NoElxyz; e++) {
    ParticlesInElement[e].insert_list(MovedParticles[e].start());
    MovedParticles[e].reset();
  }
}


void ParticleClass::CalculateConcentrationInit()
{
  // calculation of concentration for conservative particles
  // is only nessary when output is required. (btc)
  if(partrace->elements->BTCAnz==0 && partrace->elements->ElementSequenceNo==0 && decay==0 && uptake==0)
    AlwaysCalculateConcentration=false;
}


void ParticleClass::CalculateConcentration(double simtime)
{
  if(AlwaysCalculateConcentration ||
     SaveConcTime->mustHandle(simtime) ||
     SaveMomentsTime->mustHandle(simtime) ||
     SaveParticlesTime->mustHandle(simtime)) 
    CalculateConcentration2(concentration, ParticlesInElement);
}


void ParticleClass:: CalculateConcentration2(ConstData* conc, ParticleList* plist)
{
  partrace->debug<<"CalculateConcentration2    "<<partrace->RunTime.Seconds()<<"secs"<<endl;
  Particle* particle;
  double mass;
  int bufsize=partrace->elements->calculateBufferSize();
  double *buffer=conc->start_buffered_sum(bufsize);
  for(int e=0; e<NoElxyz; e++,plist++) {
    if(decay || uptake) {
      mass=0.0;
      particle=plist->start();
      while(particle) {
	mass+=particle->get_mass();
	particle=plist->next(particle);
      } // end particle loop
      *buffer=mass;
    }
    else
      *buffer=MassPerPart * plist->count();
    buffer=conc->buffered_sum(e);
  }
  conc->end_buffered_sum();
  //cout<<"time="<<partrace->RunTime.Time<<endl<<*conc<< partrace->elements->watercontent<<endl;
  conc->divide(partrace->elements->volume, partrace->elements->watercontent);
}


 void ParticleClass::SaveConcentration(double simtime)
{
  if(!SaveConcTime->mustHandle(simtime)) return;
  // generate output filename
  ostringstream sout;
  string fn;
  int output_type;
  if(partrace->PrintMass) output_type=2;
  else output_type=1;
  // write concentration to file
  sout.width(6); sout.fill('0');
  sout<<SaveConcTime->get_counter()<<ends;
  fn=fnConc + sout.str();
  ofstream out(fn.c_str());
  if(out.fail()) partrace->error(string("Can not open concentration output file: ") + fn);
  out<<"# "<<Name<<'\n';
  out<<"# "<<simtime
     <<' '<<partrace->elements->NoElx
     <<' '<<partrace->elements->NoEly
     <<' '<<partrace->elements->NoElz
     <<' '<<partrace->elements->GeometryType
     <<' '<<output_type
     <<'\n';
  out<<"#";
  if(partrace->elements->GeometryType==1)
    out<<' '<<partrace->elements->length_x
       <<' '<<partrace->elements->length_y
       <<' '<<partrace->elements->length_z;
  else 
    out<<' '<<partrace->elements->fnCoordinates;
  out<<'\n';
  if(partrace->PrintMass) // mass output
    for(long i=0;i<NoElxyz;i++){
      //cout << Decay->decayfield->get(i)<<endl;
      out<<Concentration(i)/C0<<' '<<ConcentrationSorbed(i)/C0<<' '
	 <<partrace->elements->volume->get(i)<<' '<<partrace->elements->watercontent->get(i)<<'\n';
     }
  else // no mass output
    for(long i=0;i<NoElxyz;i++)
      out<<Concentration(i)/C0<<' '<<ConcentrationSorbed(i)/C0<<'\n';
  out.close();
  SaveConcTime->set_next_time(simtime);
}


void ParticleClass::SaveParticlePositions(double simtime)
{
  if(!SaveParticlesTime->mustHandle(simtime)) return;
  ofstream out;
  if(SaveParticlesTime->get_counter()==0) {
     out.open(fnParticles.c_str(), ios::trunc);
     out.close();
  }
  out.open(fnParticles.c_str(), ios::app);
  if(out.fail()) partrace->error(string("Can not open the output file for the particle positions: ") + fnParticles);
  out<<"#time: "<<simtime<<'\n';

  Particle* particle;
  ParticleList* plist=ParticlesInElement;
  for(long e=0; e<NoElxyz; e++,plist++) {
    particle=plist->start();
    while(particle) {
      out<<particle->position[0]<<' '<<particle->position[1]<<' '<<particle->position[2]
	 <<'\n';
      particle=plist->next(particle);
    } // end particle loop
  }
  out.close();
  SaveParticlesTime->set_next_time(simtime);
}


void ParticleClass::OpenBTC(double simtime)
{
  // open single BTC file
  if(simtime==0.0) BTCout.open(fnBTC.c_str());
  else BTCout.open(fnBTC.c_str(), ios::app);
  if(BTCout.fail()) partrace->error(string("Can not open BTC-file: ")+fnBTC);
  // open sequences of BTCs, one file per sequence
  int seqno=partrace->elements->ElementSequenceNo;
  if(seqno>0) {
    partrace->SetNewHandlerText("BTCSeqout");
    BTCSeqout=new ofstream[seqno];
    string fn;
    ostringstream sout;
    int len=int(ceil(log10(double(seqno)))+0.1);
    for(int i=0;i<seqno;i++) {
      sout.str(fn);
      sout<<fnBTCSeq;
      sout.width(len); sout.fill('0');
      sout<<i;
      fn=sout.str();
      if(simtime==0.0) BTCSeqout[i].open(fn.c_str());
      else BTCSeqout[i].open(fn.c_str(), ios::app);
      if(BTCSeqout[i].fail()) partrace->error(string("Can not open BTCSeq-file: ")+fn);
    }
  }
}


void ParticleClass::WriteBTC(double simtime)
{
  long elem;
  ElementsClass *elements=partrace->elements;
  // write single BTC file
  BTCout<<simtime;
  for(int i=0;i<elements->BTCAnz;i++) {
    elem=elements->BTCListElNo[i];
    BTCout<<' '<<Concentration(elem)/C0;
  }
  BTCout<<endl;
  // write sequence of BTC files
  int seqno=elements->ElementSequenceNo;
  if(seqno>0) {
    int e, ex,ey,ez, ex_to,ey_to,ez_to, ey_old,ez_old;
    ElementSequenceClass *sequence;
    for(int i=0;i<seqno;i++) {
      BTCSeqout[i]<<simtime;
      sequence=elements->ElementSequence+i;
      e=sequence->startindex;
      ex_to=sequence->ele_x;
      ey_to=sequence->ele_y;
      ez_to=sequence->ele_z;
      for(ez=0;ez<ez_to;ez++) {
	ez_old=e;
	for(ey=0;ey<ey_to;ey++) {
	  ey_old=e;
	  for(ex=0; ex<ex_to; ex++,e++)
	    BTCSeqout[i]<<' '<<Concentration(e)/C0;
	  e=ey_old+elements->NoElx;
	}
	e=ez_old+elements->NoElxy;
      }
      BTCSeqout[i]<<endl;
    }
  }
}


void ParticleClass::WriteMoments(double simtime)
{
  if(!SaveMomentsTime->mustHandle(simtime)) return;

  // calculate moments
  // 1.index=mom 1,2,3  2.index=x/y/z
  double mom[3][3];
  mom[0][0]=mom[0][1]=mom[0][2]=0.0;
  mom[1][0]=mom[1][1]=mom[1][2]=0.0;
  mom[2][0]=mom[2][1]=mom[2][2]=0.0;
  double mom0=0.0;
  double *p;
  double mass;
  long number=0;
  Particle* particle;
  ParticleList* plist=ParticlesInElement;
  for(long e=0; e<NoElxyz; e++,plist++) {
    number+=plist->count(); // calculating mass
    particle=plist->start();
    while(particle) {
      p=particle->position.data();
      if(decay || uptake) {
	mass=particle->get_mass();
	mom0+=mass;
	mom[0][0]+=p[0]*mass;
	mom[0][1]+=p[1]*mass;
	mom[0][2]+=p[2]*mass;
	mom[1][0]+=p[0]*p[0]*mass;
	mom[1][1]+=p[1]*p[1]*mass;
	mom[1][2]+=p[2]*p[2]*mass;
	mom[2][0]+=p[0]*p[0]*p[0]*mass;
	mom[2][1]+=p[1]*p[1]*p[1]*mass;
	mom[2][2]+=p[2]*p[2]*p[2]*mass;
      }
      else {
	mom[0][0]+=p[0];
	mom[0][1]+=p[1];
	mom[0][2]+=p[2];
	mom[1][0]+=p[0]*p[0];
	mom[1][1]+=p[1]*p[1];
	mom[1][2]+=p[2]*p[2];
	mom[2][0]+=p[0]*p[0]*p[0];
	mom[2][1]+=p[1]*p[1]*p[1];
	mom[2][2]+=p[2]*p[2]*p[2];
      }
      particle=plist->next(particle);
    }
  }

  partrace->parallel->sum(&number, 1);
  if(decay || uptake) {
    partrace->parallel->sum( mom0 );
  }
  else mom0=number*MassPerPart;
  partrace->parallel->sum((double *) mom, 9);
  if(ProcNo==0) {
    int i;
    if(decay || uptake) mass=1.0/mom0;
    else mass=MassPerPart/mom0;
    if(number>0) for(i=0;i<3;i++) {
      mom[i][0]*=mass;
      mom[i][1]*=mass;
      mom[i][2]*=mass;
    }
    double mom2z[3], mom3z[3];
    // mom2z(x) = var(x) = E( (x-E(x))**2 ) = E(x**2) - E(x)**2
    // mom3z(x) = E( (x-E(x))**3 ) = E(x**3) - 3*E(x**2)*E(x) + 2*E(x)**3
    for(i=0;i<3;i++) {
      mom2z[i]=mom[1][i]-mom[0][i]*mom[0][i];
      mom3z[i]=mom[2][i]-3.0*mom[1][i]*mom[0][i]+2.0*mom[0][i]*mom[0][i]*mom[0][i];  
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
    out<<"Number of Particles : "<<setw(13)<<number<<'\n';
    out<<"Total Mass          : "<<setw(13)<<mom0<<'\n';
    out<<"Center of Mass      : "; for(i=0;i<3;i++) out<<setw(13)<<mom[0][i];  out<<'\n';
    out<<"Variance            : "; for(i=0;i<3;i++) out<<setw(13)<<mom2z[i]; out<<'\n';
    out<<"Third Moment        : "; for(i=0;i<3;i++) out<<setw(13)<<mom3z[i]; out<<'\n';
    out<<"***************************************\n";
    out<<"Time used: "<<partrace->RunTime.Seconds()<<"secs"<<endl;
    out.close();
  }
  SaveMomentsTime->set_next_time(simtime);
}


void ParticleClass::StartDirac(ulong pno, InjectionClass &inj)
{
  long ElementNo=partrace->elements->GetElementNo(inj.StartPos);
  if(ElementNo<0) 
     partrace->error("Particle(s) in start-position outside the volume");
  for(ulong i=0;i<pno;i++)
    ParticlesInElement[ElementNo].insert(new_particle(inj.StartPos));
}


void ParticleClass::StartCube(ulong pno, InjectionClass &inj)   
{
  double xvon=inj.StartPos[0]-inj.StartRadius[0];
  double xbis=inj.StartPos[0]+inj.StartRadius[0];
  double yvon=inj.StartPos[1]-inj.StartRadius[1];
  double ybis=inj.StartPos[1]+inj.StartRadius[1];
  double zvon=inj.StartPos[2]-inj.StartRadius[2];
  double zbis=inj.StartPos[2]+inj.StartRadius[2];
  Particle* p;
  long e;
  RandomClass *rand=partrace->Random;
  for(ulong i=0;i<pno;i++) {
    p=new_particle();
    if(xvon>=xbis) p->position[0]=xvon;
    else p->position[0]=rand->draw(xvon,xbis);
    if(yvon>=ybis) p->position[1]=yvon;
    else p->position[1]=rand->draw(yvon,ybis);
    if(zvon>=zbis) p->position[2]=zvon;
    else p->position[2]=rand->draw(zvon,zbis);
    if( (e=partrace->elements->GetElementNo(p->position)) < 0 )
       partrace->error("Particle(s) in start-position outside the volume");
    ParticlesInElement[e].insert(p);
  }
}


void ParticleClass::StartEllipse(ulong pno, InjectionClass &inj)
{
  double xvon, yvon, zvon;
  Particle* p;
  long e;
  RandomClass *rand=partrace->Random;
  for(ulong i=0;i<pno;i++) {
    do {
      xvon=rand->draw(-1.0,1.0);
      yvon=rand->draw(-1.0,1.0);
      zvon=rand->draw(-1.0,1.0);
    } while (xvon*xvon+yvon*yvon+zvon*zvon>1.0);
    p=new_particle();
    p->position[0]=inj.StartPos[0]+xvon*inj.StartRadius[0];
    p->position[1]=inj.StartPos[1]+yvon*inj.StartRadius[1];
    p->position[2]=inj.StartPos[2]+zvon*inj.StartRadius[2];
    if( (e=partrace->elements->GetElementNo(p->position)) < 0 )
       partrace->error("Particle(s) in start-position outside the volume");
    ParticlesInElement[e].insert(p);
  }
}


void ParticleClass::setInitialConcentration(int &pno, InjectionClass &inj)
{
  inj.ParticlesToInject=0; // inject only once.
  ParticleList* plist=ParticlesInElement;
  double conc=inj.InitConcentration;
  double conc_sorbed;
  double mass=0.0;
  double mass_sorbed=0.0;
  ConstData *volume=partrace->elements->volume;
  ConstData *watercontent=partrace->elements->watercontent;
  int n, lastpe=0;
  ifstream inp;
  const char nl='\n';
  const int linesize=1024;
  int geom, vol;

  if(inj.StartModel==6) {
    string line;
    double t;
    triple nel;
    inp.open(RestartFromConcFile.c_str());
    if(inp.fail())
      partrace->error(string("Can not open concentration input file: ")
		      + RestartFromConcFile);
    getline(inp, line);      // ignore line 1
    inp.ignore(1);           // ignore character #
    inp>>t>>nel[0]>>nel[1]>>nel[2]>>geom>>vol;
    inp.ignore(linesize,nl); // ignore rest of line 2
    if(partrace->elements->CheckNEL(nel))
      partrace->error(string("Wrong number of elements in concentration input file: ")
		      + RestartFromConcFile);
    inp.ignore(linesize,nl); // ignore line 3
  }

  pno=0;
  for(long e=0; e<NoElxyz; e++,plist++) {
    // c=m/V <=> m=c*V
    if(inj.StartModel==6) {
      inp>>conc>>conc_sorbed;
      inp.ignore(linesize,nl); // ignore rest of line
    }
    mass+=conc*volume->get(e)*watercontent->get(e);
    n=int(mass/MassPerPart+0.5);
    if(n>0) {
      mass-=double(n)*MassPerPart;
      pno+=n;
      Distribute(n, lastpe);
      partrace->elements->injectInElement(e, n, this, ParticlesInElement);
    }
    if(ParticlesInElementSorbed && inj.StartModel==6) {
      mass_sorbed+=conc_sorbed*volume->get(e)*watercontent->get(e);
      n=int(mass_sorbed/MassPerPart+0.5);
      if(n>0) {
	mass_sorbed-=double(n)*MassPerPart;
	pno+=n;
	Distribute(n, lastpe);
	partrace->elements->injectInElement(e, n, this, ParticlesInElementSorbed);
      }
    }
  }
  if(inj.StartModel==6) inp.close();
}


void ParticleClass::SaveParticles(fstream& res)
{
  if(partrace->elements->RestartType != 2) {
		res.write((char *) &injectionNo, sizeof(int));
  for(int i=0;i<injectionNo;i++)
    res.write((char *) &(injection[i].ParticlesToInject), sizeof(ulong));
	}
  SaveParticlesList(res, ParticlesInElement);
  if(ParticlesInElementSorbed) SaveParticlesList(res, ParticlesInElementSorbed);
}


void ParticleClass::SaveParticlesList(fstream& res, ParticleList* plist)
{
  Particle* particle;
  ulong no;
  for(long e=0; e<NoElxyz; e++,plist++) {
    no=plist->count();
    res.write((char *) &no, sizeof(no));
    particle=plist->start();
    while(particle) {
      particle->write(res);
      particle=plist->next(particle);
    }
  }
}


void ParticleClass::LoadParticles(fstream& res)
{
  if(res.eof()) return;
  if(partrace->elements->RestartType != 2) {
  int injno;
  res.read((char *) &injno, sizeof(int));
  if(injno!=injectionNo) partrace->error("Invalid number of injections in restart file.");
  for(int i=0;i<injectionNo;i++)
    res.read((char *) &(injection[i].ParticlesToInject), sizeof(ulong));
	}
  LoadParticlesList(res, ParticlesInElement);
  if(ParticlesInElementSorbed) LoadParticlesList(res, ParticlesInElementSorbed);
}


void ParticleClass::LoadParticlesList(fstream& res, ParticleList* plist)
{
  if(res.eof()) return;
  ulong no, n;
  Particle *particle;
  for(long e=0; e<NoElxyz; e++,plist++) {
    res.read((char *) &no, sizeof(no));
    for(n=0;n<no;n++) {
      particle=new_particle();
      particle->read(res);
      plist->insert(particle);
    }
  }
}
