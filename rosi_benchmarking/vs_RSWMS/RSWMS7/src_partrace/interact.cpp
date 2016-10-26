#include <cmath>
#include <iostream>
#include "interact.hpp"
#include "elements.hpp"
#include "random.hpp"
#include "partrace.hpp"
#include "parallel.hpp"
#include "distdata.hpp"

// pointer to first/last interaction
InteractClass* InteractClass::top=0;
InteractClass* InteractClass::bottom=0; 


InteractClass::InteractClass(std::string &name, ParticleClass* target,
			     ParticleClass* part1, ParticleClass* part2,
			     double* p, PartraceClass *pt):
  TargetParticle(target), Particle1(part1), Particle2(part2), resultsize(8),
  partrace(pt)
{
  InteractName=name;
  // add interaction at the bottom of the internal list
  if(bottom==0) top=bottom=this;
  else {
    bottom->next=this;
    bottom=this;
  }
  next=0;
  // parameters
  param[0]=p[0]; // model: Cab = p[0] * Ca * Cb
}


void InteractClass::CalculateAllInteractions()
{
  InteractClass* p=TopOfList();
  while(p) {
     p->CalculateInteraction();
     p=p->NextInList();
  }
}


// Model:
// Sa  = ka  *  Ca**na   <==>  NoSa*PPa   = ka  * (NoCa*PPa)**na
// Sb  = kb  *  Cb**nb   <==>  NoSb*PPb   = kb  * (NoCb*PPb)**nb
// Sab = kab * Cab**nab  <==>  NoSab*PPab = kab * (NoCab*PPab)**nab
//
// Noa  = NoSa + NoCa + NoSab + NoCab
// Nob  = NoSb + NoCb + NoSab + NoCab
//
// NoCab = G * NoCa * NoCb
//
//      6 equations with 6 unknows NoCa, NoCb, NoCab, NoSa, NoSb, NoSab
//
// ==>
// Noa  = ka/PPa * (NoCa*PPa)**na + NoCa + kab/PPab * (NoCab*PPab)**nab + NoCab
// Nob  = kb/PPb * (NoCb*PPb)**nb + NoCb + kab/PPab * (NoCab*PPab)**nab + NoCab
// ==>
// Noa  = ka/PPa * (NoCa*PPa)**na + NoCa + kab/PPab * ((G*NoCa*NoCb)*PPab)**nab + (G*NoCa*NoCb)
// Nob  = kb/PPb * (NoCb*PPb)**nb + NoCb + kab/PPab * ((G*NoCa*NoCb)*PPab)**nab + (G*NoCa*NoCb)
//
//      2 equations with 2 unknows NoCa and NoCb
//
// Noa  = ka/PPa*PPa**na *NoCa**na + NoCa + kab/PPab * (G*PPab)**nab *(NoCa*NoCb)**nab + G*NoCa*NoCb
// Nob  = kb/PPb*PPb**nb *NoCb**nb + NoCb + kab/PPab * (G*PPab)**nab *(NoCa*NoCb)**nab + G*NoCa*NoCb
//
// Newton: x2 = x1 - (df(x1)/dx1)**-1 * f(x1)
//
//  | a11   a12 |-1    |  a22  -a12 |
//  | a21   a22 |   =  | -a21   a11 | / det(a)
//
// old: *********************************************************************
//
//        a1=ka/PPa*PPa**na, a2=kab/PPab * (G*PPab)**nab
//        a4=kb/PPb*PPb**nb
// ==>
// Noa  = a1*NoCa**na + NoCa + a2*(NoCa*NoCb)**nab + G*NoCa*NoCb
// Nob  = a4*NoCb**nb + NoCb + a2*(NoCa*NoCb)**nab + G*NoCa*NoCb
//
// f1   = a1*NoCa**na + NoCa + a2*(NoCa*NoCb)**nab + G*NoCa*NoCb - Noa
// f2   = a4*NoCb**nb + NoCb + a2*(NoCa*NoCb)**nab + G*NoCa*NoCb - Nob
//
//  df1/dNoCa = na*a1*NoCa**(na-1) +1 + nab*a2*NoCb*(NoCa*NoCb)**(nab-1) + G*NoCb = a11
//  df2/dNoCa =                         nab*a2*NoCb*(NoCa*NoCb)**(nab-1) + G*NoCb = a21
//  df1/dNoCb =                         nab*a2*NoCa*(NoCa*NoCb)**(nab-1) + G*NoCa = a12
//  df2/dNoCb = nb*a4*NoCb**(nb-1) +1 + nab*a2*NoCa*(NoCa*NoCb)**(nab-1) + G*NoCa = a22
//
// new: *********************************************************************
//
//      x=log(NoCa),  y=log(NoCb)  <==> exp(x)=NoCa,  exp(y)=NoCb
// ==>  Noa = ka*PPa**(na-1)*exp(na*x) + exp(x) + kab/PPab * (G*PPab)**nab * exp(nab*(x+y)) + G*exp(x+y)
// ==>  Nob = kb*PPb**(nb-1)*exp(nb*y) + exp(y) + kab/PPab * (G*PPab)**nab * exp(nab*(x+y)) + G*exp(x+y)
//
//      A = ka*PPa**(na-1),   B = kb*PPb**(nb-1)
//      C = G                 D = kab * (G*PPab)**(nab-1)
//
// ==>  f1 =A*exp(na*x) + exp(x) + D*exp(nab*(x+y)) + C*exp(x+y) - Noa
//      f2 =B*exp(nb*y) + exp(y) + D*exp(nab*(x+y)) + C*exp(x+y) - Nob
//
// ==>  df1/dx = na*A*exp(na*x) + exp(x) + nab*D*exp(nab*(x+y)) + C*exp(x+y)
// ==>  df2/dx =                           nab*D*exp(nab*(x+y)) + C*exp(x+y)
// ==>  df1/dy =                           nab*D*exp(nab*(x+y)) + C*exp(x+y)
// ==>  df2/dy = nb*B*exp(nb*y) + exp(y) + nab*D*exp(nab*(x+y)) + C*exp(x+y)

void InteractClass::CalculateInteraction()
{
  long e, ele, evon, i, maxele;
  double na, nb, nab, ka, kb, kab;
  double ConcPerPartA, ConcPerPartB, ConcPerPartAB;
  double noa, nob, noca, nocb, nocab, nosa, nosb, nosab;
  double A, B, C, D;
  double a11, a12, a22, deta;
  double f1, f2, x, y, xnew, ynew, sum;
  double expx, expy, expnx, expny, expxy, expnxy;
  const int MaxIteration=20;
  int convergence, noab;
  int cpus=partrace->parallel->cpus();
  int cpu=partrace->parallel->mycpu();
  static int ProcLast=0;
  static double *sorparma, *sorparmb, *sorparmab;
  const int result_length=9;
  static double *result=0;
  static double *resultarray1=0;
  static double *resultarray=0;
  double wc, wcvol;
  const double* soilprop;

  if(resultarray==0) {
     resultarray1=new double[result_length];
     resultarray=new double[result_length*cpus];
  }
  ElementsClass* elements=partrace->elements;  // pointer to element-data
  maxele=elements->NoElxyz;
  sorparma=Particle1->adsorption->new_buffer();
  sorparmb=Particle2->adsorption->new_buffer();
  sorparmab=TargetParticle->adsorption->new_buffer();

  // loop over all elements
  for(evon=0;evon<maxele;evon+=cpus) {
     ele=evon+cpu;
     result=resultarray1;
     if(ele<maxele) {
       // distribute computation of interactions
       wc=elements->watercontent->get(ele);
       wcvol = ( elements->volume->get(ele) * wc );
       ConcPerPartA=Particle1->MassPerPart/wcvol;
       ConcPerPartB=Particle2->MassPerPart/wcvol;
       ConcPerPartAB=TargetParticle->MassPerPart/wcvol;
       Particle1->adsorption->get_line(ele, sorparma);
       Particle2->adsorption->get_line(ele, sorparmb);
       TargetParticle->adsorption->get_line(ele, sorparmab);
       soilprop=elements->soilprop->get_line(ele);
       na=sorparma[1];
       nb=sorparmb[1];
       nab=sorparmab[1];
       ka=sorparma[0]*soilprop[4]/wc;
       kb=sorparmb[0]*soilprop[4]/wc;
       kab=sorparmab[0]*soilprop[4]/wc;
       noca=Particle1->Concentration(ele)/ConcPerPartA;
       nocb=Particle2->Concentration(ele)/ConcPerPartB;
       nocab=TargetParticle->Concentration(ele)/ConcPerPartAB;
       nosa=Particle1->ConcentrationSorbed(ele)/ConcPerPartA;
       nosb=Particle2->ConcentrationSorbed(ele)/ConcPerPartB;
       nosab=TargetParticle->ConcentrationSorbed(ele)/ConcPerPartAB;
       noa=nosa+noca+nosab+nocab;
       nob=nosb+nocb+nosab+nocab;
       if(noa+nob<0.5) result[8]=0.0; // no particles in element, nothing to do
       else {
	 result[8]=1.0;
	 /* Newton iteration */
	 convergence=0;
	 A=ka*pow(ConcPerPartA,na-1.0);
	 B=kb*pow(ConcPerPartB,nb-1.0);
	 C=param[0];
	 D=kab*pow(ConcPerPartAB,nab-1.0)*pow(param[0],nab);
	 
	 // set  start values close to solution
	 if(ka>3.0) {  // for large k nosa=noca+nosa
	   if(noa==0.0) noa=0.01; // small value to avoid noa=0
	   noca=pow(noa*ConcPerPartA/ka,1/na);
	 }
	 else {  // for small k noca=noca+nosa 
	   if(noa==0.0) noa=noca=0.01; // small value to avoid noa=0
	   else noca=noa;
	 } 
	 
	 if(kb>3.0) {  // for large k nosb=nocb+nosb
	   if(nob==0.0) nob=0.01; // small value to avoid nob=0
	   nocb=pow(nob*ConcPerPartB/kb,1/nb);
	 }
	 else { // for small k nocb=nocb+nosb 
	   if(nob==0.0) nob=nocb=0.01; // small value to avoid nob=0
	   else nocb=nob;
	 }
	 
	 /* Newton iteration */
	 convergence=0;
	 x=log(noca);
	 y=log(nocb); 
	   
	 for(i=0;convergence==0 && i<MaxIteration;i++) {
	   expx=exp(x);
	   expnx=exp(na*x);
	   expny=exp(nb*y);
	   expy=exp(y);
	   expxy=exp(x+y);
	   expnxy=exp(nab*(x+y));
	   a12=C*expxy+nab*D*expnxy;       
	   a11=expx+na*A*expnx+a12;
	   a22=expy+nb*B*expny+a12;
	   deta=a11*a22-a12*a12;
	   if(deta==0.0) {
	     std::cout<<"Determinante null"<<std::endl;
	     break; // what else?
	   }
	   f1=f2=C*expxy+D*expnxy;
	   f1+=expx+A*expnx-noa; 
	   f2+=expy+B*expny-nob;
	   xnew=x - (a22*f1-a12*f2)/deta;
	   ynew=y - (a11*f2-a12*f1)/deta;
	   // convergence?
	   sum=fabs(x)+fabs(y);
	   if(sum==0.0) sum=1;
	   if( (fabs(xnew-x)+fabs(ynew-y))/sum < 1e-4) convergence=1;
	   x=xnew;
	   y=ynew;
	 }
	 
	 noca=exp(x);
	 nocb=exp(y);
	   
	 nocab=param[0]*noca*nocb;
	 nosa=ka*pow(noca*ConcPerPartA,na)/ConcPerPartA;
	   
	 nosab=noa-nosa-noca-nocab;
	 nosb=nob-nocb-nosab-nocab;
	 result[0]=nocab+nosab-(TargetParticle->Concentration(ele)
			       +TargetParticle->ConcentrationSorbed(ele))/ConcPerPartAB;
	 result[1]=nosa;  result[2]=noca;
	 result[3]=nosb;  result[4]=nocb;
	 result[5]=nosab; result[7]=nocab;
       } // no particles in element
     } // ele<maxele
     
//      if(i==MaxIteration) {
//        cout<<"element "<<ele<<endl;
//        cout<<"start noca:"<<Particle1->Concentration(ele)/Particle1->ConcPerPart;
//        cout<<" start nocb:"<<Particle2->Concentration(ele)/Particle2->ConcPerPart;
//        cout<<" start nocab:"<<TargetParticle->Concentration(ele)/TargetParticle->ConcPerPart<<endl;
//        cout<<"start nosa:"<<Particle1->ConcentrationSorbed(ele)/Particle1->ConcPerPart;
//        cout<<" start nosb:"<<Particle2->ConcentrationSorbed(ele)/Particle2->ConcPerPart;
//        cout<<" start nosab:"<<TargetParticle->ConcentrationSorbed(ele)/TargetParticle->ConcPerPart<<endl;
//        cout<<"noca,nocb,nocab,noab: "<<noca<<' '<<nocb<<' '<<nocab<<' '<<noab<<endl;
//        cout<<"nosa,nosb,nosab     : "<<nosa<<' '<<nosb<<' '<<nosab<<endl;
//      }
//     else cout<<"convergence"<<endl;

     // Broadcast results
     if(cpus>1) {
       partrace->parallel->gather(result, resultarray, result_length);
       // create/delete particles for calculated elements
       result=resultarray;
     }
     for(e=evon;e<evon+cpus && e<maxele;e++) {
       if(result[8]>0.5) {
	 if(result[0]>=0.0) noab=int(result[0]+0.5);
	 else noab=int(result[0]-0.5);
	 // distribute particles on processors
	 TargetParticle->Distribute(noab, ProcLast);
	 // create/crack molecules
	 if(noab>0) CreateMolecules(e, noab);
	 else if(noab<0) CrackMolecules(e, noab);
	 // sorption of particles
	 if( (sum=result[1]+result[2])<=0.0 ) x=0.0;
	 else x=result[1]/sum;
	 Sorption(e, Particle1, x);
	 if( (sum=result[3]+result[4])<=0.0 ) x=0.0;
	 else x=result[3]/sum;
	 Sorption(e, Particle2, x);
	 if( (sum=result[5]+result[7])<=0.0 ) x=0.0;
	 else x=result[5]/sum;
	 Sorption(e, TargetParticle, x);
       }
       result+=result_length;
     }
  } // end element loop
  // delete local buffers
  Particle1->adsorption->delete_buffer(sorparma);
  Particle2->adsorption->delete_buffer(sorparmb);
  TargetParticle->adsorption->delete_buffer(sorparmab);
}


void InteractClass::CreateMolecules(long e, int noab)
{
  // create molecule ab, delete molecule a and b
  ParticleList* plist_a=Particle1->ParticlesInElement + e;
  ParticleList* plist_b=Particle2->ParticlesInElement + e;
  ParticleList* plist_ab=TargetParticle->ParticlesInElement + e;
  Particle *prev1, *prev2;
  Particle* pa=plist_a->start(prev1);
  Particle* pb=plist_b->start(prev2);
  for(int i=0;i<noab;i++) {
     if(pa==0 || pb==0) break;
     pa->position[0]=0.5*(pa->position[0]+pb->position[0]);
     pa->position[1]=0.5*(pa->position[1]+pb->position[1]);
     pa->position[2]=0.5*(pa->position[2]+pb->position[2]);
     pa=plist_a->moveto(pa, prev1, plist_ab);
     pb=plist_b->remove(pb, prev2);
  }
}


void InteractClass::CrackMolecules(long e, int noab)
{
  // create molecule a and b, delete molecule ab
  ParticleList* plist_a=Particle1->ParticlesInElement + e;
  ParticleList* plist_b=Particle2->ParticlesInElement + e;
  ParticleList* plist_ab=TargetParticle->ParticlesInElement + e;
  Particle *prev;
  Particle *pab=plist_ab->start(prev);
  for(int i=noab;i<0;i++) {
     if(pab==0) break;
     plist_a->insert(new Particle(pab->position));
     pab=plist_ab->moveto(pab, prev, plist_b);
  }
}


void InteractClass::Sorption(long e, ParticleClass* particles, double x)
{
  ParticleList* plist=particles->ParticlesInElement + e;
  ParticleList* plist2=particles->ParticlesInElementSorbed + e;
  Particle *p2, *prev, *ptemp;
  RandomClass *rand=partrace->Random;

  // process list of not sorbed particles
  Particle *p=plist->start(prev);
  plist2->start(p2); // set p2
  while(p) {
    if( rand->fastdraw() <= x ) {
      // particle is sorbed
      ptemp=p;
      p=plist->remove_from_list(p, prev);
      p2=plist2->insert_after(ptemp,p2);
    }
    else p=plist->next(p,prev);
  }
  // process list of sorbed particles
  prev=p2;
  p=prev->next;
  while(p) {
    if( rand->fastdraw() > x ) {
      // particle is not sorbed
      ptemp=p;
      p=plist2->remove_from_list(p, prev);
      plist->insert(ptemp);
    }
    else p=plist->next(p,prev);
  }
}

