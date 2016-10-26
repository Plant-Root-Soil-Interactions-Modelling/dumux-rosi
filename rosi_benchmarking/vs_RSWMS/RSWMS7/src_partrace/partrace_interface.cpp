#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <string>
// partrace includes
#include "partrace_interface.hpp"
#include "partrace.hpp"
#include "distdata.hpp"
#include "elements.hpp"
#include "parallel.hpp"
#include "particleclass.hpp"
#include "particlelist.hpp"
#include "runtime.hpp"
#include "decay.hpp"
#include "uptake.hpp"


extern "C" void initpartrace_()
{
  
  std::string fn("Input_PARTRACE.pTraceInpV12");
  parswms2partrace=0; partrace=0; partrace_parallel=0; watercontentOnNodes=0;
  const int PartraceCaller=0;
  partrace_parallel=new Parallel;
  partrace=new PartraceClass(PartraceCaller, partrace_parallel);
  std::set_new_handler(PartraceNewHandler);
  partrace->start_debug(false);
  partrace->input(fn);
  //  if(fe.globalNumNP!=partrace->elements->NoNoxyz)
  //    Error("Number of node points differs in Parswms and Partrace!");
  //  std::set_new_handler(NewHandler);  
  //  SetNewHandlerText("setPartraceFields");
  //  parswms2partrace=new int[geometry.petscNNP]; 
  parswms2partrace=new int[partrace->elements->NoNoxyz]; 
  if(partrace->distributedData)
    watercontentOnNodes=new DistData(partrace_parallel, "watercontentOnNodes", partrace->elements->NoNoxyz, 1);
  else
    watercontentOnNodes=new LocalData(partrace_parallel, "watercontentOnNodes", partrace->elements->NoNoxyz, 1);
  // set the partrace index for the parswms nodes
  // for(int i=0;i<geometry.petscNNP;i++)
  //   parswms2partrace[i]=geometry.nodeindex.get(geometry.petsc_appindex.get(i));
}


extern "C" void closepartrace_()
{
  delete [] parswms2partrace;
  delete watercontentOnNodes;
  delete partrace_parallel;
  delete partrace;
}


extern "C" void runpartrace_(double &t, double &dt_ext)
{
  partrace->RunTime.SetRunTime(t, dt_ext); // set the target simulation time for partrace
  //setPartraceFields();
  partrace->elements->calculateMaxTimeStepSize();
  partrace->run();
}


extern "C" void setpartracefields_(int &inodes, double v1[], double v2[], double v3[], double theta[], double sinkE[],int &iElm)
{
  // set velocity and watercontent in partrace
  ConstData *velocity=partrace->elements->velocity; // velocity for each node
  int i; 
  double velo[3];
  ConstData *uptakefield;
  ParticleClass* p=partrace->particles.TopOfList();

  while(p) {  
    // Calculate Uptake
    if(uptakefield=p->Uptake->data())
      for (i=0;i<iElm;i++){
	uptakefield->set_dataline(i, sinkE+i);
	//std::cout << " i "<<i << " uptake " << sinkE[i]  <<std::endl;
      }
    p=partrace->particles.NextInList(p);
  }
  
  for (i=0;i<inodes;i++){
      
    watercontentOnNodes->set_dataline(i, theta+i);
    
    velo[0] = v1[i];
    velo[1] = v2[i];
    velo[2] = v3[i];
    //std::cout << " i "<<i << " velo " << *velo+i << std::endl;
    velocity->set_dataline(i, velo);
  }
  
  partrace->elements->calculate_watercontentOnElements(watercontentOnNodes); 
}


extern "C" void setconcentration_(int &iElm, double concPar[])
{
  // set the parswms concentration from partrace
  int i;
  ParticleClass *particle=partrace->particles.TopOfList(); // pointer to ParticleClass

  for(i=0;i<iElm;i++){
    concPar[i]=0.0;
    concPar[i]=particle->Concentration(i)/particle->Czero();
    //  std::cout << " i "<<i << " conc " << concPar[i] << std::endl;
  }

}


extern "C" void setmass_(int &iElm, double massPar[])
{
 
  int i;
  double mass;
  Particle *particle;
  ParticleList *pl;
  ParticleClass * plist ; // =pclass->firstParticleList();
  // ParticleClass *particle=partrace->particles.TopOfList();

  plist =  partrace ->particles.TopOfList();   
  while(plist){
    pl = plist->firstParticleList();
    for(i=0; i<iElm; i++ ){
      massPar[i]=0.0;
      particle = pl->start();
      while(particle){
	mass =particle->get_mass();   // here should go the mass of each element 
	massPar[i]=massPar[i]+mass;
	particle=pl->next(particle);
      }
      pl++;
      // std::cout << " massPar "<< massPar[i] <<  std::endl;
    }
    plist=partrace ->particles.NextInList(plist);
  }
  // std::cout << " massPar "<< massPar[i] <<  std::endl;
}
