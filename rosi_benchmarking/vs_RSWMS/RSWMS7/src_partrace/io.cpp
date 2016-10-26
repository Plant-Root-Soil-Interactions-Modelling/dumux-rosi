//Source File io.C
//Jul 1997
//Horst Hardelauf
//Function: Reading the Partrace input file
//          and doing some initialization.

#include "io.hpp"
#include "partrace.hpp"
#include "parallel.hpp"
#include "interact.hpp"
#include "event.hpp"
#include "random.hpp"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstring>

using namespace std;

indata::indata(): BTC_List(0), Inter_List(0), Boundaries(0)
{
  PartTypes=BTCNo=InterNo=BoundariesNo=ElementSequenceNo=0;
}


indata::~indata()
{
  delete [] BTC_List;
  delete [] Inter_List;
  delete [] Boundaries;
}


void PartraceClass::input(string& name)
{
  int me=parallel->mycpu();
  int i, j, n;
  const int linesize=1000;
  const char nl='\n';
  string line;
  indata idata;
  partdata *pdata=0;
  string InpPath, OutPath, DispField, VelField, fn_conc, fn_btc, fn_mom;
  string fn_sorption, fn_wc;
  string RestartFile, RestartFromConcFile, RestartFileOut, fnParticles, CoordFile;
  string PartName, PartAbbrev, fn1, fn2, DecayField, UptakeField;
  string InterName, Part1, Part2, TargetPart;
  int *esequence=0;

  if(me==0) {
    cout<<"Reading Partrace input file: "<<name<<endl;
    ifstream fin(name.c_str());
    if(fin.fail()) error("Open error for the Partrace input-file");
    getline(fin,line);
    if(line!=string("*** Begin of Partrace-Input File V10 ***"))
      error("Wrong file format for the partrace-input file");
    getline(fin,InpPath); Strip(InpPath);
    appendSlash(InpPath);
    getline(fin,OutPath); Strip(OutPath);
    appendSlash(OutPath);
    getline(fin,line);
    if(line.find("PrintMass")==string::npos) idata.PrintMass=0;
    else idata.PrintMass=1;
    getline(fin,VelField); Strip(VelField);
    getline(fin,fn_conc); Strip(fn_conc);
    getline(fin,fn_btc); Strip(fn_btc);
    getline(fin,fn_mom); Strip(fn_mom);

    fin>>idata.simtime; fin.ignore(linesize,nl); 
    fin>>idata.SaveDistanceConc;  fin.ignore(linesize,nl); 
    fin>>idata.NumberOfElements[0]; fin.ignore(linesize,nl); 
    fin>>idata.NumberOfElements[1]; fin.ignore(linesize,nl); 
    fin>>idata.NumberOfElements[2]; fin.ignore(linesize,nl); 
    fin>>idata.SideLength[0];  fin.ignore(linesize,nl); 
    fin>>idata.SideLength[1];  fin.ignore(linesize,nl); 
    fin>>idata.SideLength[2]; fin.ignore(linesize,nl); 

    fin>>idata.UseDispField;  fin.ignore(linesize,nl); 
    fin>>idata.DispMean[0];   fin.ignore(linesize,nl); 
    fin>>idata.DispMean[1];   fin.ignore(linesize,nl); 
    getline(fin,DispField); Strip(DispField);
    fin>>idata.DispWeight[0]; fin.ignore(linesize,nl); 
    fin>>idata.DispWeight[1]; fin.ignore(linesize,nl); 
    fin>>idata.DispWeight[2]; fin.ignore(linesize,nl); 
    fin>>idata.DispDiffusion; fin.ignore(linesize,nl); 
    fin>>idata.UseVelField;   fin.ignore(linesize,nl); 
    fin>>idata.VelMean[0];    fin.ignore(linesize,nl); 
    fin>>idata.VelMean[1];    fin.ignore(linesize,nl); 
    fin>>idata.VelMean[2];    fin.ignore(linesize,nl); 
    fin>>idata.VelVariance[0];fin.ignore(linesize,nl); 
    fin>>idata.VelVariance[1];fin.ignore(linesize,nl); 
    fin>>idata.VelVariance[2];fin.ignore(linesize,nl); 
    fin>>idata.VelTimeStepFactor;fin.ignore(linesize,nl);
    fin>>idata.MaxTime;fin.ignore(linesize,nl); //MB
    fin>>idata.DiffusionType; fin.ignore(linesize,nl); 
    fin>>idata.WCType; fin.ignore(linesize,nl); 
    fin>>idata.Porosity;      fin.ignore(linesize,nl); 
    fin>>idata.BulkDensity;   fin.ignore(linesize,nl); 

    fin.ignore(linesize,nl); 
    fin>>idata.VelType;   fin.ignore(linesize,nl); 
    fin>>idata.DistributeData; fin.ignore(linesize,nl);
    fin.ignore(linesize,nl);
    getline(fin,RestartFile); Strip(RestartFile);
    fin.ignore(linesize,nl);
    fin.ignore(linesize,nl); 
    fin.ignore(linesize,nl);
    fin>>idata.CompRand;  fin.ignore(linesize,nl); 
    fin>>idata.CompSeed;  fin.ignore(linesize,nl);
    //
    getline(fin,RestartFileOut); Strip(RestartFileOut);
    fin>>idata.SaveDistanceMom; fin.ignore(linesize,nl);
    fin>>idata.SaveDistanceRestart; fin.ignore(linesize,nl);
    fin>>idata.FirstTimeRestart; fin.ignore(linesize,nl);
    fin>>idata.SaveParticlesTimeStep; fin.ignore(linesize,nl);
    getline(fin,fnParticles); Strip(fnParticles);
    getline(fin,fn_wc); Strip(fn_wc);
    fin>>idata.GeometryType; fin.ignore(linesize,nl);
    getline(fin,CoordFile); Strip(CoordFile);
    fin>>idata.UseWC; fin.ignore(linesize,nl);
    // expand filenames
    if(!InpPath.empty()) {
      if(!DispField.empty() && DispField[0]!='/') DispField=InpPath + DispField;
      if(!VelField.empty() && VelField[0]!='/') VelField=InpPath + VelField;
      if(!fn_wc.empty() && fn_wc[0]!='/') fn_wc=InpPath + fn_wc;
      if( !RestartFile.empty() && RestartFile[0]!='/' )
	RestartFile=InpPath + RestartFile;
      if(!CoordFile.empty() && CoordFile[0]!='/') CoordFile=InpPath + CoordFile;
    }
    if(!OutPath.empty()) {
      if(!fn_conc.empty() && fn_conc[0]!='/') fn_conc=OutPath + fn_conc;
      if(!fn_btc.empty() && fn_btc[0]!='/')  fn_btc=OutPath + fn_btc;
      if(!fn_mom.empty() && fn_mom[0]!='/')  fn_mom=OutPath + fn_mom;
      if( !RestartFileOut.empty() && RestartFileOut[0]!='/' )
	RestartFileOut=OutPath + RestartFileOut;
      if(!fnParticles.empty() && fnParticles[0]!='/') fnParticles=OutPath + fnParticles;
    }
    strncpy(idata.InpPath,InpPath.c_str(),199);
    strncpy(idata.OutPath,OutPath.c_str(),199);
    strncpy(idata.DispField,DispField.c_str(),199);
    strncpy(idata.VelField,VelField.c_str(),199);
    strncpy(idata.fn_wc,fn_wc.c_str(),199);
    strncpy(idata.fn_conc,fn_conc.c_str(),199);
    strncpy(idata.fn_btc,fn_btc.c_str(),199);
    strncpy(idata.fn_mom,fn_mom.c_str(),199);
    strncpy(idata.RestartFile,RestartFile.c_str(),199);
    strncpy(idata.RestartFileOut,RestartFileOut.c_str(),199);
    strncpy(idata.fnParticles,fnParticles.c_str(),199);
    strncpy(idata.CoordFile,CoordFile.c_str(),199);

    fin>>idata.PartTypes; fin.ignore(linesize,nl); 
    SetNewHandlerText("idata.PartTypes");
    pdata=new partdata[idata.PartTypes];
    cout<<idata.PartTypes<<" molecule types"<<endl;
    for(i=0;i<idata.PartTypes;i++) {
      getline(fin,PartName);
      getline(fin,PartAbbrev);
      getline(fin,fn_sorption);
      fin>>pdata[i].PartSorpType; fin.ignore(linesize,nl); 
      fin>>pdata[i].PartUseSorpField; fin.ignore(linesize,nl); 
      fin>>pdata[i].PartMean[0]; fin.ignore(linesize,nl); 
      fin>>pdata[i].PartMean[1]; fin.ignore(linesize,nl); 
      pdata[i].PartMean[2]=1.0;
      fin>>n; fin.ignore(linesize,nl);
      pdata[i].injectionNo=n;
      pdata[i].injection=new InjectionClass[n];
      pdata[i].PartNumber=0;
      for(j=0;j<n;j++) {
	fin>>pdata[i].injection[j].StartModel; fin.ignore(linesize,nl); 
	fin>>pdata[i].injection[j].StartPos[0]; fin.ignore(linesize,nl); 
	fin>>pdata[i].injection[j].StartPos[1]; fin.ignore(linesize,nl); 
	fin>>pdata[i].injection[j].StartPos[2]; fin.ignore(linesize,nl); 
	fin>>pdata[i].injection[j].StartRadius[0]; fin.ignore(linesize,nl); 
	fin>>pdata[i].injection[j].StartRadius[1]; fin.ignore(linesize,nl); 
	fin>>pdata[i].injection[j].StartRadius[2]; fin.ignore(linesize,nl); 
	fin>>pdata[i].injection[j].Start; fin.ignore(linesize,nl); 
	fin>>pdata[i].injection[j].End; fin.ignore(linesize,nl); 
	if(pdata[i].injection[j].StartModel==5) { // constant initial concentration
	  pdata[i].injection[j].Start=pdata[i].injection[j].End=0.0;
	  pdata[i].injection[j].ParticlesToInject=1;
	  fin>>pdata[i].injection[j].InitConcentration;
	  fin.ignore(linesize,nl); 
	}
	else if(pdata[i].injection[j].StartModel==6) { // initial concentration from file
	  pdata[i].injection[j].ParticlesToInject=1;
	  getline(fin,RestartFromConcFile);
	}
	else {
	  pdata[i].injection[j].InitConcentration=0.0;
	  fin>>pdata[i].injection[j].ParticlesToInject;
	  fin.ignore(linesize,nl); 
	}
	pdata[i].PartNumber+=pdata[i].injection[j].ParticlesToInject;
      }
      fin>>pdata[i].PartMass; fin.ignore(linesize,nl); 
      fin>>pdata[i].PartC0; fin.ignore(linesize,nl); 
      fin>>pdata[i].PartKdNB[0]; fin.ignore(linesize,nl); 
      fin>>pdata[i].PartKdNB[1]; fin.ignore(linesize,nl); 
      fin>>pdata[i].PartKdNB[2]; fin.ignore(linesize,nl); 
      fin>>pdata[i].DecayOrder; fin.ignore(linesize,nl); 
      fin>>pdata[i].UseDecayField; fin.ignore(linesize,nl); 
      fin>>pdata[i].DecayRate; fin.ignore(linesize,nl); 
      getline(fin,DecayField);
      fin>>pdata[i].UptakeOrder; fin.ignore(linesize,nl); 
      fin>>pdata[i].UseUptakeField; fin.ignore(linesize,nl); 
      fin>>pdata[i].UptakeRate; fin.ignore(linesize,nl); 
      getline(fin,UptakeField);
      // expand filenames
      if(!InpPath.empty()) {
	if(fn_sorption[0]!='/') fn_sorption=InpPath+fn_sorption;
	if(DecayField[0]!='/') DecayField=InpPath+DecayField;
      }
      strncpy(pdata[i].PartName,PartName.c_str(),199);
      strncpy(pdata[i].PartAbbrev,PartAbbrev.c_str(),199);
      strncpy(pdata[i].fn_sorption,fn_sorption.c_str(),199);
      strncpy(pdata[i].DecayField,DecayField.c_str(),199);
      strncpy(pdata[i].RestartFromConcFile,RestartFromConcFile.c_str(),199);
      cout<<"molecule type: "<<i+1<<": Sorptionstype: "<<pdata[i].PartSorpType<<endl;
    }
    // read BTC-data
    fin>>idata.BTCNo; fin.ignore(linesize,nl); 
    cout<<"No. of BTCs: "<<idata.BTCNo<<endl;
    if(idata.BTCNo<0) idata.BTCNo=0;
    SetNewHandlerText("idata.BTCNo");
    idata.BTC_List=new vector3[idata.BTCNo];
    for(i=0;i<idata.BTCNo;i++) {
      fin>>idata.BTC_List[i][0]>>idata.BTC_List[i][1]>>idata.BTC_List[i][2]; 
      fin.ignore(linesize,nl); 
    }
    // read interaction-data
    fin>>idata.InterNo; fin.ignore(linesize,nl); 
    cout<<"No. of interactions: "<<idata.InterNo<<endl;
    if(idata.InterNo>0) {
      SetNewHandlerText("idata.Inter_List");
      idata.Inter_List=new interdata[idata.InterNo];
      for(i=0;i<idata.InterNo;i++) {
	getline(fin,InterName);
	getline(fin,TargetPart);
	getline(fin,Part1);
	getline(fin,Part2);
	strncpy(idata.Inter_List[i].name,InterName.c_str(),199);
	strncpy(idata.Inter_List[i].target,TargetPart.c_str(),199);
	strncpy(idata.Inter_List[i].part1,Part1.c_str(),199);
	strncpy(idata.Inter_List[i].part2,Part2.c_str(),199);
	for(j=0;j<idata.Inter_List[i].paramNo;j++) {
	  fin>>idata.Inter_List[i].param[j];
	  fin.ignore(linesize,nl); 
	}
      }
    }
    // read boundary data
    fin>>idata.BoundariesNo; fin.ignore(linesize,nl);
    if(idata.BoundariesNo>0) {
      SetNewHandlerText("idata.Boundaries");
      idata.Boundaries=new boundarydata[idata.BoundariesNo];
      for(i=0;i<idata.BoundariesNo;i++) {
	fin>>idata.Boundaries[i].type; fin.ignore(linesize,nl);
	fin>>idata.Boundaries[i].side; fin.ignore(linesize,nl);
	fin>>idata.Boundaries[i].x1; fin.ignore(linesize,nl);
	fin>>idata.Boundaries[i].x2; fin.ignore(linesize,nl);
	fin>>idata.Boundaries[i].y1; fin.ignore(linesize,nl);
	fin>>idata.Boundaries[i].y2; fin.ignore(linesize,nl);
	switch(idata.Boundaries[i].type) {
	case 1:
	  break;
	case 2:
	  getline(fin, line); Strip(line);
	  strncpy(idata.Boundaries[i].particlename, line.c_str(), 199);
	  fin>>idata.Boundaries[i].value; fin.ignore(linesize,nl);
	  fin>>idata.Boundaries[i].InjectFrom; fin.ignore(linesize,nl);
	  fin>>idata.Boundaries[i].InjectTo; fin.ignore(linesize,nl);
	  break;
	case 3:
	  getline(fin, line); Strip(line);
	  strncpy(idata.Boundaries[i].particlename, line.c_str(), 199);
	  getline(fin, line); Strip(line);
	  strncpy(idata.Boundaries[i].injectfile, line.c_str(), 199);
	  break;
	case 4:
	  getline(fin, line); Strip(line);
	  strncpy(idata.Boundaries[i].particlename, line.c_str(), 199);
	  break;
	case 5:
	  break;
	default:
	  error("Invalid boundary type.");
	  break;
	}
      }
    }
    // read elment sequences
    fin>>idata.ElementSequenceNo; fin.ignore(linesize,nl);
    cout<<"No. of element sequences: "<<idata.ElementSequenceNo<<endl;
    if(idata.ElementSequenceNo>0) {
      SetNewHandlerText("idata.ElementSequence");
      esequence=new int[6*idata.ElementSequenceNo];
      for(j=i=0; i<idata.ElementSequenceNo; i++, j+=6) {
	fin>>esequence[j  ]>>esequence[j+1]>>esequence[j+2]
	   >>esequence[j+3]>>esequence[j+4]>>esequence[j+5];
	fin.ignore(linesize,nl);
      }
    }
    // read additional data
		// default
    idata.OutputType=1;
    idata.RestartType = 1;  // 1=full 2=particles only
    idata.ConvectionMethod = 0; // Integration of velocity, Pollock, 0=off, 1=on //works only with finite volume
    idata.CoefAtBoundary = 1;  // Calculate Reflection Coefficients at boundary 0=off from at starting pos and to at boundary, fastest, 1=on 2=on old and new position "mirror"
    idata.SecDispersiveNew = 1;  // New dispersive displacement is random
    idata.DispersionMethod = 3;  // Accounting for discontinuous dispersion, 1=off, 2=reflection method, 3=interpolation method
    idata.ReflectionMethod = 3;  // ReflectionCoefficient 1=Hoteit et al. 2002, 2=Lim 2006, 3=Bechtold 2010
    idata.NoSurfaceReflection = 1;  // 1=on inhibit reflection at surface, important for upward flow
    idata.TimeSplitting = 2;  // Time Splitting 1=linear 2= nonlinear
    idata.RandomNumber= 2; // gauss
    idata.Recalculate = 1; // Recalculate Convection at interfaces 0=off 1=on
    idata.EleVxVy0=0; // useful to test how solute accumulation is affected by advective "bouncing"
    idata.sub = 1; // Subdivision/Refinement along one axes // works only for finite volume
    idata.ReflectInside=0; // read inner reflexion boudaries from file
    string text1, text2;
    while(!fin.eof()) {
      fin>>text1>>text2;
      if(text1=="***") break;
      if(text1=="OutputType") fin>>idata.OutputType;
      else if(text1=="RestartType") fin>>idata.RestartType;
      else if(text1=="ConvectionMethod") fin>>idata.ConvectionMethod; 
      else if(text1=="DispersionMethod") fin>>idata.DispersionMethod;
      else if(text1=="CoefAtBoundary") fin>>idata.CoefAtBoundary;
      else if(text1=="SecDispersiveNew") fin>>idata.SecDispersiveNew;
      else if(text1=="ReflectionMethod") fin>>idata.ReflectionMethod;
      else if(text1=="NoSurfaceReflection") fin>>idata.NoSurfaceReflection;
      else if(text1=="EleVxVy0") fin>>idata.EleVxVy0;
      else if(text1=="TimeSplitting") fin>>idata.TimeSplitting;
      else if(text1=="RandomNumber") fin>>idata.RandomNumber;
      else if(text1=="Recalculate") fin>>idata.Recalculate;
      else if(text1=="sub") fin>>idata.sub;
      else if(text1=="reflectInside") { fin>>idata.fnReflectInside; idata.ReflectInside=1; }
      else cerr<<"Invalid additional input data: "<<text1<<endl;
      fin.ignore(linesize,nl);
    }
    fin.close();
    cout<<"Closing Partrace input file"<<endl;
  } // me==0

  // distribute input data
  parallel->broadcast( (void *) &idata, sizeof(idata) );
  PrintMass=(idata.PrintMass==1);

  if(me>0) {
    SetNewHandlerText("partdata");
    pdata=new partdata[idata.PartTypes];
    SetNewHandlerText("idata.BTC_List");
    idata.BTC_List=new vector3[idata.BTCNo];
    if(idata.InterNo>0) {
      SetNewHandlerText("idata.Inter_List");
      idata.Inter_List=new interdata[idata.InterNo];
    }
    if(idata.BoundariesNo>0) {
      SetNewHandlerText("idata.Boundaries");
      idata.Boundaries=new boundarydata[idata.BoundariesNo];
    }
    if(idata.ElementSequenceNo>0) {
      SetNewHandlerText("idata.ElementSequenceNo");
      esequence=new int[6*idata.ElementSequenceNo];
    }
  }
  parallel->broadcast( (void *) pdata, (idata.PartTypes)*sizeof(*pdata) );
  parallel->broadcast( (void *) idata.BTC_List, (idata.BTCNo)*sizeof(*idata.BTC_List) );
  if(idata.InterNo>0)
    parallel->broadcast( (void *) idata.Inter_List, (idata.InterNo)*sizeof(*idata.Inter_List) );
  if(idata.BoundariesNo>0)
    parallel->broadcast( (void *) idata.Boundaries, (idata.BoundariesNo)*sizeof(*idata.Boundaries) );
  if(idata.ElementSequenceNo>0)
    parallel->broadcast( (void *) esequence, (6*idata.ElementSequenceNo)*sizeof(*esequence) );
  SetNewHandlerText("injection data");
  for(i=0;i<idata.PartTypes;i++) {
    if(me>0) pdata[i].injection=new InjectionClass[pdata[i].injectionNo];
    parallel->broadcast( (void *) pdata[i].injection, pdata[i].injectionNo*sizeof(InjectionClass) );
  }

  // init Random
  SetNewHandlerText("RandomClass");
  if(idata.CompRand==0) Random=new RandomClass(me+idata.CompSeed-1);
  else {
    time_t t;
    Random=new RandomClass((double) (time(&t) + me*me*me));
  }

  // init runtime
  if(me==0) cout<<"Initializing runtime"<<endl;
  RunTime.Init(this, idata.simtime);
  
  //init geometry
  if(me==0) cout<<"Initializing geometry"<<endl;
  if(idata.DistributeData==2) {
    distributedData=1;
    if(me==0) cout<<"Distributing arrays over PEs"<<endl;
  }
  else distributedData=0;
  SetNewHandlerText("ElementsClass");
  switch(idata.GeometryType) {
  case 1: // simple
    elements=new ElementsClass(idata.NumberOfElements,idata.SideLength, this);
    break;
  case 2: // xy equally spaced, z not
    elements=new ElementsClassZ(idata.NumberOfElements,idata.SideLength, this);
    break;
  case 3: // xyz not equally spaced
    elements=new ElementsClassXYZ(idata.NumberOfElements,idata.SideLength, this);
    break;
  case 4: // xyz rectilinear
    elements=new ElementsClassRectilinear(idata.NumberOfElements,idata.SideLength, this);
    break;
  default: error("Invalid geometry type");
    break;
  }

  fn_conc=idata.fn_conc; fn_btc=idata.fn_btc; fn_mom=idata.fn_mom;
  CoordFile=idata.CoordFile;
  fnParticles=idata.fnParticles;
  elements->InitFilenames(fn_conc, fn_btc, fn_mom, fnParticles, CoordFile);
  elements->ReadCoordinates();
  elements->InitSides();
  elements->CalculateVolumes();
  elements->OutputType=idata.OutputType;
  elements->RestartType=idata.RestartType;	

  if(me==0) cout<<"Initializing dispersion"<<endl;
  DispField=idata.DispField;
  elements->InitDispersion(DispField, idata.UseDispField, idata.DispMean,
			   idata.DispDiffusion, idata.DispWeight,
			   idata.Porosity, idata.BulkDensity, idata.DiffusionType, idata.ConvectionMethod, idata.DispersionMethod, idata.CoefAtBoundary, idata.SecDispersiveNew, idata.ReflectionMethod, idata.NoSurfaceReflection, idata.EleVxVy0, idata.TimeSplitting, idata.RandomNumber, idata.Recalculate, idata.sub);

  if(me==0) cout<<"Initializing velocity"<<endl;
  VelField=idata.VelField; // FileName of velocity field
  elements->InitVelocity(VelField, idata.UseVelField, idata.VelType,
			 idata.VelMean, idata.VelVariance, idata.VelTimeStepFactor, idata.MaxTime); //MB
  
  if(me==0) cout<<"Initializing water content"<<endl;
  fn_wc=idata.fn_wc; //Filename of water content
  elements->InitWaterContent(fn_wc, idata.UseWC, idata.WCType);

  if(me==0) cout<<"Initializing breakthrough curves"<<endl;
  elements->InitBTC(idata.BTC_List,idata.BTCNo);
  elements->InitElementSequence(esequence, idata.ElementSequenceNo);
  delete [] esequence;

  if(me==0) cout<<"Initializing internal reflections"<<endl;
  elements->InitInternalReflections(idata.ReflectInside, idata.fnReflectInside);

  if(me==0) cout<<"Initializing particles"<<endl;
  ParticleClass *pclass=0;
  long part;
  double mass;
  for(i=0;i<idata.PartTypes;i++) {
    //init Particles
    PartName=pdata[i].PartName; PartAbbrev=pdata[i].PartAbbrev;
    part=pdata[i].PartNumber;
    if(part>1) mass=pdata[i].PartMass/part;
    else mass=pdata[i].PartMass;
    SetNewHandlerText("ParticleClass");
    switch(pdata[i].PartSorpType) {
    case 1: // conservative
      pclass=new ParticleClass(this, part, mass, PartName, PartAbbrev,
			       pdata[i].injectionNo, pdata[i].injection,
			       pdata[i].PartC0);
      break;
    case 2: // linear
      pclass=new ParticleClassTPLinear(this, part, mass, PartName, PartAbbrev,
				       pdata[i].injectionNo, pdata[i].injection,
				       pdata[i].PartC0);
      break;
    case 3: // freundlich
      pclass=new ParticleClassTPFreundlich(this, part, mass, PartName, PartAbbrev,
					   pdata[i].injectionNo, pdata[i].injection,
					   pdata[i].PartC0);
      break;
    case 4: // langmuir
      pclass=new ParticleClassTPLangmuir(this, part, mass, PartName, PartAbbrev,
					 pdata[i].injectionNo, pdata[i].injection,
					 pdata[i].PartC0);
      break;
    case 5: // linear non-equilibrium
      pclass=new ParticleClassTPLinearNE(this, part, mass, PartName, PartAbbrev,
					 pdata[i].injectionNo, pdata[i].injection,
					 pdata[i].PartC0);
      break;
    case 6: // freundlich non-equilibrium
      pclass=new ParticleClassTPFreundlichNE(this, part, mass, PartName, PartAbbrev,
					     pdata[i].injectionNo, pdata[i].injection,
					     pdata[i].PartC0);
      break;
    case 7: // langmuir non-equilibrium
      pclass=new ParticleClassTPLangmuirNE(this, part, mass, PartName, PartAbbrev,
					   pdata[i].injectionNo, pdata[i].injection,
					   pdata[i].PartC0);
      break;
    case 8: // freundlich with Retardation
      pclass=new ParticleClassR(this, part, mass, PartName, PartAbbrev,
				pdata[i].injectionNo, pdata[i].injection,
				pdata[i].PartC0);
      break;
    case 9: // freundlich with interaction
      pclass=new ParticleClassInterF(this, part, mass, PartName, PartAbbrev,
				     pdata[i].injectionNo, pdata[i].injection,
				     pdata[i].PartC0);
      break;
    default: error("Invalid Sorption type"); break;
    }
    particles.AddInList(pclass);
    pclass->SaveConcTime=new Event("concentration", RunTime.Time, idata.SaveDistanceConc);
    RunTime.insert(pclass->SaveConcTime);
    pclass->SaveMomentsTime=new Event("moments", RunTime.Time, idata.SaveDistanceMom);
    RunTime.insert(pclass->SaveMomentsTime);
    pclass->SaveParticlesTime=new Event("particles", RunTime.Time, idata.SaveParticlesTimeStep);
    RunTime.insert(pclass->SaveParticlesTime);
    if(fnParticles.empty()) pclass->SaveParticlesTime->deactivate();
    fn_sorption=pdata[i].fn_sorption;
    // init sorption
    pclass->InitSorption(fn_sorption, pdata[i].PartKdNB, pdata[i].PartUseSorpField);
    DecayField=pdata[i].DecayField;
    pclass->InitDecay(DecayField, pdata[i].DecayOrder, pdata[i].UseDecayField, pdata[i].DecayRate);
    UptakeField=pdata[i].UptakeField;
    pclass->InitUptake(UptakeField, pdata[i].UptakeOrder, pdata[i].UseUptakeField, pdata[i].UptakeRate);
    pclass->CalculateConcentrationInit();
    pclass->RestartFromConcFile=pdata[i].RestartFromConcFile;
    if(me==0) {
      cout<<"Creating "<<pdata[i].PartNumber<<" Particles "<<pclass->ParticleName()
	  <<" with sorptions-typ: "<<pclass->SorptionType()<<endl;
      cout<<pclass->DecayType()<<endl;
      cout<<pclass->UptakeType()<<endl;
    }
  } // init particles

  delete [] pdata;

  // interaction initialization
  ParticleClass *target, *p1, *p2;
  InteractClass *inter;
  for(i=0;i<idata.InterNo;i++) {
    InterName=idata.Inter_List[i].name;
    TargetPart=idata.Inter_List[i].target;
    Part1=idata.Inter_List[i].part1;
    Part2=idata.Inter_List[i].part2;
    target=particles.FindParticleType(TargetPart);
    p1=particles.FindParticleType(Part1);
    p2=particles.FindParticleType(Part2);
    SetNewHandlerText("InteractClass");
    inter=new InteractClass(InterName, target, p1, p2, idata.Inter_List[i].param, this);
  }

  // restart initialization
  RestartFile=idata.RestartFile;
  RestartFileOut=idata.RestartFileOut;
  particles.SetRestart(RestartFile, RestartFileOut,
		       idata.FirstTimeRestart, idata.SaveDistanceRestart, me);
  if(particles.LoadRestart) particles.LoadRestartFile();

  if(!PartraceCaller) elements->ReadFlowFields(particles.LoadRestart);
  if(me==0) {
    particles.OpenAllBTCs(RunTime.Time);
    cout<<"Time used for input and initialization: "
	<<RunTime.Seconds()<<"secs"<<endl;
  }

  // init boundaries
  elements->InitBoundaries(idata.Boundaries, idata.BoundariesNo);

}

