#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include "elements.hpp"
#include "particleclass.hpp"
#include "random.hpp"
#include "partrace.hpp"
#include "parallel.hpp"
#include "distdata.hpp"
#include "event.hpp"
#include "runtime.hpp" //MB
#include <cstdlib>
#include <math.h>

using namespace std;


BoundaryConcentrationClass::BoundaryConcentrationClass(): node(0), mass(0)
{
}

BoundaryConcentrationClass::~BoundaryConcentrationClass()
{
  delete [] node;
  delete [] mass;
}


BoundaryInjectionClass::BoundaryInjectionClass(): mass(0), ready(false), particle(0)
{
}

BoundaryInjectionClass::~BoundaryInjectionClass()
{
  delete [] mass;
}


GlobalConcentrationClass::GlobalConcentrationClass():
  mass(0), massz(0), volcon(0), particle(0)
{
}

GlobalConcentrationClass::~GlobalConcentrationClass()
{
  delete [] mass;
  delete [] massz;
  delete [] volcon;
}


ZLevels::ZLevels(): numint(5), z(0), vol(0), ds(0)
{
  z=new double[numint];
  double dz=1.0/double(numint);
  nullmass=0.05/double(numint);
  for(int i=0; i<numint; i++) z[i]=double(i+1)*dz;
}

ZLevels::~ZLevels()
{
  delete [] z;
  delete [] vol;
  delete [] ds;
}

void ZLevels::init(int nz)
{
  vol=new double[nz];
  ds=new double[nz];
}


BoundaryClass::BoundaryClass(long ex, long ey): dimx(ex), dimy(ey), typedata(0)
{
  // new handler is set in InitBoundaries
  int ndim=dimx*dimy;
  typedata=new int[ndim];
  fill(typedata, typedata+ndim, 0);
}


BoundaryClass::~BoundaryClass() { 
  /// destructor.
  delete [] typedata;
}


InternalReflectionClass2::InternalReflectionClass2(string200 fn, int x, int y, int z):
  InternalReflectionClass(fn, x,y,z)
{
  elementtype=new int[ex*ey];
}

InternalReflectionClass2::~InternalReflectionClass2()
{
  delete elementtype;
}

int InternalReflectionClass2::ReadField(string200 fn, string& msg)
{
  ifstream fin(fn);
  if(fin.fail()) {
    msg=string(fn)+" could not be opened";
    return 1;
  }
  string line;
  getline(fin,line);
  if(line!=string("#partrace internal reflections file")) {
    msg=string("Invalid first line in file: ")+string(fn);
    return 2;
  }
  int x,y,z;
  fin>>x>>y>>z;
  if(ex!=x || ey!=y || ez!=z) {
    msg=string("Invalid dimensions in file: ")+string(fn);
    return 3;
  }
  for(int i=0;i<ex*ey && fin;i++) {
    fin>>elementtype[i];
    if(fin.eof()) {
      msg=string("Unexpected EOF in file: ")+string(fn);
      return 4;
    }
  }
  if(fin.bad()) {
    msg=string("Unexpected error in file: ")+string(fn);
    return 5;
  }
  fin.close();
  return 0;
}


bool InternalReflectionClass2::reflect(int x1, int y1, int z1,
				       int x2, int y2, int z2)
{
  if(elementtype[x1+y1*ex]==elementtype[x2+y2*ex]) return false;
  else return true;
}


ElementsClass::ElementsClass(triple noel, vector3& sidelen, PartraceClass *pt):
  transmat(0), coords_x(0), coords_y(0), coords_z(0),
  boundaryConcentrationNo(0), boundaryConcentration(0),
  boundaryInjectionNo(0), boundaryInjection(0),
  globalConcentrationNo(0), globalConcentration(0),
  SideLength_x(0), SideLength_y(0), SideLength_z(0), volume(0), InternalReflections(0),
  BTCListElNo(0), ElementSequence(0),velocity(0), watercontent(0),
  watercontentOnNodes(0), NextVTime(0), NextWCTime(0), soilprop(0), partrace(pt)
{
  GeometryType=1;
  NoElx=noel[0]; NoEly=noel[1]; NoElz=noel[2]; 
  NoElx_1=NoElx-1; NoEly_1=NoEly-1; NoElz_1=NoElz-1;
  NoNox=noel[0]+1; NoNoy=noel[1]+1; NoNoz=noel[2]+1; 
  NoElxy=NoElx*NoEly;
  NoElxz=NoElx*NoElz;
  NoElyz=NoEly*NoElz;
  NoElxyz=NoElxy*NoElz;
  NoNoxy=NoNox*NoNoy;
  NoNoxyz=NoNoxy*NoNoz;
  DomainSize=sidelen;
  Boundaries[0]=Boundaries[1]=Boundaries[2]=0.0;
  Boundaries[3]=Boundaries[0]+DomainSize[0];
  Boundaries[4]=Boundaries[1]+DomainSize[1];
  Boundaries[5]=Boundaries[2]+DomainSize[2];
  ParticlesOutside[0]=ParticlesOutside[1]=ParticlesOutside[2]=
    ParticlesOutside[3]=ParticlesOutside[4]=ParticlesOutside[5]=0;
  velocity=watercontent=0;
  fnBTC="btc";
  fnConc="savcon";
  fnMom="moments";
  fnParticles="particles";
  VFieldPos=0;
  epsilon=numeric_limits<double>::epsilon();
  vMin=sqrt(numeric_limits<double>::min());
  // coords are not distributed.
  partrace->SetNewHandlerText("ElementsClassRectilinear");
  coords_x=new double[NoNox];
  coords_y=new double[NoNoy];
  coords_z=new double[NoNoz];
  // must be overwriten by derived classes
  length_x=DomainSize[0]/NoElx;
  length_y=DomainSize[1]/NoEly;
  length_z=DomainSize[2]/NoElz;
  VFieldRestartPos=VFieldPos=0;
  processor=partrace->parallel->mycpu();
}


ElementsClass::~ElementsClass() {
  delete velocity;
  delete watercontent;
  delete watercontentOnNodes;
  delete NextVTime;
  delete NextWCTime;
  delete soilprop;
  delete [] boundaryConcentration;
  delete [] boundaryInjection;
  delete [] globalConcentration;
  delete SideLength_x;
  delete SideLength_y;
  delete SideLength_z;
  delete [] coords_x;
  delete [] coords_y;
  delete [] coords_z;
  delete volume;
  delete InternalReflections;
  delete [] BTCListElNo;
  delete [] ElementSequence;
}


void ElementsClass::CalculateVolumes()
{
  int n;
  // initialize the volume(s)
  partrace->SetNewHandlerText("volume");
  volume=new ConstData(partrace->parallel, "volume", NoElxyz, 1);
  volume->init(length_x*length_y*length_z);
  // initialize the side lengths of the elements
  SideLength_x=new ConstData(partrace->parallel, "SideLength_x", NoElxyz,1);
  SideLength_x->init(length_x);
  SideLength_y=new ConstData(partrace->parallel, "SideLength_y", NoElxyz,1);
  SideLength_y->init(length_y);
  SideLength_z=new ConstData(partrace->parallel, "SideLength_z", NoElxyz,1);
  SideLength_z->init(length_z);
  // init rectilinear coordinates
  for(n=0;n<NoNox;n++) coords_x[n]=n*length_x;
  for(n=0;n<NoNoy;n++) coords_y[n]=n*length_y;
  for(n=0;n<NoNoz;n++) coords_z[n]=n*length_z;
}


/** set the maximum size of a simulation time step.
 * Calculates the maximum delta t for all elements using the formula:
 * v=ds/dt ==> dt=ds/v. The velocity of the element is the averaged
 * value of the velocity values of the nodes of the element.
 * Every particle should stay in every element
 * for at least \a TimeStepFactor time steps.
 * \return the calculated maximum time step size.
 * \see CalculateVolumes()
 */
void ElementsClass::calculateMaxTimeStepSize()
{
  double dt=partrace->RunTime.remainingRunTime();
  int i, index;
  long nodes[8];
  long e, ex,ey,ez;
  double *v;
  vector3 pos, vel;
  //ElementsClass *elements=partrace->elements;

  switch(VelType) {
  case 1:
    // velocity defined on elements
    for(e=0; e<NoElxyz; e++) {
      v=velocity->get_line(e);
      if(v[0]!=0.0) dt=min(dt, SideLength_x->get(e) / fabs(v[0]));
      if(v[1]!=0.0) dt=min(dt, SideLength_y->get(e) / fabs(v[1]));
      if(v[2]!=0.0) dt=min(dt, SideLength_z->get(e) / fabs(v[2]));
    }
    break;
  case 2:
    // velocity defined on nodes
    for(e=ez=0; ez<NoElz; ez++) {
      for(ey=0; ey<NoEly; ey++) {
	for(ex=0; ex<NoElx; ex++,e++) {
	  GetNodesOfElement(e, ex,ey,ez, nodes);
	  for(i=0;i<8;i++) {
	    v=velocity->get_line(nodes[i]);
	    if(v[0]!=0.0) dt=min(dt, SideLength_x->get(e) / fabs(v[0]));
	    if(v[1]!=0.0) dt=min(dt, SideLength_y->get(e) / fabs(v[1]));
	    if(v[2]!=0.0) dt=min(dt, SideLength_z->get(e) / fabs(v[2]));
	  }
	}
      }
    }
    break;
  case 3:
    // velocity coming from finite volume
    for(e=ez=0; ez<NoElz; ez++) {
      for(ey=0; ey<NoEly; ey++) {
	for(ex=0; ex<NoElx; ex++,e++) {
	  GetElementsPosition(e, pos);
	  index=6*(ez+ey*NoElz+ex*NoElyz); 
	  vel.set(finitevolvelocity[index  ]*pos[0]+finitevolvelocity[index+1],
		  finitevolvelocity[index+2]*pos[1]+finitevolvelocity[index+3],
		  finitevolvelocity[index+4]*pos[2]+finitevolvelocity[index+5]);
	  vel/=watercontent->get(e);
	  if(vel[0]!=0.0) dt=min(dt, SideLength_x->get(e) / fabs(vel[0]));
	  if(vel[1]!=0.0) dt=min(dt, SideLength_y->get(e) / fabs(vel[1]));
	  if(vel[2]!=0.0) dt=min(dt, SideLength_z->get(e) / fabs(vel[2]));
	}
      }
    }
    break;
  } // switch
  // let the particle spend at least n=1/TimeStepFactor time steps in one element
  if (TimeStepFactor*dt>MaxTime) { //MB
    partrace->RunTime.setMaxTimeStep(MaxTime); }
  else {
    partrace->RunTime.setMaxTimeStep(TimeStepFactor*dt);
  }
}


bool ElementsClass::CheckNEL(triple n, int plus)
{
  if(n[0]==NoElx+plus && n[1]==NoEly+plus && n[2]==NoElz+plus) return false;
  else return true;
}


void ElementsClass::GetElementXYZ(long e, long& ex, long& ey, long& ez) 
{
  ez=e/NoElxy;
  e-=ez*NoElxy;
  ey=e/NoElx;
  ex=e-ey*NoElx;

}

void ElementsClass::GetDxyzFxyz(long ex, long ey, long ez, vector3& pos, double& Dx, double& Dy, double& Dz, double& Fx, double& Fy, double& Fz) // Get values for factor calculation for trilinear interpolation
{
  Dx=coords_x[ex+1] - coords_x[ex];
  Dy=coords_y[ey+1] - coords_y[ey];
  Dz=coords_z[ez+1] - coords_z[ez];
  Fx=(pos[0]-coords_x[ex])/Dx;
  Fy=(pos[1]-coords_y[ey])/Dy;
  Fz=(pos[2]-coords_z[ez])/Dz;
}

void ElementsClass::GetTrilinearFactors(double Fx, double Fy, double Fz, double *f) //
{
  f[0] = (1 -Fx)*(1 -Fy)*(1 -Fz);
  f[1] = (   Fx)*(1 -Fy)*(1 -Fz);
  f[2] = (1 -Fx)*(   Fy)*(1 -Fz);
  f[3] = (   Fx)*(   Fy)*(1 -Fz);
  f[4] = (1 -Fx)*(1 -Fy)*(   Fz);
  f[5] = (   Fx)*(1 -Fy)*(   Fz);
  f[6] = (1 -Fx)*(   Fy)*(   Fz);
  f[7] = (   Fx)*(   Fy)*(   Fz);
}

void ElementsClass::InterpolateVelocityTriLinear(long ex, long ey, long ez, double *f, int *ni, long *nodes, vector3& v)
{
  ElementsClass *elements=partrace->elements;	
  double *v0=velocity->new_buffer();
  double *v1=velocity->new_buffer();
  double *v2=velocity->new_buffer();
  double *v3=velocity->new_buffer();
  double *v4=velocity->new_buffer();
  double *v5=velocity->new_buffer();
  double *v6=velocity->new_buffer();
  double *v7=velocity->new_buffer();
  //double dx=coords_x[ex+1]-coords_x[ex];
  //double dy=coords_y[ey+1]-coords_y[ey];
  //double dz=coords_z[ez+1]-coords_z[ez];
  //double dxyz=dx*dy*dz;
  ni[0]=((ez  )*elements->NoNoxy + (ey  )*elements->NoNox + (ex  )); // node index
  ni[1]=((ez  )*elements->NoNoxy + (ey  )*elements->NoNox + (ex+1));
  ni[2]=((ez  )*elements->NoNoxy + (ey+1)*elements->NoNox + (ex  ));
  ni[3]=((ez  )*elements->NoNoxy + (ey+1)*elements->NoNox + (ex+1));
  ni[4]=((ez+1)*elements->NoNoxy + (ey  )*elements->NoNox + (ex  ));
  ni[5]=((ez+1)*elements->NoNoxy + (ey  )*elements->NoNox + (ex+1));
  ni[6]=((ez+1)*elements->NoNoxy + (ey+1)*elements->NoNox + (ex  ));
  ni[7]=((ez+1)*elements->NoNoxy + (ey+1)*elements->NoNox + (ex+1));
  switch(elements->VelType) {
  case 1:
    v[0]=f[0]*elements->nodevelocity[3*ni[0]+0] + f[1]*elements->nodevelocity[3*ni[1]+0] + f[2]*elements->nodevelocity[3*ni[2]+0] + f[3]*elements->nodevelocity[3*ni[3]+0]+
      f[4]*elements->nodevelocity[3*ni[4]+0] + f[5]*elements->nodevelocity[3*ni[5]+0] + f[6]*elements->nodevelocity[3*ni[6]+0] + f[7]*elements->nodevelocity[3*ni[7]+0];
    v[1]=f[0]*elements->nodevelocity[3*ni[0]+1] + f[1]*elements->nodevelocity[3*ni[1]+1] + f[2]*elements->nodevelocity[3*ni[2]+1] + f[3]*elements->nodevelocity[3*ni[3]+1]+
      f[4]*elements->nodevelocity[3*ni[4]+1] + f[5]*elements->nodevelocity[3*ni[5]+1] + f[6]*elements->nodevelocity[3*ni[6]+1] + f[7]*elements->nodevelocity[3*ni[7]+1];
    v[2]=f[0]*elements->nodevelocity[3*ni[0]+2] + f[1]*elements->nodevelocity[3*ni[1]+2] + f[2]*elements->nodevelocity[3*ni[2]+2] + f[3]*elements->nodevelocity[3*ni[3]+2]+
      f[4]*elements->nodevelocity[3*ni[4]+2] + f[5]*elements->nodevelocity[3*ni[5]+2] + f[6]*elements->nodevelocity[3*ni[6]+2] + f[7]*elements->nodevelocity[3*ni[7]+2];
    break;
  case 2:
    // get velocity at nodes
    velocity->get_line(nodes[0], v0);
    velocity->get_line(nodes[1], v1);
    velocity->get_line(nodes[2], v2);
    velocity->get_line(nodes[3], v3);
    velocity->get_line(nodes[4], v4);
    velocity->get_line(nodes[5], v5);
    velocity->get_line(nodes[6], v6);
    velocity->get_line(nodes[7], v7);
    v[0]=f[0]*v0[0] + f[1]*v1[0] + f[2]*v2[0] + f[3]*v3[0]+
      f[4]*v4[0] + f[5]*v5[0] + f[6]*v6[0] + f[7]*v7[0];
    v[1]=f[0]*v0[1] + f[1]*v1[1] + f[2]*v2[1] + f[3]*v3[1]+
      f[4]*v4[1] + f[5]*v5[1] + f[6]*v6[1] + f[7]*v7[1];
    v[2]=f[0]*v0[2] + f[1]*v1[2] + f[2]*v2[2] + f[3]*v3[2]+
      f[4]*v4[2] + f[5]*v5[2] + f[6]*v6[2] + f[7]*v7[2];
    break;
  case 3:
    v[0]=f[0]*elements->nodevelocity[3*ni[0]+0] + f[1]*elements->nodevelocity[3*ni[1]+0] + f[2]*elements->nodevelocity[3*ni[2]+0] + f[3]*elements->nodevelocity[3*ni[3]+0]+
      f[4]*elements->nodevelocity[3*ni[4]+0] + f[5]*elements->nodevelocity[3*ni[5]+0] + f[6]*elements->nodevelocity[3*ni[6]+0] + f[7]*elements->nodevelocity[3*ni[7]+0];
    v[1]=f[0]*elements->nodevelocity[3*ni[0]+1] + f[1]*elements->nodevelocity[3*ni[1]+1] + f[2]*elements->nodevelocity[3*ni[2]+1] + f[3]*elements->nodevelocity[3*ni[3]+1]+
      f[4]*elements->nodevelocity[3*ni[4]+1] + f[5]*elements->nodevelocity[3*ni[5]+1] + f[6]*elements->nodevelocity[3*ni[6]+1] + f[7]*elements->nodevelocity[3*ni[7]+1];
    v[2]=f[0]*elements->nodevelocity[3*ni[0]+2] + f[1]*elements->nodevelocity[3*ni[1]+2] + f[2]*elements->nodevelocity[3*ni[2]+2] + f[3]*elements->nodevelocity[3*ni[3]+2]+
      f[4]*elements->nodevelocity[3*ni[4]+2] + f[5]*elements->nodevelocity[3*ni[5]+2] + f[6]*elements->nodevelocity[3*ni[6]+2] + f[7]*elements->nodevelocity[3*ni[7]+2];
    break;
  }
}


void ElementsClass::InterpolateVelocityTriLinearDX(double *f, int *ni, long *nodes, vector3& v)
{
  ElementsClass *elements=partrace->elements;	
  double *v0=velocity->new_buffer();
  double *v1=velocity->new_buffer();
  double *v2=velocity->new_buffer();
  double *v3=velocity->new_buffer();
  double *v4=velocity->new_buffer();
  double *v5=velocity->new_buffer();
  double *v6=velocity->new_buffer();
  double *v7=velocity->new_buffer();
  //double dx=coords_x[ex+1]-coords_x[ex];
  //double dy=coords_y[ey+1]-coords_y[ey];
  //double dz=coords_z[ez+1]-coords_z[ez];
  //double dxyz=dx*dy*dz;
  switch(elements->VelType) {
  case 1: // elementvelocity to nodevelocity not yet implemented
    v[0]=f[0]*elements->nodevelocity[3*ni[0]+0] + f[1]*elements->nodevelocity[3*ni[1]+0] + f[2]*elements->nodevelocity[3*ni[2]+0] + f[3]*elements->nodevelocity[3*ni[3]+0]+
      f[4]*elements->nodevelocity[3*ni[4]+0] + f[5]*elements->nodevelocity[3*ni[5]+0] + f[6]*elements->nodevelocity[3*ni[6]+0] + f[7]*elements->nodevelocity[3*ni[7]+0];
    v[1]=f[0]*elements->nodevelocity[3*ni[0]+1] + f[1]*elements->nodevelocity[3*ni[1]+1] + f[2]*elements->nodevelocity[3*ni[2]+1] + f[3]*elements->nodevelocity[3*ni[3]+1]+
      f[4]*elements->nodevelocity[3*ni[4]+1] + f[5]*elements->nodevelocity[3*ni[5]+1] + f[6]*elements->nodevelocity[3*ni[6]+1] + f[7]*elements->nodevelocity[3*ni[7]+1];
    v[2]=f[0]*elements->nodevelocity[3*ni[0]+2] + f[1]*elements->nodevelocity[3*ni[1]+2] + f[2]*elements->nodevelocity[3*ni[2]+2] + f[3]*elements->nodevelocity[3*ni[3]+2]+
      f[4]*elements->nodevelocity[3*ni[4]+2] + f[5]*elements->nodevelocity[3*ni[5]+2] + f[6]*elements->nodevelocity[3*ni[6]+2] + f[7]*elements->nodevelocity[3*ni[7]+2];
    break;
  case 2: // get velocity at nodes
    velocity->get_line(nodes[0], v0);
    velocity->get_line(nodes[1], v1);
    velocity->get_line(nodes[2], v2);
    velocity->get_line(nodes[3], v3);
    velocity->get_line(nodes[4], v4);
    velocity->get_line(nodes[5], v5);
    velocity->get_line(nodes[6], v6);
    velocity->get_line(nodes[7], v7);
    v[0]=f[0]*v0[0] + f[1]*v1[0] + f[2]*v2[0] + f[3]*v3[0]+
      f[4]*v4[0] + f[5]*v5[0] + f[6]*v6[0] + f[7]*v7[0];
    v[1]=f[0]*v0[1] + f[1]*v1[1] + f[2]*v2[1] + f[3]*v3[1]+
      f[4]*v4[1] + f[5]*v5[1] + f[6]*v6[1] + f[7]*v7[1];
    v[2]=f[0]*v0[2] + f[1]*v1[2] + f[2]*v2[2] + f[3]*v3[2]+
      f[4]*v4[2] + f[5]*v5[2] + f[6]*v6[2] + f[7]*v7[2];
    break;
  case 3: // finite volume VelType
    v[0]=f[0]*elements->nodevelocity[3*ni[0]+0] + f[1]*elements->nodevelocity[3*ni[1]+0] + f[2]*elements->nodevelocity[3*ni[2]+0] + f[3]*elements->nodevelocity[3*ni[3]+0]+
      f[4]*elements->nodevelocity[3*ni[4]+0] + f[5]*elements->nodevelocity[3*ni[5]+0] + f[6]*elements->nodevelocity[3*ni[6]+0] + f[7]*elements->nodevelocity[3*ni[7]+0];
    v[1]=f[0]*elements->nodevelocity[3*ni[0]+1] + f[1]*elements->nodevelocity[3*ni[1]+1] + f[2]*elements->nodevelocity[3*ni[2]+1] + f[3]*elements->nodevelocity[3*ni[3]+1]+
      f[4]*elements->nodevelocity[3*ni[4]+1] + f[5]*elements->nodevelocity[3*ni[5]+1] + f[6]*elements->nodevelocity[3*ni[6]+1] + f[7]*elements->nodevelocity[3*ni[7]+1];
    v[2]=f[0]*elements->nodevelocity[3*ni[0]+2] + f[1]*elements->nodevelocity[3*ni[1]+2] + f[2]*elements->nodevelocity[3*ni[2]+2] + f[3]*elements->nodevelocity[3*ni[3]+2]+
      f[4]*elements->nodevelocity[3*ni[4]+2] + f[5]*elements->nodevelocity[3*ni[5]+2] + f[6]*elements->nodevelocity[3*ni[6]+2] + f[7]*elements->nodevelocity[3*ni[7]+2];
    break;
  }
}

void ElementsClass::InterpolateWCTriLinear(double *f, int *ni, long *nodes, double& wc)
{
  ElementsClass *elements=partrace->elements;	

  double *wc0=watercontent->new_buffer();
  double *wc1=watercontent->new_buffer();
  double *wc2=watercontent->new_buffer();
  double *wc3=watercontent->new_buffer();
  double *wc4=watercontent->new_buffer();
  double *wc5=watercontent->new_buffer();
  double *wc6=watercontent->new_buffer();
  double *wc7=watercontent->new_buffer();
  //double dx=coords_x[ex+1]-coords_x[ex];
  //double dy=coords_y[ey+1]-coords_y[ey];
  //double dz=coords_z[ez+1]-coords_z[ez];
  //double dxyz=dx*dy*dz;
  switch(elements->VelType){
  case 1: // type elementwatercontent get watercontent at nodes // to be implemented
    watercontent->get_line(nodes[0], wc0);
    watercontent->get_line(nodes[1], wc1);
    watercontent->get_line(nodes[2], wc2);
    watercontent->get_line(nodes[3], wc3);
    watercontent->get_line(nodes[4], wc4);
    watercontent->get_line(nodes[5], wc5);
    watercontent->get_line(nodes[6], wc6);
    watercontent->get_line(nodes[7], wc7);
    wc=f[0]* *wc0 + f[1]* *wc1 + f[2]* *wc2 + f[3]* *wc3+
      f[4]* *wc4 + f[5]* *wc5 + f[6]* *wc6 + f[7]* *wc7;
    break;
  case 2: // get watercontent at nodes
    watercontentOnNodes->get_line(nodes[0], wc0);
    watercontentOnNodes->get_line(nodes[1], wc1);
    watercontentOnNodes->get_line(nodes[2], wc2);
    watercontentOnNodes->get_line(nodes[3], wc3);
    watercontentOnNodes->get_line(nodes[4], wc4);
    watercontentOnNodes->get_line(nodes[5], wc5);
    watercontentOnNodes->get_line(nodes[6], wc6);
    watercontentOnNodes->get_line(nodes[7], wc7);
    //cout<<" InterpolateWCTriLinear "<<*wc0 << " " << *wc1 << " " << *wc2 << endl;
    wc=f[0]* *wc0 + f[1]* *wc1 + f[2]* *wc2 + f[3]* *wc3+
      f[4]* *wc4 + f[5]* *wc5 + f[6]* *wc6 + f[7]* *wc7;
    break;
  case 3: // Finite Volume VelType
    wc=f[0]*elements->nodewc[ni[0]] + f[1]*elements->nodewc[ni[1]] + f[2]*elements->nodewc[ni[2]] + f[3]*elements->nodewc[ni[3]]+
      f[4]*elements->nodewc[ni[4]] + f[5]*elements->nodewc[ni[5]] + f[6]*elements->nodewc[ni[6]] + f[7]*elements->nodewc[ni[7]];
    break;
  }
}

void ElementsClass::GetNodeXYZ(long n, long& nx, long& ny, long& nz) 
{
  nz=n/((NoElx+1)*(NoEly+1));
  n-=nz*((NoElx+1)*(NoEly+1));
  ny=n/(NoElx+1);
  nx=n-(ny*(NoElx+1));
}

long ElementsClass::GetElementNo(const vector3& pos)
{
  long x,y,z;
  if(pos[0]<Boundaries[0]) return -1;
  if(pos[1]<Boundaries[1]) return -1;
  if(pos[2]<Boundaries[2]) return -1;
  x=long(pos[0]/length_x);
  if(x>=NoElx) return -1;
  y=long(pos[1]/length_y);
  if(y>=NoEly) return -1;
  z=long(pos[2]/length_z);
  if(z>=NoElz) return -1;
  return x+y*NoElx+z*NoElxy;
}


void ElementsClass::GetStartCoordinates(long nodes[8], long ex, long ey, long ez, vector3& startcoord)
{
  // get the coordinates of the first node of the element
  startcoord[0]=coords_x[ex]; startcoord[1]=coords_y[ey]; startcoord[2]=coords_z[ez];
}


bool ElementsClass::reflect(vector3& x, double* norm)
{
  double p=2.0*(x[0]*norm[0]+x[1]*norm[1]+x[2]*norm[2]);
  x[0]-=p*norm[0];
  x[1]-=p*norm[1];
  x[2]-=p*norm[2];
  return true;
}


long ElementsClass::MoveParticle(long e, long ex, long ey, long ez, //MB
				 vector3& pos, DispersionClass& dispersionData)
{
  static vector3 newpos;
  static double x1,x2, y1,y2, z1,z2;
  static int side;
  static double diff,diff2;
  static vector3 convection, dispersion;
  dispersion=dispersionData.ds;
  double dt=partrace->RunTime.dt;
  ElementsClass *elements=partrace->elements;
  int loop=0;
  //cout.precision( 20 );
  // ******************************* CONVECTION ****************************************
  // loop over boundary passings/reflections
  static int numdt=0;
  static int numtilt=0;

  // cout<<"NoElx "<<NoElx<<" "<<"NoEly "<<NoEly<<" "<<"NoElz "<<NoElz<<endl;
  //sleep(1);

  // cout<<"e "<<e<<" "<<"ex "<<ex<<" "<<"ey "<<ey<<" "<<"ez "<<ez<<" "<<"pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
  //sleep(0.5);
  
  numdt+=1;
  if (numdt==10000000) { // useful print out to check time step in runs
    cout<<"dt "<<dt<<endl;
    numdt=0; }
  
  calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
  if ( (NoSurfaceReflection==1) && (ez>NoElz_1-EleVxVy0) ) { // 1 cm
    convection[0]=0.0; convection[1]=0.0; 
  }
  // cout<<"Convection" <<convection[0]<<" "<<convection[1]<<" "<<convection[2]<<endl;
  while(1) {
    x1=coords_x[ex]; x2=coords_x[ex+1];
    y1=coords_y[ey]; y2=coords_y[ey+1];
    z1=coords_z[ez]; z2=coords_z[ez+1];
    newpos.setsum(pos, convection);
    // cout<<"NewPos Convection" << newpos[0]<<" "<< newpos[1]<<" "<< newpos[2]<<endl;
   
    // x-axis
    if(newpos[0]<x1) {
      diff=(x1-pos[0])/(newpos[0]-pos[0]);
      side=1;
    }
    else if(newpos[0]>x2) {
      side=2;
      diff=(x2-pos[0])/(newpos[0]-pos[0]);
    }
    else {
      diff=2.0;
      side=0;
    }
    // y-axis
    if(newpos[1]<y1) {
      diff2=(y1-pos[1])/(newpos[1]-pos[1]);
      if(diff2<diff) { diff=diff2; side=3; }
    }
    else if(newpos[1]>y2) {
      diff2=(y2-pos[1])/(newpos[1]-pos[1]);
      if(diff2<diff) { diff=diff2; side=4; }
    }
    // z-axis
    if(newpos[2]<z1) {
      diff2=(z1-pos[2])/(newpos[2]-pos[2]);
      if(diff2<diff) { diff=diff2; side=5; }
    }
    else if(newpos[2]>z2) {
      diff2=(z2-pos[2])/(newpos[2]-pos[2]);
      if(diff2<diff) { diff=diff2; side=6; }
    }
    if(side==0) {
      pos=newpos;
      // cout<<"pos=newpos"<<endl;
      // cout<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<"="<<newpos[0]<<" "<<newpos[1]<<" "<<newpos[2]<<endl;
      break;
    }
    pos.addxsum(diff, convection); //move to one of the boundaries of the element,
    diff2=1.0-diff; //calculate the residual diff
    convection*=diff2; //use the residual tiff to calculate remaining movement
    dt*=diff2;
    
    switch(side) {
    case 1: // left
      if(ex==0) {
	if(BoundaryXmin->type(ey,ez)==1) {
	  reflectYZ(convection);
	}
	else if(BoundaryXmin->type(ey,ez)==5) {
	  // cout <<"pos 0 "<< pos[0] << " " <<"coords_x[NoElx] "<<coords_x[NoElx] << endl;
	  // cout <<"ex "<< ex << " " <<"NoElx "<<NoElx << endl;
	  // cout <<"e "<< e << " " <<"e neu "<<e + NoElx_1 << endl;
	  pos[0] = coords_x[NoElx];
	  ex = NoElx_1;
	  e = e + NoElx_1;
	  calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	}
	else { // particle outside
	  ParticlesOutside[0]++;
	  return -1;  // particle has left the simulation area
	}
      }
      else if(InternalReflections->reflect(ex,ey,ez, ex-1,ey,ez)) {
	reflectYZ(convection);
      }
      else {
	//if (elements->Recalculate==0) {ex--;e--;} else {
	calculateConvection(convection, e-1, pos, dt, dispersionData.particle->Retardation(e-1));
	if (convection[0]<0.0) {ex--;e--;}
	else {calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	  if (convection[0]<0.0) {pos[0]+=DomainSize[0]*10e-10; convection[0]=0.0; }}
	//} 
      }
      break; 
      
    case 2: // right
      if(ex==NoElx_1) {
	if(BoundaryXmax->type(ey,ez)==1) {
	  reflectYZ(convection);
	}
	else if(BoundaryXmax->type(ey,ez)==5) {
	  pos[0] = coords_x[0];
	  ex = 0;
	  e = e - NoElx_1;
	  calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	}
	else  { // particle outside
	  ParticlesOutside[3]++;
	  return -1;
	}
      }
      else if(InternalReflections->reflect(ex,ey,ez, ex+1,ey,ez)) {
	reflectYZ(convection);
      }
      else {
	//if (elements->Recalculate==0) {ex++;e++;} else {
	calculateConvection(convection, e+1, pos, dt, dispersionData.particle->Retardation(e+1)); 
	if (convection[0]>0.0) {ex++;e++;}
	else {calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	  if (convection[0]>0.0) {pos[0]-=DomainSize[0]*10e-10; convection[0]=0.0; }}
	//}
      }
      break; 
      
    case 3: // front
      if(ey==0) {
	if(BoundaryYmin->type(ex,ez)==1) {
	  reflectXZ(convection);
	}	
	else if(BoundaryYmin->type(ex,ez)==5) {
	  pos[1] = coords_y[NoEly];
	  ey = NoEly_1;
	  e = e + NoElx * NoEly_1;
	  calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	}
	else  { // particle outside
	  ParticlesOutside[1]++;
	  return -1;
 	}
      }
      else if(InternalReflections->reflect(ex,ey,ez, ex,ey-1,ez)) {
	reflectXZ(convection);
      }
      else {
	//if (elements->Recalculate==0) {ey--;e-=NoElx;} else {
	calculateConvection(convection, e-NoElx, pos, dt, dispersionData.particle->Retardation(e-NoElx)); 
	if (convection[1]<0.0) {ey--;e-=NoElx;}
	else {calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	  if (convection[1]<0.0) {pos[1]+=DomainSize[1]*10e-10; convection[1]=0.0; }}
	//}
      }
      break;
      
    case 4: // back
      if(ey==NoEly_1) {
	if(BoundaryYmax->type(ex,ez)==1) {
	  reflectXZ(convection);
	}
	else if(BoundaryYmax->type(ex,ez)==5) {
	  pos[1] = coords_y[0];
	  ey = 0;
	  e = e - NoElx * NoEly_1;
	  calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	}
	else  { // particle outside
	  ParticlesOutside[4]++;
	  return -1;
	}
      }
      else if(InternalReflections->reflect(ex,ey,ez, ex,ey+1,ez)) {
	reflectXZ(convection);
      }
      else {
	//if (elements->Recalculate==0) {ey++;e+=NoElx;} else {
	calculateConvection(convection, e+NoElx, pos, dt, dispersionData.particle->Retardation(e+NoElx)); 
	if (convection[1]>0.0) {ey++;e+=NoElx;}
	else {calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	  if (convection[1]>0.0) {pos[1]-=DomainSize[1]*10e-10; convection[1]=0.0; }}
	//}
      }
      break;
      
    case 5: // bottom
      if(ez==0) {
	if(BoundaryZmin->type(ex,ey)==1) {
	  //if (convection[2]-fabs(dispersion[2])<null) {ParticlesOutside[2]++; return -1;} // possible alternative
	  reflectXY(convection);
	}
	else  { // particle outside 
	  ParticlesOutside[2]++;
	  return -1;
	}
      }
      else if(InternalReflections->reflect(ex,ey,ez, ex,ey,ez-1)) {
	reflectXY(convection);
      }
      else {
	//if (elements->Recalculate==0) {ez--;e-=NoElxy;} else {
	calculateConvection(convection, e-NoElxy, pos, dt, dispersionData.particle->Retardation(e-NoElxy)); 
	if (convection[2]<0.0) {ez--;e-=NoElxy;}
	else {calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	  if (convection[2]<0.0) {pos[2]+=DomainSize[2]*10e-10; convection[2]=0.0; }}
	//}
      }
      break;
      
    case 6: // top
      if(ez==NoElz_1) {
	if(BoundaryZmax->type(ex,ey)==1) {
	  reflectXY(convection);
	  // Option only reflection of dispersion
	  if (NoSurfaceReflection==1) {
	    convection[0]=0.0; convection[1]=0.0; convection[2]*=DomainSize[0]*10e-10; }
	}
	else  { // particle outside
	  ParticlesOutside[5]++;
	  return -1;
	}
      }
      else if(InternalReflections->reflect(ex,ey,ez, ex,ey,ez+1)) {
	reflectXY(convection);
      }
      else {
	//if (elements->Recalculate==0) {ez++;e+=NoElxy;} else {
	calculateConvection(convection, e+NoElxy, pos, dt, dispersionData.particle->Retardation(e+NoElxy)); 
	if (convection[2]>0.0) {ez++;e+=NoElxy;}
	else {calculateConvection(convection, e, pos, dt, dispersionData.particle->Retardation(e));
	  if (convection[2]>0.0) {pos[2]-=DomainSize[2]*10e-10; convection[2]=0.0; }}
	//} 
      }
      break;
    } // end switch side
    loop+=1;
    if (loop>500) { // TILT because particle is trapped in a corner
      x1=coords_x[ex]; x2=coords_x[ex+1];
      y1=coords_y[ey]; y2=coords_y[ey+1];
      z1=coords_z[ez]; z2=coords_z[ez+1];
      numtilt+=1;
      if (numtilt==1e+6){
	cout<<"MORE THAN "<<numtilt<<" tilts, SOMETHING GOES WRONG" << endl; loop=0; numtilt=0; 
      }
      if (pos[0]-x1 < DomainSize[0]*10e-10) {pos[0]=x1+DomainSize[0]*10e-10;}
      if (x2-pos[0] < DomainSize[0]*10e-10) {pos[0]=x2-DomainSize[0]*10e-10;}
      if (pos[1]-y1 < DomainSize[1]*10e-10) {pos[1]=y1+DomainSize[1]*10e-10;}
      if (y2-pos[1] < DomainSize[1]*10e-10) {pos[1]=y2-DomainSize[1]*10e-10;}
      if (pos[2]-z1 < DomainSize[2]*10e-10) {pos[2]=z1+DomainSize[2]*10e-10;}
      if (z2-pos[2] < DomainSize[2]*10e-10) {pos[2]=z2-DomainSize[2]*10e-10;}
      convection.set(0.0,0.0,0.0);
    }
    
  } // while
  
  // ********************************* DISPERSION *********************************************
  
  // set type of random number
  switch(elements->RandomNumber) {
  case 1:
    dispersionData.rand[0]=-1.0+2.0*partrace->Random->fastdraw();
    dispersionData.rand[1]=-1.0+2.0*partrace->Random->fastdraw();
    dispersionData.rand[2]=-1.0+2.0*partrace->Random->fastdraw();
    break;
  case 2:
    dispersionData.rand[0]=partrace->Random->fastgauss();
    dispersionData.rand[1]=partrace->Random->fastgauss();
    dispersionData.rand[2]=partrace->Random->fastgauss();
    break;
  } // switch RandomNumber
  
  
  dt=partrace->RunTime.dt;    
  elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
  // cout<<"Dispersion" <<dispersionData.ds[0]<<" "<<dispersionData.ds[1]<<" "<<dispersionData.ds[2]<<endl;
  
  // loop over boundary passings/reflections
  while(1) {
    x1=coords_x[ex]; x2=coords_x[ex+1];
    y1=coords_y[ey]; y2=coords_y[ey+1];
    z1=coords_z[ez]; z2=coords_z[ez+1];
    // if (pos[0]<x1 || pos[0]>x2 || pos[1]<y1 || pos[1]>y2 || pos[2]<z1 || pos[2]>z2) {
    // cout << "Disp nan particle is outside " << pos[0] << " "<< pos[1]<<" " << pos[2]<<" " << x1 << " "<< x2<<" " << y1 <<" "<< y2 << " "<< z1<<" " << z2<<" "<< endl; }
    // //newpos=pos+dispersion;
    // cout<<dispersionData.ds[0]<<" "<<dispersionData.ds[1]<<" "<<dispersionData.ds[2]<<endl;
    newpos.setsum(pos, dispersionData.ds);
    // cout<<"NewPos Dispersion" << newpos[0]<<" "<< newpos[1]<<" "<< newpos[2]<<endl;
   
    // x-axis
    if(newpos[0]<x1) {
      diff=(x1-pos[0])/(newpos[0]-pos[0]);
      side=1;
    }
    else if(newpos[0]>x2) {
      side=2;
      diff=(x2-pos[0])/(newpos[0]-pos[0]);
    }
    else {
      diff=2.0;
      side=0;
    }
    // y-axis
    if(newpos[1]<y1) {
      diff2=(y1-pos[1])/(newpos[1]-pos[1]);
      if(diff2<diff) { diff=diff2; side=3; }
    }
    else if(newpos[1]>y2) {
      diff2=(y2-pos[1])/(newpos[1]-pos[1]);
      if(diff2<diff) { diff=diff2; side=4; }
    }
    // z-axis
    if(newpos[2]<z1) {
      diff2=(z1-pos[2])/(newpos[2]-pos[2]);
      if(diff2<diff) { diff=diff2; side=5; }
    }
    else if(newpos[2]>z2) {
      diff2=(z2-pos[2])/(newpos[2]-pos[2]);
      if(diff2<diff) { diff=diff2; side=6; }
    }
    if(side==0) {
      pos=newpos;
      return e;
    }
    loop+=1;
    
    if (loop>500) { // TILT because particle is trapped in a corner
      x1=coords_x[ex]; x2=coords_x[ex+1];
      y1=coords_y[ey]; y2=coords_y[ey+1];
      z1=coords_z[ez]; z2=coords_z[ez+1];
      numtilt+=1;
      if (numtilt==1e+6){
	cout<<"MORE THAN "<<numtilt<<" Dispersion tilts, SOMETHING GOES WRONG" << endl; loop=0; numtilt=0; 
      }
      if (pos[0]-x1 < DomainSize[0]*10e-10) {pos[0]=x1+DomainSize[0]*10e-10;}
      if (x2-pos[0] < DomainSize[0]*10e-10) {pos[0]=x2-DomainSize[0]*10e-10;}
      if (pos[1]-y1 < DomainSize[1]*10e-10) {pos[1]=y1+DomainSize[1]*10e-10;}
      if (y2-pos[1] < DomainSize[1]*10e-10) {pos[1]=y2-DomainSize[1]*10e-10;}
      if (pos[2]-z1 < DomainSize[2]*10e-10) {pos[2]=z1+DomainSize[2]*10e-10;}
      if (z2-pos[2] < DomainSize[2]*10e-10) {pos[2]=z2-DomainSize[2]*10e-10;}
      return e;
    }
    
    
    pos.addxsum(diff, dispersionData.ds); //move to one of the boundaries of the element,
    diff2=1.0-diff; //calculate the residual diff
    dispersionData.ds*=diff2;
    //if (elements->DispersionMethod==2) {
    //switch(elements->TimeSplitting) {
    //	case 1:
    //		dt*=diff2;
    //		break;
    //	case 2:
    dt*=diff2*diff2; // diff2*diff2 because for dispersion it is sqrt(dt) // if dt is used for recalculate Dispersion 
    //		break;
    //} // switch TimeSplitting
    //} // end if
    
    switch(elements->DispersionMethod) {
      
    case 1: // begin DispersionMethod off
      switch(side) {
      case 1: // left
	if(ex==0) reflectXYZ(dispersionData.ds);
	else if(InternalReflections->reflect(ex,ey,ez, ex-1,ey,ez))
	  reflectXYZ(dispersionData.ds);
	else { ex--;e--; }
	break;
      case 2: // right
	if(ex==NoElx_1) reflectXYZ(dispersionData.ds);
	else if(InternalReflections->reflect(ex,ey,ez, ex+1,ey,ez))
	  reflectXYZ(dispersionData.ds);
	else { ex++;e++; }
	break;
      case 3: // front
	if(ey==0) reflectXYZ(dispersionData.ds);
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey-1,ez))
	  reflectXYZ(dispersionData.ds);
	else { ey--;e-=NoElx; }
	break;
      case 4: // back
	if(ey==NoEly_1) reflectXYZ(dispersionData.ds);
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey+1,ez))
	  reflectXYZ(dispersionData.ds);
	else { ey++;e+=NoElx; }
	break;
      case 5: // bottom
	if(ez==0) reflectXYZ(dispersionData.ds);
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey,ez-1))
	  reflectXYZ(dispersionData.ds);
	else { ez--;e-=NoElxy; }
	break;
      case 6: // top
	if(ez==NoElz_1) reflectXYZ(dispersionData.ds);
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey,ez+1))
	  reflectXYZ(dispersionData.ds);
	else { ez++;e+=NoElxy; }
	break;
      } // switch side
      break; //end of DisperionMethod = 1
      
    case 2: // begin reflection method
      switch(side) {
      case 1: // left
	if(ex==0) {
	  //if (SecDispersiveNew==1) {elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	  //if (dispersionData.ds[0]<0.0) {reflectXYZ(dispersionData.ds);} }
	  //if (SecDispersiveNew==1) {reflectXYZ(dispersionData.ds);}
	  //else {reflectXYZ(dispersionData.ds);}
	  reflectXYZ(dispersionData.ds);
	  //reflectXYZ(dispersionData.rand); 
	} //
	else if(InternalReflections->reflect(ex,ey,ez, ex-1,ey,ez))
	  reflectXYZ(dispersionData.ds);
	else if (recalculateDispersion(dispersionData, e-1, e, pos, dt, dispersionData.particle->Retardation(e-1))) {     // always recalculate dispersion
	  if (reflectDispersion(dispersionData.dref_from[0], dispersionData.dref_to[0])) {
	    //	if (SecDispersiveNew==1 ) {elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	    //		if (dispersionData.ds[0]<0.0) {reflectXYZ(dispersionData.ds);} }
	    //	else {reflectXYZ(dispersionData.ds);}
	    reflectXYZ(dispersionData.ds);
	    //reflectXYZ(dispersionData.rand); 
	  }
	  else {
	    if (fabs(1-dispersionData.dref_from[0]/dispersionData.dref_to[0])>1e-8 || fabs(1-dispersionData.dref_from[1]/dispersionData.dref_to[1])>1e-8 || fabs(1-dispersionData.dref_from[2]/dispersionData.dref_to[2])>1e-8) {
	      dispersionData.ds=dispersionData.ds_to;
	      dispersionData.d_from=dispersionData.d_to;
	      dispersionData.dref_from=dispersionData.dref_to;
	      //dispersionData.rand[0]=dispersionData.rand_to[0]; dispersionData.rand[1]=dispersionData.rand_to[1]; dispersionData.rand[2]=dispersionData.rand_to[2];
	      dispersionData.v[0]=dispersionData.v_to[0]; dispersionData.v[1]=dispersionData.v_to[1]; dispersionData.v[2]=dispersionData.v_to[2];
	      dispersionData.vNorm=dispersionData.vNorm_to;
	      if (dispersionData.ds[0]>0.0) {reflectXYZ(dispersionData.ds); cout<<"reflectexception"<<endl;}
	    }
	    ex--;e--;
	  }
	}
	break;
	
      case 2: // right
	if(ex==NoElx_1) {
	  //if (SecDispersiveNew==1) {elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	  //if (dispersionData.ds[0]>0.0) {reflectXYZ(dispersionData.ds);} }
	  //if (SecDispersiveNew==1) {reflectXYZ(dispersionData.ds);}
	  //else {reflectXYZ(dispersionData.ds);}
	  reflectXYZ(dispersionData.ds);
	  //reflectXYZ(dispersionData.rand); 
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex+1,ey,ez))
	  reflectXYZ(dispersionData.ds);
	else if (recalculateDispersion(dispersionData, e+1, e, pos, dt, dispersionData.particle->Retardation(e+1))) {     // always recalculate dispersion
	  if (reflectDispersion(dispersionData.dref_from[0], dispersionData.dref_to[0])) {
	    //if (SecDispersiveNew==1 ) {elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	    //	if (dispersionData.ds[0]>0.0) {reflectXYZ(dispersionData.ds);}	}
	    //else {reflectXYZ(dispersionData.ds);}
	    reflectXYZ(dispersionData.ds);
	    //reflectXYZ(dispersionData.rand); 
	  }
	  else {
	    if (fabs(1-dispersionData.dref_from[0]/dispersionData.dref_to[0])>1e-8 || fabs(1-dispersionData.dref_from[1]/dispersionData.dref_to[1])>1e-8 || fabs(1-dispersionData.dref_from[2]/dispersionData.dref_to[2])>1e-8) {
	      dispersionData.ds=dispersionData.ds_to;
	      dispersionData.d_from=dispersionData.d_to;
	      dispersionData.dref_from=dispersionData.dref_to; 
	      //dispersionData.rand[0]=dispersionData.rand_to[0]; dispersionData.rand[1]=dispersionData.rand_to[1]; dispersionData.rand[2]=dispersionData.rand_to[2]; 
	      dispersionData.v[0]=dispersionData.v_to[0]; dispersionData.v[1]=dispersionData.v_to[1]; dispersionData.v[2]=dispersionData.v_to[2];
	      dispersionData.vNorm=dispersionData.vNorm_to;
	      if (dispersionData.ds[0]<0.0) {reflectXYZ(dispersionData.ds);cout<<"reflectexception"<<endl;}
	    }
	    ex++;e++;
	  }
	}
	break;
	
      case 3: // front
	if(ey==0) {
	  //if (SecDispersiveNew==1) {	elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	  //if (dispersionData.ds[1]<0.0) {reflectXYZ(dispersionData.ds);}	}	
	  //if (SecDispersiveNew==1) {reflectXYZ(dispersionData.ds);}
	  //else {reflectXYZ(dispersionData.ds);}
	  reflectXYZ(dispersionData.ds);
	  //reflectXYZ(dispersionData.rand); 
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey-1,ez))
	  reflectXYZ(dispersionData.ds);
	else if (recalculateDispersion(dispersionData, e-NoElx, e, pos, dt, dispersionData.particle->Retardation(e-NoElx))) {     // always recalculate dispersion
	  if (reflectDispersion(dispersionData.dref_from[1], dispersionData.dref_to[1])) {
	    //if (SecDispersiveNew==1 ) {elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	    //	if (dispersionData.ds[1]<0.0) {reflectXYZ(dispersionData.ds);}		}
	    //else {reflectXYZ(dispersionData.ds);}
	    reflectXYZ(dispersionData.ds);
	    //reflectXYZ(dispersionData.rand); 
	  }
	  else {
	    if (fabs(1-dispersionData.dref_from[0]/dispersionData.dref_to[0])>1e-8 || fabs(1-dispersionData.dref_from[1]/dispersionData.dref_to[1])>1e-8 || fabs(1-dispersionData.dref_from[2]/dispersionData.dref_to[2])>1e-8) {
	      dispersionData.ds=dispersionData.ds_to;
	      dispersionData.d_from=dispersionData.d_to;
	      dispersionData.dref_from=dispersionData.dref_to; 
	      //dispersionData.rand[0]=dispersionData.rand_to[0]; dispersionData.rand[1]=dispersionData.rand_to[1]; dispersionData.rand[2]=dispersionData.rand_to[2]; 
	      dispersionData.v[0]=dispersionData.v_to[0]; dispersionData.v[1]=dispersionData.v_to[1]; dispersionData.v[2]=dispersionData.v_to[2];
	      dispersionData.vNorm=dispersionData.vNorm_to;
	      if (dispersionData.ds[1]>0.0) {reflectXYZ(dispersionData.ds);cout<<"reflectexception"<<endl;}
	    }
	    ey--;e-=NoElx;
	  }
	}
	break;
	
      case 4: // back
	if(ey==NoEly_1) {
	  //if (SecDispersiveNew==1) {	elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	  //if (dispersionData.ds[1]>0.0) {reflectXYZ(dispersionData.ds);}		}
	  //if (SecDispersiveNew==1) {reflectXYZ(dispersionData.ds);}
	  //else {reflectXYZ(dispersionData.ds);}
	  reflectXYZ(dispersionData.ds);
	  //reflectXYZ(dispersionData.rand); 
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey+1,ez))
	  reflectXYZ(dispersionData.ds);
	else if (recalculateDispersion(dispersionData, e+NoElx, e, pos, dt, dispersionData.particle->Retardation(e+NoElx))) {     // always recalculate dispersion
	  if (reflectDispersion(dispersionData.dref_from[1], dispersionData.dref_to[1])) {
	    //if (SecDispersiveNew==1 ) {	elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	    //	if (dispersionData.ds[1]>0.0) {reflectXYZ(dispersionData.ds);}	}	
	    //else {reflectXYZ(dispersionData.ds);}
	    reflectXYZ(dispersionData.ds);
	    //reflectXYZ(dispersionData.rand); 
	  }
	  else {
	    if (fabs(1-dispersionData.dref_from[0]/dispersionData.dref_to[0])>1e-8 || fabs(1-dispersionData.dref_from[1]/dispersionData.dref_to[1])>1e-8 || fabs(1-dispersionData.dref_from[2]/dispersionData.dref_to[2])>1e-8) {
	      dispersionData.ds=dispersionData.ds_to;
	      dispersionData.d_from=dispersionData.d_to;
	      dispersionData.dref_from=dispersionData.dref_to;
	      //dispersionData.rand[0]=dispersionData.rand_to[0]; dispersionData.rand[1]=dispersionData.rand_to[1]; dispersionData.rand[2]=dispersionData.rand_to[2]; 
	      dispersionData.v[0]=dispersionData.v_to[0]; dispersionData.v[1]=dispersionData.v_to[1]; dispersionData.v[2]=dispersionData.v_to[2];
	      dispersionData.vNorm=dispersionData.vNorm_to;
	      if (dispersionData.ds[1]<0.0) {reflectXYZ(dispersionData.ds);cout<<"reflectexception"<<endl;}
	    }
	    ey++;e+=NoElx;
	  }
	}
	break;
	
      case 5: // bottom
	if(ez==0) {
	  //if (SecDispersiveNew==1) {	elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	  //if (dispersionData.ds[2]<0.0) {reflectXYZ(dispersionData.ds);}	}	
	  //if (SecDispersiveNew==1) {reflectXYZ(dispersionData.ds);}
	  //else {reflectXYZ(dispersionData.ds);}
	  reflectXYZ(dispersionData.ds);
	  //reflectXYZ(dispersionData.rand); 
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey,ez-1))
	  reflectXYZ(dispersionData.ds);
	else if (recalculateDispersion(dispersionData, e-NoElxy, e, pos, dt, dispersionData.particle->Retardation(e-NoElxy))) {     // always recalculate dispersion
	  if (reflectDispersion(dispersionData.dref_from[2], dispersionData.dref_to[2]) ) {
	    //if (SecDispersiveNew==1 ) {	elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	    //	if (dispersionData.ds[2]<0.0) {reflectXYZ(dispersionData.ds);}	}	
	    //else {reflectXYZ(dispersionData.ds); }
	    reflectXYZ(dispersionData.ds);
	    //reflectXYZ(dispersionData.rand); 
	  }
	  else {
	    if (fabs(1-dispersionData.dref_from[0]/dispersionData.dref_to[0])>1e-8 || fabs(1-dispersionData.dref_from[1]/dispersionData.dref_to[1])>1e-8 || fabs(1-dispersionData.dref_from[2]/dispersionData.dref_to[2])>1e-8) {
	      dispersionData.ds=dispersionData.ds_to;
	      dispersionData.d_from=dispersionData.d_to;
	      dispersionData.dref_from=dispersionData.dref_to; 
	      //dispersionData.rand[0]=dispersionData.rand_to[0]; dispersionData.rand[1]=dispersionData.rand_to[1]; dispersionData.rand[2]=dispersionData.rand_to[2]; 
	      dispersionData.v[0]=dispersionData.v_to[0]; dispersionData.v[1]=dispersionData.v_to[1]; dispersionData.v[2]=dispersionData.v_to[2];
	      dispersionData.vNorm=dispersionData.vNorm_to;
	      if (dispersionData.ds[2]>0.0) {reflectXYZ(dispersionData.ds);cout<<"reflectexception"<<endl;} 
	    }
	    ez--;e-=NoElxy;
	  }
	}
	break;
	
      case 6: // top
	if(ez==NoElz_1) {
	  //if (SecDispersiveNew==1) {	elements->Dispersion(dt, e, pos, dispersionData, dispersionData.particle->Retardation(e)); //MB
	  //if (dispersionData.ds[2]>0.0) {reflectXYZ(dispersionData.ds);}		}
	  //if (SecDispersiveNew==1) {reflectXYZ(dispersionData.ds);}
	  //else {reflectXYZ(dispersionData.ds);}
	  reflectXYZ(dispersionData.ds);
	  //reflectXYZ(dispersionData.rand);
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey,ez+1))
	  reflectXYZ(dispersionData.ds);
	else if (recalculateDispersion(dispersionData, e+NoElxy, e, pos, dt, dispersionData.particle->Retardation(e+NoElxy))) {     // always recalculate dispersion
	  if (reflectDispersion(dispersionData.dref_from[2], dispersionData.dref_to[2])) {
	    reflectXYZ(dispersionData.ds);
	    //reflectXYZ(dispersionData.rand); 
	  }
	  else {
	    if (fabs(1-dispersionData.dref_from[0]/dispersionData.dref_to[0])>1e-8 || fabs(1-dispersionData.dref_from[1]/dispersionData.dref_to[1])>1e-8 || fabs(1-dispersionData.dref_from[2]/dispersionData.dref_to[2])>1e-8) {
	      dispersionData.ds=dispersionData.ds_to;
	      dispersionData.d_from=dispersionData.d_to;
	      dispersionData.dref_from=dispersionData.dref_to; 
	      //dispersionData.rand[0]=dispersionData.rand_to[0]; dispersionData.rand[1]=dispersionData.rand_to[1]; dispersionData.rand[2]=dispersionData.rand_to[2]; 
	      dispersionData.v[0]=dispersionData.v_to[0]; dispersionData.v[1]=dispersionData.v_to[1]; dispersionData.v[2]=dispersionData.v_to[2];
	      dispersionData.vNorm=dispersionData.vNorm_to;
	      if (dispersionData.ds[2]<0.0) {reflectXYZ(dispersionData.ds);cout<<"reflectexception"<<endl;}
	    }
	    ez++;e+=NoElxy;
	  }
	}
	break;
      } // switch side
      break; // end of reflection method
      
    case 3: // begin of interpolation method
      switch(side) {
      case 1: // left
	if(ex==0) {
	  if(BoundaryXmin->type(ey,ez)==1) {
	    reflectXYZ(dispersionData.ds);
	  }
	  else if(BoundaryXmin->type(ey,ez)==5) {
	    pos[0] = coords_x[NoElx];
	    ex = NoElx_1;
	    e = e + NoElx_1;
	  }
	} //
	else if(InternalReflections->reflect(ex,ey,ez, ex-1,ey,ez))
	  reflectXYZ(dispersionData.ds);
	else {
	  if (SecDispersiveNew == 2) { 
	    recalculateDispersion(dispersionData, e-1, e, pos, dt, dispersionData.particle->Retardation(e-1));
	    dispersionData.ds=dispersionData.ds_to;
	    dispersionData.d_from=dispersionData.d_to; }
	  ex--;e--;
	}
	break;
	
      case 2: // right
	if(ex==NoElx_1) {
	  if(BoundaryXmax->type(ey,ez)==1) {
	    reflectXYZ(dispersionData.ds);
	  }
	  else if(BoundaryXmax->type(ey,ez)==5) {
	    pos[0] = coords_x[0];
	    ex = 0;
	    e = e - NoElx_1;
	  }
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex+1,ey,ez))
	  reflectXYZ(dispersionData.ds);
	else {
	  if (SecDispersiveNew == 2) { 
	    recalculateDispersion(dispersionData, e+1, e, pos, dt, dispersionData.particle->Retardation(e+1));
	    dispersionData.ds=dispersionData.ds_to;
	    dispersionData.d_from=dispersionData.d_to; }
	  ex++;e++;
	}
	break;
	
      case 3: // front
	if(ey==0) {
	  if(BoundaryYmin->type(ex,ez)==1) {
	    reflectXYZ(dispersionData.ds);
	  }
	  else if(BoundaryYmin->type(ex,ez)==5) {
	    pos[1] = coords_y[NoEly];
	    ey = NoEly_1;
	    e = e + NoElx * NoEly_1;
	  }
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey-1,ez))
	  reflectXYZ(dispersionData.ds);
	else {
	  if (SecDispersiveNew == 2) { 
	    recalculateDispersion(dispersionData, e-NoElx, e, pos, dt, dispersionData.particle->Retardation(e-NoElx));
	    dispersionData.ds=dispersionData.ds_to;
	    dispersionData.d_from=dispersionData.d_to; }
	  ey--;e-=NoElx;
	}
	break;
	
      case 4: // back
	if(ey==NoEly_1) {
	  if(BoundaryYmax->type(ex,ez)==1) {
	    reflectXYZ(dispersionData.ds);
	  }
	  else if(BoundaryYmax->type(ex,ez)==5) {
	    pos[1] = coords_y[0];
	    ey = 0;
	    e = e - NoElx * NoEly_1;
	  };
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey+1,ez))
	  reflectXYZ(dispersionData.ds);
	else {
	  if (SecDispersiveNew == 2) { 
	    recalculateDispersion(dispersionData, e+NoElx, e, pos, dt, dispersionData.particle->Retardation(e+NoElx));
	    dispersionData.ds=dispersionData.ds_to;
	    dispersionData.d_from=dispersionData.d_to; }
	  ey++;e+=NoElx;
	}
	break;
	
      case 5: // bottom
	if(ez==0) {
	  reflectXYZ(dispersionData.ds);
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey,ez-1))
	  reflectXYZ(dispersionData.ds);
	else {
	  if (SecDispersiveNew == 2) { 
	    recalculateDispersion(dispersionData, e-NoElxy, e, pos, dt, dispersionData.particle->Retardation(e-NoElxy));
	    dispersionData.ds=dispersionData.ds_to;
	    dispersionData.d_from=dispersionData.d_to; }
	  ez--;e-=NoElxy;
	}
	break;
	
      case 6: // top
	if(ez==NoElz_1) {
	  reflectXYZ(dispersionData.ds);
	}
	else if(InternalReflections->reflect(ex,ey,ez, ex,ey,ez+1))
	  reflectXYZ(dispersionData.ds);
	else {
	  if (SecDispersiveNew == 2) { 
	    recalculateDispersion(dispersionData, e+NoElxy, e, pos, dt, dispersionData.particle->Retardation(e+NoElxy));
	    dispersionData.ds=dispersionData.ds_to;
	    dispersionData.d_from=dispersionData.d_to; }
	  ez++;e+=NoElxy;
	}
	break;
      } // switch side
      
      break; //end of interpolation method
    } //end switch DispersionMethod
    
  } // while move particle, exit when particle is not on any boundary
}


bool ElementsClass::reflectDispersion(double d1, double d2) //MB
{
  static double p1;
  bool a=0; //answer
  //ElementsClass *elements=partrace->elements;	
  //	switch(elements->ReflectionMethod) {
  //		case 1:
  // Hoteit et al 2002
  //			p1=d1/(d1+d2); // (after Hoteit et al, 2002): p1 is the probability that particle enters region 1, (1-p1 for region 2)
  //			a = partrace->Random->fastdraw() < p1 ; // p1 is therefore the probability for reflection
  //		break;
  //		case 2:		
  //			// Lim 2006
  //			p1=d1/(d1+d2); // p1 is the probability that particle enters region 1, (1-p1 for region 2)
  //			a = partrace->Random->fastdraw() < p1 ; // p1 is therefore the probability for reflection
  //		break;
  //		case 3:
  // Michel, one-sided reflection coefficient, higher accuracy for same time step size
  p1=d2/d1; // 0.1 -- > 90 % reflection     1.0 --> 0 % reflection
  //if (p1!=1) {cout<<p1<<endl;}
  //if (p1<0.99 || p1>1.01) {cout<<"p1 "<<p1<<endl;}
  //if (p1<1.0) {p1*=0.45;}
  a = partrace->Random->fastdraw() > p1 ; // from one side p1 is always > 1.0, therefore no reflection
  //		break;
  //	} //switch
  return a;
}


void ElementsClass::calculateConvection(vector3& convection, //MB
					long e, vector3& pos, double dt, double Retardation)

{ // Two things to do in future: i) how to implement Retardation into gradient of disersion tensor and ii) interpolation of disp.parameters between elements if not constant
  static long ex,ey,ez;
  int finitevolindex=0;
  static long nodes[8];
  //static double interpol[24];
  //static double interpolWC[8];
  static vector3 xcoor1, v1, v2, v1Tri, v2Trix, v2Triy, v2Triz, posnew, pos2, pos2dx, v2TriTemp;
  static double x1,x2,y1,y2,z1,z2;
  v1[0]=0.0; v1[1]=0.0; v1[2]=0.0; 	v2[0]=0.0; v2[1]=0.0; v2[2]=0.0;
  static int ni[8]; // nodeindex
  static double wcFiniteVol;
  static double vNorm, vNorm1, vNorm2, vNorm3;
  double ddx=0.0, ddy=0.0, ddz=0.0, ddx1=0.0, ddy1=0.0, ddz1=0.0, ddx2=0.0, ddy2=0.0, ddz2=0.0;
  static double Dx, Dy, Dz, Fx, Fy, Fz;
  static double* dispersion_parameter;
  static double TriFactor1[8], TriFactor2[8], TriFactor3[8], TriFactor4[8];
  static double wc1; // , wctemp;
  static vector3 wc2;
  static double dWCx, dWCy, dWCz, Diffu1;
  dWCx=0.0; dWCy=0.0; dWCz=0.0; Diffu1=0.0;
  static vector3 Diffu2;
  Diffu2[0]=0.0; Diffu2[1]=0.0; Diffu2[2]=0.0;
  const double wc_exponent=7.0/3.0;
  const double dx=1e-7;
  static double aLT=0.0;
  static double D11, D12, D13, D21, D22, D23, D31, D32, D33;
  ElementsClass *elements=partrace->elements;

  // get element dependant information
  elements->GetElementXYZ(e,ex,ey,ez);

  // ********************************** Velocity determination ******************************************
  switch(elements->VelType) {
  case 1: // on elements
    v1.set(elements->velocity->get_line(e));
    break;
  case 2: // on nodes
    /* old: elements->GetNodesOfElement(e, ex,ey,ez, nodes); // node number of element
       elements->GetStartCoordinates(nodes, ex,ey,ez, xcoor1); // coordinates of first node
       elements->initInterpolationFunction(nodes); // skipped for rectilinear grids
       elements->getInterpolationFunction(e, ex,ey,ez, nodes, interpol); // nodes
       elements->getInterpolationFunctionWC(e, ex,ey,ez, nodes, interpolWC);
       elements->interpolateVelocity(xcoor1, pos, interpol, v1);*/
    elements->GetNodesOfElement(e, ex,ey,ez, nodes);
    elements->GetDxyzFxyz(ex, ey, ez, pos, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
    elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor1); //
    elements->InterpolateVelocityTriLinear(ex, ey, ez, TriFactor1, ni, nodes, v1);			
    break;
  case 3: // finite volume
    finitevolindex=6*(ez+ey*NoElz+ex*NoElyz);
    wcFiniteVol=elements->watercontent->get(e);
    switch(elements->ConvectionMethod) {
    case 0: // No Velocity Integration
      v1.set(elements->finitevolvelocity[finitevolindex  ]*pos[0]+elements->finitevolvelocity[finitevolindex+1],
	     elements->finitevolvelocity[finitevolindex+2]*pos[1]+elements->finitevolvelocity[finitevolindex+3],
	     elements->finitevolvelocity[finitevolindex+4]*pos[2]+elements->finitevolvelocity[finitevolindex+5]);
      v1/=wcFiniteVol;
      v2.set(elements->finitevolvelocity[finitevolindex  ]*(pos[0]+dx)+elements->finitevolvelocity[finitevolindex+1],
	     elements->finitevolvelocity[finitevolindex+2]*(pos[1]+dx)+elements->finitevolvelocity[finitevolindex+3],
	     elements->finitevolvelocity[finitevolindex+4]*(pos[2]+dx)+elements->finitevolvelocity[finitevolindex+5]);
      v2/=wcFiniteVol;
      break;
    case 1: // Velocity Integration, Pollock
      x1=coords_x[ex]; x2=coords_x[ex+1];
      y1=coords_y[ey]; y2=coords_y[ey+1];
      z1=coords_z[ez]; z2=coords_z[ez+1];
      //position pos2 after convection
      //
      // to be implemented --> error when one direction with zero velocity
      pos2.set((elements->finitevolvelocityP[finitevolindex  ]*(pos[0]-x1)+elements->finitevolvelocityP[finitevolindex+1])/elements->finitevolvelocityP[finitevolindex  ] * exp(elements->finitevolvelocityP[finitevolindex  ]*dt)-elements->finitevolvelocityP[finitevolindex+1]/elements->finitevolvelocityP[finitevolindex  ] + x1,
	       (elements->finitevolvelocityP[finitevolindex+2]*(pos[1]-y1)+elements->finitevolvelocityP[finitevolindex+3])/elements->finitevolvelocityP[finitevolindex+2] * exp(elements->finitevolvelocityP[finitevolindex+2]*dt)-elements->finitevolvelocityP[finitevolindex+3]/elements->finitevolvelocityP[finitevolindex+2] + y1,
	       (elements->finitevolvelocityP[finitevolindex+4]*(pos[2]-z1)+elements->finitevolvelocityP[finitevolindex+5])/elements->finitevolvelocityP[finitevolindex+4] * exp(elements->finitevolvelocityP[finitevolindex+4]*dt)-elements->finitevolvelocityP[finitevolindex+5]/elements->finitevolvelocityP[finitevolindex+4] + z1);
      pos2dx.set((elements->finitevolvelocityP[finitevolindex  ]*(pos[0]+dx-x1)+elements->finitevolvelocityP[finitevolindex+1])/elements->finitevolvelocityP[finitevolindex  ] * exp(elements->finitevolvelocityP[finitevolindex  ]*dt)-elements->finitevolvelocityP[finitevolindex+1]/elements->finitevolvelocityP[finitevolindex  ] + x1,
		 (elements->finitevolvelocityP[finitevolindex+2]*(pos[1]+dx-y1)+elements->finitevolvelocityP[finitevolindex+3])/elements->finitevolvelocityP[finitevolindex+2] * exp(elements->finitevolvelocityP[finitevolindex+2]*dt)-elements->finitevolvelocityP[finitevolindex+3]/elements->finitevolvelocityP[finitevolindex+2] + y1,
		 (elements->finitevolvelocityP[finitevolindex+4]*(pos[2]+dx-z1)+elements->finitevolvelocityP[finitevolindex+5])/elements->finitevolvelocityP[finitevolindex+4] * exp(elements->finitevolvelocityP[finitevolindex+4]*dt)-elements->finitevolvelocityP[finitevolindex+5]/elements->finitevolvelocityP[finitevolindex+4] + z1);
      v1.set((pos2[0]-pos[0])/dt,
	     (pos2[1]-pos[1])/dt,
	     (pos2[2]-pos[2])/dt);
      v2.set((pos2dx[0]-(pos[0]+dx))/dt,
	     (pos2dx[1]-(pos[1]+dx))/dt,
	     (pos2dx[2]-(pos[2]+dx))/dt);
      //Comment, finitevolvelocityP already devided by watercontent
      break;
    } // switch ConvectionMethod
    break; // case 3: finite volume
  } // end switch VelType

  // ********************************** Drift determination ******************************************
  // divergence of D, component of convective displacement, e.g. Kitanidis et al. 1994
  // D after Lichtner et al. 2002
  dispersion_parameter=elements->soilprop->get_line(e);
  aLT = fabs(dispersion_parameter[0]-dispersion_parameter[1]);
  switch(elements->DispersionMethod) {
  case 1: //
    ddx1=0.0; ddy1=0.0; ddz1=0.0;
    ddx2=0.0; ddy2=0.0; ddz2=0.0;
    ddx=0.0;
    ddy=0.0;
    ddz=0.0;
    break;
  case 2: // reflection method ///
    // calculation of divergence using differences
    // v 1
    vNorm=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
    if (vNorm<=vMin) {
      ddx=0.0; ddy=0.0; ddz=0.0;
    }
    else {
      ddx1=(dispersion_parameter[1]*vNorm + aLT*v1[0]*v1[0]/vNorm    + aLT*v1[0]*v1[1]/vNorm   + aLT*v1[0]*v1[2]/vNorm) * DispWeight[0];
      ddy1=(dispersion_parameter[1]*vNorm + aLT*v1[1]*v1[1]/vNorm    + aLT*v1[1]*v1[0]/vNorm   + aLT*v1[1]*v1[2]/vNorm) * DispWeight[1];
      ddz1=(dispersion_parameter[1]*vNorm + aLT*v1[2]*v1[2]/vNorm    + aLT*v1[2]*v1[0]/vNorm   + aLT*v1[2]*v1[1]/vNorm) * DispWeight[2];
      // v 2
      vNorm1=sqrt(v2[0]*v2[0]+v1[1]*v1[1]+v1[2]*v1[2]); // vNorm dx in x direction
      vNorm2=sqrt(v1[0]*v1[0]+v2[1]*v2[1]+v1[2]*v1[2]); // vNorm dx in y direction
      vNorm3=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v2[2]*v2[2]); // vNorm dx in z direction
      ddx2=(dispersion_parameter[1]*vNorm1 + aLT*v2[0]*v2[0]/vNorm1    + aLT*v1[0]*v2[1]/vNorm2   + aLT*v1[0]*v2[2]/vNorm3) * DispWeight[0];
      ddy2=(dispersion_parameter[1]*vNorm2 + aLT*v2[1]*v2[1]/vNorm2    + aLT*v1[1]*v2[0]/vNorm1   + aLT*v1[1]*v2[2]/vNorm3) * DispWeight[1];
      ddz2=(dispersion_parameter[1]*vNorm3 + aLT*v2[2]*v2[2]/vNorm3    + aLT*v1[2]*v2[0]/vNorm1   + aLT*v1[2]*v2[1]/vNorm2) * DispWeight[2];
      if (fabs(ddx2-ddx1)<1e-80) {
	ddx=0.0; ddy=0.0; ddz=0.0; 
      }
      else {
	ddx=(ddx2-ddx1)/dx;
	ddy=(ddy2-ddy1)/dx;
	ddz=(ddz2-ddz1)/dx;
      }
    }
    break;
  case 3: // interpolation method //
    ////////////////////////////////// trilinear interpolation of velocity
    elements->GetNodesOfElement(e, ex,ey,ez, nodes);
    elements->GetDxyzFxyz(ex, ey, ez, pos, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
    elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor1); //
    elements->InterpolateVelocityTriLinear(ex, ey, ez, TriFactor1, ni, nodes, v1Tri);			
    //cout<<"F and TriFactor "<<Fx<<" "<<Fy<<" "<<Fz<<" "<<TriFactor1[0]<<" "<<TriFactor1[1]<<" "<<TriFactor1[2]<<" "<<TriFactor1[3]<<" "<<TriFactor1[4]<<" "<<TriFactor1[5]<<" "<<TriFactor1[6]<<" "<<TriFactor1[7]<<endl;
    // v2Tri[0]
    posnew=pos;
    posnew[0]=pos[0]+dx;
    elements->GetDxyzFxyz(ex, ey, ez, posnew, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
    elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor2); //
    elements->InterpolateVelocityTriLinearDX(TriFactor2, ni, nodes, v2Trix);
    // v2Tri[1]
    posnew=pos;
    posnew[1]=pos[1]+dx;
    elements->GetDxyzFxyz(ex, ey, ez, posnew, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
    elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor3); //
    elements->InterpolateVelocityTriLinearDX(TriFactor3, ni, nodes, v2Triy);
    // v2Tri[2]
    posnew=pos;
    posnew[2]=pos[2]+dx;
    elements->GetDxyzFxyz(ex, ey, ez, posnew, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
    elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor4); //
    elements->InterpolateVelocityTriLinearDX(TriFactor4, ni, nodes, v2Triz);
    //cout<<v1Tri[0]<<" "<<v1Tri[1]<<" "<<v1Tri[2]<<" "<<v2Triz[0]<<" "<<v2Triz[1]<<" "<<v2Triz[2]<<endl;
    //cout<<v2Trix[0]<<" "<<v2Trix[1]<<" "<<v2Trix[2]<<" "<<v2Triy[0]<<" "<<v2Triy[1]<<" "<<v2Triy[2]<<endl;
						
    // Gradient of water content / porosity 
    //////////////////////////////////////// trilinear interpolation of water content for water content and diffusion drift
    elements->InterpolateWCTriLinear(TriFactor1, ni, nodes, wc1);			
    // wc2[0]
    elements->InterpolateWCTriLinear(TriFactor2, ni, nodes, wc2[0]);
    // wc2[1]
    elements->InterpolateWCTriLinear(TriFactor3, ni, nodes, wc2[1]);
    // wc2[2]
    elements->InterpolateWCTriLinear(TriFactor4, ni, nodes, wc2[2]);
    //cout<<"wc1 "<<wc1<<endl;

    //////////////////////////////// calculate diffusion for diffusion drift
    switch(elements->DiffusionType) {
    case 1: // simple
      Diffu1=dispersion_parameter[2];
      Diffu2[0]=Diffu1; Diffu2[1]=Diffu1; Diffu2[2]=Diffu1;
      break;
    case 2: // water content based
      // D = Di * wc**(7/3) / porosity**2
      Diffu1=dispersion_parameter[2] * pow(wc1, wc_exponent) /
	(dispersion_parameter[3] * dispersion_parameter[3]);
      Diffu2[0]=dispersion_parameter[2] * pow(wc2[0], wc_exponent) /
	(dispersion_parameter[3] * dispersion_parameter[3]);
      Diffu2[1]=dispersion_parameter[2] * pow(wc2[1], wc_exponent) /
	(dispersion_parameter[3] * dispersion_parameter[3]);
      Diffu2[2]=dispersion_parameter[2] * pow(wc2[2], wc_exponent) /
	(dispersion_parameter[3] * dispersion_parameter[3]);
      break;
    }
    //cout<<"Diffu1 "<<Diffu1<<" "<<Diffu2[0]<<" "<<endl;
    // calculation of divergence using differences
    // v 1
    vNorm=sqrt(v1Tri[0]*v1Tri[0]+v1Tri[1]*v1Tri[1]+v1Tri[2]*v1Tri[2]); 
    if (vNorm<=vMin) {
      ddx=0.0; ddy=0.0; ddz=0.0;
      D11=Diffu1;
      D12=0.0 ;
      D13=0.0 ;
      D21=0.0 ;
      D22=Diffu1 ;
      D23=0.0 ;
      D31=0.0 ;
      D32=0.0 ;
      D33=Diffu1 ;
    }
    else {
      D11=Diffu1 + dispersion_parameter[1]*vNorm + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[0]*v1Tri[0]/vNorm;
      D12=fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[0]*v1Tri[1]/vNorm ;
      D13=fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[0]*v1Tri[2]/vNorm;
      D21=fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[1]*v1Tri[0]/vNorm;
      D22=Diffu1 + dispersion_parameter[1]*vNorm + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[1]*v1Tri[1]/vNorm;
      D23=fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[1]*v1Tri[2]/vNorm;
      D31=fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[2]*v1Tri[0]/vNorm;
      D32=fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[2]*v1Tri[1]/vNorm;
      D33=Diffu1 + dispersion_parameter[1]*vNorm + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v1Tri[2]*v1Tri[2]/vNorm;
    }
    ddx1=(D11 + D12 + D13) * DispWeight[0];
    ddy1=(D22 + D21 + D23) * DispWeight[1];
    ddz1=(D33 + D31 + D32) * DispWeight[2];
    //cout<<Diffu1<<" "<<vNorm<<" "<<v1Tri[0]<<" "<<endl;
    // v2Trix is v2 in x direction at a dx step
    vNorm1=sqrt(v2Trix[0]*v2Trix[0]+v2Trix[1]*v2Trix[1]+v2Trix[2]*v2Trix[2]); // vNorm dx in x direction
    vNorm2=sqrt(v2Triy[0]*v2Triy[0]+v2Triy[1]*v2Triy[1]+v2Triy[2]*v2Triy[2]); // vNorm dx in y direction
    vNorm3=sqrt(v2Triz[0]*v2Triz[0]+v2Triz[1]*v2Triz[1]+v2Triz[2]*v2Triz[2]); // vNorm dx in z direction
    // ddx2 Dispersion 
    if (vNorm1>vMin) {
      ddx2=(Diffu2[0] + dispersion_parameter[1]*vNorm1 + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Trix[0]*v2Trix[0]/vNorm1    + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Triy[0]*v2Triy[1]/vNorm2   + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Triz[0]*v2Triz[2]/vNorm3) * DispWeight[0];
    }
    else {
      ddx2=Diffu2[0];
    }
    if (vNorm2>vMin) {
      ddy2=(Diffu2[1] + dispersion_parameter[1]*vNorm2 + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Triy[1]*v2Triy[1]/vNorm2    + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Trix[1]*v2Trix[0]/vNorm1   + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Triz[1]*v2Triz[2]/vNorm3) * DispWeight[1];
    }
    else {
      ddy2=Diffu2[1];
    }
    if (vNorm3>vMin) {
      ddz2=(Diffu2[2] + dispersion_parameter[1]*vNorm3 + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Triz[2]*v2Triz[2]/vNorm3    + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Trix[2]*v2Trix[0]/vNorm1   + fabs(dispersion_parameter[0]-dispersion_parameter[1])*v2Triy[2]*v2Triy[1]/vNorm2) * DispWeight[2];
    }
    else {
      ddz2=Diffu2[2];
    }
    //cout<<"ddx2 "<<ddx2<<endl;
    if (fabs(ddx2-ddx1)<1e-80) {
      ddx=0.0; ddy=0.0; ddz=0.0; 
    }
    else {
      ddx=(ddx2-ddx1)/dx;
      ddy=(ddy2-ddy1)/dx;
      ddz=(ddz2-ddz1)/dx;
      //cout<<"ddx2 ddx1 "<<ddx2<<" "<<ddx1<<endl;
    }
    // calculate water content drift
    // add drift of water content
    if (wc1>vMin) {
      dWCx = 1/wc1*(D11*(wc2[0]-wc1) + D12*(wc2[1]-wc1) + D13*(wc2[2]-wc1))/dx;
      dWCy = 1/wc1*(D21*(wc2[0]-wc1) + D22*(wc2[1]-wc1) + D23*(wc2[2]-wc1))/dx;
      dWCz = 1/wc1*(D31*(wc2[0]-wc1) + D32*(wc2[1]-wc1) + D33*(wc2[2]-wc1))/dx;
      ddx+= dWCx;
      ddy+= dWCy;
      ddz+= dWCz; 
    }

    break; // end Interpolation method, very long part as drift correction for dispersion, diffusion and water content
  } // switch DispersionMethod

    // convection
  //cout<<"v1 and ddx"<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<" "<<ddx<<" "<<endl;
  convection[0]=(v1[0]+ddx)*Retardation*dt;
  convection[1]=(v1[1]+ddy)*Retardation*dt;
  convection[2]=(v1[2]+ddz)*Retardation*dt;
}

bool ElementsClass::recalculateDispersion(DispersionClass& disp, //MB
					  long e_to, long e_from, vector3& pos, double dt, double Retardation)

{
  // compare the dispersion from the starting location with the dispersion
  // in the new element e for computing the probability for the reflection
  // of the dispersive vector.
  //static double* dispersion_parameter_to;
  double* dispersion_parameter_to;
  double* dispersion_parameter_from;
  const double wc_exponent=7.0/3.0;
  static long ex,ey,ez, ex_from,ey_from,ez_from, ex_to, ey_to, ez_to;
  int finitevolindex=0;
  static long nodes[8];
  static double interpol[24];
  static vector3 xcoor1, posit;
  static double Dx, Dy, Dz, Fx, Fy, Fz;
  static int ni[8]; 
  static double TriFactor[8];
  static double aNorm, c1, c2, c3, d1, d2, d3; 
  ElementsClass *elements=partrace->elements;
  static vector3 e;
  //static double vx, vy, vz, vNorm;

  //try
  //static vector3 B_from_Nc_v;
  //static vector3 B_from_Nr_v;
  //static vector3 B_to_Nr_v;
  //static vector3 B_from_Nc_1, B_from_Nc_2, B_from_Nc_3;
  //static vector3 B_from_Nr_1, B_from_Nr_2, B_from_Nr_3;
  //static vector3 B_to_Nr_1, B_to_Nr_2, B_to_Nr_3;
  //static vector3 C1_1, C1_2, C1_3;
  //static vector3 C2_1, C2_2, C2_3;
  //static vector3 C3_1, C3_2, C3_3;
  static double a1_to, a2_to, a1_from, a2_from;
  static double Df11, Df12, Df13;
  static double Df21, Df22, Df23;
  static double Df31, Df32, Df33;
  static double Dt11, Dt12, Dt13;
  static double Dt21, Dt22, Dt23;
  static double Dt31, Dt32, Dt33;
  //static vector3 Ct_1, Ct_2, Ct_3;
  //static double Z1, Z2, Z3;
  //static int ZTransform;
  int side=1;
  static double x1_from, x2_from, x3_from;
  //int i=0;
  static double ConditionedVariance_to_1, ConditionedVariance_to_2, ConditionedVariance_to_3;
  static double ConditionedVariance_from_1, ConditionedVariance_from_2, ConditionedVariance_from_3;
  //static double varia=0;
  //static double variai=0;
  //static double varia_from=0;
  //static double variai_from=0;
  //static double uncvaria=0;
  //static double cvaria=0;
  //double fac, meanUNC;

  //////////////////////// e_to //////////////////////////////////
	
  // e_to // get element dependant information
  dispersion_parameter_to=elements->soilprop->get_line(e_to);
	
  elements->GetElementXYZ(e_to,ex,ey,ez);
  disp.wc_to=elements->watercontent->get(e_to);
	
  //////////////////////////////// Velocity Determination for Recalculate Dispersion method
  switch(DispersionMethod) { 
  case 1: //no correction, same as reflection method
    switch(elements->VelType) {
    case 1: // velocity defined on elements
      disp.v_to.set(elements->velocity->get_line(e_to));
      break;
    case 2: // interpolate between nodes
      elements->GetNodesOfElement(e_to, ex,ey,ez, nodes);
      elements->GetStartCoordinates(nodes, ex,ey,ez, xcoor1);
      elements->initInterpolationFunction(nodes);
      elements->getInterpolationFunction(e_to, ex,ey,ez, nodes, interpol);
      elements->interpolateVelocity(xcoor1, pos, interpol, disp.v_to);
      break;
    case 3: // finite volume  
      finitevolindex=6*(ez+ey*elements->NoElz+ex*elements->NoElyz);
      disp.v_to.set(elements->finitevolvelocity[finitevolindex  ]*pos[0]+elements->finitevolvelocity[finitevolindex+1],
		    elements->finitevolvelocity[finitevolindex+2]*pos[1]+elements->finitevolvelocity[finitevolindex+3],
		    elements->finitevolvelocity[finitevolindex+4]*pos[2]+elements->finitevolvelocity[finitevolindex+5]);
      disp.v_to/=disp.wc_to;
      break;
    }
    break;
  case 2:  /////////////// Reflection Method ////////////////
    switch(elements->VelType) {
    case 1: // velocity defined on elements
      disp.v_to.set(elements->velocity->get_line(e_to));
      break;
    case 2: // interpolate between nodes
      elements->GetNodesOfElement(e_to, ex,ey,ez, nodes);
      elements->GetStartCoordinates(nodes, ex,ey,ez, xcoor1);
      elements->initInterpolationFunction(nodes);
      elements->getInterpolationFunction(e_to, ex,ey,ez, nodes, interpol);
      elements->interpolateVelocity(xcoor1, pos, interpol, disp.v_to);
      break;
    case 3: // finite volume  
      finitevolindex=6*(ez+ey*elements->NoElz+ex*elements->NoElyz);
      disp.v_to.set(elements->finitevolvelocity[finitevolindex  ]*pos[0]+elements->finitevolvelocity[finitevolindex+1],
		    elements->finitevolvelocity[finitevolindex+2]*pos[1]+elements->finitevolvelocity[finitevolindex+3],
		    elements->finitevolvelocity[finitevolindex+4]*pos[2]+elements->finitevolvelocity[finitevolindex+5]);
      disp.v_to/=disp.wc_to;
      if (fabs(disp.v_to[0])<1e-25) disp.v_to[0]=1e-25; // evtl vMin, better value, needs to be checked
      if (fabs(disp.v_to[1])<1e-25) disp.v_to[1]=1e-25;
      if (fabs(disp.v_to[2])<1e-25) disp.v_to[2]=1e-25;
      break;
    }
    break;
  case 3: /////////////// Interpolation Method /////////////////
    elements->GetNodesOfElement(e_to, ex,ey,ez, nodes);
    elements->GetDxyzFxyz(ex, ey, ez, pos, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
    elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor); //
    elements->InterpolateVelocityTriLinear(ex, ey, ez, TriFactor, ni, nodes, disp.v_to);
    elements->InterpolateWCTriLinear(TriFactor, ni, nodes, disp.wc_to);							
    break; // end of Interpolation Method /////////////////
  } //end switch Dispersion Methods
	
  // calculate diffusion
  switch(elements->DiffusionType) {
  case 1: // simple
    disp.dmc_to=dispersion_parameter_to[2];
    break;
  case 2: // water content based
	  // D = Di * wc**(7/3) / porosity**2
    disp.dmc_to=dispersion_parameter_to[2] * pow(disp.wc_to, wc_exponent) /
      (dispersion_parameter_to[3] * dispersion_parameter_to[3]);
    break;
  }
	
  disp.vNorm_to=sqrt(disp.v_to[2]*disp.v_to[2]+disp.v_to[0]*disp.v_to[0]+disp.v_to[1]*disp.v_to[1]); 
	
  if(disp.vNorm_to<vMin  || (dispersion_parameter_to[0]<vMin && dispersion_parameter_to[1]<vMin) ) {	
    Dt11=disp.dmc_to;
    Dt22=Dt11;
    Dt33=Dt11;
    disp.d_to[0]=Retardation*sqrt(disp.dmc_to) * DispWeight[0];
    disp.d_to[1]=Retardation*sqrt(disp.dmc_to) * DispWeight[1];
    disp.d_to[2]=Retardation*sqrt(disp.dmc_to) * DispWeight[2];
	  
    // for Reflection coefficient
    switch(elements->ReflectionMethod) {
    case 1: //Hoteit et al. 2002
      disp.dref_to[0]=disp.d_to[0];
      disp.dref_to[1]=disp.d_to[1];
      disp.dref_to[2]=disp.d_to[2];
      break;
    case 2: //Lim 2006
      disp.dref_to[0]=disp.d_to[0]* disp.wc_to;
      disp.dref_to[1]=disp.d_to[1]* disp.wc_to;
      disp.dref_to[2]=disp.d_to[2]* disp.wc_to;
      break;
    case 3: //Bechtold et al. 2010
      disp.dref_to[0]=disp.d_to[0]* disp.wc_to;
      disp.dref_to[1]=disp.d_to[1]* disp.wc_to;
      disp.dref_to[2]=disp.d_to[2]* disp.wc_to;
      break;
    } //switch
	  
  }
  else {
    // Dispersion coefficient Lichtner et al. 2002, for scaling of second dispersive displacement only use of variance, i.e. mean displacement in x, y and z
    a1_to = dispersion_parameter_to[1];
    a2_to = fabs(dispersion_parameter_to[0]-dispersion_parameter_to[1]);
    // only diagonal as mean displacement
    Dt11=disp.dmc_to + a1_to*disp.vNorm_to + a2_to*disp.v_to[0]*disp.v_to[0]/disp.vNorm_to ;
    Dt22=disp.dmc_to + a1_to*disp.vNorm_to + a2_to*disp.v_to[1]*disp.v_to[1]/disp.vNorm_to ;
    Dt33=disp.dmc_to + a1_to*disp.vNorm_to + a2_to*disp.v_to[2]*disp.v_to[2]/disp.vNorm_to ;
    disp.d_to[0]=Retardation*sqrt(Dt11 ) * DispWeight[0];
    disp.d_to[1]=Retardation*sqrt(Dt22 ) * DispWeight[1];
    disp.d_to[2]=Retardation*sqrt(Dt33 ) * DispWeight[2];
	  
    // for Reflection coefficient
    switch(elements->ReflectionMethod) {
    case 1: //Hoteit et al. 2002
      disp.dref_to[0]=disp.d_to[0];
      disp.dref_to[1]=disp.d_to[1];
      disp.dref_to[2]=disp.d_to[2];
      break;
    case 2: //Lim 2006
      disp.dref_to[0]=disp.d_to[0]* disp.wc_to;
      disp.dref_to[1]=disp.d_to[1]* disp.wc_to;
      disp.dref_to[2]=disp.d_to[2]* disp.wc_to;
      break;
    case 3: //Bechtold et al. 2010
      disp.dref_to[0]=disp.d_to[0]* disp.wc_to;
      disp.dref_to[1]=disp.d_to[1]* disp.wc_to;
      disp.dref_to[2]=disp.d_to[2]* disp.wc_to;
      break;
    } //switch
  }
	
  if (CoefAtBoundary==1 && DispersionMethod==2) {
    //////////////////////// CoefAtBoundary, Recalculate dref_from //////////////////////////////////
    // index _to used here for from, to save variables
    // e_from // get element dependant information
    dispersion_parameter_from=elements->soilprop->get_line(e_from);
    elements->GetElementXYZ(e_from,ex_from,ey_from,ez_from);
    disp.wc_from=elements->watercontent->get(e_from);
    switch(elements->VelType) { 
    case 1: // velocity defined on elements
      disp.v_from.set(elements->velocity->get_line(e_from));
      break;
    case 2:
      elements->GetNodesOfElement(e_from, ex_from,ey_from,ez_from, nodes);
      elements->GetStartCoordinates(nodes, ex_from,ey_from,ez_from, xcoor1);
      elements->initInterpolationFunction(nodes);
      elements->getInterpolationFunction(e_from, ex_from,ey_from,ez_from, nodes, interpol);
      break;
    case 3:
      finitevolindex=6*(ez_from+ey_from*elements->NoElz+ex_from*elements->NoElyz);
      break;
    }
    switch(elements->VelType) {
    case 2: // interpolate between nodes
      elements->interpolateVelocity(xcoor1, pos, interpol, disp.v_from);
      break;
    case 3: // finite volume  
      disp.v_from.set(elements->finitevolvelocity[finitevolindex  ]*pos[0]+elements->finitevolvelocity[finitevolindex+1],
		      elements->finitevolvelocity[finitevolindex+2]*pos[1]+elements->finitevolvelocity[finitevolindex+3],
		      elements->finitevolvelocity[finitevolindex+4]*pos[2]+elements->finitevolvelocity[finitevolindex+5]);
      disp.v_from/=disp.wc_from;
      break;
    }
	  
    // calculate diffusion
    switch(elements->DiffusionType) {
    case 1: // simple
      disp.dmc_from=dispersion_parameter_from[2];
      break;
    case 2: // water content based
	    // D = Di * wc**(7/3) / porosity**2
      disp.dmc_from=dispersion_parameter_from[2] * pow(disp.wc_from, wc_exponent) /
	(dispersion_parameter_from[3] * dispersion_parameter_from[3]);
      break;
    }
	  
    disp.vNorm_from=sqrt(disp.v_from[2]*disp.v_from[2]+disp.v_from[0]*disp.v_from[0]+disp.v_from[1]*disp.v_from[1]); 
	  
    if(disp.vNorm_from<vMin  || (dispersion_parameter_from[0]<vMin && dispersion_parameter_from[1]<vMin) ) {
      Df11=disp.dmc_from;
      Df22=Df11;
      Df33=Df11;
      disp.d_from[0]=Retardation*sqrt(disp.dmc_from) * DispWeight[0];
      disp.d_from[1]=Retardation*sqrt(disp.dmc_from) * DispWeight[1];
      disp.d_from[2]=Retardation*sqrt(disp.dmc_from) * DispWeight[2];
	    
      // for Reflection coefficient
      switch(elements->ReflectionMethod) {
      case 1: //Hoteit et al. 2002
	disp.dref_from[0]=disp.d_from[0];
	disp.dref_from[1]=disp.d_from[1];
	disp.dref_from[2]=disp.d_from[2];
	break;
      case 2: //Lim 2006
	disp.dref_from[0]=disp.d_from[0]* disp.wc_from;
	disp.dref_from[1]=disp.d_from[1]* disp.wc_from;
	disp.dref_from[2]=disp.d_from[2]* disp.wc_from;
	break;
      case 3: //Bechtold et al. 2010
	disp.dref_from[0]=disp.d_from[0]* disp.wc_from;
	disp.dref_from[1]=disp.d_from[1]* disp.wc_from;
	disp.dref_from[2]=disp.d_from[2]* disp.wc_from;
	break;
      } //switch
    }
    else {
      // Dispersion coefficient Lichtner et al. 2002, for scaling of second dispersive displacement only use of variance, i.e. mean displacement in x, y and z
      a1_from = dispersion_parameter_from[1];
      a2_from = fabs(dispersion_parameter_from[0]-dispersion_parameter_from[1]);
      // only diagonal as mean displacement
      Df11=disp.dmc_from + a1_from*disp.vNorm_from + a2_from*disp.v_from[0]*disp.v_from[0]/disp.vNorm_from ;
      Df22=disp.dmc_from + a1_from*disp.vNorm_from + a2_from*disp.v_from[1]*disp.v_from[1]/disp.vNorm_from ;
      Df33=disp.dmc_from + a1_from*disp.vNorm_from + a2_from*disp.v_from[2]*disp.v_from[2]/disp.vNorm_from ;
      disp.d_from[0]=Retardation*sqrt( Df11 ) * DispWeight[0];
      disp.d_from[1]=Retardation*sqrt( Df22 ) * DispWeight[1];
      disp.d_from[2]=Retardation*sqrt( Df33 ) * DispWeight[2];
	    
      // for Reflection coefficient
      switch(elements->ReflectionMethod) {
      case 1: //Hoteit et al. 2002
	disp.dref_from[0]=disp.d_from[0];
	disp.dref_from[1]=disp.d_from[1];
	disp.dref_from[2]=disp.d_from[2];
	break;
      case 2: //Lim 2006
	disp.dref_from[0]=disp.d_from[0]* disp.wc_from;
	disp.dref_from[1]=disp.d_from[1]* disp.wc_from;
	disp.dref_from[2]=disp.d_from[2]* disp.wc_from;
	break;
      case 3: //Bechtold et al. 2010
	disp.dref_from[0]=disp.d_from[0]* disp.wc_from;
	disp.dref_from[1]=disp.d_from[1]* disp.wc_from;
	disp.dref_from[2]=disp.d_from[2]* disp.wc_from;
	break;
      } //switch
    }
	  
  } // end if CoefAtBoundary
	
  //////////////////////// SecDispersiveNew = second dispersive displacement with transformed Z //////////////////////////////////
  if (SecDispersiveNew==1 && DispersionMethod==2){
	  
	  
    //disp.rand[0]=-1.0+2.0*partrace->Random->fastdraw();
    //disp.rand[1]=-1.0+2.0*partrace->Random->fastdraw();
    //disp.rand[2]=-1.0+2.0*partrace->Random->fastdraw();
    // calculate diffusion disp.dmc_to
    if(disp.vNorm_to<vMin) { //  || (dispersion_parameter_to[0]<vMin && dispersion_parameter_to[1]<vMin) ) {   
      // velocity==0, only diffusion
      c1=sqrt( DispFak*Retardation*dt*disp.dmc_to );
      disp.ds_to[0]= disp.rand[0] * c1 * DispWeight[0] ;
      disp.ds_to[1]= disp.rand[1] * c1 * DispWeight[1] ;
      disp.ds_to[2]= disp.rand[2] * c1 * DispWeight[2] ;
      //disp.rand_to[0]=disp.rand[0];
      //disp.rand_to[1]=disp.rand[1];
      //disp.rand_to[2]=disp.rand[2];
    }
    else {
	    
	    
      //ZTransform=1;
      if (disp.vNorm_to>vMin) {   //&& (ZTransform==1)) {
	////////////////////////////// B_from ///////////////////
	//vx=disp.v[0];
	//vy=disp.v[1];
	//vz=disp.v[2];
	/*		Df12=a2_from*disp.v_from[0]*disp.v_from[1]/disp.vNorm_from;
			Df13=a2_from*disp.v_from[0]*disp.v_from[2]/disp.vNorm_from;
			Df21=a2_from*disp.v_from[1]*disp.v_from[0]/disp.vNorm_from;
			Df23=a2_from*disp.v_from[1]*disp.v_from[2]/disp.vNorm_from;
			Df31=a2_from*disp.v_from[2]*disp.v_from[0]/disp.vNorm_from;
			Df32=a2_from*disp.v_from[2]*disp.v_from[1]/disp.vNorm_from;*/
	      
	// override because use of starting position, dmc_from a1_from and a2_from stay same
	Df11=disp.dmc_from + a1_from*disp.vNorm + a2_from*disp.v[0]*disp.v[0]/disp.vNorm ;
	Df22=disp.dmc_from + a1_from*disp.vNorm + a2_from*disp.v[1]*disp.v[1]/disp.vNorm ;
	Df33=disp.dmc_from + a1_from*disp.vNorm + a2_from*disp.v[2]*disp.v[2]/disp.vNorm ;
	Df12=a2_from*disp.v[0]*disp.v[1]/disp.vNorm;
	Df13=a2_from*disp.v[0]*disp.v[2]/disp.vNorm;
	Df21=a2_from*disp.v[1]*disp.v[0]/disp.vNorm;
	Df23=a2_from*disp.v[1]*disp.v[2]/disp.vNorm;
	Df31=a2_from*disp.v[2]*disp.v[0]/disp.vNorm;
	Df32=a2_from*disp.v[2]*disp.v[1]/disp.vNorm;
	Dt12=a2_to*disp.v_to[0]*disp.v_to[1]/disp.vNorm_to;
	Dt13=a2_to*disp.v_to[0]*disp.v_to[2]/disp.vNorm_to;
	Dt21=a2_to*disp.v_to[1]*disp.v_to[0]/disp.vNorm_to;
	Dt23=a2_to*disp.v_to[1]*disp.v_to[2]/disp.vNorm_to;
	Dt31=a2_to*disp.v_to[2]*disp.v_to[0]/disp.vNorm_to;
	Dt32=a2_to*disp.v_to[2]*disp.v_to[1]/disp.vNorm_to;
	      
	x1_from=disp.ds[0];///sqrt(dt); remaining displacement
	x2_from=disp.ds[1];///sqrt(dt);
	x3_from=disp.ds[2];///sqrt(dt);
	      
	elements->GetElementXYZ(e_from,ex_from,ey_from,ez_from);
	elements->GetElementXYZ(e_to,ex_to,ey_to,ez_to);
	if (ex_to!=ex_from) side=1;
	else if (ey_to!=ey_from) side=2;
	else if (ez_to!=ez_from) side=3;
	      
	if (side==1) {
	  disp.ds_to_correlated[0]=sqrt(Dt11)/sqrt(Df11)*x1_from;
	  disp.ds_to_correlated[1]=Dt21/sqrt(Dt11)/sqrt(Df11)*x1_from;
	  disp.ds_to_correlated[2]=Dt31/sqrt(Dt11)/sqrt(Df11)*x1_from;
	  //disp.ds_from_correlated[0]=Df11/Df11*x1_from;
	  disp.ds_from_correlated[1]=Df21/Df11*x1_from;
	  disp.ds_from_correlated[2]=Df31/Df11*x1_from;
	  ConditionedVariance_to_2= sqrt( Dt22 * (1- pow(((Dt21)/sqrt(Dt11)/sqrt(Dt22)),2)));
	  ConditionedVariance_to_3= sqrt( Dt33 * (1- pow(((Dt31)/sqrt(Dt11)/sqrt(Dt33)),2)));
	  ConditionedVariance_from_2= sqrt( Df22 * (1- pow(((Df21)/sqrt(Df11)/sqrt(Df22)),2)));
	  ConditionedVariance_from_3= sqrt( Df33 * (1- pow(((Df31)/sqrt(Df11)/sqrt(Df33)),2)));
	  disp.ds_to_uncorrelated[0]=0.0;
	  disp.ds_to_uncorrelated[1]=(x2_from-disp.ds_from_correlated[1])*ConditionedVariance_to_2/ConditionedVariance_from_2;
	  disp.ds_to_uncorrelated[2]=(x3_from-disp.ds_from_correlated[2])*ConditionedVariance_to_3/ConditionedVariance_from_3;
	}
	else if (side==2) {
	  disp.ds_to_correlated[0]=Dt12/sqrt(Dt22)/sqrt(Df22)*x2_from;
	  disp.ds_to_correlated[1]=sqrt(Dt22)/sqrt(Df22)*x2_from;
	  disp.ds_to_correlated[2]=Dt32/sqrt(Dt22)/sqrt(Df22)*x2_from;
	  disp.ds_from_correlated[0]=Df12/Df22*x2_from;
	  //disp.ds_from_correlated[1]=Df22/Df22*x2_from;
	  disp.ds_from_correlated[2]=Df32/Df22*x2_from;
	  ConditionedVariance_to_1= sqrt( Dt11 * (1- pow(((Dt12)/sqrt(Dt22)/sqrt(Dt11)),2)));
	  ConditionedVariance_to_3= sqrt( Dt33 * (1- pow(((Dt32)/sqrt(Dt22)/sqrt(Dt33)),2)));
	  ConditionedVariance_from_1= sqrt( Df11 * (1- pow(((Df12)/sqrt(Df22)/sqrt(Df11)),2)));
	  ConditionedVariance_from_3= sqrt( Df33 * (1- pow(((Df32)/sqrt(Df22)/sqrt(Df33)),2)));
	  disp.ds_to_uncorrelated[0]=(x1_from-disp.ds_from_correlated[0])*ConditionedVariance_to_1/ConditionedVariance_from_1;
	  disp.ds_to_uncorrelated[1]=0.0;
	  disp.ds_to_uncorrelated[2]=(x3_from-disp.ds_from_correlated[2])*ConditionedVariance_to_3/ConditionedVariance_from_3;
	}
	else if (side==3) {
	  disp.ds_to_correlated[0]=Dt13/sqrt(Dt33)/sqrt(Df33)*x3_from;
	  disp.ds_to_correlated[1]=Dt23/sqrt(Dt33)/sqrt(Df33)*x3_from;
	  disp.ds_to_correlated[2]=sqrt(Dt33)/sqrt(Df33)*x3_from;
	  disp.ds_from_correlated[0]=Df13/Df33*x3_from;
	  disp.ds_from_correlated[1]=Df23/Df33*x3_from;
	  //disp.ds_from_correlated[2]=Df33/Df33*x3_from;
	  ConditionedVariance_to_1= sqrt( Dt11 * (1- pow(((Dt13)/sqrt(Dt33)/sqrt(Dt11)),2)));
	  ConditionedVariance_to_2= sqrt( Dt22 * (1- pow(((Dt23)/sqrt(Dt33)/sqrt(Dt22)),2)));
	  ConditionedVariance_from_1= sqrt( Df11 * (1- pow(((Df13)/sqrt(Df33)/sqrt(Df11)),2)));
	  ConditionedVariance_from_2= sqrt( Df22 * (1- pow(((Df23)/sqrt(Df33)/sqrt(Df22)),2)));
	  disp.ds_to_uncorrelated[0]=(x1_from-disp.ds_from_correlated[0])*ConditionedVariance_to_1/ConditionedVariance_from_1;
	  disp.ds_to_uncorrelated[1]=(x2_from-disp.ds_from_correlated[1])*ConditionedVariance_to_2/ConditionedVariance_from_2;
	  disp.ds_to_uncorrelated[2]=0.0;
	}
	      
	disp.ds_to[0]= disp.ds_to_correlated[0] + disp.ds_to_uncorrelated[0] ;
	disp.ds_to[1]= disp.ds_to_correlated[1] + disp.ds_to_uncorrelated[1] ;
	disp.ds_to[2]= disp.ds_to_correlated[2] + disp.ds_to_uncorrelated[2] ;
	/*		cout<<disp.ds_to_correlated[0]<< " " << disp.ds_to_correlated[1]<< " " << disp.ds_to_correlated[2]<< " " << endl;
			cout<<disp.ds_from_correlated[0]<< " " << disp.ds_from_correlated[1]<< " " << disp.ds_from_correlated[2]<< " " << endl;		
			cout<<ConditionedVariance_to_1 << " " <<ConditionedVariance_to_2 << " " <<ConditionedVariance_from_1 << " " <<ConditionedVariance_from_2 << endl;
			cout<<disp.ds_to_uncorrelated[0]<< " " << disp.ds_to_uncorrelated[1]<< " " << disp.ds_to_uncorrelated[2]<< " " << endl;
			cout<<Df11<<" "<<Df12<<" "<<Df13<<" "<<Df21<<" "<<Df22<<" "<<Df23<<" "<<Df31<<" "<<Df32<<" "<<Df33<<endl;
			cout<<Dt11<<" "<<Dt12<<" "<<Dt13<<" "<<Dt21<<" "<<Dt22<<" "<<Dt23<<" "<<Dt31<<" "<<Dt32<<" "<<Dt33<<endl;
			cout<<"ds_to "<< disp.ds_to[0] <<" " << disp.ds_to[1] << " "<< disp.ds_to[2] <<endl;
	*/
	/*if (ez_to==0) {
	//cout << "Conditioned variance " << ConditionedVariance_from_1 << endl;
	fstream filestr;
	filestr.open ("uncor_from.txt", fstream::in | fstream::out | fstream::app);
	filestr<<(x1_from-disp.ds_from_correlated[0])<<" "<<(x3_from-disp.ds_from_correlated[2])<<endl;
	filestr.close();
	fstream filestr2;
	filestr2.open ("cor_from.txt", fstream::in | fstream::out | fstream::app);
	filestr2<<disp.ds_from_correlated[0]<<" "<<disp.ds_from_correlated[2]<<endl;
	filestr2.close();
	fstream filestr3;
	filestr3.open ("total_from.txt", fstream::in | fstream::out | fstream::app);
	filestr3<<x1_from<<" "<<x3_from<<endl;
	filestr3.close();
	}*/
	/*		if (ez_to==0) {
			varia_from+=x1_from*x1_from;
			variai_from+=1;
			cout<<"x1_from Var: "<<B_from_norm_1<<"Var(performance): "<<sqrt(varia_from/variai_from)<<endl;
			} 
	*/
      }
      else { // No Z Transform
	// generate 2 vectors a,b orthographical to v
	if(fabs(disp.v_to[0])<vMin && fabs(disp.v_to[1])<vMin) {
	  // v parallel to z-axis
	  disp.v1[0]=1.0; disp.v1[1]=0.0; disp.v1[2]=0.0;
	  disp.v2[0]=0.0; disp.v2[1]=1.0; disp.v2[2]=0.0;
	}
	else if (fabs(disp.v_to[0])<vMin && fabs(disp.v_to[2])<vMin) {
	  // v parallel to y-axis
	  disp.v1[0]=0.0; disp.v1[1]=0.0; disp.v1[2]=1.0;
	  disp.v2[0]=1.0; disp.v2[1]=0.0; disp.v2[2]=0.0;
	}
	else if (fabs(disp.v_to[1])<vMin && fabs(disp.v_to[2])<vMin) {
	  // v parallel to x-axis
	  disp.v1[0]=0.0; disp.v1[1]=1.0; disp.v1[2]=0.0;
	  disp.v2[0]=0.0; disp.v2[1]=0.0; disp.v2[2]=1.0;
	}
	else {
	  aNorm=disp.v_to[0]*disp.v_to[0]+disp.v_to[1]*disp.v_to[1];
	  aNorm=sqrt(aNorm); //MB
	  disp.v1[0]= - (disp.v_to[0]*disp.v_to[2]) / (disp.vNorm_to * aNorm);
	  disp.v1[1]= - (disp.v_to[1]*disp.v_to[2]) / (disp.vNorm_to * aNorm);
	  disp.v1[2]=	aNorm/disp.vNorm_to;												
	  // 
	  disp.v2[0]= - disp.v_to[1]/aNorm;
	  disp.v2[1]= disp.v_to[0]/aNorm;
	  disp.v2[2]= 0.0;
	}
	d1=dispersion_parameter_to[0]*disp.vNorm_to+disp.dmc_to;
	// dispersion coefficient D orthographical to v
	d2=dispersion_parameter_to[1]*disp.vNorm_to+disp.dmc_to;
	d3=d2;
	// dispersion in v-direction
	c1=disp.rand[0] * sqrt( DispFak*Retardation*dt*(d1) );
	// dispersion orthographical to v
	c2=sqrt( DispFak*Retardation*dt*(d2) );
	c3=disp.rand[2]*c2;
	c2*=disp.rand[1];
	disp.ds_to[0]=( c1*disp.v_to[0]/disp.vNorm_to + c2*disp.v1[0] + c3*disp.v2[0] ) * DispWeight[0];
	disp.ds_to[1]=( c1*disp.v_to[1]/disp.vNorm_to + c2*disp.v1[1] + c3*disp.v2[1] ) * DispWeight[1];
	disp.ds_to[2]=( c1*disp.v_to[2]/disp.vNorm_to + c2*disp.v1[2] + c3*disp.v2[2] ) * DispWeight[2];
	//disp.rand_to[0]=disp.rand[0];
	//disp.rand_to[1]=disp.rand[1];
	//disp.rand_to[2]=disp.rand[2];
      } // end no Z transform
    } // end of else (not only diffusion with vx vy vz =0 for which no transformation is needed)
  } // end of SecDispersiveNew=1

  return true;
}


void ElementsClass::Dispersion(double dt, long e, vector3& pos, DispersionClass& disp, double Retardation)
{ 
  static double aNorm;
  static double c1, c2, c3, d1, d2, d3;
  const double wc_exponent=7.0/3.0;
  static long ex,ey,ez;
  int finitevolindex=0;
  static long nodes[8];
  //static double interpol[24]; 
  //static double interpolWC[8];
  //static vector3 xcoor1;
  static double Dx, Dy, Dz, Fx, Fy, Fz;
  static int ni[8]; 
  static double TriFactor[8];
  ElementsClass *elements=partrace->elements;
	
  //static double varia1=0;
  //static double varia1i=0;
  //double meanDisp;
  //double a1, a2;
  //////////////////////////////// Velocity Determination for Dispersion coefficient ////////////////
  elements->soilprop->get_line(e, disp.parameter);  // for interpolation method even the disp.parameter must be interpolated, to be implemented if needed
  elements->GetElementXYZ(e,ex,ey,ez);
  disp.wc=elements->watercontent->get(e);
  switch(elements->DispersionMethod) {
  case 1: // DisperionMethod off
    switch(elements->VelType) {
    case 1: // velocity defined on elements
      disp.v.set(elements->velocity->get_line(e));
      break;
    case 2: // interpolate between nodes
      /*elements->GetNodesOfElement(e, ex,ey,ez, nodes);
	elements->GetStartCoordinates(nodes, ex,ey,ez, xcoor1);
	elements->initInterpolationFunction(nodes);
	elements->getInterpolationFunction(e, ex,ey,ez, nodes, interpol);
	elements->interpolateVelocity(xcoor1, pos, interpol, disp.v);*/
      elements->GetNodesOfElement(e, ex,ey,ez, nodes);
      elements->GetDxyzFxyz(ex, ey, ez, pos, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
      elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor); //
      elements->InterpolateVelocityTriLinear(ex, ey, ez, TriFactor, ni, nodes, disp.v);
      elements->InterpolateWCTriLinear(TriFactor, ni, nodes, disp.wc);					
      break;
    case 3: // finite volume  
      finitevolindex=6*(ez+ey*elements->NoElz+ex*elements->NoElyz);
      disp.v.set(elements->finitevolvelocity[finitevolindex  ]*pos[0]+elements->finitevolvelocity[finitevolindex+1],
		 elements->finitevolvelocity[finitevolindex+2]*pos[1]+elements->finitevolvelocity[finitevolindex+3],
		 elements->finitevolvelocity[finitevolindex+4]*pos[2]+elements->finitevolvelocity[finitevolindex+5]);
      disp.v/=disp.wc;
      break;
    }
    break; // DisperionMethod off
  case 2: //// Reflection Method /////////////////
    switch(elements->VelType) {
    case 1: // velocity defined on elements
      disp.v.set(elements->velocity->get_line(e));
      break;
    case 2: // interpolate between nodes
      /*elements->GetNodesOfElement(e, ex,ey,ez, nodes);
	elements->GetStartCoordinates(nodes, ex,ey,ez, xcoor1);
	elements->initInterpolationFunction(nodes);
	elements->getInterpolationFunction(e, ex,ey,ez, nodes, interpol);
	elements->interpolateVelocity(xcoor1, pos, interpol, disp.v);*/
      elements->GetNodesOfElement(e, ex,ey,ez, nodes);
      elements->GetDxyzFxyz(ex, ey, ez, pos, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
      elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor); //
      elements->InterpolateVelocityTriLinear(ex, ey, ez, TriFactor, ni, nodes, disp.v);
      elements->InterpolateWCTriLinear(TriFactor, ni, nodes, disp.wc);					
      break;
    case 3: // finite volume  
      finitevolindex=6*(ez+ey*elements->NoElz+ex*elements->NoElyz);
      disp.v.set(elements->finitevolvelocity[finitevolindex  ]*pos[0]+elements->finitevolvelocity[finitevolindex+1],
		 elements->finitevolvelocity[finitevolindex+2]*pos[1]+elements->finitevolvelocity[finitevolindex+3],
		 elements->finitevolvelocity[finitevolindex+4]*pos[2]+elements->finitevolvelocity[finitevolindex+5]);
      disp.v/=disp.wc;
      if (fabs(disp.v[0])<1e-25) disp.v[0]=1e-25;
      if (fabs(disp.v[1])<1e-25) disp.v[1]=1e-25;
      if (fabs(disp.v[2])<1e-25) disp.v[2]=1e-25;
      break;
    }
    break; // end of Reflection Method //
  case 3: //// Interpolation Method /////////////////
    switch(elements->VelType) {
    case 1: //  velocity defined on elements
      //if (numdjy==3545) {cout<<"velocity defined on elements"<<endl; cout<<"elements->NoNoxy "<<elements->NoNoxy<< " "<<elements->NoNox<<endl;}
      elements->GetNodesOfElement(e, ex,ey,ez, nodes);
      elements->GetDxyzFxyz(ex, ey, ez, pos, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
      elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor); //
      elements->InterpolateVelocityTriLinear(ex, ey, ez, TriFactor, ni, nodes, disp.v);
      elements->InterpolateWCTriLinear(TriFactor, ni, nodes, disp.wc);													
      break;
    case 2: // interpolate between nodes
      /* Horst elements->GetNodesOfElement(e, ex,ey,ez, nodes); // node number of element
	 elements->GetStartCoordinates(nodes, ex,ey,ez, xcoor1); // coordinates of first node
	 elements->initInterpolationFunction(nodes); // skipped for rectilinear grids
	 elements->getInterpolationFunction(e, ex,ey,ez, nodes, interpol); // nodes
	 elements->getInterpolationFunctionWC(e, ex,ey,ez, nodes, interpolWC);
	 elements->interpolateVelocity(xcoor1, pos, interpol, disp.v);
	 elements->interpolateWaterContent(xcoor1, pos, interpolWC, disp.wc); */
      elements->GetNodesOfElement(e, ex,ey,ez, nodes);
      elements->GetDxyzFxyz(ex, ey, ez, pos, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
      elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor); //
      elements->InterpolateVelocityTriLinear(ex, ey, ez, TriFactor, ni, nodes, disp.v);
      elements->InterpolateWCTriLinear(TriFactor, ni, nodes, disp.wc);													
      break;
    case 3: // finite volume
      elements->GetNodesOfElement(e, ex,ey,ez, nodes);					
      elements->GetDxyzFxyz(ex, ey, ez, pos, Dx, Dy, Dz, Fx, Fy, Fz); // Get values for factor calculation for trilinear interpolation
      elements->GetTrilinearFactors(Fx, Fy, Fz, TriFactor); //
      elements->InterpolateVelocityTriLinear(ex, ey, ez, TriFactor, ni, nodes, disp.v);			
      elements->InterpolateWCTriLinear(TriFactor, ni, nodes, disp.wc);													
      break;
    }
    break; // end of Interpolation Method /////////////////
  } // end of switch Dispersion Method
  //
  //if (numdjy==3545) {
  //cout<<"disp.v "<<disp.v[0]<<" "<<disp.v[1]<<" "<<disp.v[2]<<" "<<endl;//}
  aNorm=disp.v[0]*disp.v[0]+disp.v[1]*disp.v[1];
  disp.vNorm=sqrt(aNorm+disp.v[2]*disp.v[2]);
  aNorm=sqrt(aNorm); //MB
  // calculate diffusion
  switch(elements->DiffusionType) {
  case 1: // simple
    disp.dmc=disp.parameter[2];
    break;
  case 2: // water content based
    // D = Di * wc**(7/3) / porosity**2
    disp.dmc=disp.parameter[2] * pow(disp.wc, wc_exponent) /
      (disp.parameter[3]*disp.parameter[3]);
    break;
  }
  //if(disp.vNorm<diffusion*1e-07) {
  if(disp.vNorm<vMin  || (disp.parameter[0]<vMin && disp.parameter[1]<vMin) ) {   
    // velocity==0, only diffusion
    c1=sqrt( DispFak*Retardation*dt*disp.dmc );
    disp.ds[0]= disp.rand[0] * c1 * DispWeight[0] ;
    disp.ds[1]= disp.rand[1] * c1 * DispWeight[1] ;
    disp.ds[2]= disp.rand[2] * c1 * DispWeight[2] ;
  }
  else {
    // generate 2 vectors a,b orthographical to v
    if(fabs(disp.v[0])<vMin && fabs(disp.v[1])<vMin) {
      // v parallel to z-axis
      disp.v1[0]=1.0; disp.v1[1]=0.0; disp.v1[2]=0.0;
      disp.v2[0]=0.0; disp.v2[1]=1.0; disp.v2[2]=0.0;
    }
    else if (fabs(disp.v[0])<vMin && fabs(disp.v[2])<vMin) {
      // v parallel to y-axis
      disp.v1[0]=0.0; disp.v1[1]=0.0; disp.v1[2]=1.0;
      disp.v2[0]=1.0; disp.v2[1]=0.0; disp.v2[2]=0.0;
    }
    else if (fabs(disp.v[1])<vMin && fabs(disp.v[2])<vMin) {
      // v parallel to x-axis
      disp.v1[0]=0.0; disp.v1[1]=1.0; disp.v1[2]=0.0;
      disp.v2[0]=0.0; disp.v2[1]=0.0; disp.v2[2]=1.0;
    }
    else {
      disp.v1[0]= - (disp.v[0]*disp.v[2]) / (disp.vNorm * aNorm);
      disp.v1[1]= - (disp.v[1]*disp.v[2]) / (disp.vNorm * aNorm);
      disp.v1[2]=	aNorm/disp.vNorm;												
      // 
      disp.v2[0]= - disp.v[1]/aNorm;
      disp.v2[1]= disp.v[0]/aNorm;
      disp.v2[2]= 0.0;
    }

    d1=disp.parameter[0]*disp.vNorm+disp.dmc;
    // dispersion coefficient D orthographical to v
    d2=disp.parameter[1]*disp.vNorm+disp.dmc;
    d3=d2;

    // dispersion in v-direction
    //cout<<"disprand "<<disp.rand[0]<<" "<<disp.rand[1]<<" "<<disp.rand[2]<<" "<<DispFak<< " " <<Retardation<< " "<<dt<<" "<< d1 << " "<< d2 << " " << disp.parameter[0] << " " << disp.parameter[1] << " " << disp.dmc << endl;
    c1=disp.rand[0] * sqrt( DispFak*Retardation*dt*(d1) );
    // dispersion orthographical to v
    c2=sqrt( DispFak*Retardation*dt*(d2) );
    c3=disp.rand[2]*c2;
    c2*=disp.rand[1];
    //////////////////////////////////// Dispersion step ///////////////////////////////////
    disp.ds[0]=( c1*disp.v[0]/disp.vNorm + c2*disp.v1[0] + c3*disp.v2[0] ) * DispWeight[0];
    disp.ds[1]=( c1*disp.v[1]/disp.vNorm + c2*disp.v1[1] + c3*disp.v2[1] ) * DispWeight[1];
    disp.ds[2]=( c1*disp.v[2]/disp.vNorm + c2*disp.v1[2] + c3*disp.v2[2] ) * DispWeight[2];
	
    //cout<<disp.ds[0]<<" "<<disp.ds[1]<<" "<<disp.ds[2]<<endl;
    /*		if (fabs(disp.v[0])<1e-12 && pos[2]<0.15) {
		varia1+=disp.ds[0]*disp.ds[0]/dt;
		varia1i+=1;		
		a1 = disp.parameter[1];
		a2 = fabs(disp.parameter[0]-disp.parameter[1]);
		meanDisp = sqrt( 2* (a1*disp.vNorm + a2*disp.v[0]*disp.v[0]/disp.vNorm) );
		cout<<"Lower Variance: "<< meanDisp <<" Variance Performance: "<<sqrt(varia1/varia1i) <<endl;
		} */


  }
}

void ElementsClass::CountParticlesOutside()
{
  partrace->parallel->sum(ParticlesOutside, 6);
  if(processor==0) { 
    cout<<'\n';
    if(ParticlesOutside[0]>0) cout<<ParticlesOutside[0]<<" particles outside because x<xmin\n";
    if(ParticlesOutside[1]>0) cout<<ParticlesOutside[1]<<" particles outside because y<ymin\n";
    if(ParticlesOutside[2]>0) cout<<ParticlesOutside[2]<<" particles outside because z<zmin\n";
    if(ParticlesOutside[3]>0) cout<<ParticlesOutside[3]<<" particles outside because x>xmax\n";
    if(ParticlesOutside[4]>0) cout<<ParticlesOutside[4]<<" particles outside because y>ymax\n";
    if(ParticlesOutside[5]>0) cout<<ParticlesOutside[5]<<" particles outside because z>zmax\n";
    cout<<endl;
  }
}


void ElementsClass::GetElementsPosition(long i, vector3& pos) 
{
  if(i<0 || i>NoElxyz) partrace->error("Invalid arguement for GetElementsPosition");
  long n;
  n=i/NoElxy;
  pos[2]=0.5*(coords_z[n]+coords_z[n+1]);
  i%=NoElxy;
  n=i/NoElx;
  pos[1]=0.5*(coords_y[n]+coords_y[n+1]);
  n=i%NoElx;
  pos[0]=0.5*(coords_x[n]+coords_x[n+1]);
}


void ElementsClass::GetNodeCoordinates(long n, vector3& pos) 
{
  if(n<0 || n>NoNoxyz) partrace->error("Invalid node number for GetNodeCoordinates");
  pos[2]=coords_z[n/NoNoxy];
  n%=NoNoxy;
  pos[1]=coords_y[n/NoNox];
  pos[0]=coords_x[n%NoNox];
}


void ElementsClass::GetNodesOfElement(long el, long ex, long ey, long ez, long nodes[8]) const
{
  nodes[0]=ex+ey*NoNox+ez*NoNoxy;
  nodes[1]=nodes[0]+1;
  nodes[2]=nodes[0]+NoNox;
  nodes[3]=nodes[2]+1;
  nodes[4]=nodes[0]+NoNoxy;
  nodes[5]=nodes[1]+NoNoxy;
  nodes[6]=nodes[2]+NoNoxy;
  nodes[7]=nodes[3]+NoNoxy;
}


void ElementsClass::GetNoNo(long x, long y, long z, long nodes[8]) 
{
  nodes[0]= x+y*NoNox+z*NoNoxy;
  nodes[1]= nodes[0]+1;
  nodes[2]= nodes[0]+NoNox;
  nodes[3]= nodes[2]+1;
  nodes[4]= nodes[0]+NoNoxy;
  nodes[5]= nodes[1]+NoNoxy;
  nodes[6]= nodes[2]+NoNoxy;
  nodes[7]= nodes[3]+NoNoxy;
}


bool ElementsClass::GetNoNo(const vector3& pos, long nodes[8]) 
{
  long x,y,z;
  if(pos[0]<Boundaries[0] || pos[1]<Boundaries[1] || pos[2]<Boundaries[2]) return false;
  x=(long) (pos[0]/length_x);
  if(x>=NoElx) return false;
  y=(long) (pos[1]/length_y);
  if(y>=NoEly) return false;
  z=(long) (pos[2]/length_z);
  if(z>=NoElz) return false;

  nodes[0]= x+y*NoNox+z*NoNoxy;
  nodes[1]= nodes[0]+1;
  nodes[2]= nodes[0]+NoNox;
  nodes[3]= nodes[2]+1;
  nodes[4]= nodes[0]+NoNoxy;
  nodes[5]= nodes[1]+NoNoxy;
  nodes[6]= nodes[2]+NoNoxy;
  nodes[7]= nodes[3]+NoNoxy;
  return true;
}


int ElementsClass::calculateBufferSize(int minsize)
{
  if(NoElx>minsize) return NoElx;
  if(NoEly>minsize) return NoEly;
  if(NoElz>minsize) return NoElz;
  if(NoElxy>minsize) return NoElxy;
  if(NoElxz>minsize) return NoElxz;
  if(NoElyz>minsize) return NoElyz;
  return NoElxyz;
}


void ElementsClass::calculate_watercontentOnElements(ConstData *wc)
{
  long e, ex,ey,ez;
  long n[8];
  long from=watercontent->begin();
  long to=watercontent->end();

  // calculate water content for the element
  for(e=from;e<to;e++) {
    GetElementXYZ(e, ex,ey,ez);
    GetNodesOfElement(e, ex,ey,ez, n);
    watercontent->set_data(e, 0, 0.125*
			   ( wc->get(n[0])+wc->get(n[1])+wc->get(n[2])+wc->get(n[3])+
			     wc->get(n[4])+wc->get(n[5])+wc->get(n[6])+wc->get(n[7]) ) );
  }
}


// not used at the moment
void ElementsClass::GetFieldNo(int pro, int procs, int fields, int &von, int &bis)
{
  int i;
  double d, dpro;

  von=1; bis=0;
  if(procs==1) bis=fields;
  else if(fields==1) {
    if(pro==1) bis=1;
  }
  else if(procs==fields) {
    von=bis=pro;
  }
  else if(procs>fields) {
    d=(double) (procs-1) / (double) (fields-1);
    for(i=0;i<fields;i++) {
      von=(int) (d*(double) i) +1;
      if(von==pro) {
	bis=von=i+1;
	return;
      }
    }
  }
  else { // procs<fields
    d=(double) fields / (double) procs;
    dpro=(double) pro;
    bis=(int) (d*dpro);
    von=(int) (d*(dpro-1.0)) + 1;
  }
}


void ElementsClass::ReadField(const string& fn, const char *text, const char *head,
			      int ibuffer[4], ConstData* &iv)
{
  ifstream fin;
  string line;
  const char nl='\n';
  const int linesize=1024;
  long e, ex,ey,ez;
  int i, nx,ny;

  if(processor==0) {
    fin.open(fn.c_str());
    if(fin.fail()) partrace->error(string(text)+": "+fn+" could not be opened");
    getline(fin,line);
    if(line!=string(head))
      partrace->error(string("Invalid first line in ")+text+string(": ")+fn);

    while(fin.peek()=='#') getline(fin,line);
    fin>>ibuffer[0]; fin.ignore(linesize,nl);

    while(fin.peek()=='#') getline(fin,line);
    fin>>ex>>ey>>ez; fin.ignore(linesize,nl);
    if(ex!=NoElx || ey!=NoEly || ez!=NoElz)
      partrace->error(string("Invalid number of elements in ")+text+string(": ")+fn);

    while(fin.peek()=='#') getline(fin,line);
    fin>>ibuffer[1]; fin.ignore(linesize,nl);

    while(fin.peek()=='#') getline(fin,line);
    fin>>ibuffer[2]; fin.ignore(linesize,nl);
    if(ibuffer[2]<1)
      partrace->error(string("Invalid number of parameters in ")+text+string(": ")+fn);

    while(fin.peek()=='#') getline(fin,line);
    fin>>ibuffer[3]; fin.ignore(linesize,nl);
    if(ibuffer[3]<0)
      partrace->error(string("Invalid number of parameter sets in ")+text+string(": ")+fn);
    while(fin.peek()=='#') getline(fin,line);
  }
  partrace->parallel->broadcast(ibuffer, 4);
  nx=ibuffer[3]; ny=ibuffer[2];
  partrace->SetNewHandlerText(text);
  if(nx==0) {
    // no indexing
    if(partrace->distributedData)
      iv=new DistData(partrace->parallel, text, NoElxyz, ny);
    else
      iv=new LocalData(partrace->parallel, text, NoElxyz, ny);
  }
  else {
    if(partrace->distributedData)
      iv=new IndexedDistData(partrace->parallel, text, NoElxyz, ny, nx);
    else
      iv=new IndexedData(partrace->parallel, text, NoElxyz, ny, nx);
  }
  // read parameter sets
  double* param=new double[ny];
  for(i=0;i<ny;i++) param[i]=0.0;
  long pto=iv->global_datalines();
  for(e=0;e<pto;e++) {
    if(processor==0) {
      for(i=0;i<ny;i++) fin>>param[i];
      if(fin.eof()) partrace->error(string("Unexpected EOF in file: ")+fn);
    }
    iv->set_dataline_from_one_PE(e, param);
  }
  delete [] param;
  // read index
  if(nx==0) return; // no indexing
  int indexfile;
  string fnindex;
  ConstDataIndex *ivi;
  if(processor==0) {
    fin.ignore(linesize,nl);
    while(fin.peek()=='#') getline(fin,line);
    fin>>indexfile; fin.ignore(linesize,nl);
    if(fin.eof()) partrace->error(string("Unexpected EOF in file: ")+fn);
    if(indexfile<1 || indexfile>2 )
      partrace->error(string("Unexpected index mode in file: ")+fn);
    while(fin.peek()=='#') getline(fin,line);
    getline(fin,fnindex);
  }
  partrace->parallel->broadcast( &indexfile, 1);
  partrace->parallel->broadcast(fnindex);
  if(indexfile==1) {
    // read embedded index
    ivi=iv->alloc_index();
    if(processor==0) while(fin.peek()=='#') getline(fin,line);
    for(e=0;e<NoElxyz;e++) {
      if(processor==0) {
	fin>>ex;
	ex--;
	if(ex<0 || ex>=pto)
	  partrace->error(string("Invalid index in index file: ")+fnindex);
	if(fin.eof()) partrace->error(string("Unexpected EOF in file: ")+fnindex);
      }
      iv->set_index_data_from_one_PE(e, ex);
    }    
    if(processor==0) fin.close();
  }
  else if(indexfile==2) {
    // read index from separate file
    ivi=IndexListe[fnindex];
    if(ivi) {
      // index is already there
      iv->set_index(ivi);
    }
    else {
      // read index from file
      if(processor==0) {
        fin.clear(); // gcc bug
        fin.close();
        fin.open(fnindex.c_str());
        if(fin.fail()) partrace->error(string("Can not open index file: ")+fnindex);
	cout<<"Reading index file "<<fnindex<<endl;
        getline(fin,line);
        if(line!="#partrace index file") partrace->error(string("Wrong index file: ")+fnindex);
        while(fin.peek()=='#') getline(fin,line);
        fin>>e; fin.ignore(linesize,nl);
        if(e!=NoElxyz) partrace->error(string("Wrong number of elements in index file: ")+fnindex);
        while(fin.peek()=='#') getline(fin,line);
      }
      ivi=iv->alloc_index();
      for(e=0;e<NoElxyz;e++) {
	if(processor==0) {
	  fin>>ex;
	  ex--;
	  if(ex<0 || ex>=pto)
	    partrace->error(string("Invalid index in index file: ")+fnindex);
	  if(fin.eof()) partrace->error(string("Unexpected EOF in file: ")+fnindex);
	}
	iv->set_index_data_from_one_PE(e, ex);
      }
      if(processor==0) fin.close();
      IndexListe[fnindex]=ivi;
    }
  }
}


void ElementsClass::ReadFlowFields(bool restart)
{
  if(NextVTime->isActivated() && (restart || NextVTime->mustHandle(partrace->RunTime.Time))) {
    ReadVField(restart);
    calculateMaxTimeStepSize();
  }
  if(NextWCTime->isActivated() && (restart || NextWCTime->mustHandle(partrace->RunTime.Time)))
    ReadWCField(restart);
}


void ElementsClass::FiniteVolToNodes()
{	
  int index=0;
  int ex1, ex2, ey1, ey2, ez1, ez2, Num;
  long n, nx, ny, nz;
  vector3 vel1, vel2, vel3, vel4, vel5, vel6, vel7, vel8;
  double wc1, wc2, wc3, wc4, wc5, wc6, wc7, wc8;
  vector3 x; // position of node
  ElementsClass *elements=partrace->elements;	
  elements->nodevelocity.reserve(3*elements->NoNoxyz);
  elements->nodewc.reserve(elements->NoNoxyz);
	
  for(n=0;n<(elements->NoNoxyz);n++) {
    // XYZ  coordinate
    //elements->GetElementXYZ(n,nx,ny,nz);
    elements->GetNodeXYZ(n,nx,ny,nz);
    ex1=nx-1;
    ex2=nx;
    ey1=ny-1;
    ey2=ny;   
    ez1=nz-1;
    ez2=nz;
    GetNodeCoordinates(n, x);
    Num=8;
    if (n==0) cout<<"position: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<endl;
    // 1
    if (ex1>=0 && ey1>=0 && ez1>=0 && ex1<elements->NoElx && ey1<elements->NoEly && ez1<elements->NoElz) {
      index=6*(ez1+ey1*elements->NoElz+ex1*elements->NoElyz);
      vel1.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	       finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	       finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
      vel1/=watercontent->get(ex1+ey1*elements->NoElx+ez1*elements->NoElxy);
    }
    else {	
      vel1.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 2
    if (ex2>=0 && ey1>=0 && ez1>=0 && ex2<elements->NoElx && ey1<elements->NoEly && ez1<elements->NoElz) {
      index=6*(ez1+ey1*elements->NoElz+ex2*elements->NoElyz);
      vel2.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	       finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	       finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
      vel2/=watercontent->get(ex2+ey1*elements->NoElx+ez1*elements->NoElxy);
    }
    else {
      vel2.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 3
    if (ex1>=0 && ey2>=0 && ez1>=0 && ex1<elements->NoElx && ey2<elements->NoEly && ez1<elements->NoElz) {
      index=6*(ez1+ey2*elements->NoElz+ex1*elements->NoElyz);
      vel3.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	       finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	       finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
      vel3/=watercontent->get(ex1+ey2*elements->NoElx+ez1*elements->NoElxy);
    }
    else {
      vel3.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 4
    if (ex2>=0 && ey2>=0 && ez1>=0 && ex2<elements->NoElx && ey2<elements->NoEly && ez1<elements->NoElz) {
      index=6*(ez1+ey2*elements->NoElz+ex2*elements->NoElyz);
      vel4.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	       finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	       finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
      vel4/=watercontent->get(ex2+ey2*elements->NoElx+ez1*elements->NoElxy);
    }
    else {
      vel4.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 5
    if (ex1>=0 && ey1>=0 && ez2>=0 && ex1<elements->NoElx && ey1<elements->NoEly && ez2<elements->NoElz) {
      index=6*(ez2+ey1*elements->NoElz+ex1*elements->NoElyz);
      vel5.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	       finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	       finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
      vel5/=watercontent->get(ex1+ey1*elements->NoElx+ez2*elements->NoElxy);
    }
    else {
      vel5.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 6
    if (ex2>=0 && ey1>=0 && ez2>=0 && ex2<elements->NoElx && ey1<elements->NoEly && ez2<elements->NoElz) {
      index=6*(ez2+ey1*elements->NoElz+ex2*elements->NoElyz);
      vel6.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	       finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	       finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
      vel6/=watercontent->get(ex2+ey1*elements->NoElx+ez2*elements->NoElxy);
    }
    else {
      vel6.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 7
    if (ex1>=0 && ey2>=0 && ez2>=0 && ex1<elements->NoElx && ey2<elements->NoEly && ez2<elements->NoElz) {
      index=6*(ez2+ey2*elements->NoElz+ex1*elements->NoElyz);
      vel7.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	       finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	       finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
      vel7/=watercontent->get(ex1+ey2*elements->NoElx+ez2*elements->NoElxy);
    }
    else {
      vel7.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 8
    if (ex2>=0 && ey2>=0 && ez2>=0 && ex2<elements->NoElx && ey2<elements->NoEly && ez2<elements->NoElz) {
      index=6*(ez2+ey2*elements->NoElz+ex2*elements->NoElyz);
      vel8.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	       finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	       finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
      vel8/=watercontent->get(ex2+ey2*elements->NoElx+ez2*elements->NoElxy);
    }	
    else {
      vel8.set(0.0, 0.0, 0.0);
      Num-=1; }
    elements->nodevelocity[n*3+0] = (vel1[0]+vel2[0]+vel3[0]+vel4[0]+vel5[0]+vel6[0]+vel7[0]+vel8[0] ) /Num;
    elements->nodevelocity[n*3+1] = (vel1[1]+vel2[1]+vel3[1]+vel4[1]+vel5[1]+vel6[1]+vel7[1]+vel8[1] ) /Num;
    elements->nodevelocity[n*3+2] = (vel1[2]+vel2[2]+vel3[2]+vel4[2]+vel5[2]+vel6[2]+vel7[2]+vel8[2] ) /Num;
	
  } // end of for loop
  ex1=0; ey1=0; ez1=0; n=0;
  GetNodeCoordinates(n, x);
  index=6*(ez1+ey1*elements->NoElz+ex1*elements->NoElyz);
  vel1.set(finitevolvelocity[index  ]*x[0]+finitevolvelocity[index+1],
	   finitevolvelocity[index+2]*x[1]+finitevolvelocity[index+3],
	   finitevolvelocity[index+4]*x[2]+finitevolvelocity[index+5]);
  vel1/=watercontent->get(ex1+ey1*elements->NoElx+ez1*elements->NoElxy);
  cout<<"Successfully transferred finitevolvelocity to nodevelocity, check first node"<<endl;
	
  // Water content interpolation
  for(n=0;n<(elements->NoNoxyz);n++) {
    // XYZ  coordinate
    elements->GetNodeXYZ(n,nx,ny,nz);
    ex1=nx-1;
    ex2=nx;
    ey1=ny-1;
    ey2=ny;   
    ez1=nz-1;
    ez2=nz;
    GetNodeCoordinates(n, x);
    Num=8;
    //if (n==0) cout<<"position: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<endl;
    // 1
    if (ex1>=0 && ey1>=0 && ez1>=0 && ex1<elements->NoElx && ey1<elements->NoEly && ez1<elements->NoElz) {
      index=(ex1+ey1*elements->NoElx+ez1*elements->NoElxy);
      wc1=watercontent->get(index);
    }
    else {
      wc1=0.0;
      Num-=1; }
    // 2
    if (ex2>=0 && ey1>=0 && ez1>=0 && ex2<elements->NoElx && ey1<elements->NoEly && ez1<elements->NoElz) {
      index=(ex2+ey1*elements->NoElx+ez1*elements->NoElxy);
      wc2=watercontent->get(index);
    }
    else {
      wc2=0.0;
      Num-=1; }
    // 3
    if (ex1>=0 && ey2>=0 && ez1>=0 && ex1<elements->NoElx && ey2<elements->NoEly && ez1<elements->NoElz) {
      index=(ex1+ey2*elements->NoElx+ez1*elements->NoElxy);
      wc3=watercontent->get(index);
    }
    else {
      wc3=0.0;
      Num-=1; }
    // 4
    if (ex2>=0 && ey2>=0 && ez1>=0 && ex2<elements->NoElx && ey2<elements->NoEly && ez1<elements->NoElz) {
      index=(ex2+ey2*elements->NoElx+ez1*elements->NoElxy);
      wc4=watercontent->get(index);
    }
    else {
      wc4=0.0;
      Num-=1; }
    // 5
    if (ex1>=0 && ey1>=0 && ez2>=0 && ex1<elements->NoElx && ey1<elements->NoEly && ez2<elements->NoElz) {
      index=(ex1+ey1*elements->NoElx+ez2*elements->NoElxy);
      wc5=watercontent->get(index);
    }
    else {
      wc5=0.0;
      Num-=1; }
    // 6
    if (ex2>=0 && ey1>=0 && ez2>=0 && ex2<elements->NoElx && ey1<elements->NoEly && ez2<elements->NoElz) {
      index=(ex2+ey1*elements->NoElx+ez2*elements->NoElxy);
      wc6=watercontent->get(index);
    }
    else {
      wc6=0.0;
      Num-=1; }
    // 7
    if (ex1>=0 && ey2>=0 && ez2>=0 && ex1<elements->NoElx && ey2<elements->NoEly && ez2<elements->NoElz) {
      index=(ex1+ey2*elements->NoElx+ez2*elements->NoElxy);
      wc7=watercontent->get(index);
    }
    else {
      wc7=0.0;
      Num-=1; }
    // 8
    if (ex2>=0 && ey2>=0 && ez2>=0 && ex2<elements->NoElx && ey2<elements->NoEly && ez2<elements->NoElz) {
      index=(ex2+ey2*elements->NoElx+ez2*elements->NoElxy);
      wc8=watercontent->get(index);
    }
    else {
      wc8=0.0;
      Num-=1; }

    elements->nodewc[n] = ( wc1+wc2+wc3+wc4+wc5+wc6+wc7+wc8 ) /Num;
  }

}

void ElementsClass::ElementsToNodes()
{	
  int index=0;
  int ex1, ex2, ey1, ey2, ez1, ez2, Num;
  long n, nx, ny, nz;
  vector3 vel1, vel2, vel3, vel4, vel5, vel6, vel7, vel8;
  double wc1, wc2, wc3, wc4, wc5, wc6, wc7, wc8;
  vector3 x; // position of node
  ElementsClass *elements=partrace->elements;	
  elements->nodevelocity.reserve(3*elements->NoNoxyz);
  elements->nodewc.reserve(elements->NoNoxyz);

  //cout<<"elements to nodes"<<endl;
  for(n=0;n<(elements->NoNoxyz);n++) {
    // XYZ  coordinate
    //elements->GetElementXYZ(n,nx,ny,nz);
    elements->GetNodeXYZ(n,nx,ny,nz);
    ex1=nx-1;
    ex2=nx;
    ey1=ny-1;
    ey2=ny;   
    ez1=nz-1;
    ez2=nz;
    GetNodeCoordinates(n, x);
    Num=8;
    //if (n==0) cout<<"position: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<endl;
    // 1
    if (ex1>=0 && ey1>=0 && ez1>=0 && ex1<elements->NoElx && ey1<elements->NoEly && ez1<elements->NoElz) {
      index=ex1+ey1*elements->NoElx+ez1*elements->NoElxy;
      vel1.set(elements->velocity->get_line(index));
    }
    else {
      vel1.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 2
    if (ex2>=0 && ey1>=0 && ez1>=0 && ex2<elements->NoElx && ey1<elements->NoEly && ez1<elements->NoElz) {
      index=ex2+ey1*elements->NoElx+ez1*elements->NoElxy;
      vel2.set(elements->velocity->get_line(index));
    }
    else {
      vel2.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 3
    if (ex1>=0 && ey2>=0 && ez1>=0 && ex1<elements->NoElx && ey2<elements->NoEly && ez1<elements->NoElz) {
      index=ex1+ey2*elements->NoElx+ez1*elements->NoElxy;
      vel3.set(elements->velocity->get_line(index));
    }
    else {
      vel3.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 4
    if (ex2>=0 && ey2>=0 && ez1>=0 && ex2<elements->NoElx && ey2<elements->NoEly && ez1<elements->NoElz) {
      index=ex2+ey2*elements->NoElx+ez1*elements->NoElxy;
      vel4.set(elements->velocity->get_line(index));
    }
    else {
      vel4.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 5
    if (ex1>=0 && ey1>=0 && ez2>=0 && ex1<elements->NoElx && ey1<elements->NoEly && ez2<elements->NoElz) {
      index=ex1+ey1*elements->NoElx+ez2*elements->NoElxy;
      vel5.set(elements->velocity->get_line(index));
    }
    else {
      vel5.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 6
    if (ex2>=0 && ey1>=0 && ez2>=0 && ex2<elements->NoElx && ey1<elements->NoEly && ez2<elements->NoElz) {
      index=ex2+ey1*elements->NoElx+ez2*elements->NoElxy;
      vel6.set(elements->velocity->get_line(index));
    }
    else {
      vel6.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 7
    if (ex1>=0 && ey2>=0 && ez2>=0 && ex1<elements->NoElx && ey2<elements->NoEly && ez2<elements->NoElz) {
      index=ex1+ey2*elements->NoElx+ez2*elements->NoElxy;
      vel7.set(elements->velocity->get_line(index));
    }
    else {
      vel7.set(0.0, 0.0, 0.0);
      Num-=1; }
    // 8
    if (ex2>=0 && ey2>=0 && ez2>=0 && ex2<elements->NoElx && ey2<elements->NoEly && ez2<elements->NoElz) {
      index=ex2+ey2*elements->NoElx+ez2*elements->NoElxy;
      vel8.set(elements->velocity->get_line(index));
    }
    else {
      vel8.set(0.0, 0.0, 0.0);
      Num-=1; }

    elements->nodevelocity[n*3+0] = (vel1[0]+vel2[0]+vel3[0]+vel4[0]+vel5[0]+vel6[0]+vel7[0]+vel8[0] ) /Num;
    elements->nodevelocity[n*3+1] = (vel1[1]+vel2[1]+vel3[1]+vel4[1]+vel5[1]+vel6[1]+vel7[1]+vel8[1] ) /Num;
    elements->nodevelocity[n*3+2] = (vel1[2]+vel2[2]+vel3[2]+vel4[2]+vel5[2]+vel6[2]+vel7[2]+vel8[2] ) /Num;
  } // end of for loop
  ex1=0; ey1=0; ez1=0; n=0;
  GetNodeCoordinates(n, x);
  //cout<<"Successfully transferred elementvelocity to nodevelocity"<<endl;
  /*cout<<"x-velocity "<<nodevelocity[0]<<endl;
    cout<<"z-velocity "<<nodevelocity[2]<<endl;

    cout<<"x-velocity "<<nodevelocity[909]<<endl;
    cout<<"z-velocity "<<nodevelocity[911]<<endl;


    cout<<"x-velocity "<<nodevelocity[9333]<<endl;
    cout<<"z-velocity "<<nodevelocity[9335]<<endl;

    cout<<"x-velocity "<<nodevelocity[12120]<<endl;
    cout<<"z-velocity "<<nodevelocity[12122]<<endl;*/

  // Water content interpolation
  for(n=0;n<(elements->NoNoxyz);n++) {
    // XYZ  coordinate
    elements->GetNodeXYZ(n,nx,ny,nz);
    ex1=nx-1;
    ex2=nx;
    ey1=ny-1;
    ey2=ny;   
    ez1=nz-1;
    ez2=nz;
    GetNodeCoordinates(n, x);
    Num=8;
    //if (n==0) cout<<"position: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<endl;
    // 1
    if (ex1>=0 && ey1>=0 && ez1>=0 && ex1<elements->NoElx && ey1<elements->NoEly && ez1<elements->NoElz) {
      index=(ex1+ey1*elements->NoElx+ez1*elements->NoElxy);
      wc1=watercontent->get(index);
    }
    else {
      wc1=0.0;
      Num-=1; }
    // 2
    if (ex2>=0 && ey1>=0 && ez1>=0 && ex2<elements->NoElx && ey1<elements->NoEly && ez1<elements->NoElz) {
      index=(ex2+ey1*elements->NoElx+ez1*elements->NoElxy);
      wc2=watercontent->get(index);
    }
    else {
      wc2=0.0;
      Num-=1; }
    // 3
    if (ex1>=0 && ey2>=0 && ez1>=0 && ex1<elements->NoElx && ey2<elements->NoEly && ez1<elements->NoElz) {
      index=(ex1+ey2*elements->NoElx+ez1*elements->NoElxy);
      wc3=watercontent->get(index);
    }
    else {
      wc3=0.0;
      Num-=1; }
    // 4
    if (ex2>=0 && ey2>=0 && ez1>=0 && ex2<elements->NoElx && ey2<elements->NoEly && ez1<elements->NoElz) {
      index=(ex2+ey2*elements->NoElx+ez1*elements->NoElxy);
      wc4=watercontent->get(index);
    }
    else {
      wc4=0.0;
      Num-=1; }
    // 5
    if (ex1>=0 && ey1>=0 && ez2>=0 && ex1<elements->NoElx && ey1<elements->NoEly && ez2<elements->NoElz) {
      index=(ex1+ey1*elements->NoElx+ez2*elements->NoElxy);
      wc5=watercontent->get(index);
    }
    else {
      wc5=0.0;
      Num-=1; }
    // 6
    if (ex2>=0 && ey1>=0 && ez2>=0 && ex2<elements->NoElx && ey1<elements->NoEly && ez2<elements->NoElz) {
      index=(ex2+ey1*elements->NoElx+ez2*elements->NoElxy);
      wc6=watercontent->get(index);
    }
    else {
      wc6=0.0;
      Num-=1; }
    // 7
    if (ex1>=0 && ey2>=0 && ez2>=0 && ex1<elements->NoElx && ey2<elements->NoEly && ez2<elements->NoElz) {
      index=(ex1+ey2*elements->NoElx+ez2*elements->NoElxy);
      wc7=watercontent->get(index);
    }
    else {
      wc7=0.0;
      Num-=1; }
    // 8
    if (ex2>=0 && ey2>=0 && ez2>=0 && ex2<elements->NoElx && ey2<elements->NoEly && ez2<elements->NoElz) {
      index=(ex2+ey2*elements->NoElx+ez2*elements->NoElxy);
      wc8=watercontent->get(index);
    }
    else {
      wc8=0.0;
      Num-=1; }

    elements->nodewc[n] = ( wc1+wc2+wc3+wc4+wc5+wc6+wc7+wc8 ) /Num;
  }
  //cout<<"elements to nodes end"<<endl;	
}

void ElementsClass::findFieldPos(ifstream &inp, double simtime)
{
  // find the corresponding field in stream.
  string line;
  double newtime;
  double oldtime=-1.0;
  while(oldtime<simtime) {
    getline(inp,line);
    if(inp.eof()) {
      inp.clear();
      break;
    }
    if(line.find("#simulation time:")==0) {
      istringstream s(line.substr(line.find(':')+1));
      s>>newtime;
      if(oldtime==-1.0) oldtime=newtime;
      if(newtime<simtime) oldtime=newtime;
      else break;
    }
  }
  newtime=-1.0;
  inp.seekg(0, ios::beg);
  while(newtime<oldtime) {
    getline(inp,line);
    if(line.find("#simulation time:")==0) {
      istringstream s(line.substr(line.find(':')+1));
      s>>newtime;
    }
  }
}


void ElementsClass::ReadVField(bool restart)
{
  int i, n;
  static bool firstread=true;
  double *v;
  string line;
  const char nl='\n';
  string text("Velocity file");

  if(processor==0) {
    if(firstread) {
      firstread=false;
      velinp.open(fnVelo.c_str());
      if(velinp.fail()) partrace->error(text+" "+fnVelo+" could not be opened");
      getline(velinp,line);
      if(line!=string("#velocity"))
	partrace->error(string("Invalid first line in ")+text+string(": ")+fnVelo);
      findFieldPos(velinp, partrace->RunTime.Time);
    }
  } // processor==0

  int to=velocity->global_datalines();
  v=velocity->start_buffered_distribute();
  for(n=0;n<to;n++) {
    if(processor==0) {
      velinp>>i>>v[0]>>v[1]>>v[2]; velinp.ignore(1000,nl);
      if(velinp.eof()) partrace->error(string("Unexpected EOF in file: ")+fnVelo);
      if(i-1!=n) {
	ostringstream s;
	s<<"Invalid element/nodenumber: "<<i<<" while reading "<<fnVelo;
	partrace->error(s.str());
      }
    }
    v=velocity->buffered_distribute(n);
  }
  velocity->end_buffered_distribute();

  // read time for next v-field
  double nexttime;
  if(processor==0) {
    getline(velinp,line);
    if(velinp.eof()) nexttime=NextVTime->deactivate();
    else {
      if(line.find("#simulation time:")!=0)
	partrace->error(fnVelo+" has wrong file format: simulation time missing\n"+line);
      istringstream s(line.substr(line.find(':')+1));
      s>>nexttime;
    }
    cout<<"Finished reading "<<text<<' '<<fnVelo<<" on processor "<<processor
	<<",  next field at time: "<<nexttime<<endl;
  }
  // distribute time
  partrace->parallel->broadcast(&nexttime, 1);
  NextVTime->set_next_time(partrace->RunTime.Time, nexttime);
  //  cout<<velocity;
}


void ElementsClass::ReadWCField(bool restart)
{
  int i, n;
  static bool firstread=true;
  string line;
  const char nl='\n';
  string text("Water Content file");

  if(processor==0) {
    if(firstread) {
      firstread=false;
      wcinp.open(fnWC.c_str());
      if(wcinp.fail()) partrace->error(text+" "+fnWC+" could not be opened");
      getline(wcinp,line);
      if(line!=string("#water content"))
	partrace->error(string("Invalid first line in ")+text+string(": ")+fnWC);
      findFieldPos(wcinp, partrace->RunTime.Time);
    }
  } // processor==0

  ConstData *wcdata;
  if(watercontentOnNodes) wcdata=watercontentOnNodes;
  else wcdata=watercontent;
  int to=wcdata->global_datalines();
  double* wc=wcdata->start_buffered_distribute();
  for(n=0;n<to;n++) {
    if(processor==0) {
      wcinp>>i>>*wc;
      if(wcinp.eof()) partrace->error(string("Unexpected EOF in file: ")+fnWC);
      if(i-1!=n) {
	ostringstream s;
	s<<"Invalid element/nodenumber: "<<i<<" while reading "<<fnWC;
	partrace->error(s.str());
      }
    }
    wc=wcdata->buffered_distribute(n);
  }
  wcdata->end_buffered_distribute();

  // read time for next field
  double nexttime;
  if(processor==0) {
    WCFieldRestartPos=WCFieldPos;
    wcinp.ignore(1000,nl);
    getline(wcinp,line);
    if(wcinp.eof()) nexttime=NextWCTime->deactivate();
    else {
      if(line.find("#simulation time:")!=0)
	partrace->error(fnWC+" has wrong file format: simulation time missing");
      istringstream s(line.substr(line.find(':')+1));
      s>>nexttime;
    }
    cout<<"Finished reading "<<text<<' '<<fnWC
	<<",  next field at time: "<<nexttime<<endl;
  }

  // distribute time
  partrace->parallel->broadcast( &nexttime, 1);
  NextWCTime->set_next_time(partrace->RunTime.Time, nexttime);
  // calculate the water content for the elements.
  if(watercontentOnNodes) calculate_watercontentOnElements(watercontentOnNodes);
}

// void ElementsClass::ReadDecayField(bool restart)
// {
//   int i, n;
//   static bool firstread=true;
//   string line;
//   const char nl='\n';
//   string text("Decay Input file");

//   if(processor==0) {
//     if(firstread) {
//       firstread=false;
//       decayinp.open(DecayField.c_str());
//       if(decayinp.fail()) partrace->error(text+" "+DecayField+" could not be opened");
//       getline(decayinp,line);
//       if(line!=string("#water content"))
// 	partrace->error(string("Invalid first line in ")+text+string(": ")+DecayField);
//       findFieldPos(decayinp, partrace->RunTime.Time);
//     }
//   } // processor==0

//   ConstData *wcdata;
//   if(watercontentOnNodes) wcdata=watercontentOnNodes;
//   else wcdata=watercontent;
//   int to=wcdata->global_datalines();
//   double* wc=wcdata->start_buffered_distribute();
//   for(n=0;n<to;n++) {
//     if(processor==0) {
//       wcinp>>i>>*wc;
//       if(wcinp.eof()) partrace->error(string("Unexpected EOF in file: ")+fnWC);
//       if(i-1!=n) {
// 	ostringstream s;
// 	s<<"Invalid element/nodenumber: "<<i<<" while reading "<<fnWC;
// 	partrace->error(s.str());
//       }
//     }
//     wc=wcdata->buffered_distribute(n);
//   }
//   wcdata->end_buffered_distribute();

//   // read time for next field
//   double nexttime;
//   if(processor==0) {
//     WCFieldRestartPos=WCFieldPos;
//     wcinp.ignore(1000,nl);
//     getline(wcinp,line);
//     if(wcinp.eof()) nexttime=NextWCTime->deactivate();
//     else {
//       if(line.find("#simulation time:")!=0)
// 	partrace->error(fnWC+" has wrong file format: simulation time missing");
//       istringstream s(line.substr(line.find(':')+1));
//       s>>nexttime;
//     }
//     cout<<"Finished reading "<<text<<' '<<fnWC
// 	<<",  next field at time: "<<nexttime<<endl;
//   }

//   // distribute time
//   partrace->parallel->broadcast( &nexttime, 1);
//   NextWCTime->set_next_time(partrace->RunTime.Time, nexttime);
//   // calculate the water content for the elements.
//   if(watercontentOnNodes) calculate_watercontentOnElements(watercontentOnNodes);
// }


void ElementsClass::InitBoundaries(boundarydata* bdata, int no){

  partrace->SetNewHandlerText("BoundaryClass");
  BoundaryXmin=new BoundaryClass(NoEly,NoElz);
  BoundaryXmax=new BoundaryClass(NoEly,NoElz);
  BoundaryYmin=new BoundaryClass(NoElx,NoElz);
  BoundaryYmax=new BoundaryClass(NoElx,NoElz);
  BoundaryZmin=new BoundaryClass(NoElx,NoEly);
  BoundaryZmax=new BoundaryClass(NoElx,NoEly);
  BoundaryClass* boundary=0;
  long x,y, yfrom,yto, xfrom,xto;
  for(int i=0;i<no;i++) {
    switch(bdata[i].type) {
    case 1: // reflection boundary
      switch(bdata[i].side) {
      case 1: boundary=BoundaryZmin; break;     
      case 2: boundary=BoundaryZmax; break;     
      case 3: boundary=BoundaryYmin; break;     
      case 4: boundary=BoundaryYmax; break;     
      case 5: boundary=BoundaryXmin; break;     
      case 6: boundary=BoundaryXmax; break;     
      default: partrace->error("Invalid boundary side."); break;
      }
      yfrom=bdata[i].y1-1;
      yto=bdata[i].y2;
      xfrom=bdata[i].x1-1;
      xto=bdata[i].x2;
      for(y=yfrom;y<yto;y++)
	for(x=xfrom;x<xto;x++)
	  boundary->settype(x,y,bdata[i].type);
      break;
    case 2: // concentration based injection boundary
      boundaryConcentrationNo++;  // calculate the number of boundary injections.
      break;
    case 3: // injection boundary profiles
      boundaryInjectionNo++;
      break;
    case 4: // global concentration boundary
      globalConcentrationNo++;
      break;
    case 5: // reflection boundary
      switch(bdata[i].side) {
      case 1: boundary=BoundaryZmin; break;     
      case 2: boundary=BoundaryZmax; break;     
      case 3: boundary=BoundaryYmin; break;     
      case 4: boundary=BoundaryYmax; break;     
      case 5: boundary=BoundaryXmin; break;     
      case 6: boundary=BoundaryXmax; break;     
      default: partrace->error("Invalid boundary side."); break;
      }
      yfrom=bdata[i].y1-1;
      yto=bdata[i].y2;
      xfrom=bdata[i].x1-1;
      xto=bdata[i].x2;
      for(y=yfrom;y<yto;y++)
	for(x=xfrom;x<xto;x++)
	  boundary->settype(x,y,bdata[i].type);
      break;
    default: partrace->error("Invalid boundary type.");
      break;
    }
  }
  InitBoundaryConcentration(bdata, no);
  InitBoundaryInjections(bdata, no);
  InitGlobalConcentration(bdata, no);
}


void ElementsClass::InitBoundaryConcentration(boundarydata* bdata, int no)
{
  int i,num;
  string name;
  if(boundaryConcentrationNo==0) return;
  partrace->SetNewHandlerText("InitBoundaryConcentration");
  boundaryConcentration=new BoundaryConcentrationClass[boundaryConcentrationNo];
  BoundaryConcentrationClass *injection=boundaryConcentration;
  Event *event;
  for(i=0; i<no; i++) {
    if(bdata[i].type==2 ) {
      injection->type=bdata[i].type;
      injection->side=bdata[i].side;
      injection->x1=bdata[i].x1-1;
      injection->x2=bdata[i].x2-1;
      injection->y1=bdata[i].y1-1;
      injection->y2=bdata[i].y2-1;
      injection->value=bdata[i].value;
      injection->time1=bdata[i].InjectFrom;
      injection->time2=bdata[i].InjectTo;
      injection->dirac=(injection->time1==injection->time2);
      name=bdata[i].particlename;
      injection->particle=partrace->particles.FindParticleType(name);
      injection->nnpx=injection->x2-injection->x1+2;
      injection->nnpy=injection->y2-injection->y1+2;
      injection->nelx=injection->nnpx-1;
      injection->nely=injection->nnpy-1;
      injection->nnp=injection->nnpx*injection->nnpy;
      injection->nel=injection->nelx*injection->nely;
      if(VelType==2) num=injection->nnp;
      else num=1;
      injection->node=new int[num];
      injection->mass=new double[injection->nel];
      if(!injection->dirac)
	fill(injection->mass, injection->mass+injection->nel, 0.0);
      setBoundaryConcentration(injection);
      event=new Event("boundaryConcentrationStart", injection->time1, 1);
      partrace->RunTime.insert(event);
      if( !injection->dirac ) {
	event=new Event("boundaryConcentrationStop", injection->time2, 1);
	partrace->RunTime.insert(event);
      }
      injection++;
    }
  }
  if(processor==0)
    cout<<"Initializing "<<boundaryConcentrationNo<<" injection boundaries."<<endl;
}


void ElementsClass::InitBoundaryInjections(boundarydata* bdata, int no)
{
  if(boundaryInjectionNo==0) return;
  double simtime=partrace->RunTime.Time;
  int i,j;
  int e1,e2, e12;
  string line, name;
  double x;
  partrace->SetNewHandlerText("InitBoundaryInjection");
  boundaryInjection=new BoundaryInjectionClass[boundaryInjectionNo];
  BoundaryInjectionClass *injection=boundaryInjection;
  Event *event;
  bool found;
  int ibuffer[2];
  for(i=0; i<no; i++) {
    if(bdata[i].type!=3 ) continue;
    injection->fn=bdata[i].injectfile;
    if(processor==0) {
      injection->inp.open(bdata[i].injectfile);
      if( injection->inp.fail() ) 
	partrace->error(string("Could not open file: ")+injection->fn);
      getline(injection->inp, line);
      if(line!=string("boundary injection")) 
	partrace->error(string("Invalid first line in file: ")+injection->fn+"  line="+line);
      injection->inp>>ibuffer[0]>>ibuffer[1]; injection->inp.ignore(1024, '\n');
    }
    partrace->parallel->broadcast( ibuffer, 2);
    e1=ibuffer[0]; e2=ibuffer[1]; 
    injection->exy=e12=e1*e2;
    injection->side=bdata[i].side;
    switch(injection->side) {
    case 1:
      injection->ex=NoElx; injection->ey=NoEly;
      injection->e1=0;
      injection->dx=1; injection->dy=NoElx;
      break;
    case 2: 
      injection->ex=NoElx; injection->ey=NoEly;
      injection->e1=NoElxyz-NoElxy;
      injection->dx=1; injection->dy=NoElx;
      break;
    case 3:
      injection->ex=NoElx; injection->ey=NoElz;
      injection->e1=0;
      injection->dx=1; injection->dy=NoElxy;
      break;
    case 4: 
      injection->ex=NoElx; injection->ey=NoElz;
      injection->e1=NoElxy-NoElx;
      injection->dx=1; injection->dy=NoElxy;
      break;
    case 5:
      injection->ex=NoEly; injection->ey=NoElz;
      injection->e1=0;
      injection->dx=NoElx; injection->dy=NoElxy;
      break;
    case 6: 
      injection->ex=NoEly; injection->ey=NoElz;
      injection->e1=NoElx-1;
      injection->dx=NoElx; injection->dy=NoElxy;
      break;
    }
    if( processor==0 && (injection->inp.fail() || injection->ex!=e1 || injection->ey!=e2) )
      partrace->error(string("Invalid second line in file: ")+injection->fn);

    name=bdata[i].particlename;
    injection->particle=partrace->particles.FindParticleType(name);
    injection->mass=new double[injection->exy];
    do {
      found=readBoundaryInjection(injection);
    } while(found && injection->time2<simtime);
    injection->dirac=(injection->time1==injection->time2);
    if(found) {
      if(injection->time1<simtime && simtime<=injection->time2) {
	// reduce mass if partrace has been restarted.
	x=(injection->time2-simtime) / (injection->time2-injection->time1);
	for(j=0;j<e12;j++) injection->mass[j]*=x;
      }
      event=new Event("boundaryInjectionStart", injection->time1, 1);
      partrace->RunTime.insert(event);
      event=new Event("boundaryInjectionStop", injection->time2, 1);
      partrace->RunTime.insert(event);
    }
    injection++;
  }
  if(processor==0)
    cout<<"Initializing "<<boundaryInjectionNo<<" injection boundary profiles."<<endl;
}


void ElementsClass::InitGlobalConcentration(boundarydata* bdata, int no)
{
  int i;
  string name;
  if(globalConcentrationNo==0) return;
  partrace->SetNewHandlerText("InitGlobalConcentration");
  globalConcentration=new GlobalConcentrationClass[globalConcentrationNo];
  zLevels.init(NoElz);
  GlobalConcentrationClass *injection=globalConcentration;
  for(i=0; i<no; i++) {
    if(bdata[i].type==4 ) {
      injection->mass=new double[NoElz];
      injection->massz=new double[zLevels.numint];
      injection->volcon=new double[zLevels.numint];
      name=bdata[i].particlename;
      injection->particle=partrace->particles.FindParticleType(name);
      injection->particle->AlwaysCalculateConcentration=true;
      injection++;
    }
  }
  if(processor==0)
    cout<<"Initializing "<<globalConcentrationNo<<" global concentration boundaries."<<endl;
}


void ElementsClass::setAllGlobalConcentration(RunTimeClass &RunTime)
{
  if(globalConcentrationNo<=0) return;
  double dt=RunTime.dt;
  int e,n, i,m;

  // reset the masses for the z-intervals
  for(m=0; m<globalConcentrationNo; m++) {
    // reset outflow for all z-levels
    for(i=0;i<zLevels.numint;i++)
      globalConcentration[m].massz[i]=globalConcentration[m].volcon[i]=0.0;
  }

  // calculate outflowing masses.
  // front side
  //  cout<<"front"<<endl;
  e=n=0;
  calculateOutFlow(dt, e, n, 1,1, NoElx, 1, false, 3, 0,0,1,0);
  // back side
  //  cout<<"back"<<endl;
  e=NoElxy-NoElx;
  n=NoNoxy-NoNox;
  calculateOutFlow(dt, e, n, 1,1, NoElx, 1, true, 4, 0,NoEly-1,1,0);
  // left side
  //  cout<<"left"<<endl;
  e=n=0;
  calculateOutFlow(dt, e, n, NoElx,NoNox, NoEly, 0, false, 5, 0,0,0,1);
  // right side
  //  cout<<"right"<<endl;
  e=NoElx-1;
  n=NoNox-1;
  calculateOutFlow(dt, e, n, NoElx,NoNox, NoEly, 0, true, 6, NoElx-1,0,0,1);

  // calculate concentrations for each interval
  for(m=0; m<globalConcentrationNo; m++) {
    // reset outflow for all z-levels
    for(i=0;i<zLevels.numint;i++) {
      if(globalConcentration[m].volcon[i]>0.0)
	globalConcentration[m].volcon[i]=globalConcentration[m].massz[i]/globalConcentration[m].volcon[i];
    }
  }

  // inject particles at element sides with inflow.
  // front side
  //  cout<<"inflow front"<<endl;
  e=n=0;
  calculateInFlow(dt, e, n, 1,1, NoElx, 1, true, 3, 0,0,1,0);
  // back side
  //  cout<<"inflow back"<<endl;
  e=NoElxy-NoElx;
  n=NoNoxy-NoNox;
  calculateInFlow(dt, e, n, 1,1, NoElx, 1, false, 4, 0,NoEly-1,1,0);
  // left side
  //  cout<<"inflow left"<<endl;
  e=n=0;
  calculateInFlow(dt, e, n, NoElx,NoNox, NoEly, 0, true, 5, 0,0,0,1);
  // right side
  //  cout<<"inflow right"<<endl;
  e=NoElx-1;
  n=NoNox-1;
  calculateInFlow(dt, e, n, NoElx,NoNox, NoEly, 0, false, 6, NoElx-1,0,0,1);
}


void ElementsClass::calculateOutFlow(double dt, int e, int n, int e_inc, int n_inc,
				     int i2_to, int vcomp, bool plus, int side,
				     int ex, int ey, int ex_inc, int ey_inc)
{
  int ez,i2, ibegin,iend;
  int m,iz;
  double mass, v=0.0, vol, m0, area;
  double z,z1,z2, dz;
  int e1, n1, index, xyind;
  double outflow=0;
  double *pos=0;
  ElementsClass *elements=partrace->elements;
  GlobalConcentrationClass *globalconc;
  if(vcomp==0) pos=coords_x;
  else if(vcomp==1) pos=coords_y;
  //  cout<<"time="<<partrace->RunTime.Time<<endl;
  //cout<<"****outflow="<<outflow<<endl;
  for(i2=0;i2<i2_to;i2++) {
    e1=e;
    n1=n;
    //    cout<<"  element="<<e<<" node="<<n<<endl;
    for(m=0; m<globalConcentrationNo; m++) globalConcentration[m].summass=0.0;
    for(ez=0;ez<NoElz;ez++) {
      switch(elements->VelType) {
      case 1: // velocity defined on elements
	v=velocity->get(e1,vcomp);
	break;
      case 2: // velocity defined on nodes
	v=velocity->get(n1,vcomp);
	v+=velocity->get(n1+n_inc,vcomp);
	v+=velocity->get(n1+NoNoxy,vcomp);
	v+=velocity->get(n1+NoNoxy+n_inc,vcomp);
	v*=0.25;
	break;
      case 3: // velocity defined on finite volumes
	index=6*(ez+ey*NoElz+ex*NoElyz)+2*vcomp;
	if(vcomp==0) xyind=ex;
	else xyind=ey;
	v=(finitevolvelocity[index]*pos[xyind]
	   +finitevolvelocity[index+1]);
	break;
      }
      if( (plus && v>0.0) || (!plus && v<0.0) ) {
	area=areaOfElementSide(e1, ex,ey,ez, side);
	vol=area * fabs(v) * watercontent->get(e1) * dt;
	for(m=0; m<globalConcentrationNo; m++) {
	  mass=globalConcentration[m].particle->Concentration(e1) * vol;
	  //	  cout<<"mass="<<mass<<" conc="<<globalConcentration[m].particle->Concentration(e1)
	  //	      <<" vol="<<vol<<" area="<<area<<" v="<<v<<" dt="<<dt<<endl;
	  globalConcentration[m].mass[ez]=mass;
	  globalConcentration[m].summass+=mass;
	  zLevels.vol[ez]=vol;
	}
      }
      else {
	for(m=0; m<globalConcentrationNo; m++)
	  globalConcentration[m].mass[ez]=zLevels.vol[ez]=0.0;
      }
      e1+=NoElxy; n1+=NoNoxy;
    }
    //    cout<<"  summass="<<globalConcentration[0].summass<<endl;
    outflow+=globalConcentration[0].summass;
    // calculate the area with outflow
    globalconc=globalConcentration;
    for(m=0; m<globalConcentrationNo; m++,++globalconc) {
      m0=globalconc->summass*zLevels.nullmass;
      for(ez=0;ez<NoElz;ez++) {
	if(globalconc->mass[ez]>m0) break;
      }
      ibegin=ez;
      for(ez=NoElz-1;ez>=0;ez--) {
	if(globalconc->mass[ez]>m0) break;
      }
      iend=ez;
      //      cout<<"  ibegin="<<ibegin<<" iend="<<iend<<endl;
      if(iend>=ibegin) {
	z1=ibegin*length_z;
	z2=(iend+1)*length_z;
	//	cout<<"    z1="<<z1<<"   z2="<<z2<<endl;
	dz=z2-z1;
	iz=0;
	for(ez=ibegin;ez<=iend;ez++) {
	  z=( (double(ez)+0.5)*length_z-z1 ) / dz;
	  while(iz<zLevels.numint-1) {
	    if(z<zLevels.z[iz]) break;
	    iz++;
	  }
	  //	  cout<<"   iz="<<iz<<" ez="<<ez<<" mass="<<globalconc->mass[ez]<<endl;
	  globalconc->massz[iz]+=globalconc->mass[ez];
	  globalconc->volcon[iz]+=zLevels.vol[ez];
	}
      }
    }
    // go to next z-column.
    e+=e_inc; n+=n_inc;
    ex+=ex_inc; ey+=ey_inc;
  }
  //   cout<<"massz= "; m0=0;
  //   for(ez=0;ez<zLevels.numint;ez++) {
  //     cout<<' '<<globalConcentration[0].massz[ez];
  //     m0+=globalConcentration[0].massz[ez];
  //   }
  //   cout<<" sum="<<m0<<endl;
}


void ElementsClass::calculateInFlow(double dt, int e, int n, int e_inc, int n_inc,
				    int i2_to, int vcomp, bool plus, int side,
				    int ex, int ey, int ex_inc, int ey_inc)
{
  int ez, i2, ibegin,iend;
  int m,iz, npart;
  double mass, massperpart, massrest, v=0.0, m0, area;
  double z,z1,z2, dz;
  int e1, n1, index, xyind;
  double inflow=0;
  static int sumnpart=0;
  static int PElast=0;
  GlobalConcentrationClass *globalconc;
  ParticleClass *particle;
  ElementsClass *elements=partrace->elements;
  double *pos=0;
  if(vcomp==0) pos=coords_x;
  else if(vcomp==1) pos=coords_y;

  //   cout<<"time="<<partrace->RunTime.Time<<endl;
  //   cout<<"****inflow="<<inflow<<endl;
  for(i2=0;i2<i2_to;i2++) {
    e1=e;
    n1=n;
    //     cout<<"  element="<<e<<" node="<<n<<endl;
    for(m=0; m<globalConcentrationNo; m++) globalConcentration[m].summass=0.0;
    for(ez=0;ez<NoElz;ez++) {
      switch(elements->VelType) {
      case 1: // velocity defined on elements
	v=velocity->get(e1,vcomp);
	break;
      case 2: // velocity defined on nodes
	v=velocity->get(n1,vcomp);
	v+=velocity->get(n1+n_inc,vcomp);
	v+=velocity->get(n1+NoNoxy,vcomp);
	v+=velocity->get(n1+NoNoxy+n_inc,vcomp);
	v*=0.25;
	break;
      case 3: // velocity defined on finite volumes
	index=6*(ez+ey*NoElz+ex*NoElyz)+2*vcomp;
	if(vcomp==0) xyind=ex;
	else xyind=ey;
	v=(finitevolvelocity[index]*pos[xyind]
	   +finitevolvelocity[index+1]);
	break;
      }
      if( (plus && v>0.0) || (!plus && v<0.0) ) {
	area=areaOfElementSide(e1, ex,ey,ez, side);
	zLevels.ds[ez]=fabs(v) * dt;
	zLevels.vol[ez]=area * watercontent->get(e1) * zLevels.ds[ez];
	for(m=0; m<globalConcentrationNo; m++) {
	  mass=globalConcentration[m].particle->Concentration(e1) * zLevels.vol[ez];
	  // 	  cout<<"mass="<<mass<<" conc="<<globalConcentration[m].particle->Concentration(e1)
 	  //    <<" vol="<<zLevels.vol[ez]<<" area="<<area<<" v="<<v<<" dt="<<dt<<endl;
	  globalConcentration[m].mass[ez]=mass;
	  globalConcentration[m].summass+=mass;
	}
      }
      else {
	for(m=0; m<globalConcentrationNo; m++) globalConcentration[m].mass[ez]=0.0;
      }
      e1+=NoElxy; n1+=NoNoxy;
    }
    //    cout<<"  summass="<<globalConcentration[0].summass<<endl;
    inflow+=globalConcentration[0].summass;
    // calculate the area with inflow
    globalconc=globalConcentration;
    for(m=0; m<globalConcentrationNo; m++,++globalconc) {
      particle=globalconc->particle;
      massperpart=particle->massPerPart();
      m0=globalconc->summass*zLevels.nullmass;
      for(ez=0;ez<NoElz;ez++) {
	if(globalconc->mass[ez]>m0) break;
      }
      ibegin=ez;
      for(ez=NoElz-1;ez>=0;ez--) {
	if(globalconc->mass[ez]>m0) break;
      }
      iend=ez;
      //       cout<<"  ibegin="<<ibegin<<" iend="<<iend<<endl;
      if(iend>=ibegin) {
	z1=ibegin*length_z;
	z2=(iend+1)*length_z;
	// 	cout<<"    z1="<<z1<<"   z2="<<z2<<endl;
	dz=z2-z1;
	iz=0;
	massrest=0.0;
	for(ez=ibegin;ez<=iend;ez++) {
	  z=( (double(ez)+0.5)*length_z-z1 ) / dz;
	  while(iz<zLevels.numint-1) {
	    if(z<zLevels.z[iz]) break;
	    iz++;
	  }
	  mass=globalconc->volcon[iz]*zLevels.vol[ez]+massrest;
	  npart=int( mass / massperpart + 0.5);
	  massrest=mass-double(npart)*massperpart;
	  particle->Distribute(npart, PElast);
	  // 	  cout<<"   iz="<<iz<<" ez="<<ez<<" npart="<<npart<<" volcon="<<globalconc->volcon[iz]<<endl;
	  injectInElement(e+ez*NoElxy, ex,ey,ez, side, zLevels.ds[ez], npart, particle);
	  sumnpart+=npart;
	}
      } // loop over globalConcentration
    }
    // go to next z-column.
    e+=e_inc; n+=n_inc;
    ex+=ex_inc; ey+=ey_inc;
  }
  //   cout<<"sumnpart="<<sumnpart<<endl;
}


void ElementsClass::setBoundaryConcentration(BoundaryConcentrationClass *injection)
{
  int i,j, n,n1, nx,ny, dx,dy;
  nx=injection->nnpx;
  ny=injection->nnpy;
  if(injection->side<3) {
    // bottom=1, top=2
    n1=injection->x1+NoNox*injection->y1;
    injection->e1=injection->x1+NoElx*injection->y1;
    if(injection->side==2) {
      n1+=NoNoxy*NoElz;
      injection->e1+=NoElxy*NoElz_1;
    }
    dx=1; dy=NoNox;
    injection->ex=1; injection->ey=NoElx;
    injection->component=2;
  }
  else if(injection->side<5) {
    // front=3, back=4
    n1=injection->x1+NoNoxy*injection->y1;
    injection->e1=injection->x1+NoElxy*injection->y1;
    if(injection->side==4) {
      n1+=NoNox*NoEly;
      injection->e1+=NoElx*NoEly_1;
    }
    dx=1; dy=NoNoxy;
    injection->ex=1; injection->ey=NoElxy;
    injection->component=1;
  }
  else {
    // left=5, right=6
    n1=NoNox*injection->x1+NoNoxy*injection->y1;
    injection->e1=NoElx*injection->x1+NoElxy*injection->y1;
    if(injection->side==6) {
      n1+=NoElx;
      injection->e1+=NoElx_1;
    }
    dx=NoNox; dy=NoNoxy;
    injection->ex=NoElx; injection->ey=NoElxy;
    injection->component=0;
  }
  injection->value /= partrace->parallel->cpus();
  if(injection->dirac) injection->value/=double(injection->nel);
  if(VelType==2) {
    // calculate the nodes for the injections.
    int inode=0;
    for(j=0;j<ny;j++) {
      n=n1+j*dy;
      for(i=0;i<nx;i++) {
	injection->node[inode]=n;
	inode++;
	n+=dx;
      }
    }
  }
  else injection->node[0]=n1;
}


bool ElementsClass::setAllBoundaryConcentration(RunTimeClass &RunTime, int mode)
{
  if(boundaryConcentrationNo<=0) return false;
  long e,e0,e1,e2, ex,ey,ez, n, nx;
  int side, index;
  double dt=RunTime.dt;
  double area, ds, v=0.0, wc;
  int injected=0;
  BoundaryConcentrationClass *injection=boundaryConcentration;
  ElementsClass *elements=partrace->elements;	
  vector3 pos;
  for(int i=0; i<boundaryConcentrationNo; i++, injection++) {
    if( ( ( injection->time1<RunTime.Time && RunTime.Time<=injection->time2 ) ||
	  ( injection->dirac && injection->time1==RunTime.oldTime )) &&
	( (injection->dirac && mode==1) || (!injection->dirac && mode==2) ) ) {
      injected=1;
      side=injection->side%2;
      nx=injection->nnpx;
      //       cout<<"time="<<RunTime.Time
      // 	  <<" conc="<<injection->value<<endl;
      // inject the particles.
      GetNodeCoordinates(injection->node[0], pos);
      for(n=e0=e2=0; e2<injection->nely; e2++,n++) {
	e=injection->e1+e2*injection->ey;
	for(e1=0; e1<injection->nelx; e1++,e0++,n++) {
	  GetElementXYZ(e, ex, ey, ez);
	  if(injection->dirac) {
	    ds=0.0;
	    injection->mass[e0]=injection->value;
	  }
	  else {
	    wc=watercontent->get(e);
	    // vol = area * v * dt * wc
	    // conc = mass / vol
	    // mass = conc * vol = conc * area * v * dt * wc
	    switch(elements->VelType) {
	    case 1:
	      v=velocity->get(e, injection->component);
	      break;
	    case 2:
	      v= velocity->get(injection->node[n], injection->component);
	      v+=velocity->get(injection->node[n+1], injection->component);
	      v+=velocity->get(injection->node[n+nx], injection->component);
	      v+=velocity->get(injection->node[n+nx+1], injection->component);
	      v*=0.25;
	      break;
	    case 3:
	      index=6*(ez+ey*NoElz+ex*NoElyz)+2*injection->component;
	      v=(finitevolvelocity[index]*pos[injection->component]
		 +finitevolvelocity[index+1])/wc;
	      break;
	    }
	    if(v>=0.0) {
	      if(side==0) v=0.0;
	    }
	    else {
	      if(side==0) v=-v;
	      else v=0.0;
	    }
	    ds=dt*v;
	    area=areaOfElementSide(e, ex,ey,ez, injection->side);
	    injection->mass[e0]+=injection->value * area * ds * wc;
	  }
	  // 	  cout<<"element="<<e
	  // 	      <<" pe="<<processor
	  //  	      <<" nel="<<injection->nel
	  //  	      <<" ds="<<ds
	  //  	      <<" area="<<area
	  //  	      <<" mass="<<injection->mass[e0]
	  //  	      <<" wc="<<watercontent->get(e)
	  //  	      <<" v="<<v
	  //  	      <<" dt="<<dt
	  // 	      <<endl;
	  // 	      <<" n="<<n<<' '<<n+1<<' '<<n+nx<<' '<<n+nx+1
	  //  	      <<" nodes="<<injection->node[n]<<' '<<injection->node[n+1]<<' '<<injection->node[n+nx]<<' '<<injection->node[n+nx+1]
	  //  	      <<" v="<<velocity->get(injection->node[n], injection->component)<<' '<<velocity->get(injection->node[n+1], injection->component)<<' '<<velocity->get(injection->node[n+nx], injection->component)<<' '<<velocity->get(injection->node[n+nx+1], injection->component)
	  //  	      <<endl;
	  injectArea(injection->mass[e0], e, ex,ey,ez, injection->side, ds, injection->particle);
	  e+=injection->ex;
	}
      }
    }
  }
  partrace->parallel->sum(injected);
  return injected;
}


bool ElementsClass::doAllBoundaryInjections(RunTimeClass &RunTime)
{
  if(boundaryInjectionNo<=0) return false;
  double simtime=RunTime.Time;
  double oldsimtime=RunTime.oldTime;
  double dt=RunTime.dt;
  double m, mpp;
  double *mass;
  int e,e0,e1,e2, numpar;
  int PElast=0;
  int injected=0;
  BoundaryInjectionClass *injection=boundaryInjection;

  for(int i=0; i<boundaryInjectionNo; i++, injection++) {
    if( (injection->time1<simtime && simtime<=injection->time2) ||
	( injection->dirac && injection->time1==oldsimtime ) ) {
      // inject the particles.
      mpp=injection->particle->massPerPart();
      e0=injection->e1;
      mass=injection->mass;
      for(e2=0; e2<injection->ey; e2++) {
	e=e0;
	for(e1=0; e1<injection->ex; e1++) {
	  if(injection->dirac) m=*mass;
	  else m=dt/(injection->time2-oldsimtime) * *mass;
	  numpar=int(m/mpp+0.5);
	  if(numpar>0) {
	    injected=1;
	    *mass-=double(numpar)*mpp;
	    if(*mass<0.0) *mass=0.0;
	    injection->particle->Distribute(numpar, PElast);
	    injectAtElementSide(e, injection->side, numpar, injection->particle);
	    //	    cout<<"injecting at time="<<simtime<<" mass="<<m<<" no particles="<<numpar
	    //		<<" element="<<e<<" side="<<injection->side<<endl;
	  }
	  e+=injection->dx;
	  mass++;
	}
	e0+=injection->dy;
      }
    }
    if(simtime>=injection->time2) readBoundaryInjection(injection);
  }
  partrace->parallel->sum(injected);
  return injected;
}


bool ElementsClass::readBoundaryInjection(BoundaryInjectionClass *injection)
{ 
  int ok=1;
  double buffer[2];
  if(processor==0) {
    if(injection->inp.eof()) ok=0;
    else {
      injection->inp>>buffer[0]>>buffer[1];
      if(injection->inp.eof()) ok=0;
    }
  }
  partrace->parallel->broadcast( &ok, 1);
  if(ok==0) return false;
  partrace->parallel->broadcast( buffer, 2);
  injection->time1=buffer[0];
  injection->time2=buffer[1];
  if(processor==0) {
    for(int j=0; j<injection->exy; j++) {
      injection->inp>>injection->mass[j];
      if(injection->inp.eof()) 
	partrace->error(string("EOF in file: ")+injection->fn);
    }
  }
  partrace->parallel->broadcast( injection->mass, injection->exy);
  return true;
}


void ElementsClass::injectAtElementSide(int e, int side, int npart, ParticleClass* particle)
{
  if(npart<=0) return;
  long ex,ey,ez;
  GetElementXYZ(e, ex,ey,ez);
  double x1=coords_x[ex];
  double x2=coords_x[ex+1];
  double y1=coords_y[ey];
  double y2=coords_y[ey+1];
  double z1=coords_z[ez];
  double z2=coords_z[ez+1];
  RandomClass *rand=partrace->Random;
  Particle* p;
  ParticleList* plist=particle->firstParticleList()+e;
  for(int i=0;i<npart;i++) {
    p=particle->new_particle();

    if(side==5) p->position[0]=x1;
    else if(side==6) p->position[0]=x2;
    else p->position[0]=rand->draw(x1,x2);

    if(side==3) p->position[1]=y1;
    else if(side==4) p->position[1]=y2;
    else p->position[1]=rand->draw(y1,y2);

    if(side==1) p->position[2]=z1;
    else if(side==2) p->position[2]=z2;
    else p->position[2]=rand->draw(z1,z2);

    //    cout<<"element="<<e<<" pos="<<p->position[0]<<' '<<p->position[1]<<' '<<p->position[2]<<endl;
    plist->insert(p);
  }
}


void ElementsClass::injectArea(double &mass, int e, int ex, int ey, int ez, int side,
			       double ds, ParticleClass* particle)
{
  int npart=int(mass/particle->massPerPart()+0.5);
  if(npart<=0) return;
  mass-=double(npart)*particle->massPerPart();
  double x1=coords_x[ex];
  double x2=coords_x[ex+1];
  double y1=coords_y[ey];
  double y2=coords_y[ey+1];
  double z1=coords_z[ez];
  double z2=coords_z[ez+1];
  switch(side) {
  case 1: z2=min(z2, z1+ds); break;
  case 2: z1=max(z1, z2-ds); break;
  case 3: y2=min(y2, y1+ds); break;
  case 4: y1=max(y1, y2-ds); break;
  case 5: x2=min(x2, x1+ds); break;
  case 6: x1=max(x1, x2-ds); break;
  }
  Particle* p;
  ParticleList* plist=particle->firstParticleList()+e;
  RandomClass *rand=partrace->Random;
  //    cout<<"injecting "<<npart<<" particles in element "<<e
  //        <<" at pos "<<pos[0]<<' '<<pos[1]<<' '<<pos[2]
  //        <<" radius "<<radius[0]<<' '<<radius[1]<<' '<<radius[2]<<endl;
  for(int i=0;i<npart;i++) {
    p=particle->new_particle();
    if(x1==x2) p->position[0]=x1;
    else p->position[0]=rand->draw(x1,x2);
    if(y1==y2) p->position[1]=y1;
    else p->position[1]=rand->draw(y1,y2);
    if(z1==z2) p->position[2]=z1;
    else p->position[2]=rand->draw(z1,z2);
    plist->insert(p);
  }
}


void ElementsClass::injectInElement(long e, int npart,
				    ParticleClass* particle,
				    ParticleList* particlelist)
{
  if(npart<=0) return;
  long ex,ey,ez;
  GetElementXYZ(e, ex, ey, ez);
  double x1=coords_x[ex];
  double x2=coords_x[ex+1];
  double y1=coords_y[ey];
  double y2=coords_y[ey+1];
  double z1=coords_z[ez];
  double z2=coords_z[ez+1];
  RandomClass *rand=partrace->Random;
  Particle* p;
  ParticleList* plist=particlelist+e;
  for(int i=0;i<npart;i++) {
    p=particle->new_particle();
    p->position[0]=rand->draw(x1,x2);
    p->position[1]=rand->draw(y1,y2);
    p->position[2]=rand->draw(z1,z2);
    plist->insert(p);
  }
}


void ElementsClass::injectInElement(int e, int ex, int ey, int ez, int side,
				    double ds, int npart, ParticleClass* particle)
{
  if(npart<=0) return;
  double x1=coords_x[ex];
  double x2=coords_x[ex+1];
  double y1=coords_y[ey];
  double y2=coords_y[ey+1];
  double z1=coords_z[ez];
  double z2=coords_z[ez+1];
  switch(side) {
  case 1: z2=min(z2, z1+ds); break;
  case 2: z1=max(z1, z2-ds); break;
  case 3: y2=min(y2, y1+ds); break;
  case 4: y1=max(y1, y2-ds); break;
  case 5: x2=min(x2, x1+ds); break;
  case 6: x1=max(x1, x2-ds); break;
  }
  RandomClass *rand=partrace->Random;
  Particle* p;
  ParticleList* plist=particle->firstParticleList()+e;
  for(int i=0;i<npart;i++) {
    p=particle->new_particle();
    if(x1==x2) p->position[0]=x1;
    else p->position[0]=rand->draw(x1,x2);
    if(y1==y2) p->position[1]=y1;
    else p->position[1]=rand->draw(y1,y2);
    if(z1==z2) p->position[2]=z1;
    else p->position[2]=rand->draw(z1,z2);
    plist->insert(p);
  }
}


double ElementsClass::areaOfElementSide(int e, int ex, int ey, int ez, int side)
{
  switch(side) {
  case 1: 
  case 2: return length_x*length_y;
  case 3: 
  case 4: return length_x*length_z;
  case 5:
  case 6: return length_y*length_z;
  default: return 0;
  }  
}


void ElementsClass::InitDispersion(const string& fn, int usefield, vector3& disp,
				   double diffusion, vector3& dispweight, double poros, double bd, int difftype, int convectionmethod, int dispersionmethod, int coefatboundary, int secdispersivenew, int reflectionmethod, int nosurfacereflection, int elevxvy0, int timesplitting, int randomnumber, int recalculate, int SUB)
{
  partrace->SetNewHandlerText("soilprop");
  if(usefield==2) {
    // heterogeneous soil properties
    int ibuffer[4];
    ReadField(fn,"soil properties file",
	      "#partrace soil properties file",ibuffer,soilprop);
    partrace->Report("heterogeneous soil properties");
  }
  else {
    // homogeneous soil properties
    double data[5];
    soilprop=new ConstData(partrace->parallel, "soilprop", NoElxyz, 5);
    data[0]=disp[0]; data[1]=disp[1]; data[2]=diffusion; // durch 3?
    data[3]=poros; data[4]=bd;
    soilprop->set_dataline(0,data);
    partrace->Report("homogeneous soil properties");
  }
  DispWeight[0]=dispweight[0];
  DispWeight[1]=dispweight[1];
  DispWeight[2]=dispweight[2];
  DiffusionType=difftype;
  ConvectionMethod=convectionmethod;
  DispersionMethod=dispersionmethod;
  CoefAtBoundary=coefatboundary;
  SecDispersiveNew=secdispersivenew;
  ReflectionMethod=reflectionmethod;
  NoSurfaceReflection=nosurfacereflection;
  EleVxVy0=elevxvy0;
  TimeSplitting=timesplitting;
  RandomNumber=randomnumber;
  Recalculate=recalculate;
  sub=SUB;
  switch(RandomNumber) {
  case 1:
    DispFak=6.0;
    cout<<"You draw uniformly distributed random number"<<endl;
    break;
  case 2:
    DispFak=2.0;
    cout<<"You draw normally distributed random number"<<endl;
    break;
  } // switch RandomNumber
}


void ElementsClass::InitVelocity(const string& fn, int usefield, int VType,
				 vector3& mean, vector3& variance, double timestepfactor, double maxtime)
{
  switch(partrace->PartraceCaller) {
  case 0: VelType=VType; break;
  case 1: VelType=2; break;
  case 2: VelType=3; break;
  }
  if(VelType==1) NoVxyz=NoElxyz; // velocity defined on elements
  else NoVxyz=NoNoxyz;           // velocity defined in nodes
  TimeStepFactor=timestepfactor;
  MaxTime=maxtime;  

  partrace->SetNewHandlerText("velocity");
  NextVTime=new Event("velocity", partrace->RunTime.Time);
  if((partrace->PartraceCaller==0 && usefield==2) || partrace->PartraceCaller==1) {
    // velocity from file or from parswms
    if(partrace->PartraceCaller==0) partrace->RunTime.insert(NextVTime);
    if(partrace->distributedData)
      velocity=new  DistData(partrace->parallel, "velocity", NoVxyz, 3);
    else
      velocity=new LocalData(partrace->parallel, "velocity", NoVxyz, 3);
    fnVelo=fn;
  }
  else if(partrace->PartraceCaller==0 && (variance[0]!=0.0 || variance[1]!=0.0 || variance[2]!=0.0)) {
    // velocity generated by random generator
    if(partrace->distributedData)
      velocity=new  DistData(partrace->parallel, "velocity", NoVxyz, 3);
    else
      velocity=new LocalData(partrace->parallel, "velocity", NoVxyz, 3);
    VelocityFieldGenerator(mean,variance);
    // velocity not time dependant
    NextVTime->deactivate();
    calculateMaxTimeStepSize();
  }
  else if(partrace->PartraceCaller==0) {
    // constant velocity
    velocity=new ConstData(partrace->parallel, "velocity", NoVxyz, 3);
    velocity->set_dataline(0, mean.data());
    NextVTime->deactivate();
    calculateMaxTimeStepSize();
  }
}


void ElementsClass::InitWaterContent(const string& fn, int predefined, int wctype)
{
  double *sprop;
  WCType=wctype;
  partrace->SetNewHandlerText("watercontent");
  NextWCTime=new Event("watercontent", partrace->RunTime.Time);
  partrace->RunTime.insert(NextWCTime);
  if(partrace->PartraceCaller) {
    if(partrace->distributedData)
      watercontent=new DistData(partrace->parallel, "watercontent", NoElxyz, 1);
    else
      watercontent=new LocalData(partrace->parallel, "watercontent", NoElxyz, 1);
    NextWCTime->deactivate(); // water content not time dependant
  }
  else if(predefined==1) {
    // saturated wc
    if(soilprop->constant()) {
      // soil properties homogenous
      watercontent=new ConstData(partrace->parallel, "watercontent", NoElxyz, 1);
      sprop=soilprop->get_line(0);
      // set watercontent constant to (constant) porosity
      watercontent->set_data(0, 0, sprop[3]);
      watercontentOnNodes=new ConstData(partrace->parallel, "watercontentOnNodes", NoNoxyz, 1);
      watercontentOnNodes->set_data(0, 0, sprop[3]);
    }
    else {
      // soil properties not constant
      if(soilprop->indexed()) {
	if(soilprop->distributed())
	  watercontent=new IndexedDistData(partrace->parallel, "watercontent",
					   NoElxyz, 1, soilprop->global_datalines());
	else
	  watercontent=new IndexedData(partrace->parallel, "watercontent",
				       NoElxyz, 1, soilprop->global_datalines());
	watercontent->set_index(soilprop->get_index());
      }
      else {
	if(soilprop->distributed())
	  watercontent=new DistData(partrace->parallel, "watercontent", NoElxyz, 1);
	else
	  watercontent=new LocalData(partrace->parallel, "watercontent", NoElxyz, 1);
      }
      // copy the porosity to the watercontent
      watercontent->copy(soilprop, 3);
    }
    NextWCTime->deactivate(); // water content not time dependant
  }
  else {
    // read water content from file
    if(partrace->distributedData)
      watercontent=new DistData(partrace->parallel, "watercontent", NoElxyz, 1);
    else
      watercontent=new LocalData(partrace->parallel, "watercontent", NoElxyz, 1);
    if(WCType==2) {
      if(partrace->distributedData)
	watercontentOnNodes=new DistData(partrace->parallel, "watercontentOnNodes", NoNoxyz, 1);
      else
	watercontentOnNodes=new LocalData(partrace->parallel, "watercontentOnNodes", NoNoxyz, 1);
    }
    fnWC=fn;
  }
}


void ElementsClass::VelocityFieldGenerator(vector3& mean, vector3& variance)
{
  int i;
  int from=velocity->begin();
  int to=velocity->end();
  double data[3];

  for(i=from;i<to;i++)  {
    data[0]=partrace->Random->gauss(mean[0], variance[0]);
    data[1]=partrace->Random->gauss(mean[1], variance[1]);
    data[2]=partrace->Random->gauss(mean[2], variance[2]);
    velocity->set_dataline(i, data);
  }
}


void ElementsClass::InitBTC(vector3 *btclist, int btcanz)
{
  BTCAnz=btcanz;
  if(BTCAnz>0) {
    partrace->SetNewHandlerText("InitBTC");
    BTCListElNo=new long[BTCAnz];
    // Saving the ElementNumber of each BTC in array BTCListElNo
    for(int i=0;i<BTCAnz;i++) {
      if( (BTCListElNo[i]=GetElementNo(btclist[i]))<0 ) 
	partrace->error("Position of BTC is outside the volume");
    }
  }
}


void ElementsClass::InitInternalReflections(int t, string200 fn)
{
  if(t) {
    string msg;
    InternalReflections=new InternalReflectionClass2(fn, NoElx, NoEly, NoElz);
    if(InternalReflections->ReadField(fn, msg))
      partrace->error(msg);
    else
      partrace->Report("internal reflections used");
  }
  else {
    InternalReflections=new InternalReflectionClass(fn, NoElx, NoEly, NoElz);
    partrace->Report("no internal reflections");
  }
}


void ElementsClass::InitElementSequence(int *esequence, int no)
{
  ElementSequenceNo=no;
  if(ElementSequenceNo>0) {
    partrace->SetNewHandlerText("InitSequence");
    ElementSequence=new ElementSequenceClass[ElementSequenceNo];
    int i,j;
    for(j=i=0; i<ElementSequenceNo; i++,j+=6) {
      if(esequence[j]<1 || esequence[j+2]<1 || esequence[j+4]<1 ||
	 esequence[j+1]>NoElx || esequence[j+3]>NoEly || esequence[j+5]>NoElz ||
	 esequence[j]>esequence[j+1] || esequence[j+2]>esequence[j+3] || esequence[j+4]>esequence[j+5])
	partrace->error("Invalid element number in element sequence.");
      ElementSequence[i].ele_x=esequence[j+1]-esequence[j]+1;
      ElementSequence[i].ele_y=esequence[j+3]-esequence[j+2]+1;
      ElementSequence[i].ele_z=esequence[j+5]-esequence[j+4]+1;
      ElementSequence[i].startindex=(esequence[j+4]-1)*NoElxy
	+(esequence[j+2]-1)*NoElx
	+esequence[j]-1;
    }
  }
}

void ElementsClass::InitFilenames(string& fn_conc, string& fn_btc, string& fn_mom,
				  string& fn_pos, string& fn_coord)
{
  fnConc=fn_conc; fnBTC=fn_btc; fnMom=fn_mom; fnParticles=fn_pos;
  fnCoordinates=fn_coord;
}
 

void ElementsClass::getInterpolationFunction(long e, long ex, long ey, long ez, long nodes[], double *f)
{
  // regular hexahedron (cuboid)
  // new coordinate system with origin in node0
  // x-axis 0-1, y-axis 0-2, z-axis 0-4
  // node numbering: 
  //         lower side 2 3  upper side 6 7
  //                    0 1             4 5
  // wanted trilinear function:
  // f(x,y,z)=f0*x*y*z + f1*x*y + f2*x*z + f3*y*z + f4*x + f5*y + f6*z + f7
  // calculating coefficients f0..f7:
  // f( 0, 0, 0)=f7=data0
  // f( 0, 0,dz)=f6*dz+f7=data4 ==> f6=(data4-data0)/dz
  // f( 0,dy, 0)=f5*dy+f7=data2 ==> f5=(data2-data0)/dy
  // f(dx, 0, 0)=f4*dx+f7=data1 ==> f4=(data1-data0)/dx
  // f( 0,dy,dz)=f3*dy*dz+f5*dy+f6*dz+f7=data6 ==> f3=(data6+data0-data2-data4)/(dy*dz)
  // f(dx, 0,dz)=f2*dx*dz+f4*dx+f6*dz+f7=data5 ==> f2=(data5+data0-data1-data4)/(dx*dz)
  // f(dx,dy, 0)=f1*dx*dy+f4*dx+f5*dy+f7=data3 ==> f1=(data3+data0-data1-data2)/(dx*dy)
  // f(dx,dy,dz)=f0*dx*dy*dz + f1*dx*dy + f2*dx*dz + f3*dy*dz + f4*dx + f5*dy + f6*dz + f7=data7
  //         ==> f0=(data7-f1*dx*dy-f2*dx*dz-f3*dy*dz-f4*dx-f5*dy-f6*dz-f7)/(dx*dy*dz)
  //               =(data1+data2+data4+data7-data0-data3-data5-data6)/(dx*dy*dz)
  double *v0=velocity->new_buffer();
  double *v1=velocity->new_buffer();
  double *v2=velocity->new_buffer();
  double *v3=velocity->new_buffer();
  double *v4=velocity->new_buffer();
  double *v5=velocity->new_buffer();
  double *v6=velocity->new_buffer();
  double *v7=velocity->new_buffer();
  double dx=coords_x[ex+1]-coords_x[ex];
  double dy=coords_y[ey+1]-coords_y[ey];
  double dz=coords_z[ez+1]-coords_z[ez];
  double dxyz=dx*dy*dz;
  // get velocity at nodes
  velocity->get_line(nodes[0], v0);
  velocity->get_line(nodes[1], v1);
  velocity->get_line(nodes[2], v2);
  velocity->get_line(nodes[3], v3);
  velocity->get_line(nodes[4], v4);
  velocity->get_line(nodes[5], v5);
  velocity->get_line(nodes[6], v6);
  velocity->get_line(nodes[7], v7);
  // x - component
  f[ 7]= v0[0];
  f[ 6]=(v4[0]-v0[0])/dz;
  f[ 5]=(v2[0]-v0[0])/dy;
  f[ 4]=(v1[0]-v0[0])/dx;
  f[ 3]=(v6[0]+v0[0]-v2[0]-v4[0])/(dy*dz);
  f[ 2]=(v5[0]+v0[0]-v1[0]-v4[0])/(dx*dz);
  f[ 1]=(v3[0]+v0[0]-v1[0]-v2[0])/(dx*dy);
  f[ 0]=(v1[0]+v2[0]+v4[0]+v7[0]-v0[0]-v3[0]-v5[0]-v6[0])/dxyz;
  // y - component
  f[15]= v0[1];
  f[14]=(v4[1]-v0[1])/dz;
  f[13]=(v2[1]-v0[1])/dy;
  f[12]=(v1[1]-v0[1])/dx;
  f[11]=(v6[1]+v0[1]-v2[1]-v4[1])/(dy*dz);
  f[10]=(v5[1]+v0[1]-v1[1]-v4[1])/(dx*dz);
  f[ 9]=(v3[1]+v0[1]-v1[1]-v2[1])/(dx*dy);
  f[ 8]=(v1[1]+v2[1]+v4[1]+v7[1]-v0[1]-v3[1]-v5[1]-v6[1])/dxyz;
  // z - component
  f[23]= v0[2];
  f[22]=(v4[2]-v0[2])/dz;
  f[21]=(v2[2]-v0[2])/dy;
  f[20]=(v1[2]-v0[2])/dx;
  f[19]=(v6[2]+v0[2]-v2[2]-v4[2])/(dy*dz);
  f[18]=(v5[2]+v0[2]-v1[2]-v4[2])/(dx*dz);
  f[17]=(v3[2]+v0[2]-v1[2]-v2[2])/(dx*dy);
  f[16]=(v1[2]+v2[2]+v4[2]+v7[2]-v0[2]-v3[2]-v5[2]-v6[2])/dxyz;
  // deleting buffers
  velocity->delete_buffer(v0);
  velocity->delete_buffer(v1);
  velocity->delete_buffer(v2);
  velocity->delete_buffer(v3);
  velocity->delete_buffer(v4);
  velocity->delete_buffer(v5);
  velocity->delete_buffer(v6);
  velocity->delete_buffer(v7);
}

void ElementsClass::getInterpolationFunctionWC(long e, long ex, long ey, long ez, long nodes[], double *f)
{
  // modify for water content 
  // regular hexahedron (cuboid)
  // new coordinate system with origin in node0
  // x-axis 0-1, y-axis 0-2, z-axis 0-4
  // node numbering: 
  //         lower side 2 3  upper side 6 7
  //                    0 1             4 5
  // wanted trilinear function:
  // f(x,y,z)=f0*x*y*z + f1*x*y + f2*x*z + f3*y*z + f4*x + f5*y + f6*z + f7
  // calculating coefficients f0..f7:
  // f( 0, 0, 0)=f7=data0
  // f( 0, 0,dz)=f6*dz+f7=data4 ==> f6=(data4-data0)/dz
  // f( 0,dy, 0)=f5*dy+f7=data2 ==> f5=(data2-data0)/dy
  // f(dx, 0, 0)=f4*dx+f7=data1 ==> f4=(data1-data0)/dx
  // f( 0,dy,dz)=f3*dy*dz+f5*dy+f6*dz+f7=data6 ==> f3=(data6+data0-data2-data4)/(dy*dz)
  // f(dx, 0,dz)=f2*dx*dz+f4*dx+f6*dz+f7=data5 ==> f2=(data5+data0-data1-data4)/(dx*dz)
  // f(dx,dy, 0)=f1*dx*dy+f4*dx+f5*dy+f7=data3 ==> f1=(data3+data0-data1-data2)/(dx*dy)
  // f(dx,dy,dz)=f0*dx*dy*dz + f1*dx*dy + f2*dx*dz + f3*dy*dz + f4*dx + f5*dy + f6*dz + f7=data7
  //         ==> f0=(data7-f1*dx*dy-f2*dx*dz-f3*dy*dz-f4*dx-f5*dy-f6*dz-f7)/(dx*dy*dz)
  //               =(data1+data2+data4+data7-data0-data3-data5-data6)/(dx*dy*dz)
  double *wc0=watercontent->new_buffer();
  double *wc1=watercontent->new_buffer();
  double *wc2=watercontent->new_buffer();
  double *wc3=watercontent->new_buffer();
  double *wc4=watercontent->new_buffer();
  double *wc5=watercontent->new_buffer();
  double *wc6=watercontent->new_buffer();
  double *wc7=watercontent->new_buffer();
  double dx=coords_x[ex+1]-coords_x[ex];
  double dy=coords_y[ey+1]-coords_y[ey];
  double dz=coords_z[ez+1]-coords_z[ez];
  double dxyz=dx*dy*dz;
  // get watercontent at nodes
  watercontent->get_line(nodes[0], wc0);
  watercontent->get_line(nodes[1], wc1);
  watercontent->get_line(nodes[2], wc2);
  watercontent->get_line(nodes[3], wc3);
  watercontent->get_line(nodes[4], wc4);
  watercontent->get_line(nodes[5], wc5);
  watercontent->get_line(nodes[6], wc6);
  watercontent->get_line(nodes[7], wc7);
  // x - component
  f[ 7]= wc0[0];
  f[ 6]=(wc4[0]-wc0[0])/dz;
  f[ 5]=(wc2[0]-wc0[0])/dy;
  f[ 4]=(wc1[0]-wc0[0])/dx;
  f[ 3]=(wc6[0]+wc0[0]-wc2[0]-wc4[0])/(dy*dz);
  f[ 2]=(wc5[0]+wc0[0]-wc1[0]-wc4[0])/(dx*dz);
  f[ 1]=(wc3[0]+wc0[0]-wc1[0]-wc2[0])/(dx*dy);
  f[ 0]=(wc1[0]+wc2[0]+wc4[0]+wc7[0]-wc0[0]-wc3[0]-wc5[0]-wc6[0])/dxyz;
  // y - component
  f[15]= wc0[1];
  f[14]=(wc4[1]-wc0[1])/dz;
  f[13]=(wc2[1]-wc0[1])/dy;
  f[12]=(wc1[1]-wc0[1])/dx;
  f[11]=(wc6[1]+wc0[1]-wc2[1]-wc4[1])/(dy*dz);
  f[10]=(wc5[1]+wc0[1]-wc1[1]-wc4[1])/(dx*dz);
  f[ 9]=(wc3[1]+wc0[1]-wc1[1]-wc2[1])/(dx*dy);
  f[ 8]=(wc1[1]+wc2[1]+wc4[1]+wc7[1]-wc0[1]-wc3[1]-wc5[1]-wc6[1])/dxyz;
  // z - component
  f[23]= wc0[2];
  f[22]=(wc4[2]-wc0[2])/dz;
  f[21]=(wc2[2]-wc0[2])/dy;
  f[20]=(wc1[2]-wc0[2])/dx;
  f[19]=(wc6[2]+wc0[2]-wc2[2]-wc4[2])/(dy*dz);
  f[18]=(wc5[2]+wc0[2]-wc1[2]-wc4[2])/(dx*dz);
  f[17]=(wc3[2]+wc0[2]-wc1[2]-wc2[2])/(dx*dy);
  f[16]=(wc1[2]+wc2[2]+wc4[2]+wc7[2]-wc0[2]-wc3[2]-wc5[2]-wc6[2])/dxyz;
  // deleting buffers
  watercontent->delete_buffer(wc0);
  watercontent->delete_buffer(wc1);
  watercontent->delete_buffer(wc2);
  watercontent->delete_buffer(wc3);
  watercontent->delete_buffer(wc4);
  watercontent->delete_buffer(wc5);
  watercontent->delete_buffer(wc6);
  watercontent->delete_buffer(wc7);
}


void ElementsClass::interpolateVelocity(vector3& x0, vector3& pos, double* f, vector3& v)
{
  static double p0,p1,p2;
  // transform global coordinates to local
  p0=pos[0]-x0[0];
  p1=pos[1]-x0[1];
  p2=pos[2]-x0[2];
  // interpolate =  f[0]*p0*p1*p2 + f[1]*p0*p1 + f[2]*p0*p2 + f[3]*p1*p2
  //              + f[4]*p0 + f[5]*p1 + f[6]*p2 + f[7]
  v[0]= p0 * ( p1 * ( f[ 0]*p2 + f[ 1] ) + f[ 2]*p2 + f[ 4] )
    + p1 * ( f[ 3]*p2 + f[ 5] ) + f[ 6]*p2 + f[ 7];
  v[1]= p0 * ( p1 * ( f[ 8]*p2 + f[ 9] ) + f[10]*p2 + f[12] )
    + p1 * ( f[11]*p2 + f[13] ) + f[14]*p2 + f[15];
  v[2]= p0 * ( p1 * ( f[16]*p2 + f[17] ) + f[18]*p2 + f[20] )
    + p1 * ( f[19]*p2 + f[21] ) + f[22]*p2 + f[23];
}


void ElementsClass::interpolateWaterContent(vector3& x0, vector3& pos, double* f, double& wc)
{
  static double p0,p1,p2;
  // transform global coordinates to local
  p0=pos[0]-x0[0];
  p1=pos[1]-x0[1];
  p2=pos[2]-x0[2];
  // interpolate =  f[0]*p0*p1*p2 + f[1]*p0*p1 + f[2]*p0*p2 + f[3]*p1*p2
  //              + f[4]*p0 + f[5]*p1 + f[6]*p2 + f[7]
  wc= p0 * ( p1 * ( f[ 0]*p2 + f[ 1] ) + f[ 2]*p2 + f[ 4] )
    + p1 * ( f[ 3]*p2 + f[ 5] ) + f[ 6]*p2 + f[ 7];
}


void ElementsClass::gaujor(double** m) {
  double divi;

  //cout<<m[0][0]<<" "<<m[0][1]<<" "<<m[0][2]<<" "<<m[0][3]<<endl;
  //cout<<m[1][0]<<" "<<m[1][1]<<" "<<m[1][2]<<" "<<m[1][3]<<endl;
  //cout<<m[2][0]<<" "<<m[2][1]<<" "<<m[2][2]<<" "<<m[2][3]<<endl;
  // swaping because very low numbers are bad for division
  // swap row 1 and row2 or row3
  if (m[0][0]<m[1][0]) {
    swap(m[0][0],m[1][0]); swap(m[0][1],m[1][1]); swap(m[0][2],m[1][2]); swap(m[0][3],m[1][3]); }
  if (m[0][0]<m[2][0]) {
    swap(m[0][0],m[2][0]); swap(m[0][1],m[2][1]); swap(m[0][2],m[2][2]); swap(m[0][3],m[2][3]); }

  divi=m[0][0];
  m[0][0]/=divi;
  m[0][1]/=divi;
  m[0][2]/=divi;
  m[0][3]/=divi;

  divi=m[1][0];
  m[1][0]-=divi*m[0][0];
  m[1][1]-=divi*m[0][1];
  m[1][2]-=divi*m[0][2];
  m[1][3]-=divi*m[0][3];

  divi=m[2][0];
  m[2][0]-=divi*m[0][0];
  m[2][1]-=divi*m[0][1];
  m[2][2]-=divi*m[0][2];
  m[2][3]-=divi*m[0][3];

  // swap row 2 and row1 or row3
  if (m[1][1]<m[2][1]) {
    swap(m[1][0],m[2][0]); swap(m[1][1],m[2][1]); swap(m[1][2],m[2][2]); swap(m[1][3],m[2][3]); }

  divi=m[1][1];
  m[1][0]/=divi;
  m[1][1]/=divi;
  m[1][2]/=divi;
  m[1][3]/=divi;

  divi=m[2][1];
  m[2][0]-=divi*m[1][0];
  m[2][1]-=divi*m[1][1];
  m[2][2]-=divi*m[1][2];
  m[2][3]-=divi*m[1][3];

  divi=m[2][2];
  m[2][0]/=divi;
  m[2][1]/=divi;
  m[2][2]/=divi;
  m[2][3]/=divi;
  // matrix is now in row echelon form

  divi=m[0][2];
  m[0][0]-=divi*m[2][0];
  m[0][1]-=divi*m[2][1];
  m[0][2]-=divi*m[2][2];
  m[0][3]-=divi*m[2][3];

  divi=m[1][2];
  m[1][0]-=divi*m[2][0];
  m[1][1]-=divi*m[2][1];
  m[1][2]-=divi*m[2][2];
  m[1][3]-=divi*m[2][3];

  divi=m[0][1];
  m[0][0]-=divi*m[1][0];
  m[0][1]-=divi*m[1][1];
  m[0][2]-=divi*m[1][2];
  m[0][3]-=divi*m[1][3];
}
