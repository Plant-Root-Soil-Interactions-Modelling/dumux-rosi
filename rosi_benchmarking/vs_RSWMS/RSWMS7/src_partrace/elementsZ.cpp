#include <iostream>
#include <fstream>
#include "elements.hpp"
#include "particleclass.hpp"
#include "random.hpp"
#include "partrace.hpp"
#include "parallel.hpp"
#include "distdata.hpp"
#include "runtime.hpp" //MB


using namespace std;

ElementsClassZ::ElementsClassZ(triple noel, vector3& sidelen, PartraceClass *pt):
  ElementsClass(noel, sidelen, pt), coordinates(0), xysides(0)
{
  GeometryType=2;
  partrace->SetNewHandlerText("ElementsClassZ");
  if(partrace->distributedData) {
    coordinates=new DistData(partrace->parallel, "coordinates", NoNoxyz, 3);
    volume=new DistData(partrace->parallel, "volume", NoElxyz, 1);
    xysides=new DistData(partrace->parallel, "xysides", NoElxyz+NoElxy, 4);
  }
  else {
    coordinates=new LocalData(partrace->parallel, "coordinates", NoNoxyz, 3);
    volume=new LocalData(partrace->parallel, "volume", NoElxyz, 1);
    xysides=new LocalData(partrace->parallel, "xysides", NoElxyz+NoElxy, 4);
  }
  transmat=new double[9];
}


ElementsClassZ::~ElementsClassZ()
{
  delete coordinates;
  delete xysides;
  delete transmat;
}


long ElementsClassZ::GetElementNo(const vector3& pos)
{ 
  static long ex,ey,ez, el;
  static int above;
  static double *pos_p;

  if(pos[0]<Boundaries[0]) return -1;
  if(pos[1]<Boundaries[1]) return -1;
  if(pos[2]<Boundaries[2]) return -1;
  ex=long( (pos[0]-Boundaries[0])/length_x );
  if(ex>=NoElx) return -1;
  ey=long( (pos[1]-Boundaries[1])/length_y );
  if(ey>=NoEly) return -1;
  ez=long( (pos[2]-Boundaries[2])/length_z );
  if(ez>=NoElz) return -1;

  el=ex+ey*NoElx+ez*NoElxy;
  if(el<0) return -1; // because of nan pos

  // search if pos is above top-surface or below bottom surface
  above=0;
  pos_p=const_cast<double*>(pos.data());
  while(Above(pos_p, el)) {
    above=1;
    el+=NoElxy;
    if(el>=NoElxyz) return -1;
  }
  if(above) return el;
  while(Below(pos_p, el)) {
    el-=NoElxy;
    if(el<0) return -1;
  }
  // return element number
  return el;
}


bool ElementsClassZ::Above(const double* pos, long element)
{
  double *area=xysides->get_line(element+NoElxy);
  return pos[0]*area[0]+pos[1]*area[1]+pos[2]*area[2]>area[3];
}


bool ElementsClassZ::Below(const double* pos, long element)
{
  double *area=xysides->get_line(element);
  return pos[0]*area[0]+pos[1]*area[1]+pos[2]*area[2]<area[3];
}


void ElementsClassZ::ReadCoordinates()
{
  int type;
  triple nodes;
  int n;
  double xmin=1e20, xmax=-xmin;
  double ymin=xmin, ymax=-xmin;
  double zmin=xmin, zmax=-xmin;
  ifstream fin;
  string line;

  if(processor==0) {
    fin.open(fnCoordinates.c_str());
    if(fin.fail())
      partrace->error(string("Error: The coordinates file: ")+fnCoordinates+" could not be opened");
    getline(fin,line);
    if(line!=string("trace geometry file"))
      partrace->error(string("Invalid first line in coordinates file: ")+fnCoordinates);
    fin>>type>>nodes[0]>>nodes[1]>>nodes[2];  fin.ignore(1000,'\n');
    if(CheckNEL(nodes,1))
      partrace->error(string("Invalid number of nodepoints in coordinates file: ")+fnCoordinates);
  }

  int to=coordinates->global_datalines();
  double *x=coordinates->start_buffered_distribute();
  for(n=0;n<to;n++) {
    if(processor==0) {
      fin>>x[0]>>x[1]>>x[2];
      if(fin.eof()) partrace->error(string("Unexpected EOF in coordinates file: ")+fnCoordinates);
      if(fin.fail()) partrace->error(string("Invalid input in coordinates file: ")+fnCoordinates);
      if(x[0]<xmin) xmin=x[0];
      if(x[0]>xmax) xmax=x[0];
      if(x[1]<ymin) ymin=x[1];
      if(x[1]>ymax) ymax=x[1];
      if(x[2]<zmin) zmin=x[2];
      if(x[2]>zmax) zmax=x[2];
    }
    x=coordinates->buffered_distribute(n);
  }
  coordinates->end_buffered_distribute();
  if(processor==0) fin.close();

  if(processor==0) {
    Boundaries[0]=xmin;
    Boundaries[3]=xmax;
    Boundaries[1]=ymin;
    Boundaries[4]=ymax;
    Boundaries[2]=zmin;
    Boundaries[5]=zmax;
  }
  partrace->parallel->broadcast(Boundaries, 6);
  if(processor==0) {
    cout<<"min x="<<Boundaries[0]<<" max x="<<Boundaries[3]<<endl;
    cout<<"min y="<<Boundaries[1]<<" max y="<<Boundaries[4]<<endl;
    cout<<"min z="<<Boundaries[2]<<" max z="<<Boundaries[5]<<endl;
  }
  DomainSize[0]=Boundaries[3]-Boundaries[0];
  DomainSize[1]=Boundaries[4]-Boundaries[1];
  DomainSize[2]=Boundaries[5]-Boundaries[2];
  length_x=DomainSize[0]/NoElx;
  length_y=DomainSize[1]/NoEly;
  length_z=DomainSize[2]/NoElz;
}


void ElementsClassZ::InitSides()
{
  setSides(xysides, 1, 2, 4, NoElx, NoEly, NoElz, 1L, NoElx, NoElxy);
}


void ElementsClassZ::CalculateVolumes()
{
  // initialize the side lengths of the elements
  SideLength_x=new ConstData(partrace->parallel, "SideLength_x",NoElxyz,1);
  SideLength_x->init(length_x);
  SideLength_y=new ConstData(partrace->parallel, "SideLength_y",NoElxyz,1);
  SideLength_y->init(length_y);
  if(partrace->distributedData)
    SideLength_z=new DistData(partrace->parallel, "SideLength_z",NoElxyz,1);
  else
    SideLength_z=new LocalData(partrace->parallel, "SideLength_z",NoElxyz,1);
  double delta_z;
  long nodes[8];
  double dxdy=length_x*length_y;
  long e, ex,ey,ez, e_from,e_to;
  e_from=volume->begin();
  e_to=volume->end();
  for(e=e_from; e<e_to; e++) {
    GetElementXYZ(e, ex,ey,ez);
    GetNodesOfElement(e, ex,ey,ez, nodes);
    delta_z=0.25*( coordinates->get(nodes[4],2)+coordinates->get(nodes[5],2)+
		   coordinates->get(nodes[6],2)+coordinates->get(nodes[7],2)-
		   coordinates->get(nodes[0],2)-coordinates->get(nodes[1],2)-
		   coordinates->get(nodes[2],2)-coordinates->get(nodes[3],2));
    SideLength_z->set_data(e,0,delta_z);
    volume->set_data(e,0, dxdy * delta_z);
  }
}


void ElementsClassZ::GetNodeCoordinates(long n, vector3& pos) 
{
  if(n<0 || n>NoNoxyz) partrace->error("Invalid node number for GetNodeCoordinates");
  double* coord=coordinates->get_line(n);
  pos[0]=coord[0]; pos[1]=coord[1]; pos[2]=coord[2];
}


void ElementsClassZ::GetStartCoordinates(long nodes[8], long ex, long ey, long ez, vector3& pos)
{
  // get the coordinates of the first node of the element
  double* coord=coordinates->get_line(nodes[0]);
  pos[0]=coord[0]; pos[1]=coord[1]; pos[2]=coord[2];
}


long ElementsClassZ::MoveParticle(long e, long ex, long ey, long ez, vector3& pos,
				  vector3& convection, DispersionClass& dispersionData)
{
  const double null=0.0;
  static vector3 newpos;
  static double x1,x2, y1,y2;
  static int side;
  static double diff,diff2, newdprod, dprod, disp1,disp2;
  static double *area;
  static vector3 dispersion, ds;
  dispersion=dispersionData.ds;

	//MB
	double dt=partrace->RunTime.dt;

  while(1) { // loop over boundary passings/reflections
    x1=ex*length_x; x2=x1+length_x;
    y1=ey*length_y; y2=y1+length_y;
    // newpos=pos+convection+dispersion
    newpos.setsum(pos, convection, dispersion);
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
    // position of point of intersection x_vec:
    // x_vec*n_vec=d  and  pos_vec+(newpos_vec-pos_vec)*a=x_vec
    //  ==> ( pos_vec+(newpos_vec-pos_vec)*a ) * n_vec = d
    // <==> a=(d-pos_vec*n_vec)/(newpos_vec*n_vec-pos_vec*n_vec)
    area=xysides->get_line(e);
    // newdprod=newpos * area
    newdprod=newpos.dotprod(area);
    if(newdprod<area[3]) {
      dprod=pos.dotprod(area);
      diff2=(area[3]-dprod)/(newdprod-dprod);
      if(diff2<diff) { diff=diff2; side=5; }
    }
    else {
      area=xysides->get_line(e+NoElxy);
      newdprod=newpos.dotprod(area);
      if(newdprod>area[3]) {
	dprod=pos.dotprod(area);
	diff2=(area[3]-dprod)/(newdprod-dprod);
	if(diff2<diff) { diff=diff2; side=6; }
      }
    }
    if(side==0) {
      pos=newpos;
      return e;
    }
    // pos+=diff*(convection+dispersion)
    pos.addxsum(diff, convection, dispersion);
    diff2=1.0-diff;
    convection*=diff2;
    dispersion*=diff2;
    switch(side) {
    case 1: // left
      if(ex==0) {
	if(BoundaryXmin->type(ey,ez)==1) {
	  reflectYZ(convection);
	  reflectYZ(dispersion);
	}
	else { // particle outside
	  ParticlesOutside[0]++;
	  return -1;
	}
      }
      else if(recalculateDispersion(dispersionData, e-1, e, ds, dt, dispersionData.particle->Retardation(e-1))) {
	if(dispersion[0]<null && reflectDispersion(dispersionData.ds[0], ds[0]))
	  reflectYZ(dispersion);
	else if(fabs(dispersionData.ds[0])>epsilon)
	  dispersion*=fabs(ds[0]/dispersionData.ds[0]);
	if(convection[0]+dispersion[0] <  0.0) {
	  dispersionData.ds=ds;
	  ex--;e--;
	}
      }
      else { ex--;e--; }
      break;
    case 2: // right
      if(ex==NoElx_1) {
	if(BoundaryXmax->type(ey,ez)==1) {
	  reflectYZ(convection);
	  reflectYZ(dispersion);
	}
	else  { // particle outside
	  ParticlesOutside[3]++;
	  return -1;
	}
      }
      else if(recalculateDispersion(dispersionData, e+1, e, ds, dt, dispersionData.particle->Retardation(e+1))) {
	if(dispersion[0]>null && reflectDispersion(dispersionData.ds[0], ds[0]))
	  reflectYZ(dispersion);
	else if(fabs(dispersionData.ds[0])>epsilon)
	  dispersion*=fabs(ds[0]/dispersionData.ds[0]);
	if(convection[0]+dispersion[0] >  null) {
	  dispersionData.ds=ds;
	  ex++;e++;
	}
      }
      else { ex++;e++; }
      break;
    case 3: // front
      if(ey==0) {
	if(BoundaryYmin->type(ex,ez)==1) {
	  reflectXZ(convection);
	  reflectXZ(dispersion);
	}
	else  { // particle outside
	  ParticlesOutside[1]++;
	  return -1;
 	}
      }
      else if(recalculateDispersion(dispersionData, e-NoElx, e, ds, dt, dispersionData.particle->Retardation(e-NoElx))) {
	if(dispersion[1]<null && reflectDispersion(dispersionData.ds[1], ds[1]))
	  reflectXZ(dispersion);
	else if(fabs(dispersionData.ds[1])>epsilon)
	  dispersion*=fabs(ds[1]/dispersionData.ds[1]);
	if(convection[1]+dispersion[1] <  null) {
	  dispersionData.ds=ds;
	  ey--;e-=NoElx;
	}
      }
      else { ey--;e-=NoElx; }
      break;
    case 4: // back
      if(ey==NoEly_1) {
	if(BoundaryYmax->type(ex,ez)==1) {
	  reflectXZ(convection);
	  reflectXZ(dispersion);
	}
	else  { // particle outside
	  ParticlesOutside[4]++;
	  return -1;
	}
      }
      else if(recalculateDispersion(dispersionData, e+NoElx, e, ds, dt, dispersionData.particle->Retardation(e+NoElx))) {
	if(dispersion[1]>null && reflectDispersion(dispersionData.ds[1], ds[1]))
	  reflectXZ(dispersion);
	else if(fabs(dispersionData.ds[1])>epsilon)
	  dispersion*=fabs(ds[1]/dispersionData.ds[1]);
	if(convection[1]+dispersion[1] >  null) {
	  dispersionData.ds=ds;
	  ey++;e+=NoElx;
	}
      }
      else { ey++;e+=NoElx; }
      break;
    case 5: // bottom
      if(ez==0) {
	if(BoundaryZmin->type(ex,ey)==1) {
	  reflect(convection, area);
	  reflect(dispersion, area);
	}
	else  { // particle outside
	  ParticlesOutside[2]++;
	  return -1;
	}
      }
      else if(recalculateDispersion(dispersionData, e-NoElxy, e, ds, dt, dispersionData.particle->Retardation(e-NoElxy))) {
	// (pos+dispersion)*n=pos*n+dispersion*n=d+dispersion*n
	disp1=dispersionData.ds.dotprod(area);
	disp2=ds.dotprod(area);
	if(reflectDispersion(disp1, disp2))
          reflect(dispersion, area);
	else if(fabs(disp1)>epsilon)
	  dispersion*=fabs(disp2/disp1);
	if(convection.dotprod(area)+dispersion.dotprod(area) <  null) {
	  dispersionData.ds=ds;
	  ez--;e-=NoElxy;
	}
      }
      else { ez--;e-=NoElxy; }
      break;
    case 6: // top
      if(ez==NoElz_1) {
	if(BoundaryZmax->type(ex,ey)==1) {
	  reflect(convection, area);
	  reflect(dispersion, area);
	}
	else  { // particle outside
	  ParticlesOutside[5]++;
	  return -1;
	}
      }
      else if(recalculateDispersion(dispersionData, e+NoElxy, e, ds, dt, dispersionData.particle->Retardation(e-NoElxy))) {
	disp1=dispersionData.ds.dotprod(area);
	disp2=ds.dotprod(area);
	if(reflectDispersion(disp1, disp2))
	  reflect(dispersion, area);
	else if(fabs(disp1)>epsilon)
	  dispersion*=fabs(disp2/disp1);
	if(convection.dotprod(area)+dispersion.dotprod(area) >  null) {
	  dispersionData.ds=ds;
	  ez++;e+=NoElxy;
	}
      }
      else { ez++;e+=NoElxy; }
      break;
    } // switch
  } // while
  return e;
}


void ElementsClassZ::initInterpolationFunction(long nodes[])
{
  // calculate transformation matrix for coordinate transformation within element
  // new coordinate system with origin in node0
  // x-axis 0-1, y-axis 0-2, z-axis 0-4:
  //   ex=x1-x0,   ey=x2-x0,   ez=x4-x0
  // node numbering: 
  //         lower side 2 3  upper side 6 7
  //                    0 1             4 5
  double x0,x1,x2, y0,y1,y2, z0,z1,z2;
  double* coord0=coordinates->new_buffer();
  coordinates->get_line(nodes[0], coord0);
  double* coord1=coordinates->get_line(nodes[1]);
  x0=coord1[0]-coord0[0];
  y0=coord1[1]-coord0[1];
  z0=coord1[2]-coord0[2];
  double* coord2=coordinates->get_line(nodes[2]);
  x1=coord2[0]-coord0[0];
  y1=coord2[1]-coord0[1];
  z1=coord2[2]-coord0[2];
  double* coord4=coordinates->get_line(nodes[4]);
  x2=coord4[0]-coord0[0];
  y2=coord4[1]-coord0[1];
  z2=coord4[2]-coord0[2];
  //       ( x0 x1 x2 )
  //  A =  ( y0 y1 y2 )
  //       ( z0 z1 z2 )
  double det=x0*y1*z2 + x1*y2*z0 + x2*y0*z1
           - x2*y1*z0 - x0*y2*z1 - x1*y0*z2;

  // calculate inverse matrix
  //   -1   ( Ad00 Ad01 Ad02 )T
  //  A   = ( Ad10 Ad11 Ad12 )  / det(A)
  //        ( Ad20 Ad21 Ad22 )
  transmat[0]= (y1*z2 - y2*z1) / det;
  transmat[1]=-(x1*z2 - x2*z1) / det;
  transmat[2]= (x1*y2 - x2*y1) / det;
  transmat[3]=-(y0*z2 - y2*z0) / det;
  transmat[4]= (x0*z2 - x2*z0) / det;
  transmat[5]=-(x0*y2 - x2*y0) / det;
  transmat[6]= (y0*z1 - y1*z0) / det;
  transmat[7]=-(x0*z1 - x1*z0) / det;
  transmat[8]= (x0*y1 - x1*y0) / det;
  // delete local buffer
  coordinates->delete_buffer(coord0);
}


void ElementsClassZ::getInterpolationFunction(long e, long ex, long ey, long ez, long nodes[], double *f)
{
  // hexahedron regular in x,y direction, not in z
  // new coordinate system with origin in node0
  // x-axis 0-1, y-axis 0-2, z-axis 0-4:
  //   ex=x1-x0,   ey=x2-x0,   ez=x4-x0
  // node numbering: 
  //         lower side 2 3  upper side 6 7
  //                    0 1             4 5
  // wanted trilinear function:
  // f(x,y,z)=f0*x*y*z + f1*x*y + f2*x*z + f3*y*z + f4*x + f5*y + f6*z + f7
  // calculating coefficients f0..f7:
  // f(0,0,0)=f7=data0
  // f(0,0,1)=f6+f7=data4 ==> f6=data4-data0
  // f(0,1,0)=f5+f7=data2 ==> f5=data2-data0
  // f(1,0,0)=f4+f7=data1 ==> f4=data1-data0
  // f(0,1,1)=f3*y6*z6+f5*y6+f6*z6+f7=data6 ==> f3=(data6-data0-f5*y6-f6*z6)/(y6*z6)
  // f(1,0,1)=f2*x5*z5+f4*x5+f6*z5+f7=data5 ==> f2=(data5-data0-f4*x5-f6*z5)/(x5*z5)
  // f(1,1,0)=f1*x3*y3+f4*x3+f5*y3+f7=data3 ==> f1=(data3-data0-f4*x3-f5*y3)/(x3*y3)
  // f(1,1,1)=f0*x*y*z+f1*x*y+f2*x*z+f3*y*z+f4*x+f5*y+f6*z+f7=data7
  //         ==> f0=(data7-data0-f1*x*y-f2*x*z-f3*y*z-f4*x-f5*y-f6*z)/(x*y*z)
  //         
  static double *v0, *v1, *v2, *v3, *v4, *v5, *v6, *v7;
  static double *coord0, *coord3, *coord5, *coord6, *coord7;
  static double x3,y3, x5,z5, y6,z6, x7,y7,z7;
  // allocate local buffers
  v0=velocity->new_buffer();
  v1=velocity->new_buffer();
  v2=velocity->new_buffer();
  v3=velocity->new_buffer();
  v4=velocity->new_buffer();
  v5=velocity->new_buffer();
  v6=velocity->new_buffer();
  v7=velocity->new_buffer();
  coord0=coordinates->new_buffer();
  coord3=coordinates->new_buffer();
  coord5=coordinates->new_buffer();
  coord6=coordinates->new_buffer();
  coord7=coordinates->new_buffer();
  // get velocity at nodes
  velocity->get_line(nodes[0], v0);
  velocity->get_line(nodes[1], v1);
  velocity->get_line(nodes[2], v2);
  velocity->get_line(nodes[3], v3);
  velocity->get_line(nodes[4], v4);
  velocity->get_line(nodes[5], v5);
  velocity->get_line(nodes[6], v6);
  velocity->get_line(nodes[7], v7);
  coordinates->get_line(nodes[0], coord0);
  coordinates->get_line(nodes[3], coord3);
  coordinates->get_line(nodes[5], coord5);
  coordinates->get_line(nodes[6], coord6);
  coordinates->get_line(nodes[7], coord7);
  // transform global coordinates to local
  x3=transmat[0]*(coord3[0]-coord0[0]) + transmat[1]*(coord3[1]-coord0[1]) + transmat[2]*(coord3[2]-coord0[2]);
  y3=transmat[3]*(coord3[0]-coord0[0]) + transmat[4]*(coord3[1]-coord0[1]) + transmat[5]*(coord3[2]-coord0[2]);
  x5=transmat[0]*(coord5[0]-coord0[0]) + transmat[1]*(coord5[1]-coord0[1]) + transmat[2]*(coord5[2]-coord0[2]);
  z5=transmat[6]*(coord5[0]-coord0[0]) + transmat[7]*(coord5[1]-coord0[1]) + transmat[8]*(coord5[2]-coord0[2]);
  y6=transmat[3]*(coord6[0]-coord0[0]) + transmat[4]*(coord6[1]-coord0[1]) + transmat[5]*(coord6[2]-coord0[2]);
  z6=transmat[6]*(coord6[0]-coord0[0]) + transmat[7]*(coord6[1]-coord0[1]) + transmat[8]*(coord6[2]-coord0[2]);
  x7=transmat[0]*(coord7[0]-coord0[0]) + transmat[1]*(coord7[1]-coord0[1]) + transmat[2]*(coord7[2]-coord0[2]);
  y7=transmat[3]*(coord7[0]-coord0[0]) + transmat[4]*(coord7[1]-coord0[1]) + transmat[5]*(coord7[2]-coord0[2]);
  z7=transmat[6]*(coord7[0]-coord0[0]) + transmat[7]*(coord7[1]-coord0[1]) + transmat[8]*(coord7[2]-coord0[2]);
  // x - component
  f[ 7]=v0[0];
  f[ 6]=v4[0]-v0[0];
  f[ 5]=v2[0]-v0[0];
  f[ 4]=v1[0]-v0[0];
  f[ 3]=(v6[0]-v0[0]-f[5]*y6-f[6]*z6) / (y6*z6);
  f[ 2]=(v5[0]-v0[0]-f[4]*x5-f[6]*z5) / (x5*z5);
  f[ 1]=(v3[0]-v0[0]-f[4]*x3-f[5]*y3) / (x3*y3);
  f[ 0]=(v7[0]-v0[0]-f[1]*x7*y7-f[2]*x7*z7-f[3]*y7*z7
	 -f[4]*x7-f[5]*y7-f[6]*z7) / (x7*y7*z7); 
  // y - component
  f[15]=v0[1];
  f[14]=v4[1]-v0[1];
  f[13]=v2[1]-v0[1];
  f[12]=v1[1]-v0[1];
  f[11]=(v6[1]-v0[1]-f[13]*y6-f[14]*z6) / (y6*z6);
  f[10]=(v5[1]-v0[1]-f[12]*x5-f[14]*z5) / (x5*z5);
  f[ 9]=(v3[1]-v0[1]-f[12]*x3-f[13]*y3) / (x3*y3);
  f[ 8]=(v7[1]-v0[1]-f[ 9]*x7*y7-f[10]*x7*z7-f[11]*y7*z7
	 -f[12]*x7-f[13]*y7-f[14]*z7) / (x7*y7*z7); 
  // z - component
  f[23]=v0[2];
  f[22]=v4[2]-v0[2];
  f[21]=v2[2]-v0[2];
  f[20]=v1[2]-v0[2];
  f[19]=(v6[2]-v0[2]-f[21]*y6-f[22]*z6) / (y6*z6);
  f[18]=(v5[2]-v0[2]-f[20]*x5-f[22]*z5) / (x5*z5);
  f[17]=(v3[2]-v0[2]-f[20]*x3-f[21]*y3) / (x3*y3);
  f[16]=(v7[2]-v0[2]-f[17]*x7*y7-f[18]*x7*z7-f[19]*y7*z7
	 -f[20]*x7-f[21]*y7-f[22]*z7) / (x7*y7*z7); 
  // delete local buffers
  velocity->delete_buffer(v0);
  velocity->delete_buffer(v1);
  velocity->delete_buffer(v2);
  velocity->delete_buffer(v3);
  velocity->delete_buffer(v4);
  velocity->delete_buffer(v5);
  velocity->delete_buffer(v6);
  velocity->delete_buffer(v7);
  coordinates->delete_buffer(coord0);
  coordinates->delete_buffer(coord3);
  coordinates->delete_buffer(coord5);
  coordinates->delete_buffer(coord6);
  coordinates->delete_buffer(coord7);
}


void ElementsClassZ::interpolateVelocity(vector3& x0, vector3& pos, double *f, vector3& v)
{
  static double p0,p1,p2, delta_x,delta_y,delta_z;
  // transform global coordinates to local
  delta_x=pos[0]-x0[0]; delta_y=pos[1]-x0[1]; delta_z=pos[2]-x0[2];
  p0=transmat[0]*delta_x + transmat[1]*delta_y + transmat[2]*delta_z;
  p1=transmat[3]*delta_x + transmat[4]*delta_y + transmat[5]*delta_z;
  p2=transmat[6]*delta_x + transmat[7]*delta_y + transmat[8]*delta_z;
  // interpolate =  f[0]*p0*p1*p2 + f[1]*p0*p1 + f[2]*p0*p2 + f[3]*p1*p2
  //              + f[4]*p0 + f[5]*p1 + f[6]*p2 + f[7]
  v[0]= p0 * ( p1 * ( f[ 0]*p2 + f[ 1] ) + f[ 2]*p2 + f[ 4] )
             + p1 * ( f[ 3]*p2 + f[ 5] ) + f[ 6]*p2 + f[ 7];
  v[1]= p0 * ( p1 * ( f[ 8]*p2 + f[ 9] ) + f[10]*p2 + f[12] )
             + p1 * ( f[11]*p2 + f[13] ) + f[14]*p2 + f[15];
  v[2]= p0 * ( p1 * ( f[16]*p2 + f[17] ) + f[18]*p2 + f[20] )
             + p1 * ( f[19]*p2 + f[21] ) + f[22]*p2 + f[23];
}


void ElementsClassZ::injectAtElementSide(int e, int side, int npart, ParticleClass* particle)
{
  if(npart<=0) return;
  long ex,ey,ez;
  long nodes[8];
  GetElementXYZ(e, ex,ey,ez);
  GetNodesOfElement(e, ex,ey,ez, nodes);
  //Horst Hardelauf 5.9.2008
  double* coord;
  double x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3;
  double a,b,c,d, x,z, h0,h1;
  int n0,n1,n2,n3;
  switch(side) {
  case 1: // bottom
    n0=0; n1=1; n2=2; n3=3;
    break;
  case 2: // top
    n0=4; n1=5; n2=6; n3=7;
    break;
  case 3: // front
    n0=0; n1=1; n2=4; n3=5;
    break;
  case 4: // back
    n0=2; n1=3; n2=6; n3=7;
    break;
  case 5: // left
    n0=0; n1=2; n2=4; n3=6;
    break;
  case 6: // right
    n0=1; n1=3; n2=5; n3=7;
    break;
  default: // should never occur
    n0=n1=n2=n3=0;
    break;
  }
  coord=coordinates->get_line(nodes[n0]);
  x0=coord[0]; y0=coord[1]; z0=coord[2];
  coord=coordinates->get_line(nodes[n1]);
  x1=coord[0]; y1=coord[1]; z1=coord[2];
  coord=coordinates->get_line(nodes[n2]);
  x2=coord[0]; y2=coord[1]; z2=coord[2];
  coord=coordinates->get_line(nodes[n3]);
  x3=coord[0]; y3=coord[1]; z3=coord[2];
  if(side<3) {
    // z(x,y) = a*x +b*y + c
    // z(x0,y0) = z0 = a*x0 + b*y0 + c
    // z(x1,y0) = z1 = a*x1 + b*y0 + c  => z1-z0 = a*(x1-x0)
    //                                 <=> a = (z1-z0) / (x1-x0)
    // z(x0,y2) = z2 = a*x0 + b*y2 + c  => z2-z0 = b*(y2-y0)
    //                                 <=> b = (z2-z0) / (y2-y0)
    //                                     c = z0 - a*x0 - b*y0
    a=(z1-z0)/(x1-x0);
    b=(z2-z0)/(y2-y0);
    c=z0-a*x0-b*y0;
    d=0;
  }
  else if(side<5) {
    // lower line
    // z(x) = a*x + b
    // z(x0) = a*x0 + b = z0
    // z(x1) = a*x1 + b = z1 => z1-z0=a*(x1-x0)
    // <=> a=(z1-z0)/(x1-x0) &&  b=z0-a*x0
    a=(z1-z0)/(x1-x0);
    b=z0-a*x0;
    // upper line
    c=(z3-z2)/(x3-x2);
    d=z2-c*x2;
  }
  else {
    a=(z1-z0)/(y1-y0);
    b=z0-a*y0;
    // upper line
    c=(z3-z2)/(y3-y2);
    d=z2-c*y2;
  }
  if(z0>z1) h0=z1;
  else h0=z0;
  if(z2>z3) h1=z2;
  else h1=z3;

  RandomClass *rand=partrace->Random;
  Particle* p;
  ParticleList* plist=particle->firstParticleList()+e;
  int i,n;
  for(i=0;i<npart;i++) {
    p=particle->new_particle();
    if(side<3) {
      p->position[0]=rand->draw(x0,x1);
      p->position[1]=rand->draw(y0,y2);
      p->position[2]=a*p->position[0] +b*p->position[1] + c;
    }
    else {
      if(side<5) {
	x=p->position[0]=rand->draw(x0,x1);
	p->position[1]=y0;
      }
      else {
	p->position[0]=x0;
	x=p->position[1]=rand->draw(y0,y1);
      }
      n=0;
      do {
	if(n>100) partrace->error("No convergence in ElementsClassZ::injectAtElementSide");
	z=rand->draw(h0,h1);
      } while(z<a*x+b || z>c*x+d);
      p->position[2]=z;
    }
    plist->insert(p);
  }
}


void ElementsClassZ::injectArea(double &mass, int e, vector3 &pos, vector3 &radius,
				ParticleClass* particle)
{
  partrace->error("injectArea not yet implemented in class ElementsClassZ.");
}


void ElementsClassZ::injectInElement(long e, int npart, ParticleClass* particle, ParticleList* particlelist)
{
  partrace->error("injectInElement not yet implemented in class ElementsClassZ.");
}


void ElementsClassZ::injectInElement(int e, int ex, int ey, int ez, int side,
				     double ds, int npart, ParticleClass* particle)
{
  if(npart<=0) return;
  static long nodes[8];
  GetNodesOfElement(e, ex,ey,ez, nodes);
  double x1=ex*length_x;
  double x2=x1+length_x;
  double y1=ey*length_y;
  double y2=y1+length_y;
  double z1=0;
  double z2=0;
  double a,z;
  double h1=coordinates->get(nodes[0], 2);
  double h2=coordinates->get(nodes[1], 2);
  double h3=coordinates->get(nodes[2], 2);
  double h4=coordinates->get(nodes[3], 2);
  double h5=coordinates->get(nodes[4], 2);
  double h6=coordinates->get(nodes[5], 2);
  double h7=coordinates->get(nodes[6], 2);
  double h8=coordinates->get(nodes[7], 2);

  switch(side) {
  case 1: // bottom 
  case 2: // top
    partrace->error("Element sides 1/2 not yet implemented in ElementsClassZ::injectInElement.");
    break;
  case 3: // front
    // minimum z value from  bottom
    if(h1<h2) z1=h1;
    else z1=h2;
    a=ds/(y2-y1);
    z=h1+a*(h3-h1);
    if(z<z1) z1=z;
    z=h2+a*(h4-h2);
    if(z<z1) z1=z;
    // maximum z value from top
    if(h5>h6) z2=h5;
    else z2=h6;
    z=h5+a*(h7-h5);
    if(z>z2) z2=z;
    z=h6+a*(h8-h6);
    if(z>z2) z2=z;
    y2=y1+ds;
    break;
  case 4: // back
    // minimum z value from bottom
    if(h3<h4) z1=h3;
    else z1=h4;
    a=ds/(y2-y1);
    z=h3-a*(h3-h1);
    if(z<z1) z1=z;
    z=h4-a*(h4-h2);
    if(z<z1) z1=z;
    // maximum z value from top
    if(h7>h8) z2=h7;
    else z2=h8;
    z=h7-a*(h7-h5);
    if(z>z2) z2=z;
    z=h8-a*(h8-h6);
    if(z>z2) z2=z;
    y1=y2-ds;
    break;
  case 5: // left
    // minimum z value from  bottom
    if(h1<h3) z1=h1;
    else z1=h3;
    a=ds/(x2-x1);
    z=h1+a*(h2-h1);
    if(z<z1) z1=z;
    z=h3+a*(h4-h3);
    if(z<z1) z1=z;
    // maximum z value from top
    if(h5>h7) z2=h5;
    else z2=h7;
    z=h5+a*(h6-h5);
    if(z>z2) z2=z;
    z=h7+a*(h8-h7);
    if(z>z2) z2=z;
    x2=x1+ds;
    break;
  case 6: // right
    // minimum z value from  bottom
    if(h2<h4) z1=h2;
    else z1=h4;
    a=ds/(x2-x1);
    z=h2-a*(h2-h1);
    if(z<z1) z1=z;
    z=h4-a*(h4-h3);
    if(z<z1) z1=z;
    // maximum z value from top
    if(h6>h8) z2=h6;
    else z2=h8;
    z=h6-a*(h6-h5);
    if(z>z2) z2=z;
    z=h8-a*(h8-h7);
    if(z>z2) z2=z;
    x1=x2-ds;
    break;
  }
  RandomClass *rand=partrace->Random;
  Particle* p;
  ParticleList* plist=particle->firstParticleList()+e;
  for(int i=0;i<npart;i++) {
    p=particle->new_particle();
    p->position[0]=rand->draw(x1,x2);
    p->position[1]=rand->draw(y1,y2);
    do {
      p->position[2]=rand->draw(z1,z2);
    } while(Above(p->position, e) || Below(p->position, e));
    plist->insert(p);
  }
}


double ElementsClassZ::areaOfElementSide(int e, int ex, int ey, int ez, int side)
{
  static long nodes[8];
  static double* coord;
  static double x0,y0,z0, x1,y1,z1, y2,z2, z3;

  GetNodesOfElement(e, ex,ey,ez, nodes);
  switch(side) {
  case 1: // bottom
    coord=coordinates->get_line(nodes[0]);
    x0=coord[0]; y0=coord[1];
    coord=coordinates->get_line(nodes[1]);
    x1=coord[0];
    coord=coordinates->get_line(nodes[2]);
                 y2=coord[1];
    return (x1-x0)*(y2-y0);
    break;
  case 2: // top
    coord=coordinates->get_line(nodes[4]);
    x0=coord[0]; y0=coord[1];
    coord=coordinates->get_line(nodes[5]);
    x1=coord[0];
    coord=coordinates->get_line(nodes[6]);
                 y2=coord[1];
    return (x1-x0)*(y2-y0);
    break;
  case 3: // front
    coord=coordinates->get_line(nodes[0]);
    x0=coord[0]; z0=coord[2];
    coord=coordinates->get_line(nodes[1]);
    x1=coord[0]; z1=coord[2];
    coord=coordinates->get_line(nodes[4]);
                 z2=coord[2];
    coord=coordinates->get_line(nodes[5]);
                 z3=coord[2];
    return (x1-x0)*0.5*(z3+z2-z1-z0);
    break;
  case 4: // back
    coord=coordinates->get_line(nodes[2]);
    x0=coord[0]; z0=coord[2];
    coord=coordinates->get_line(nodes[3]);
    x1=coord[0]; z1=coord[2];
    coord=coordinates->get_line(nodes[6]);
                 z2=coord[2];
    coord=coordinates->get_line(nodes[7]);
                 z3=coord[2];
    return (x1-x0)*0.5*(z3+z2-z1-z0);
    break;
  case 5: // left
    coord=coordinates->get_line(nodes[0]);
    y0=coord[1]; z0=coord[2];
    coord=coordinates->get_line(nodes[2]);
    y1=coord[1]; z1=coord[2];
    coord=coordinates->get_line(nodes[4]);
                 z2=coord[2];
    coord=coordinates->get_line(nodes[6]);
                 z3=coord[2];
    return (y1-y0)*0.5*(z3+z2-z1-z0);
    break;
  case 6: // right
    coord=coordinates->get_line(nodes[1]);
    y0=coord[1]; z0=coord[2];
    coord=coordinates->get_line(nodes[3]);
    y1=coord[1]; z1=coord[2];
    coord=coordinates->get_line(nodes[5]);
                 z2=coord[2];
    coord=coordinates->get_line(nodes[7]);
                 z3=coord[2];
    return (y1-y0)*0.5*(z3+z2-z1-z0);
    break;
  default: // should never occur
    return 0;
    break;
  }
}


void ElementsClassZ::setSides(ConstData *sides, int n2inc, int n3inc, int node0,
			      long e1_to,  long e2_to,  long e3_to, 
			      long e1_inc, long e2_inc, long e3_inc)
{
  double *x0=coordinates->new_buffer();
  double *x1=coordinates->new_buffer();
  double *x2=coordinates->new_buffer();
  vector3 a, b, n;
  long  e, e1,e2,e3, esave, n1,n2,n3, ex,ey,ez;
  double dprod;
  double res[4];
  long nodes[8];

  for(esave=e3=0;e3<e3_to;e3++) {
    for(e2=0;e2<e2_to;e2++) {
      e=e3*e3_inc+e2*e2_inc;
      for(e1=0;e1<e1_to;e1++) {
	GetElementXYZ(e, ex,ey,ez);
    	GetNodesOfElement(e, ex,ey,ez, nodes);
    	n1=nodes[0];
    	n2=nodes[n2inc];
    	n3=nodes[n3inc];
    	// one point in the plane
    	coordinates->get_line(n1, x0);
    	// first vector of the plane
    	coordinates->get_line(n2, x1);
    	a[0]=x1[0]-x0[0];
    	a[1]=x1[1]-x0[1];
    	a[2]=x1[2]-x0[2];
    	// second vector of the plane
    	coordinates->get_line(n3, x2);
    	b[0]=x2[0]-x0[0];
    	b[1]=x2[1]-x0[1];
    	b[2]=x2[2]-x0[2];
    	// vector orthogonal to plane (perpendicular)
    	n[0]=a[1]*b[2]-a[2]*b[1];
    	n[1]=a[2]*b[0]-a[0]*b[2];
    	n[2]=a[0]*b[1]-a[1]*b[0];
    	dprod=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  	if(dprod==0.0) partrace->error("Error in InitSides: Invalid element side(s).");
    	res[0]=n[0]/dprod;
    	res[1]=n[1]/dprod;
    	res[2]=n[2]/dprod;
    	// distance d from point x to plane: d=Dot_Product(x,n)
    	res[3]=x0[0]*res[0]+x0[1]*res[1]+x0[2]*res[2];
// 	cout<<"e1,e2,e2="<<e1<<' '<<e2<<' '<<e3<<" res="<<res[3]<<endl;
    	sides->set_dataline(esave++, res);
	e+=e1_inc;
      }
    }
  } 
  // handle last side
  e3=e3_to-1;
  for(e2=0;e2<e2_to;e2++) {
    e=e3*e3_inc+e2*e2_inc;
    for(e1=0;e1<e1_to;e1++) {
  	GetElementXYZ(e, ex,ey,ez);
  	GetNodesOfElement(e, ex,ey,ez, nodes);
  	n1=nodes[node0];
  	n2=nodes[node0+n2inc];
  	n3=nodes[node0+n3inc];
  	// one point in the plane
  	coordinates->get_line(n1, x0);
  	// first vector of the plane
  	coordinates->get_line(n2, x1);
  	a[0]=x1[0]-x0[0];
  	a[1]=x1[1]-x0[1];
  	a[2]=x1[2]-x0[2];
  	// second vector of the plane
  	coordinates->get_line(n3, x2);
  	b[0]=x2[0]-x0[0];
  	b[1]=x2[1]-x0[1];
  	b[2]=x2[2]-x0[2];
  	// vector orthogonal to plane (perpendicular)
  	n[0]=a[1]*b[2]-a[2]*b[1];
  	n[1]=a[2]*b[0]-a[0]*b[2];
  	n[2]=a[0]*b[1]-a[1]*b[0];
  	dprod=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  	res[0]=n[0]/dprod;
  	res[1]=n[1]/dprod;
  	res[2]=n[2]/dprod;
  	if(dprod==0.0) partrace->error("Error in InitSides: Invalid element side(s).");
  	// distance d from point x to plane: d=Dot_Product(x,n)
  	res[3]=x0[0]*res[0]+x0[1]*res[1]+x0[2]*res[2];
// 	cout<<"e1,e2,e2="<<e1<<' '<<e2<<' '<<e3<<" res="<<res[3]<<endl;
  	sides->set_dataline(esave++, res);
  	e+=e1_inc;
    }
  }
  // deleting buffers
  coordinates->delete_buffer(x0);
  coordinates->delete_buffer(x1);
  coordinates->delete_buffer(x2);
}
