#include <iostream>
#include <fstream>
#include <algorithm>
#include "elements.hpp"
#include "particleclass.hpp"
#include "random.hpp"
#include "partrace.hpp"
#include "parallel.hpp"
#include "distdata.hpp"

using namespace std;

ElementsClassRectilinear::ElementsClassRectilinear(triple noel, vector3& sidelen, PartraceClass *pt):
  ElementsClass(noel, sidelen, pt)
{
  GeometryType=4;
  partrace->SetNewHandlerText("ElementsClassRectilinear");
  if(partrace->distributedData)
    volume=new DistData(partrace->parallel, "volume", NoElxyz, 1);
  else
    volume=new LocalData(partrace->parallel, "volume", NoElxyz, 1);
}


ElementsClassRectilinear::~ElementsClassRectilinear()
{
}


long ElementsClassRectilinear::GetElementNo(const vector3& pos)
{ 
  static int ex,ey,ez;
  if(pos[0]<Boundaries[0] || pos[0]>Boundaries[3] || 
     pos[1]<Boundaries[1] || pos[1]>Boundaries[4] ||
     pos[2]<Boundaries[2] || pos[2]>Boundaries[5]) return -1;
  ex=int(lower_bound(coords_x, coords_x+NoNox, pos[0])-coords_x);
  if(ex>0) ex--;
  ey=int(lower_bound(coords_y, coords_y+NoNoy, pos[1])-coords_y);
  if(ey>0) ey--;
  ez=int(lower_bound(coords_z, coords_z+NoNoz, pos[2])-coords_z);
  if(ez>0) ez--;
  // return element number
  return ex+ey*NoElx+ez*NoElxy;
}


void ElementsClassRectilinear::ReadCoordinates()
{
  int type;
  triple nodes;
  int n;
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

  // read x coordinates
  if(processor==0) {
    for(n=0;n<NoNox;n++) {
      fin>>coords_x[n];
      if(fin.eof()) partrace->error(string("Unexpected EOF in coordinates file: ")+fnCoordinates);
      if(fin.fail()) partrace->error(string("Invalid input in coordinates file: ")+fnCoordinates);
      if(n>0 && coords_x[n]<=coords_x[n-1])
	partrace->error(string("Coordinates not in ascending order in coordinates file: ")+fnCoordinates);
    }
  }
  partrace->parallel->broadcast(coords_x, NoNox);
  Boundaries[0]=coords_x[0];
  Boundaries[3]=coords_x[NoNox-1];
  // read y coordinates
  if(processor==0) {
    for(n=0;n<NoNoy;n++) {
      fin>>coords_y[n];
      if(fin.eof()) partrace->error(string("Unexpected EOF in coordinates file: ")+fnCoordinates);
      if(fin.fail()) partrace->error(string("Invalid input in coordinates file: ")+fnCoordinates);
      if(n>0 && coords_y[n]<=coords_y[n-1])
	partrace->error(string("Coordinates not in ascending order in coordinates file: ")+fnCoordinates);
    }
  }
  partrace->parallel->broadcast(coords_y, NoNoy);
  Boundaries[1]=coords_y[0];
  Boundaries[4]=coords_y[NoNoy-1];
  // read z coordinates
  if(processor==0) {
    for(n=0;n<NoNoz;n++) {
      fin>>coords_z[n];
      if(fin.eof()) partrace->error(string("Unexpected EOF in coordinates file: ")+fnCoordinates);
      if(fin.fail()) partrace->error(string("Invalid input in coordinates file: ")+fnCoordinates);
      if(n>0 && coords_z[n]<=coords_z[n-1])
	partrace->error(string("Coordinates not in ascending order in coordinates file: ")+fnCoordinates);
    }
  }
  partrace->parallel->broadcast(coords_z, NoNoz);
  Boundaries[2]=coords_z[0];
  Boundaries[5]=coords_z[NoNoz-1];
  
  if(processor==0) fin.close();

  DomainSize[0]=Boundaries[3]-Boundaries[0];
  DomainSize[1]=Boundaries[4]-Boundaries[1];
  DomainSize[2]=Boundaries[5]-Boundaries[2];
  length_x=DomainSize[0]/NoElx;
  length_y=DomainSize[1]/NoEly;
  length_z=DomainSize[2]/NoElz;
}


void ElementsClassRectilinear::CalculateVolumes()
{
  // initialize the side lengths of the elements
  SideLength_x=new ConstData(partrace->parallel, "SideLength_x",NoElxyz,1);
  SideLength_x->init(length_x);
  SideLength_y=new ConstData(partrace->parallel, "SideLength_y",NoElxyz,1);
  SideLength_y->init(length_y);
  SideLength_z=new ConstData(partrace->parallel, "SideLength_z",NoElxyz,1); //MB
  SideLength_z->init(length_z);
  if(partrace->distributedData) {
    SideLength_x=new DistData(partrace->parallel, "SideLength_x",NoElxyz,1);
    SideLength_y=new DistData(partrace->parallel, "SideLength_y",NoElxyz,1);
    SideLength_z=new DistData(partrace->parallel, "SideLength_z",NoElxyz,1);
  }
  else {
    SideLength_x=new LocalData(partrace->parallel, "SideLength_x",NoElxyz,1);
    SideLength_y=new LocalData(partrace->parallel, "SideLength_y",NoElxyz,1);
    SideLength_z=new LocalData(partrace->parallel, "SideLength_z",NoElxyz,1);
  }
  double dx,dy,dz;
  long e, ex,ey,ez, e_from,e_to;
  e_from=volume->begin();
  e_to=volume->end();
  for(e=e_from; e<e_to; e++) {
    GetElementXYZ(e, ex,ey,ez);
    dx=coords_x[ex+1]-coords_x[ex];
    SideLength_x->set_data(e,0, dx);
    dy=coords_y[ey+1]-coords_y[ey];
    SideLength_y->set_data(e,0, dy);
    dz=coords_z[ez+1]-coords_z[ez];
    SideLength_z->set_data(e,0, dz);
    volume->set_data(e,0,  dx*dy*dz);
  }
}


double ElementsClassRectilinear::areaOfElementSide(int e, int ex, int ey, int ez, int side)
{
  switch(side) {
  case 1: 
  case 2: return ( coords_x[ex+1]-coords_x[ex] ) * ( coords_y[ey+1]-coords_y[ey] );
  case 3: 
  case 4: return ( coords_x[ex+1]-coords_x[ex] ) * ( coords_z[ez+1]-coords_z[ez] );
  case 5:
  case 6: return ( coords_y[ey+1]-coords_y[ey] ) * ( coords_z[ez+1]-coords_z[ez] );
  default: return 0;
  }  
}
