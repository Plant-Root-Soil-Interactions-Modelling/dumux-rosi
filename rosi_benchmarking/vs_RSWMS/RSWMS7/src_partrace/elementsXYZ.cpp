#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "elements.hpp"
#include "particleclass.hpp"
#include "random.hpp"
#include "partrace.hpp"
#include "distdata.hpp"
#include "runtime.hpp" //MB

using namespace std;

ElementsClassXYZ::ElementsClassXYZ(triple noel, vector3& sidelen, PartraceClass *pt):
  ElementsClassZ(noel, sidelen, pt), xzsides(0), yzsides(0)
{
  GeometryType=3;
  partrace->SetNewHandlerText("ElementsClassXYZ");
  if(partrace->distributedData) {
    xzsides=new DistData(partrace->parallel, "xzsides", NoElxyz+NoElxz, 4);
    yzsides=new DistData(partrace->parallel, "yzsides", NoElxyz+NoElyz, 4);
  }
  else {
    xzsides=new LocalData(partrace->parallel, "xzsides", NoElxyz+NoElxz, 4);
    yzsides=new LocalData(partrace->parallel, "yzsides", NoElxyz+NoElyz, 4);
  }
}


ElementsClassXYZ::~ElementsClassXYZ()
{
  delete [] xzsides;
  delete [] yzsides;
}


long ElementsClassXYZ::GetElementNo(const vector3& pos)
{ 
  static long ex,ey,ez, el,e;
  static double prod1, prod2, diff, x;
  static int side;
 
  // calculate element number approximately 
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
  if(el<0) return -1; // position outside the (possibly extended) volume

  // search for the element which contains pos
  // element sides: 1=xmin 2=xmax  3=ymin 4=ymax  5=zmin 6=zmax
  //  cout<<"element="<<el<<" pos="<<pos[0]<<' '<<pos[1]<<' '<<pos[2]<<endl;
  double *area1=xysides->new_buffer();
  double *area2=xysides->new_buffer();
  while(el>=0) {
    diff=-1;
    // x-direction
    e=ex*NoElyz+ez*NoEly+ey;
    area1=yzsides->get_line(e);
    prod1=pos[0]*area1[0]+pos[1]*area1[1]+pos[2]*area1[2];
    area2=yzsides->get_line(e+NoElyz);
    prod2=pos[0]*area2[0]+pos[1]*area2[1]+pos[2]*area2[2];
    if(prod1<area1[3]) {
      x=area1[3]-prod1;
      if(x>diff) { diff=x; side=1; }
    }
    if(prod2>area2[3]) {
      x=prod2-area2[3];
      if(x>diff) { diff=x; side=2; }
    }
    // y-direction
    e=ey*NoElxz+ez*NoElx+ex;
    area1=xzsides->get_line(e);
    prod1=pos[0]*area1[0]+pos[1]*area1[1]+pos[2]*area1[2];
    area2=xzsides->get_line(e+NoElxz);
    prod2=pos[0]*area2[0]+pos[1]*area2[1]+pos[2]*area2[2];
    if(prod1<area1[3]) {
      x=area1[3]-prod1;
      if(x>diff) { diff=x; side=3; }
    }
    if(prod2>area2[3]) {
      x=prod2-area2[3];
      if(x>diff) { diff=x; side=4; }
    }
    // z-direction
    e=ez*NoElxy+ey*NoElx+ex;
    area1=xysides->get_line(e);
    prod1=pos[0]*area1[0]+pos[1]*area1[1]+pos[2]*area1[2];
    area2=xysides->get_line(e+NoElxy);
    prod2=pos[0]*area2[0]+pos[1]*area2[1]+pos[2]*area2[2];
    if(prod1<area1[3]) {
      x=area1[3]-prod1;
      if(x>diff) { diff=x; side=5; }
    }
    if(prod2>area2[3]) {
      x=prod2-area2[3];
      if(x>diff) { diff=x; side=6; }
    }

    if(diff<0.0) break;
    switch(side) {
    case 1:
      if(ex==0) el=-1;
      else { ex--; el--; }
      break;
    case 2: 
      if(ex==NoElx_1) el=-1;
      else { ex++; el++; }
      break;
    case 3:
      if(ey==0) el=-1;
      else { ey--; el-=NoElx; }
      break;
    case 4: 
      if(ey==NoEly_1) el=-1;
      else { ey++; el+=NoElx; }
      break;
    case 5:
      if(ez==0) el=-1;
      else { ez--; el-=NoElxy; }
      break;
    case 6: 
      if(ez==NoElz_1) el=-1;
      else { ez++; el+=NoElxy; }
      break;
    }
    //cout<<"side="<<side<<" new element="<<el<<endl;
  }
  // return element number
  xysides->delete_buffer(area1);
  xysides->delete_buffer(area2);
  return el;
}


void ElementsClassXYZ::CalculateVolumes()
{
  int i;
  long nodes[8];
  long e, ex,ey,ez, e_from,e_to;
  double vol, delta_s;
  double* coord[8];
  for(i=0;i<8;i++) coord[i]=coordinates->new_buffer();

  // initialize the side lengths of the elements
  if(partrace->distributedData) {
    SideLength_x=new DistData(partrace->parallel, "SideLength_x", NoElxyz, 1);
    SideLength_y=new DistData(partrace->parallel, "SideLength_y", NoElxyz, 1);
    SideLength_z=new DistData(partrace->parallel, "SideLength_z", NoElxyz, 1);
  }
  else {
    SideLength_x=new LocalData(partrace->parallel, "SideLength_x", NoElxyz, 1);
    SideLength_y=new LocalData(partrace->parallel, "SideLength_y", NoElxyz, 1);
    SideLength_z=new LocalData(partrace->parallel, "SideLength_z", NoElxyz, 1);
  }
  e_from=volume->begin();
  e_to=volume->end();
  for(e=e_from; e<e_to; e++) {
    GetElementXYZ(e, ex,ey,ez);
    GetNodesOfElement(e, ex,ey,ez, nodes);
    // decomposition of the hexahedron in 5 tetrahedrons
    vol = ( VolumeOfTetrahedron(nodes[0],nodes[1],nodes[2],nodes[4])
	  + VolumeOfTetrahedron(nodes[3],nodes[1],nodes[2],nodes[7])
	  + VolumeOfTetrahedron(nodes[5],nodes[1],nodes[4],nodes[7])
	  + VolumeOfTetrahedron(nodes[6],nodes[2],nodes[4],nodes[7])
	  + VolumeOfTetrahedron(nodes[1],nodes[2],nodes[4],nodes[7])) / 6.0;
    //cout<<"element="<<e<<" volume="<<vol<<endl;
    volume->set_data(e,0, vol);
    // calculate the side lengths of all elements
    // node numbering: 23 67
    //                 01 45
    for(i=0;i<8;i++) coord[i]=coordinates->get_line(nodes[i]);
    delta_s=(coord[1][0]+coord[3][0]+coord[5][0]+coord[7][0])
      -(coord[0][0]+coord[2][0]+coord[4][0]+coord[6][0]);
    SideLength_x->set_data(e,0, 0.25*delta_s);
    delta_s=(coord[2][1]+coord[3][1]+coord[6][1]+coord[7][1])
      -(coord[0][1]+coord[1][1]+coord[4][1]+coord[5][1]);
    SideLength_y->set_data(e,0, 0.25*delta_s);
    delta_s=(coord[4][2]+coord[5][2]+coord[6][2]+coord[7][2])
      -(coord[0][2]+coord[1][2]+coord[2][2]+coord[3][2]);
    SideLength_z->set_data(e,0, 0.25*delta_s);
  }
  // delete local buffers
  for(i=0;i<8;i++) coordinates->delete_buffer(coord[i]);
}


void ElementsClassXYZ::InitSides()
{
  cout<<"xy-sides"<<endl;
  setSides(xysides, 1, 2, 4, NoElx, NoEly, NoElz, 1L, NoElx, NoElxy);
  cout<<"xz-sides"<<endl;
  setSides(xzsides, 4, 1, 2, NoElx, NoElz, NoEly, 1L, NoElxy, NoElx);
  cout<<"yz-sides"<<endl;
  setSides(yzsides, 2, 4, 1, NoEly, NoElz, NoElx, NoElx, NoElxy, 1L);
}


long ElementsClassXYZ::MoveParticle(long el, long ex, long ey, long ez, vector3& pos,
				    vector3& convection, DispersionClass& dispersionData)
{
  // problem: Does the particle on the way from pos_vec to newpos_vec
  //          pass the plane area?
  // position of point of intersection x_vec:
  // x_vec*n_vec=d  and  pos_vec+(newpos_vec-pos_vec)*a=x_vec
  //  ==> ( pos_vec+(newpos_vec-pos_vec)*a ) * n_vec = d
  // <==> a=(d-pos_vec*n_vec)/(newpos_vec*n_vec-pos_vec*n_vec)
  //
  // algorithm: 
  //   1) calculate a for each element side
  //   2) 0<a<1 means particle leaves the element,
  //      smallest a in this range belongs to the first passed side.
  //   3) if a is equal for 2 or 3 sides, then the particle leaves the
  //      element through edge or corner, the new element must be found.
  //
  //   if(particle starts from boundary) {
  //     reflect;
  //     find_new_element;
  //   }
  //   if(particle leaves element) {
  //   }
  //   else set_position;
  const double null=0.0;
  static long neighbours[8];
  static vector3 newpos;
  static double diff, diff2, disp_old,disp_new;
  static long elem, elem_x,elem_y,elem_z, el_diff,ex_diff,ey_diff,ez_diff;
  static bool inside, reflected;
  static BoundaryPosition boundpos1[3], boundpos2[3], boundpos3[3];
  static BoundaryPosition *boundstart, *boundend, *boundact, *boundtemp;
  static int n_from, n_to, n_act, i,n;
  static vector3 dispersion, ds;

  dispersion=dispersionData.ds;

  //MB
  double dt=partrace->RunTime.dt;

  boundstart=boundpos1;
  boundend=boundpos2;
  boundact=boundpos3;
  //cout<<"element="<<el<<endl;

  // look if particle leaves the element
  // newpos=pos+convection+dispersion
  newpos.setsum(pos, convection, dispersion);
  //cout<<"startpos="<<pos[0]<<' '<<pos[1]<<' '<<pos[2]<<endl;
  //cout<<"newpos="<<newpos[0]<<' '<<newpos[1]<<' '<<newpos[2]<<endl;
  inside=checkBoundaries(el, ex, ey, ez, pos, newpos, true, boundstart, boundend);
  if(inside) {
    //cout<<"particle inside"<<endl;
    // particle does not leave element
    pos=newpos;
    return el;
  }
  n_from=sortBoundaries(boundstart);
  n_to=sortBoundaries(boundend);
  if(n_from==0) {  // particle does not start on boundary
    // set particle position to boundary
    diff=boundend[boundend[0].seq].a;
    // pos+=diff*(convection+dispersion)
    pos.addxsum(diff, convection, dispersion);
    diff2=1.0-diff;
    convection*=diff2;
    dispersion*=diff2;
    n_from=n_to;
    boundtemp=boundstart; boundstart=boundend; boundend=boundtemp;
    //cout<<"actpos="<<pos[0]<<' '<<pos[1]<<' '<<pos[2]<<endl;
  }

  // particle is on boundary now
  while(1) { // loop over boundary passings/reflections
    // do all required reflections
    //cout<<"n_from="<<n_from<<" seq="<<boundstart[0].seq<<" "<<boundstart[1].seq<<endl;
    reflected=false;
    el_diff=ex_diff=ey_diff=ez_diff=0;
    for(i=0;i<n_from;i++) { // reflect particle
      switch(boundstart[i].seq) {
      case 0: // x-axis
	if(boundstart[0].side==1) { // left global boundary
	  if(ex==0) {
	    //cout<<"reflect left"<<endl;
	    if(BoundaryXmin->type(ey,ez)==1) {
	      reflected=reflect(convection, boundstart[0].area);
	      reflected=reflect(dispersion, boundstart[0].area);
	    }
	    else { // particle outside
	      ParticlesOutside[0]++;
	      return -1;
	    }
	  }
	  else if(recalculateDispersion(dispersionData, el-1, el, ds, dt, dispersionData.particle->Retardation(el-1))) {
	    // particle goes into element with different dispersion
	    disp_old=dispersionData.ds.dotprod(boundstart[0].area);
	    disp_new=ds.dotprod(boundstart[0].area);
	    if(reflectDispersion(disp_old, disp_new))
	      reflect(dispersion, boundstart[0].area);
	    else if(fabs(disp_old)>epsilon)
	      dispersion*=fabs(disp_new/disp_old);
	    if(convection.dotprod(boundstart[0].area)+
	       dispersion.dotprod(boundstart[0].area) <  null) {
	      dispersionData.ds=ds;
	      ex_diff--;el_diff--;
	    }
	  }
	  else {
	    // dispersion does not change
	    ex_diff--;el_diff--;
	  }
	}
	else { // right
	  if(ex==NoElx_1) {
	    //cout<<"reflect right"<<endl;
	    if(BoundaryXmax->type(ey,ez)==1) {
	      reflected=reflect(convection, boundstart[0].area);
	      reflected=reflect(dispersion, boundstart[0].area);
	    }
	    else  { // particle outside
	      ParticlesOutside[3]++;
	      return -1;
	    }
	  }
	  else { ex_diff++;el_diff++; }
	}
	break;
      case 1: // y-axis
	if(boundstart[1].side==1) { // front
	  if(ey==0) {
	    //cout<<"reflect front"<<endl;
	    if(BoundaryYmin->type(ex,ez)==1) {
	      reflected=reflect(convection, boundstart[1].area);
	      reflected=reflect(dispersion, boundstart[1].area);
	    }
	    else  { // particle outside
	      ParticlesOutside[1]++;
	      return -1;
	    }
	  }
	  else { ey_diff--;el_diff-=NoElx; }
	}
	else { // back
	  if(ey==NoEly_1) {
	    //cout<<"reflect back"<<endl;
	    if(BoundaryYmax->type(ex,ez)==1) {
	      reflected=reflect(convection, boundstart[1].area);
	      reflected=reflect(dispersion, boundstart[1].area);
	    }
	    else  { // particle outside
	      ParticlesOutside[4]++;
	      return -1;
	    }
	  }
	  else { ey_diff++;el_diff+=NoElx; }
	}
	break;
      case 2: // z-axis
	if(boundstart[2].side==1) { // bottom
	  if(ez==0) {
	    //cout<<"reflect bottom"<<endl;
	    if(BoundaryZmin->type(ex,ey)==1) {
	      reflected=reflect(convection, boundstart[2].area);
	      reflected=reflect(dispersion, boundstart[2].area);
	    }
	    else  { // particle outside
	      ParticlesOutside[2]++;
	      return -1;
	    }
	  }
	  else { ez_diff--;el_diff-=NoElxy; }
	}
	else { // top
	  if(ez==NoElz_1) {
	    //cout<<"reflect top"<<endl;
	    if(BoundaryZmax->type(ex,ey)==1) {
	      reflected=reflect(convection, boundstart[2].area);
	      reflected=reflect(dispersion, boundstart[2].area);
	    }
	    else  { // particle outside
	      ParticlesOutside[5]++;
	      return -1;
	    }
	  }
	  else { ez_diff++;el_diff+=NoElxy; }
	}
	break;
      } // switch
    } // for
    if(reflected) {
      // particle was reflected, change newpos
      newpos.setsum(pos, convection, dispersion);
      //cout<<"after reflect newpos="<<newpos[0]<<' '<<newpos[1]<<' '<<newpos[2]<<endl;
    }
    if(n_from==1) {
      // element passes through element side
      //cout<<"passing element side"<<endl;
      el+=el_diff; ex+=ex_diff; ey+=ey_diff; ez+=ez_diff;
      // particle passes element side
      inside=checkBoundaries(el, ex, ey, ez, pos, newpos,
			     false, boundstart, boundend);
      if(inside) {
	// particle inside element el
	pos=newpos;
	return el;
      }
      n_to=sortBoundaries(boundend);
    }
    else {
      // particle passes through element edge or corner
      cout<<"time="<<partrace->RunTime.Time
	  <<" passing edge or corner, n_from="<<n_from<<endl;
      // search for the new element is required
      n=getNeighbourElements(ex, ey, ez, boundstart, n_from, neighbours);
      for(i=0;i<n;i++) {
	elem=neighbours[i];
	GetElementXYZ(elem, elem_x, elem_y, elem_z);
	inside=checkBoundaries(elem, elem_x, elem_y, elem_z, pos, newpos,
			       false, boundstart, boundact);
	if(inside) el=elem;
	else {
	  n_act=sortBoundaries(boundact);
	  if(boundact[boundact[0].seq].a<boundend[boundend[0].seq].a) {
	    boundtemp=boundend; boundend=boundact; boundact=boundtemp;
	    n_to=n_act;
	    el=elem;
	  }
	}
      }
      //cout<<"particle not in element"<<endl;
    }
    diff=boundend[boundend[0].seq].a;
    // pos+=diff*(convection+dispersion)
    pos.addxsum(diff, convection, dispersion);
    diff2=1.0-diff;
    convection*=diff2;
    dispersion*=diff2;
    n_from=n_to;
    boundtemp=boundstart; boundstart=boundend; boundend=boundtemp;

    //cout<<"diff="<<setprecision(15)<<diff<<" ex="<<ex<<" ey="<<ey<<" ez="<<ez
    //	<<" n_to="<<n_to<<"  side=";
    //for(i=0;i<n_to;i++) cout<<' '<<boundend[i].seq;
    // set particle position to boundary
    //cout<<endl<<"after move pos="<<pos[0]<<' '<<pos[1]<<' '<<pos[2]<<endl;

  } // loop over boundary passings/reflections
  return el;
}


bool ElementsClassXYZ::checkBoundaries(long el, long ex, long ey, long ez,
				       vector3& pos, vector3& newpos, bool first,
				       BoundaryPosition *boundstart,
				       BoundaryPosition *boundend)
{
  static long e;
  static double *area;
  static double a, prod, newprod;
  static bool inside;
  const double eps=epsilon*10.0;
  const double epsm=-eps;
  const double eps1=1.0-eps;

  inside=true;
  if(first) boundstart[0].a=boundstart[1].a=boundstart[2].a=2.0;
  boundend[0].a=boundend[1].a=boundend[2].a=2.0;
  // x-axis
  e=ex*NoElyz+ez*NoEly+ey;
  area=yzsides->get_line(e);
  prod=pos.dotprod(area);
  newprod=newpos.dotprod(area);
  if(fabs(prod-newprod)>eps) {
    a=(area[3]-prod)/(newprod-prod);
    if(first && epsm<a && a<eps ) {
      boundstart[0].a=a;
      boundstart[0].area=area;
      boundstart[0].side=1;
    }
    else if(eps<a && a<eps1) {
      inside=false;
      boundend[0].a=a;
      boundend[0].area=area;
      boundend[0].side=1;
    }
  }
  area=yzsides->get_line(e+NoElyz);
  prod=pos.dotprod(area);
  newprod=newpos.dotprod(area);
  if(fabs(prod-newprod)>eps) {
    a=(area[3]-prod)/(newprod-prod);
    if(first && epsm<a && a<eps) {
      boundstart[0].a=a;
      boundstart[0].area=area;
      boundstart[0].side=2;
    }
    else if(eps<a && a<eps1) {
      inside=false;
      boundend[0].a=a;
      boundend[0].area=area;
      boundend[0].side=2;
    }
  }
  // y-axis
  e=ey*NoElxz+ez*NoElx+ex;
  area=xzsides->get_line(e);
  prod=pos.dotprod(area);
  newprod=newpos.dotprod(area);
  if(fabs(prod-newprod)>eps) {
    a=(area[3]-prod)/(newprod-prod);
    if(first && epsm<a && a<eps ) {
      boundstart[1].a=a;
      boundstart[1].area=area;
      boundstart[1].side=1;
    }
    else if(eps<a && a<eps1) {
      inside=false;
      boundend[1].a=a;
      boundend[1].area=area;
      boundend[1].side=1;
    }
  }
  area=xzsides->get_line(e+NoElxz);
  prod=pos.dotprod(area);
  newprod=newpos.dotprod(area);
  if(fabs(prod-newprod)>eps) {
    a=(area[3]-prod)/(newprod-prod);
    if(first && epsm<a && a<eps ) {
      boundstart[1].a=a;
      boundstart[1].area=area;
      boundstart[1].side=2;
    }
    else if(eps<a && a<eps1) {
      inside=false;
      boundend[1].a=a;
      boundend[1].area=area;
      boundend[1].side=2;
    }
  }
  // z-axis
  e=ez*NoElxy+ey*NoElx+ex;
  area=xysides->get_line(e);
  prod=pos.dotprod(area);
  newprod=newpos.dotprod(area);
  if(fabs(prod-newprod)>eps) {
    a=(area[3]-prod)/(newprod-prod);
    if(first && epsm<a && a<eps ) {
      boundstart[2].a=a;
      boundstart[2].area=area;
      boundstart[2].side=1;
    }
    else if(eps<a && a<eps1) {
      inside=false;
      boundend[2].a=a;
      boundend[2].area=area;
      boundend[2].side=1;
    }
  }
  area=xysides->get_line(e+NoElxy);
  prod=pos.dotprod(area);
  newprod=newpos.dotprod(area);
  if(fabs(prod-newprod)>eps) {
    a=(area[3]-prod)/(newprod-prod);
    if(first && epsm<a && a<eps ) { 
      boundstart[2].a=a;
      boundstart[2].area=area;
      boundstart[2].side=2;
    }
    else if(eps<a && a<eps1) {
      inside=false;
      boundend[2].a=a;
      boundend[2].area=area;
      boundend[2].side=2;
    }
  }
  return inside;
}


int ElementsClassXYZ::sortBoundaries(BoundaryPosition *boundpos)
{
  const double eps=epsilon*10.0;
  // sort points of intersection
  if(boundpos[0].a<=boundpos[1].a) {
    if(boundpos[1].a<=boundpos[2].a) {
      boundpos[0].seq=0; boundpos[1].seq=1; boundpos[2].seq=2;
    }
    else if(boundpos[2].a<=boundpos[0].a) {
      boundpos[0].seq=2; boundpos[1].seq=0; boundpos[2].seq=1;
    }
    else {
      boundpos[0].seq=0; boundpos[1].seq=2; boundpos[2].seq=1;
    }
  }
  else {
    if(boundpos[0].a<=boundpos[2].a) {
      boundpos[0].seq=1; boundpos[1].seq=0; boundpos[2].seq=2;
    }
    else if(boundpos[2].a<=boundpos[1].a) {
      boundpos[0].seq=2; boundpos[1].seq=1; boundpos[2].seq=0;
    }
    else {
      boundpos[0].seq=1; boundpos[1].seq=2; boundpos[2].seq=0;
    }
  }
  // more than one point of intersection?
  if(boundpos[boundpos[0].seq].a==2.0) return 0; 
  if(fabs(boundpos[boundpos[0].seq].a-boundpos[boundpos[1].seq].a)<eps) 
    if(fabs(boundpos[boundpos[1].seq].a-boundpos[boundpos[2].seq].a)<eps) return 3;
    else return 2;
  else return 1;
}


int ElementsClassXYZ::getNeighbourElements(int ex, int ey, int ez,
					   BoundaryPosition *bound,
					   int nbound, long* neighbours)
{
  static int axis1,axis2, temp, side1,side2;
  static int x,y,z, n;
  static long ex1,ex2, ey1,ey2, ez1,ez2, ix,iy,iz;
  if(nbound==2) {
    axis1=bound[0].seq;
    axis2=bound[1].seq;
    if(axis1>axis2) {
      temp=axis1;
      axis1=axis2;
      axis2=temp;
    }
    side1=bound[axis1].side;
    side2=bound[axis2].side;
    if(axis1==0 && axis2==1) {
      if(side1==1) x=-1;
      else x=1;
      if(side2==1) y=-1;
      else y=1;
      z=0;
    }
    else if(axis1==0 && axis2==2) {
      if(side1==1) x=-1;
      else x=1;
      y=0;
      if(side2==1) z=-1;
      else z=1;
    }
    else {
      x=0;
      if(side1==1) y=-1;
      else y=1;
      if(side2==1) z=-1;
      else z=1;
    }
  }
  else { // nbound==3
    if(bound[0].side==1) x=-1;
    else x=1;
    if(bound[1].side==1) y=-1;
    else y=1;
    if(bound[2].side==1) z=-1;
    else z=1;
  }
  if(ex==0 && x==-1) x=0;
  if(ex==NoElx_1 && x==1) x=0;
  if(ey==0 && y==-1) y=0;
  if(ey==NoEly_1 && y==1) y=0;
  if(ez==0 && z==-1) z=0;
  if(ez==NoElz_1 && z==1) z=0;

  switch(x) {
  case -1: ex1=ex-1; ex2=ex; break;
  case  1: ex1=ex; ex2=ex+1; break;
  default: ex1=ex2=ex; break;
  }
  switch(y) {
  case -1: ey1=ey-1; ey2=ey; break;
  case  1: ey1=ey; ey2=ey+1; break;
  default: ey1=ey2=ey; break;
  }
  switch(z) {
  case -1: ez1=ez-1; ez2=ez; break;
  case  1: ez1=ez; ez2=ez+1; break;
  default: ez1=ez2=ez; break;
  }
  n=0;
  for(iz=ez1;iz<=ez2;iz++)
    for(iy=ey1;iy<=ey2;iy++)
      for(ix=ex1;ix<=ex2;ix++)
	neighbours[n++]=ix+iy*NoElx+iz*NoElxy;
  return n;
}


bool ElementsClassXYZ::isParticleInElement(vector3& pos, long ex, long ey, long ez)
{
  // element sides: 1=xmin 2=xmax  3=ymin 4=ymax  5=zmin 6=zmax
  static double *area;
  static long e;
  // x-direction
  e=ex*NoElyz+ez*NoEly+ey;
  area=yzsides->get_line(e);
  if(pos[0]*area[0]+pos[1]*area[1]+pos[2]*area[2]<area[3]) return false;
  area=yzsides->get_line(e+NoElyz);
  if(pos[0]*area[0]+pos[1]*area[1]+pos[2]*area[2]>area[3]) return false;
  // y-direction
  e=ey*NoElxz+ez*NoElx+ex;
  area=xzsides->get_line(e);
  if(pos[0]*area[0]+pos[1]*area[1]+pos[2]*area[2]<area[3]) return false;
  area=xzsides->get_line(e+NoElxz);
  if(pos[0]*area[0]+pos[1]*area[1]+pos[2]*area[2]>area[3]) return false;
  // z-direction
  e=ez*NoElxy+ey*NoElx+ex;
  area=xysides->get_line(e);
  if(pos[0]*area[0]+pos[1]*area[1]+pos[2]*area[2]<area[3]) return false;
  area=xysides->get_line(e+NoElxy);
  if(pos[0]*area[0]+pos[1]*area[1]+pos[2]*area[2]>area[3]) return false;
  return true;
}


void ElementsClassXYZ::injectAtElementSide(int e, int side, int npart, ParticleClass* particle)
{
  if(npart<0) return;
  partrace->error("Method injectAtElementSide not yet available for class ElementsClassXYZ!");
}


double ElementsClassXYZ::VolumeOfTetrahedron(int n0, int n1, int n2, int n3)
{
  // volume(tetrahedron)=x_vec * (y_vec X z_vec) / 6
  // division is not done, happens later
  static double *coord0, *coord1;
  static double x1,y1,z1, x2,y2,z2, x3,y3,z3;
  // allocate local buffer for coord0
  coord0=coordinates->new_buffer();
  coordinates->get_line(n0, coord0);
  coord1=coordinates->get_line(n1);
  x1=coord1[0]-coord0[0];
  y1=coord1[1]-coord0[1];
  z1=coord1[2]-coord0[2];
  coord1=coordinates->get_line(n2);
  x2=coord1[0]-coord0[0];
  y2=coord1[1]-coord0[1];
  z2=coord1[2]-coord0[2];
  coord1=coordinates->get_line(n3);
  x3=coord1[0]-coord0[0];
  y3=coord1[1]-coord0[1];
  z3=coord1[2]-coord0[2];
  coordinates->delete_buffer(coord0);
  return fabs(x1*y2*z3+x3*y1*z2+x2*y3*z1-x3*y2*z1-x1*y3*z2-x2*y1*z3);
}
