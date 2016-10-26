#ifndef ELEMENTSCLASS
#define ELEMENTSCLASS
#include "definitions.hpp"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include "dispersion.hpp"
#include "decay.hpp"
#include "uptake.hpp"


/** The partrace namespace.
 *  This namespace is used for all partrace objects. It was mainly created
 *  for the integration of partrace into parswms.
 */
class ConstDataIndex;
class ConstData;
class Event;
class ParticleClass;
class PartraceClass;
class RunTimeClass;
class ParticleList;

typedef std::map<std::string,ConstDataIndex*> indexliste; 

/**
 * Class for storing a sequence of elements with BTCs.
 * One sequence of elements if defined by a start element number,
 * a number of elements in x-, y- and z-direction. Every element
 * in the sequence is marked for having a breakthrough curve (BTC).
 */
class ElementSequenceClass {
public:
  int startindex; ///< start element number.
  int ele_x; ///<number of elements in x-direction
  int ele_y; ///<number of elements in y-direction
  int ele_z; ///<number of elements in z-direction
  ElementSequenceClass() {}  ///< constructor
  ~ElementSequenceClass() {} ///< destructor
};

/// Class for storing the input boundary values.
class boundarydata {
public:
  int type; ///< type of the boundary
  int side; ///< side for the boundary
  int x1;   ///< start element in x-direction
  int x2;   ///< end  element in x-direction
  int y1;   ///< start element in y-direction
  int y2;   ///< end element in y-direction
  double value;   ///< value for the concentration for boundary type 2.
  double InjectFrom, InjectTo; ///< injection time for boundary type 2.
  string200 particlename; ///< name of the corresponding particle.
  string200 injectfile;   ///< name of the file with the injection profiles.
};


/// Class for a injection boundary condition.
class BoundaryConcentrationClass {
public:
  int type; ///< type of the boundary
  int side; ///< side for the boundary
  int x1;   ///< start element in x-direction
  int x2;   ///< end  element in x-direction
  int y1;   ///< start element in y-direction
  int y2;   ///< end element in y-direction
  int e1;   ///< element number of the first element.
  int ex;   ///< element number increment in x-direction.
  int ey;   ///< element number increment in y-direction.
  int nnp;  ///< number of node points.
  int nel;  ///< number of elements.
  int nnpx; ///< number of node points in 1, direction.
  int nnpy; ///< number of node points in 2. direction.
  int nelx; ///< number of elements in 1, direction.
  int nely; ///< number of elements in 2. direction.
  int component; ///< use x,y or z-component (0,1,2) of the velocity vector.
  double value;  ///< value for the concentration for boundary type 2.
  double time1;  ///< injection start time.
  double time2;  ///< injection end time.
  bool dirac;    ///< true if start time == end time.
  int *node;     ///< array with the node numbers.
  double *mass;  ///< mass of each particle injection.
  ParticleClass *particle; ///< particle type.
  BoundaryConcentrationClass();
  virtual ~BoundaryConcentrationClass();
};


/// Class for time dependent injecting of particles at boundaries.
class BoundaryInjectionClass {
public:
  int side;  ///< side of the element, 1=bot, 2=top, 3=front, 4=back, 5=left, 6=right
  int e1;    ///< starting element.
  int ex,ey; ///< number of elements in first and second direction.
  int exy; 
  int dx,dy;                  ///< increment in first and second direction.
  double* mass;               ///< mass to inject for each time step.
  double time1,time2;         ///< start and end time of injection.
  bool dirac;                 ///< true if injection is dirac (all mass at once).
  bool ready;                 ///< all injections have been made.
  ParticleClass *particle;    ///< particle type.
  std::string fn;             ///< filename of the boundary injection file.
  std::ifstream inp;          ///< input file stream for the profile.
  BoundaryInjectionClass();
  virtual ~BoundaryInjectionClass();
};


/// Class for time dependent injecting of particles at boundaries.
class GlobalConcentrationClass {
public:
  double* mass;               ///< mass for each z-element  (dimension NoElz).
  double* massz;              ///< mass for each z-interval (dimension ZLevels.numint).
  double* volcon;             ///< volume,concentration for each z-interval (dimension ZLevels.numint).
  double summass;             ///< sum of mass for all z-elements.
  ParticleClass *particle;    ///< particle type.
  GlobalConcentrationClass();
  virtual ~GlobalConcentrationClass();
};


/// Class for time dependent injecting of particles at boundaries.
class ZLevels {
public:
  int numint;      ///< number of z-intervals.
  double* z;       ///< defines the z-intervals.
  double* vol;     ///< volume for each z-element  (dimension NoElz).
  double* ds;      ///< inflow distance for each z-element  (dimension NoElz).
  double nullmass; ///< defines the mass that is ignored.
  ZLevels();
  virtual ~ZLevels();
  void init(int nz);
};


/**
 * Class containing the boundary conditions for one plane.
 * One object of type BoundaryClass describes the boundary conditions
 * for one plane consisting of dimx X dimy element sides.
 */
class BoundaryClass {
protected:
  long dimx;         ///< number of element sides in x-direction.
  long dimy;         ///< number of element sides in y-direction.
  int* typedata;     ///< contains the boundary type for every element side.
  BoundaryClass() {} ///< default constructor is not used.
public:
  BoundaryClass(long ex, long ey); ///< constructor.
  virtual ~BoundaryClass();        ///< destructor.
  /// Access function returning the boundary type for element side (ex,ey).
  inline int type(long ex, long ey) const { return typedata[ex+ey*dimx]; }
  /// Access function setting the boundary type for element side (ex,ey).
  inline void settype(long ex, long ey, int t) { typedata[ex+ey*dimx]=t; }
};

/**
 * This class describes the position of a particle on an element side.
 * It is used in the case of irregular elements.
 */
class BoundaryPosition {
public:
  double a;     ///< travelled way of the particle (0<a<=1).
  double *area; ///< pointer to the area description.
  int side;     ///< the side of the element (1<=side<=6),
  int seq;      ///< used for ordering the passed element sides.
  BoundaryPosition() {}  ///< constructor.
  ~BoundaryPosition() {} ///< destructor.
};


/** Class for internal reflections-
 *  This class handles internal reflections.
 */
class InternalReflectionClass {
protected:
  int ex;
  int ey;
  int ez;
  int *elementtype;                 ///<  type of the elements
public:
  InternalReflectionClass(string200 fn, int x, int y, int z):
    ex(x), ey(y), ez(z), elementtype(0)
  { ///< constructor.
  }
  virtual ~InternalReflectionClass() {}           ///< destructor.
  virtual bool reflect(int x1, int y1, int z1, int x2, int y2, int z2) { return false; }
  virtual int ReadField(string200 fn, std::string& msg) { return 0; }
};

class InternalReflectionClass2: public InternalReflectionClass {
private:
public:
  InternalReflectionClass2(string200 fn, int x, int y, int z); ///< constructor.
  virtual ~InternalReflectionClass2();
  virtual bool reflect(int x1, int y1, int z1, int x2, int y2, int z2);
  virtual int ReadField(string200 fn, std::string& msg);
};
 

/** Class for regular geometries.
 *  This class descibes the a completely regular geometry.
 */
class ElementsClass {
public:
  long NoElx, NoEly, NoElz;       // Number of Elements
  long NoElx_1, NoEly_1, NoElz_1; // Number of Elements - 1
  long NoNox, NoNoy, NoNoz;       // Number of Nodes
  long NoElxyz, NoElxy, NoElxz, NoElyz;
  long NoNoxyz, NoNoxy;
  long NoVxyz;
  vector3 DomainSize;              // Sidelength of whole domain
  double Boundaries[6];           // xmin,ymin,zmin,xmax,ymax,zmax
  double *transmat;  // transformation matrix for element coordinates
  double *coords_x; // vector for the x-coordinates.
  double *coords_y; // vector for the y-coordinates.
  double *coords_z; // vector for the z-coordinates.
  BoundaryClass* BoundaryXmin;
  BoundaryClass* BoundaryXmax;
  BoundaryClass* BoundaryYmin;
  BoundaryClass* BoundaryYmax;
  BoundaryClass* BoundaryZmin;
  BoundaryClass* BoundaryZmax;
  int boundaryConcentrationNo;                       ///< number of concentration boundaries.
  BoundaryConcentrationClass *boundaryConcentration; ///< array of boundary concentrations.
  int boundaryInjectionNo;                           ///< number of boundary injections.
  BoundaryInjectionClass *boundaryInjection;         ///< array of boundary injections.
  int globalConcentrationNo;                         ///< number of global concentration boundary.
  GlobalConcentrationClass *globalConcentration;     ///< array of global concentration boundaries.
  ZLevels zLevels;                                   ///< defines the z-dependant concentrations.
  double length_x,length_y,length_z; ///< side length of one element
  ConstData *SideLength_x;           ///< x side lengths for all elements
  ConstData *SideLength_y;           ///< y side lengths for all elements
  ConstData *SideLength_z;           ///< z side lengths for all elements
  double epsilon;                    ///< epsilon for double
  double vMin;                       ///< minimum v greater than zero
  long ParticlesOutside[6];          ///< no. of particles leaving the volume
  ConstData *volume;                 ///< vector for volume of the elements
  InternalReflectionClass* InternalReflections;   ///<  type of the elements
  int processor;
  // BTC's
  int BTCAnz;
  long *BTCListElNo;
  int ElementSequenceNo;
  ElementSequenceClass *ElementSequence;
  std::string fnBTC, fnConc, fnMom, fnParticles, fnCoordinates;
  std::string fnVelo, fnWC;       // flow files
  // geometry
  int GeometryType;
  // Velocity
  std::ifstream velinp;         ///< input file stream for the velocity.
  std::ifstream wcinp;          ///< input file stream for the water content.
  ConstData *velocity;          // velocity for each element
  ConstData *watercontent;      // water content for each element
  ConstData *watercontentOnNodes; // water content for each mode.
  std::vector<double> finitevolvelocity; ///< velocity coefficients from finite volumes.
  std::vector<double> finitevolvelocityP; ///< velocity coefficients from finite volumes, Pollock type.
  std::vector<double> nodevelocity;
  std::vector<double> nodewc;
  Event *NextVTime;    ///< next time event for reading Velo files.
  Event *NextWCTime;   ///< next time event for reading WC files.
  int VFieldPos, VFieldRestartPos;    // file position for v-field
  int WCFieldPos, WCFieldRestartPos;  // file position for wc-field
  int VelType;                    // velocity on elements (=1) or nodes (=2)
  int WCType;                     // watercontent on elements (=1) or nodes (=2)
  int OutputType;                 ///< output type for the concentration.
  int RestartType;                ///< type of restart: 1=full  2=particles only.
  double TimeStepFactor;          // factor for the calculation of dt.
  double MaxTime; // manual setting of maximum time step, in case of diffusion dominated process //MB
  // soil properties
  vector3 DispWeight; // weights for the dispersion
  int DiffusionType; // Diffusion model number  1=simple 2=water content based
  int ConvectionMethod; // 1: no integration 2: integration
  int DispersionMethod; // Accounting for discontinuous dispersion, 1=off, 2=reflection method, 3=interpolation method
  int CoefAtBoundary; // 0: off 1: on
  int SecDispersiveNew; // 0: off 1: on
  int ReflectionMethod; // ReflectionCoefficient 1=Hoteit et al. 2002, 2=Lim 2006, 3=Bechtold 2010
  int NoSurfaceReflection;
  int EleVxVy0;
  int TimeSplitting; // Time Splitting 1=linear 2= nonlinear
  int RandomNumber; // RandomNumber 1=uniform,fast 2=gauss,slow,but superaccurate
  int Recalculate; // Recalculate Convection at interfaces 0=off, 1=on
  int sub;
  double DispFak;
  // 0=longitudal dispersivity (alphaL)
  // 1=transversal dispersivity (alphaT)
  // 2=diffusion coefficient
  // 3=porosity
  // 4=bulk density
  ConstData* soilprop;
  indexliste IndexListe; // map of indexes for ConstData's
  PartraceClass *partrace;

  // constructor / destructor
  ElementsClass(triple noel, vector3& sidelen, PartraceClass *partrace);
  virtual ~ElementsClass();

  // not virtual methods
  bool reflect(vector3& x, double* norm);
  inline void reflectXY(vector3& x) { x[2]=-x[2]; }
  inline void reflectXZ(vector3& x) { x[1]=-x[1]; }
  inline void reflectYZ(vector3& x) { x[0]=-x[0]; }
  inline void reflectXYZ(vector3& x) { x[0]=-x[0]; x[1]=-x[1]; x[2]=-x[2]; }
  void GetNodesOfElement(long el, long ex, long ey, long ez, long nodes[8]) const;
  void GetElementXYZ(long e, long& ex, long& ey, long& ez);
  void GetDxyzFxyz(long ex, long ey, long ez, vector3& pos, double& Dx, double& Dy, double& Dz, double& Fx, double& Fy, double& Fz);
  void GetTrilinearFactors(double Fx, double Fy, double Fz, double *f);
  void InterpolateVelocityTriLinear(long ex, long ey, long ez, double *f, int *ni, long *nodes, vector3& v);
  void InterpolateVelocityTriLinearDX(double *f, int *ni, long *nodes, vector3& v);
  void InterpolateWCTriLinear(double *f, int *ni, long *nodes, double& wc);
  void GetNodeXYZ(long n, long& nx, long& ny, long& nz); 
  void GetNoNo(long x, long y, long z, long nodes[8]);
  bool GetNoNo(const vector3& pos, long nodes[8]);
  inline long TotalNumber() const { return NoElxyz; }
  inline bool globalConcentrationBoundary() const { return globalConcentrationNo>0; }
  bool CheckNEL(triple n, int plus=0);
  void InitDispersion(const std::string& fn, int predefined, vector3& dispersion,
		      double diffusion, vector3& dispweight,
		      double poros, double bd, int difftype,  int convectionmethod, int dispersionmethod, int coefatboundary, int secdispersivenew, int reflectionmethod, int nosurfacereflection, int elevxvy0, int timesplitting, int randomnumber, int recalculate, int SUB);
  void InitVelocity(const std::string& fn, int predefined, int VelType,
		    vector3& mean, vector3& variance, double timestepfactor, double maxtime);
  void InitWaterContent(const std::string& fn, int predefined, int wctype);
  void InitFilenames(std::string& fn_conc, std::string& fn_btc, std::string& fn_mom, 
		     std::string& fn_pos, std::string& fn_coord);
  void InitBTC(vector3 *btclist, int BTCAnz);
  void InitElementSequence(int *esequence, int no);
  void InitBoundaries(boundarydata* bdata, int no);
  void VelocityFieldGenerator(vector3& mean, vector3& variance);
  void ReadField(const std::string& fn, const char *text, const char *head,
		 int ibuffer[4], ConstData* &iv);
  void GetFieldNo(int pro, int procs, int fields, int &von, int &bis);
  void ReadFlowFields(bool restart=false);
  void FiniteVolToNodes();
  void ElementsToNodes();
  void findFieldPos(std::ifstream &inp, double simtime);
  void ReadVField(bool restart=false);
  void ReadWCField(bool restart=false);
  void CountParticlesOutside();
  void calculateMaxTimeStepSize();
  void Dispersion(double dt, long e, vector3& pos, DispersionClass& dispersionData, double Retardation);
  void calculateConvection(vector3& convection, long e, vector3& pos, double dt, double Retardation);
  bool recalculateDispersion(DispersionClass& dispersionData,
			     long e_to, long e_from, vector3& pos, double dt, double Retardation); //MB
  bool reflectDispersion(double d1, double d2);
  int calculateBufferSize(int minsize=1024);
  void calculate_watercontentOnElements(ConstData *watercontentOnNodes);
  void InitBoundaryConcentration(boundarydata* bdata, int no);
  void InitBoundaryInjections(boundarydata* bdata, int no);
  void InitGlobalConcentration(boundarydata* bdata, int no);
  bool setAllBoundaryConcentration(RunTimeClass &RunTime, int mode);
  bool doAllBoundaryInjections(RunTimeClass &RunTime);
  void setAllGlobalConcentration(RunTimeClass &RunTime);
  void calculateOutFlow(double dt, int e, int n, int e_inc, int n_inc,
			int i2_to, int vcomp, bool plus, int side,
			int ex, int ey, int ex_int, int ey_int);
  void calculateInFlow( double dt, int e, int n, int e_inc, int n_inc,
		        int i2_to, int vcomp, bool plus, int side,
			int ex, int ey, int ex_int, int ey_int);
  bool readBoundaryInjection(BoundaryInjectionClass *injection);
  void gaujor(double** m);
  void InitInternalReflections(int t, string200 fn);
  // virtual methods
  virtual void GetElementsPosition(long i, vector3& pos);
  virtual bool Above(const vector3& pos, long element) { return true; } // not needed here
  virtual bool Below(const vector3& pos, long element) { return true; }
  virtual void GetStartCoordinates(long nodes[8], long ex, long ey, long ez, vector3& startcoord);
  virtual long GetElementNo(const vector3& pos);
  virtual void ReadCoordinates() {}
  virtual void CalculateVolumes();
  virtual void InitSides() {}
  virtual long MoveParticle(long e, long ex, long ey, long ez, vector3& pos,
			    DispersionClass& dispersion);
  virtual void initInterpolationFunction(long nodes[]) {} // not required for regular grids
  virtual void getInterpolationFunction(long e, long ex, long ey, long ez, long nodes[], double *interpol);
  virtual void getInterpolationFunctionWC(long e, long ex, long ey, long ez, long nodes[], double *interpolWC);
  virtual void interpolateVelocity(vector3& x0, vector3& pos, double* interpol, vector3& v);
  virtual void interpolateWaterContent(vector3& x0, vector3& pos, double* interpolWC, double& WC);
  virtual void setBoundaryConcentration(BoundaryConcentrationClass *injection);
  virtual void GetNodeCoordinates(long n, vector3& pos);
  virtual void injectArea(double &mass, int e, int ex, int ey, int ez, int side,
			  double ds, ParticleClass* particle);
  virtual void injectAtElementSide(int e, int side, int npart, ParticleClass* particle);
  virtual void injectInElement(int e, int ex, int ey, int ez, int side,
			       double ds, int npart, ParticleClass* particle);
  virtual void injectInElement(long e, int npart, ParticleClass* particle, ParticleList* particlelist);
  virtual double areaOfElementSide(int e, int ex, int ey, int ez, int side);
};


/** Class for rectilinear coordinates.
 *  This class descibes a geometry, which is rectilinear. That means that
 *  the position of the elements/nodes is defined by one vector for x-coordinates,
 *  on for the y-coordinates and one for the z-coordinates.
 */
class ElementsClassRectilinear: public ElementsClass
{
protected:
public:
  ElementsClassRectilinear(triple noel, vector3& sidelen, PartraceClass *partrace);
  virtual ~ElementsClassRectilinear();

  virtual long GetElementNo(const vector3& pos);
  virtual void ReadCoordinates();
  virtual void CalculateVolumes();
  virtual double areaOfElementSide(int e, int ex, int ey, int ez, int side);
};


/** Class for irregular z-values.
 *  This class descibes a geometry, which is regular in x- and y-direction,
 *  but irregular in z-direction.
 */
class ElementsClassZ: public ElementsClass
{
protected:
  ConstData *coordinates; // vector for coordinates
  ConstData *xysides;     // vector for the xy sides of the elements
public:
  ElementsClassZ(triple noel, vector3& sidelen, PartraceClass *partrace);
  virtual ~ElementsClassZ();

  void setSides(ConstData *sides, int n2inc, int n3inc, int node0,
		long e1_to,  long e2_to,  long e3_to, 
		long e1_inc, long e2_inc, long e3_inc);
  virtual bool Above(const double* pos, long element);
  virtual bool Below(const double* pos, long element);
  virtual void GetStartCoordinates(long nodes[8], long ex, long ey, long ez, vector3& pos);
  virtual long GetElementNo(const vector3& pos);
  virtual void ReadCoordinates();
  virtual void CalculateVolumes();
  virtual void InitSides();
  virtual long MoveParticle(long e, long ex, long ey, long ez, vector3& pos,
			    vector3& convection, DispersionClass& dispersion);
  virtual void initInterpolationFunction(long nodes[]);
  virtual void getInterpolationFunction(long e, long ex, long ey, long ez, long nodes[], double *interpol);
  virtual void interpolateVelocity(vector3& x0, vector3& pos, double* interpol, vector3& v);
  virtual void injectArea(double &mass, int e, vector3 &pos, vector3 &radius, ParticleClass* particle);
  virtual void injectAtElementSide(int e, int side, int npart, ParticleClass* particle);
  virtual void injectInElement(long e, int npart, ParticleClass* particle, ParticleList* particlelist);
  virtual void injectInElement(int e, int ex, int ey, int ez, int side,
			       double ds, int npart, ParticleClass* particle);
  virtual double areaOfElementSide(int e, int ex, int ey, int ez, int side);
  virtual void GetNodeCoordinates(long n, vector3& pos);
};


/** Class for irregular geometries.
 *  This class descibes a geometry, which is irregular in x-, y-
 *  and z-direction.
 */
class ElementsClassXYZ: public ElementsClassZ
{
protected:
  ConstData *xzsides; // vector for the xz sides of the elements
  ConstData *yzsides; // vector for the yz sides of the elements
public:
  ElementsClassXYZ(triple noel, vector3& sidelen, PartraceClass *partrace);
  virtual ~ElementsClassXYZ();

  double VolumeOfTetrahedron(int n0, int n1, int n2, int n3);
  virtual long GetElementNo(const vector3& pos);
  virtual void CalculateVolumes();
  virtual void InitSides();
  virtual long MoveParticle(long e, long ex, long ey, long ez, vector3& pos,
			    vector3& convection, DispersionClass& dispersion);
  virtual bool checkBoundaries(long el, long ex, long ey, long ez,
			       vector3& pos, vector3& newpos, bool first,
			       BoundaryPosition *boundstart,
			       BoundaryPosition *boundend);
  virtual int sortBoundaries(BoundaryPosition *boundpos);
  virtual void injectAtElementSide(int e, int side, int npart, ParticleClass* particle);
  int getNeighbourElements(int ex, int ey, int ez,
			   BoundaryPosition *bound,
			   int nbound, long* neighbours);
  bool isParticleInElement(vector3& pos, long ex, long ey, long ez);
};

#endif
