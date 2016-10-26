#ifndef DEFINITIONS
#define DEFINITIONS

typedef unsigned int uint;
typedef unsigned long ulong;
typedef int    triple[3];
typedef char string200[200];

class vector3 {
public:
  double Data[3];
  inline vector3() { }
  inline ~vector3() { }
  inline double* data() const { return const_cast<double*>(Data); }
  inline double& operator[](int i) { return Data[i]; }
  inline vector3& operator=(vector3& v) {
    Data[0]=v.Data[0]; Data[1]=v.Data[1]; Data[2]=v.Data[2];
    return *this;
  }
  inline vector3& operator+=(vector3 &v) {
    Data[0]+=v[0]; Data[1]+=v[1]; Data[2]+=v[2]; //MB
    return *this;
  }
  inline vector3& operator*=(double x) {
    Data[0]*=x; Data[1]*=x; Data[2]*=x;
    return *this;
  }
  inline vector3& operator/=(double x) {
    Data[0]/=x; Data[1]/=x; Data[2]/=x;
    return *this;
  }
  inline void set(double x, double y, double z) {
    Data[0]=x; Data[1]=y; Data[2]=z;
  }
  inline void set(double *x) {
    Data[0]=x[0]; Data[1]=x[1]; Data[2]=x[2];
  }
  inline void setproduct(double x, vector3 &v) {
    Data[0]=x*v[0]; Data[1]=x*v[1]; Data[2]=x*v[2];
  }
  inline void setsum(vector3 &v1, vector3 &v2) {
    Data[0]=v1.Data[0]+v2.Data[0];
    Data[1]=v1.Data[1]+v2.Data[1];
    Data[2]=v1.Data[2]+v2.Data[2];
  }
  inline void setsum(vector3 &v1, vector3 &v2, vector3 &v3) {
    Data[0]=v1.Data[0]+v2.Data[0]+v3.Data[0];
    Data[1]=v1.Data[1]+v2.Data[1]+v3.Data[1];
    Data[2]=v1.Data[2]+v2.Data[2]+v3.Data[2];
  }
  inline void setxsum(double x1, vector3 &v1, double x2, vector3 &v2) {
    Data[0]=x1*v1.Data[0] + x2*v2.Data[0];
    Data[1]=x1*v1.Data[1] + x2*v2.Data[1];
    Data[2]=x1*v1.Data[2] + x2*v2.Data[2];
  }
  inline void addxsum(double x, vector3 &v1) {
    Data[0]+=x*(v1.Data[0]);
    Data[1]+=x*(v1.Data[1]);
    Data[2]+=x*(v1.Data[2]);
  }
  inline void addxsum(double x, vector3 &v1, vector3 &v2) {
    Data[0]+=x*(v1.Data[0]+v2.Data[0]);
    Data[1]+=x*(v1.Data[1]+v2.Data[1]);
    Data[2]+=x*(v1.Data[2]+v2.Data[2]);
  }
  inline double dotprod(const vector3& v) const {
    return Data[0]*v.Data[0]+Data[1]*v.Data[1]+Data[2]*v.Data[2];
  }
  inline double dotprod(const double* v) const {
    return Data[0]*v[0]+Data[1]*v[1]+Data[2]*v[2];
  }
  inline operator const double*() const { return Data; } // convert to const double*
};

#endif
