#ifndef PARTICLE_H
#define PARTICLE_H
#include "definitions.hpp"
#include <fstream>

/** Class for one particle.
 *  This class contains all properties for one particle.
 */
class Particle {
public:
  vector3 position;
  Particle* next; // pointer to next particle
  Particle() {}
  Particle(vector3& p) { position=p; }
  virtual void reduce_mass(double factor) {}
  virtual double get_mass() { return 1.0; }
  virtual void read(std::fstream &res);
  virtual void write(std::fstream &res);
  virtual const char* type() { return "simple"; }
  virtual ~Particle() {}
};


class ParticleWithMass: public Particle {
protected:
  double mass;
  ParticleWithMass() {}
public:
  ParticleWithMass(double m): mass(m) {}
  ParticleWithMass(vector3& p, double m): Particle(p), mass(m) { }
  virtual void reduce_mass(double factor) { mass*=factor; }
  virtual double get_mass() { return mass; }
  virtual void read(std::fstream &res);
  virtual void write(std::fstream &res);
  virtual const char* type() { return "withMass"; }
  virtual ~ParticleWithMass() {}
};


/** Class for a simple list of Particle(s).
 *  This class is a simple, but fast implementation for a list of Particle(s).
 *  It is used as a temporary list for Particle(s) which move from one element
 *  to another.
 */
class ParticleSimpleList {

protected:
  Particle* top;  // first particle

public:
  ParticleSimpleList(): top(0) {}
  virtual ~ParticleSimpleList() {}

  inline Particle* start() const { return top; }
  inline Particle* next(Particle* p) const { return p->next; }
  inline void insert(Particle* p) { p->next=top; top=p; }
  inline void reset() { top=0; }
};


/** Main class for keeping Particle(s) in lists.
 *  Every element has a list of the particles it contains.
 */
class ParticleList {

protected:
  Particle* top;  // first 'virtual' particle
  ulong number;   // number of particles

public:
  ParticleList();
  virtual ~ParticleList();

  inline ulong count() const { return number; }
  // iterating without removing
  inline Particle* start() const { return top->next; }
  inline Particle* next(Particle* p) const { return p->next; }

  // iterating with removing
  inline Particle* start(Particle* &prev) const { return (prev=top)->next; }
  inline Particle* next(Particle* p, Particle* &prev) const { return (prev=p)->next; }

  inline void insert(Particle* p) {
    // insert Particle p after top
    number++;
    // insert after top
    p->next=top->next;
    top->next=p;
  }

  inline Particle* insert_after(Particle* p, Particle* prev) {
    // insert Particle p after particle prev
    number++;
    // insert after prev
    p->next=prev->next;
    return prev->next=p;
  }

  inline Particle* remove(Particle* p, Particle* prev) {
    // remove particle p from list
    number--;
    prev->next=p->next;
    delete p;
    return prev->next;
  }

  inline Particle* remove_from_list(Particle* p, Particle* prev) {
    // remove particle p from list
    number--;
    return prev->next=p->next;
  }

  inline Particle* moveto(Particle* p, Particle* prev, ParticleList* plist) {
    // move particle p to list plist
    number--;
    prev->next=p->next;
    plist->insert(p);
    return prev->next;
  }

  inline void clear() {
    // clear list without deleting particles
    number=0;
    top->next=0;
  }

  void insert_list(Particle* p);  // insert simple list of Particles
  void join(ParticleList& plist); // join particle lists
};

#endif
