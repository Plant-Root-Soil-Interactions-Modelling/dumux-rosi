#include "particle.hpp"

void Particle::read(std::fstream &res)
{
  res.read((char *) position.data(), 3*sizeof(double));
}


void ParticleWithMass::read(std::fstream &res)
{
  res.read((char *) position.data(), 3*sizeof(double));
  res.read((char *) &mass, sizeof(double));
}


void Particle::write(std::fstream &res)
{
  res.write((char *) position.data(), 3*sizeof(double));
}


void ParticleWithMass::write(std::fstream &res)
{
  res.write((char *) position.data(), 3*sizeof(double));
  res.write((char *) &mass, sizeof(double));
}


ParticleList::ParticleList(): number(0)
{
  top=new Particle;
  top->next=0;
}


ParticleList::~ParticleList()
{
  Particle* p;
  Particle* prev=top;
  while(prev) {
    p=prev->next;
    delete prev;
    prev=p;
  }
}


void ParticleList::insert_list(Particle* p) {
  // insert simple list of Particles starting with p after top
  if(p==0) return;
  Particle* head=top->next;
  Particle* prev;
  top->next=p;
  while(p) {
    number++;
    prev=p;
    p=p->next;
  }
  prev->next=head;
}


void ParticleList::join(ParticleList& plist)
{
  // join 2 particles lists, list plist is cleared 
  Particle *p, *prev;
  p=start(prev);
  while(p) p=next(p,prev);
  prev->next=plist.start();
  plist.clear();
}
