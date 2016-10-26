#include "runtime.hpp"
#include "event.hpp"
#include <limits>
#include <algorithm>
#include <iostream>
#include "partrace.hpp"
#include "parallel.hpp"


RunTimeClass::RunTimeClass():
  dt(1.0), Time(0)
{
}


RunTimeClass::~RunTimeClass()
{
}


double RunTimeClass::Seconds() const
{
  return partrace->parallel->getTime()-SystemTime;
}


void RunTimeClass::Init(PartraceClass *pt, double simtime) 
{ 
  partrace=pt;
  time_eps=100.0*std::numeric_limits<double>::epsilon();
  SimulationTime=nextEventTime=dt=dt_max=simtime;
  Time=oldTime=0.0;
  TimeStep=0;
  dt_external=-1.0;
  ProgressStep=SimulationTime/20.0;
  ProgressTime=ProgressStep;
  Progress=0;
  SystemTime=partrace->parallel->getTime();    
}


void RunTimeClass::resetEvents()
{
  // initialize the events depending on the simulation time.
  for(eventlist::iterator it=EventList.begin();it!=EventList.end();++it)
    (*it)->init(Time);
}


void RunTimeClass::Increment()
{
  double next;
  nextEventTime=Event::get_never();
  eventlist::iterator it=EventList.begin();
  if(dt_external>0.0) {
    dt_max=dt_external;
  }
  else if(Time==0.0) {
    // look if the first event happens before the first initial time step.
    while(it!=EventList.end()) {
      next=(*it)->get_increment();
      if(next>Time && next<nextEventTime) nextEventTime=next;
      next=(*it)->get_next_time();
      if(next>Time && next<nextEventTime) nextEventTime=next;
      // partrace->debug<<(*it)->name<<" nextEventTime0="<<nextEventTime<<std::endl;
      ++it;
    }
  }
  else {
    // An event has been handled and a new event time was set.
    // Now the time for the next event must be found.
    Event::resetEventHappend();
    eventlist::iterator it2;
    while(it!=EventList.end()) {
      //  std::cout<<"oldtime="<<oldTime<< "time="<<Time<<" event="<<(*it)->name<<std::endl;
      if( (*it)->isDeactivated() ) {
		std::cout<<"deactivated"<<std::endl;
	// remove deactivated events from the list.
	it2=it;
	++it;
	EventList.erase(it2);
      }
      else {
	//	std::cout<<"activated"<<std::endl;
	next=(*it)->get_next_time();
	std::cout<<"event="<<(*it)->name<<' '<<next<<std::endl;
	// get the minimum time for the next event.
	if(next>Time+time_eps && next<nextEventTime) nextEventTime=next;
	// partrace->debug<<" nextEventTime="<<next<<std::endl;
	// remove and delete autodelete events from the list after handling.
	if( (*it)->autoDelete() && next<=oldTime) {
	  it2=it;
	  ++it;
	  EventList.erase(it2);
	  delete *it2;
	}
	else ++it;
      }
    }
  }
  // Now calculate the next simulation time.
  // partrace->debug<<" lastEventTime="<<nextEventTime<<std::endl;
  oldTime=Time;
  Time+=dt_max;
  if(nextEventTime>oldTime && nextEventTime<Time+0.02*dt_max) Time=nextEventTime;
  dt=Time-oldTime; 
}


bool RunTimeClass::progress(int &pro)
{
  if(Time>=ProgressTime-time_eps) {
    Progress+=5;
    ProgressTime+=ProgressStep;
    pro=Progress;
    return true;
  }
  else return false;
}


void RunTimeClass::insert(Event *e)
{
  EventList.push_back(e);
}
