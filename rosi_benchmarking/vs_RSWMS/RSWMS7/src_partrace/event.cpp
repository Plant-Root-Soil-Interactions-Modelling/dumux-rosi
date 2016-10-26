#include "event.hpp"
#include <limits>

double Event::epsilon=100.0 * std::numeric_limits<double>::epsilon();
double Event::never=std::numeric_limits<double>::max();
int Event::event_happend=0;


Event::Event(const char *eventname, double time, int autodel)
{
  name=eventname;
  next=time;
  increment=0.0;
  counter=0;
  autodelete=autodel;
}


Event::Event(const char *eventname, double time, double inc)
{
  name=eventname;
  increment=inc;
  init(time);
  autodelete=0;
}


void Event::init(double time)
{
  if(increment>0.0) {
    if(time==0.0) {
      next=0;
      counter=0;
    }
    else {
      counter=int(time/increment+epsilon);
      next=(++counter)*increment;
    }
    ++event_happend;
  }
}


void Event::set_next_time(double simtime)
{
  if(mustHandle(simtime)) {
    while(next<=simtime+epsilon)
      next=(++counter)*increment;
    ++event_happend;
  }
}


void Event::set_next_time(double simtime, double nexttime)
{
  if(mustHandle(simtime)) {
    next=nexttime;
    ++counter;
    ++event_happend;
  }
}
