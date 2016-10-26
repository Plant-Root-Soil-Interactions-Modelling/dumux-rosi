#ifndef EVENT_HPP
#define EVENT_HPP

#include <string>

/**
 * Class for handling simulation time events.
 */
class Event {
protected:
  double next;      ///< time for the next event.
  double increment; ///< increment for regular events.
  int counter;      ///< counter for the events.
  int autodelete;   ///< 0 or 1.
  static double epsilon;
  static double never;
  static int event_happend;
  Event() {} ///< protect the default constructor.
public:
  std::string name; ///< name of the event.
  /// constructor for irregular time events.
  Event(const char *eventname, double time, int autodel=0);
  /// constructor for regular time events.
  Event(const char *eventname, double time, double inc);
  ~Event() {}
  /// returns true if time simtime is greater as or equal to the next event time.
  inline bool mustHandle(double simtime) const { return simtime>next-epsilon; }
  /// returns the time stepping of the event.
  inline double get_increment() const { return increment; }
  /// initialize next and counter.
  void init(double time);
  /// returns the time for the next occurance of the event.
  inline double get_next_time() const { return next; }
  /// deactivates the event.
  inline double deactivate() { return next=never; }
  /// set the next time for regular events.
  void set_next_time(double simtime);
  /// set the next time for irregular events.
  void set_next_time(double simtime, double nexttime);
  /// returns the number of events happend since the last reset.
  static inline int EventHappend() { return event_happend; }
  /// resets the happened counter.
  static inline void resetEventHappend() { event_happend=0; }
  /// returns the number of event occurances.
  inline int get_counter() const { return counter; }
  /// returns true if the event is deactivated.
  inline bool isDeactivated() const { return next==never; }
  /// returns true if the event is activated.
  inline bool isActivated() const { return next!=never; }
  /// return time that never happens.
  inline int autoDelete() { return autodelete; }
  /// return time that never happens.
  static inline double get_never() { return never; }
};

#endif
