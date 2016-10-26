// Partrace main program.
// AUTHOR: Horst Hardelauf (h.hardelauf@fz-juelich.de)
// VERSION: August 2006
//       Agrosphere Institute
//       ICG IV
//       Forschungszentrum Juelich GmbH
//       52425 Juelich
//       (Building 16.6)
//       Tel  0049 2461 61 6392
//       Fax  0049 2461 61 2518
//


#include "partrace.hpp"
#include "parallel.hpp"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{ 
  const int PartraceCaller=0;
  Parallel parallel(argc, argv);
  PartraceClass partrace(PartraceCaller, &parallel);
  set_new_handler(PartraceNewHandler);
  int me=parallel.mycpu();

  bool debug;
  if(argc<3) debug=false;
  else debug=true;
  partrace.start_debug(debug);

  if(me==0) {
    cout<<"Calculation of Solute Transport by Particle Tracking."<<endl;
  }
  cout.flush();
  parallel.sync();

  partrace.input(argv[1]);

  // run time loop
  partrace.run();

  partrace.elements->CountParticlesOutside();
  partrace.particles.WriteRestartFile(true);
  if(me==0) cout<<"Partrace finished"<<endl;  

  return 0;
}
