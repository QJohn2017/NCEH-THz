#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_gaussian Pulse;
#include "../sfa.h"

#include "setup.h"
int main (int argc, char* argv[]) {

  set_global_params ();
  global_tag = 3; // distinguish file name
  setup_params ();

  //Pulse pulse (w0, E0, tc, fwhm, w1, E1, phi1);
  Pulse pulse (w0, E0, tc, fwhm);
  std::ofstream of_res;
 
  SFA_HHG s;
  s.prepare_pulse (pulse);
  
  return 0;
}
