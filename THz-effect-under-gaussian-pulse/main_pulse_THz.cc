#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_THz Pulse;
#include "../sfa.h"

#include "setup.h"
int main (int argc, char* argv[]) {

  set_global_params ();
  global_tag = 2; // distinguish file name
  setup_params ();

  // Pulse pulse (w0, E0, tc, fwhm, w1, E1, phi1);
  Pulse pulse (w1, E1, phi1);
  std::ofstream of_res;
 
  SFA_HHG s;
  s.prepare_pulse (pulse);
  
  return 0;
}