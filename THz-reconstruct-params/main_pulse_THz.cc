#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_gaussian_and_THz_practical Pulse;
#include "../sfa.h"

#include "setup.h"
int main (int argc, char* argv[]) {

  set_global_params ();
  global_tag = 101; // distinguish file name
  E0 = 0.;
  E1 = 2e-5;
  setup_params ();

  Pulse pulse1 (w0, E0, tc, fwhm, w1, E1, phi1);
  // Pulse pulse (w1, E1, phi1);
 
  SFA_HHG s;
  s.prepare_pulse (pulse1);



  set_global_params ();
  global_tag = 102; // distinguish file name
  E0 = 0.06;
  E1 = 0.;
  setup_params ();

  Pulse pulse2 (w0, E0, tc, fwhm, w1, E1, phi1);
  // Pulse pulse (w1, E1, phi1);
 
  s.prepare_pulse (pulse2);  
  
  return 0;
}