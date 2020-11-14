#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_mono_THz Pulse;
#include "../sfa.h"

#include "setup.h"
int main (int argc, char* argv[]) {

  set_global_params ();
  global_tag = 0; // distinguish file name
  E1 = 0.;
  setup_params ();

  std::string filename_para, filename_data;
  filename_para = "para_0.dat";
  filename_data = "data_0.dat";
  calculate <SFA_HHG> (filename_para, filename_data);

  return 0;
}
