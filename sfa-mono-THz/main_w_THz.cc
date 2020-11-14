#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_mono_THz Pulse;
#include "../sfa.h"

#include "setup.h"
int main (int argc, char* argv[]) {

  set_global_params ();
  global_tag = 1; // distinguish file name
  E1 = 2e-5;
  setup_params ();

  std::string filename_para, filename_data;
  filename_para = "para_1.dat";
  filename_data = "data_1.dat";
  calculate <SFA_HHG> (filename_para, filename_data);

  return 0;
}
