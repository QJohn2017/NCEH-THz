#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_gaussian_and_THz_practical Pulse;
#include "../sfa.h"

#include "setup.h"
int main (int argc, char* argv[]) {


  set_global_params ();
  int n_delay = 100;
  double delay_0 = 0.;
  double delay_1 = 6283.; //94274.;
  double d_delay = (delay_1 - delay_0) / (n_delay - 1.);

  for (int i_delay = 0; i_delay < n_delay; i_delay ++) {
    tc = delay_0 + i_delay * d_delay;
    t0 = tc - half_duration;
    t1 = tc + half_duration;
    global_tag = i_delay; // distinguish file name
    setup_params ();
    std::cout << "scanning at " << std::to_string (i_delay) << " of " 
              << std::to_string (n_delay) << " pts..." << std::endl;
    std::string filename_para = "res/para_" + std::to_string (i_delay) + ".dat";
    std::string filename_data = "res/data_" + std::to_string (i_delay) + ".dat";
    calculate <SFA_HHG> (filename_para, filename_data);
  }


  // calculate the HHG for tc at the peak of the THz field
/*
  set_global_params ();
  E1 = 2e-5;
  tc = 1750; //25718.;
  t0 = tc - half_duration;
  t1 = tc + half_duration;
  global_tag = 1; // distinguish file name
  setup_params ();

  std::string filename_para, filename_data;
  filename_para = "para_1.dat";
  filename_data = "data_1.dat";
  calculate <SFA_HHG> (filename_para, filename_data);
*/

  return 0;
}
