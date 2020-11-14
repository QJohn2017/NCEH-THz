#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_mono Pulse;
#include "../sfa.h"

int main (int argc, char* argv[]) {

  // parameters of time grid for integration
  nt = 11001;
  t0 = 0;
  t1 = 2.*M_PI/w0 * 50;
  tp0 = 0;
  tp1 = t1;

  global_tag = 0; // distinguish file name

  setup_params ();

  Pulse pulse (w0, E0);
  SFA_HHG s;
  s.prepare_pulse (pulse);
  s.calc_hhg_dx ();

  std::ofstream of_para ("para_0.dat");
  of_para << Ip << " " << w0 << " " << E0 << std::endl;
  of_para << t0 << " " << t1 << " " << dt << std::endl;
  of_para.close ();

  std::ofstream of_res ("data_0.dat");
  for (int it = 0; it < nt; it ++)
    of_res << s.time[it] << " " << s.Ef[it] << " " << s.Af[it] << " "
           << s.intA1[it] << " "<< real (s.d[it]) << " " << imag (s.d[it])
           << std::endl;
  of_res.close ();

  return 0;
}
