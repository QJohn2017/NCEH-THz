#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_gaussian Pulse;
#include "../sfa.h"

int main (int argc, char* argv[]) {

  // parameters of fs pulse
  fwhm  = 100 * 41.34;
  tc    = 6000.;
  half_duration = fwhm * 1.478;
  sigma = fwhm / (2. * sqrt (2. * log (2.)));
  
  // parameters of time grid for integration
  nt  = int (half_duration * 2.27);
  t0  = tc - half_duration;
  t1  = tc + half_duration;
  tp0 = -1400;
  tp1 = 14000;

  global_tag = 1; // distinguish file name

  setup_params ();

  std::ofstream of_para ("para_1.dat");
  of_para << Ip << " " << w0 << " " << E0 << std::endl;
  of_para << t0 << " " << t1 << " " << dt << std::endl;
  of_para.close ();

  Pulse pulse (w0, E0, tc, fwhm);
  std::ofstream of_res;
 
  SFA_HHG s;
  s.prepare_pulse (pulse);
  s.calc_hhg_dx ();
  of_res.open ("data_1.dat");
  for (int it = 0; it < nt; it ++)
    of_res << s.time[it] << " " << s.Ef[it] << " " << s.Af[it] << " "
           << s.intA1[it] << " "<< real (s.d[it]) << " " << imag (s.d[it])
           << std::endl;
  of_res.close ();

  return 0;
}
