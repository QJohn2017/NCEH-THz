#include "pulse.h"
#include <iostream>
#include <fstream>

// typedef Pulse_sin2 Pulse;
// int main (int argc, char* argv[]) {
//   const double w = 0.057;
//   const double E0 = 0.05;
//   const double nc = 10.;
//   const double tc = 500.;
//   Pulse pulse (w, E0, tc, nc);

//   const double t0 = - 2. * M_PI / w * (nc + 4.) / 2. + tc;
//   const double t1 =   2. * M_PI / w * (nc + 4.) / 2. + tc;
//   const int nt = 200;
//   const double dt = (t1 - t0) / (nt - 1);
//   std::ofstream outfile ("res_sin2.dat");
//   outfile << std::scientific;
//   for (int it = 0; it < nt; it ++) {
//     double t = t0 + it * dt;
//     outfile << t << " " << pulse.A (t) << " " << pulse.E (t) << std::endl;
//   }
//   return 0;
// }

typedef Pulse_gaussian Pulse;
int main (int argc, char* argv[]) {
  const double w = 0.057;
  const double E0 = 0.05;
  const double fwhm = 440.;
  const double tc = 500.;
  Pulse pulse (w, E0, tc, fwhm);

  const double t0 = -1000. + tc;
  const double t1 = 1000. + tc;
  const int nt = 200;
  const double dt = (t1 - t0) / (nt - 1);
  std::ofstream outfile ("res_gaus.dat");
  outfile << std::scientific;
  for (int it = 0; it < nt; it ++) {
    double t = t0 + it * dt;
    outfile << t << " " << pulse.A (t) << " " << pulse.E (t) << std::endl;
  }
  return 0;
}
