#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_gaussian Pulse;
#include "../sfa.h"

class SFA_HHG_GAUSS : public SFA_HHG {
public:

  double action_S (double ps, int itau, int it, int itp) {

    double tau = itau * dt;

    double phi_tau = w0 * itau * dt;
    double phi_t = w0 * (t0 + it * dt - tc);
    double phi_sig = w0 * sigma;

    const double sqrt_pi = 2. / M_2_SQRTPI; // sqrt (M_PI);
    double x = phi_t / phi_sig * M_SQRT1_2;
    double y = phi_sig * M_SQRT1_2;
    double eta = phi_tau / phi_sig * M_SQRT1_2;

    auto G = [sqrt_pi] (double x, double y, double e) {
      auto g = [] (complex z, double e) {
        return w_of_z (I * M_SQRT2 * z) - exp (4.*e*(z-e/2.)) * w_of_z (I * M_SQRT2 * (z - e));
      };
      complex z = complex (x, -y);
      double im = imag (exp (2.*I*x*y) * g (z * M_SQRT1_2, e * M_SQRT1_2));
      return real (exp (4.*I*x*y) * g (z, e) - g (x, e)) - M_SQRT2 * sqrt_pi / e * im * im;
    };

    return Ip * tau + Up * exp (-2*x*x) / M_2_SQRTPI * sigma * G (x, y, eta);
  }

/*
  double action_S (double ps, int itau, int it, int itp) {

    double tau = itau * dt;

    double phi_tau = w0 * itau * dt;
    double phi_t = w0 * (it * dt - tc);
    double phi_sig = w0 * sigma;

    const double sqrt_pi = 2. / M_2_SQRTPI; // sqrt (M_PI);
    double x = phi_t / phi_sig * M_SQRT1_2;
    double y = phi_sig * M_SQRT1_2;
    double eta = phi_tau / phi_sig * M_SQRT1_2;

    auto s = sin (0.5 * phi_tau);
    auto C0 = sin (phi_tau) - 4. * s * s / phi_tau;
    return Ip * tau + Up / w0 * exp (-2*x*x) * (phi_tau - 2./phi_tau*(1-cos (phi_tau)) - C0 * cos (2*phi_t-phi_tau));
  }
*/
};

int main (int argc, char* argv[]) {
  std::ofstream of_para ("para.dat");
  of_para << Ip << " " << w0 << " " << E0 << std::endl;
  of_para << t0 << " " << t1 << " " << dt << std::endl;
  of_para.close ();

  Pulse pulse (w0, E0, tc, fwhm);
  std::ofstream of_res;
 

  // // full calculation
  SFA_HHG s;
  s.prepare_pulse (pulse);
  s.calc_hhg_dx ();
  of_res.open ("data.dat");
  for (int it = 0; it < nt; it ++)
    of_res << s.time[it] << " " << s.Ef[it] << " " << s.Af[it] << " "
           << s.intA1[it] << " "<< real (s.d[it]) << " " << imag (s.d[it])
           << std::endl;
  of_res.close ();

/*
  // analytical results   
  SFA_HHG_GAUSS s1;
  s1.prepare_pulse (pulse);
  s1.calc_hhg_dx ();
  of_res.open ("data_gauss.dat");
  for (int it = 0; it < nt; it ++)
    of_res << s1.time[it] << " " << s1.Ef[it] << " " << s1.Af[it] << " "
           << s1.intA1[it] << " "<< real (s1.d[it]) << " " << imag (s1.d[it])
           << std::endl;
  of_res.close ();
*/
  return 0;
}
