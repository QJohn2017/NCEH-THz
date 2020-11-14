#include <iostream>
#include <fstream>
#include "../pulse.h"
typedef Pulse_gaussian_and_THz Pulse;
#include "../sfa.h"
#include "../erf.h"

class SFA_HHG_GAUSS : public SFA_HHG {
 public:
  double* Ef0;
  double* Af0;
  double* intA10;
  double* intA20;
  double* Ef1;
  double* Af1;
  double* intA11;
  double* intA21;
  double* int_A0_times_A1;

  double eps;

  SFA_HHG_GAUSS () {
    Ef0 = new double [nt];
    Af0 = new double [nt];
    intA10 = new double [nt];
    intA20 = new double [nt];
    Ef1 = new double [nt];
    Af1 = new double [nt];
    intA11 = new double [nt];
    intA21 = new double [nt];
    int_A0_times_A1 = new double [nt];
    eps = w1 / w0;
  }
  ~SFA_HHG_GAUSS () {
    delete[] Ef0;
    delete[] Af0;
    delete[] intA10;
    delete[] intA20;
    delete[] Ef1;
    delete[] Af1;
    delete[] intA11;
    delete[] intA21;
    delete[] int_A0_times_A1;
  }

  void prepare_pulse (Pulse& pulse) {
    double t_cur;
    double Ef0_cur, Af0_cur, intA10_cur, intA20_cur;  
    double Ef0_old, Af0_old, intA10_old, intA20_old;
    double Ef1_cur, Af1_cur, intA11_cur, intA21_cur, int_A0_times_A1_cur;  
    double Ef1_old, Af1_old, intA11_old, intA21_old, int_A0_times_A1_old;

    if (t0 < tp0) {
      std::cerr << "invalid t0 < tp0 when prepare_pulse!" << std::endl;
      std::cerr << "recommended value of tp0 < " << t0 << std::endl; // based on the Gaussian pulse
    }
    if (t1 > tp1) {
      std::cerr << "invalid t1 > tp1 when prepare_pulse!" << std::endl;
      std::cerr << "recommended value of tp1 > " << t1 << std::endl; // based on the Gaussian pulse
    }
    int it_valid_0 = int (ceil ((t0 - tp0) / dt));


    double tp0_adjusted = t0 - it_valid_0 * dt;
    double tp1_adjusted = t1 + int (ceil ((tp1 - t1) / dt)) * dt;

/*
    std::cout << "it_valid_0: " << it_valid_0 << std::endl;
    std::cout << "t0: " << t0 << " tp0_adjusted: " << tp0_adjusted << std::endl;
    std::cout << "t1: " << t1 << " tp1_adjusted: " << tp1_adjusted << std::endl;
    std::cout << "ntp: " << (tp1_adjusted - tp0_adjusted) / dt << std::endl;
    std::cout << "nt: " << nt << std::endl; 
*/

    int ntp = (tp1_adjusted - tp0_adjusted) / dt + 1;

    // std::ostringstream filename;
    // filename << "field_" << global_tag << ".dat";
    // std::ofstream of_field (filename.str ());

    Pulse_gaussian pulse0 (w0, E0, tc, fwhm);
    Pulse_THz pulse1 (w1, E1, phi1);

    for (int itp = 0; itp < ntp; itp ++) {
      /// start
      t_cur = tp0_adjusted + itp * dt;
      Ef0_cur = pulse0.E (t_cur);
      Af0_cur = pulse0.A (t_cur);
      Ef1_cur = pulse1.E (t_cur);
      Af1_cur = pulse1.A (t_cur);

      if (itp == 0) {
        intA10_cur = 0.;
        intA20_cur = 0.;
        intA11_cur = 0.;
        intA21_cur = 0.;
        int_A0_times_A1_cur = 0.;
      } else {
        intA10_cur = intA10_old + Af0_old * dt;
        intA20_cur = intA20_old + Af0_old * Af0_old * dt;
        intA11_cur = intA11_old + Af1_old * dt;
        intA21_cur = intA21_old + Af1_old * Af1_old * dt;
        int_A0_times_A1_cur = int_A0_times_A1_old + Af0_old * Af1_old * dt;
      }

      // of_field << t_cur << " " << Ef_cur << " " << Af_cur << std::endl;

      int it = itp - it_valid_0;
      if (it >= 0 && it < nt) {
        time[it] = t_cur;
        Ef0[it] = Ef0_cur;
        Af0[it] = Af0_cur;
        intA10[it] = intA10_cur;
        intA20[it] = intA20_cur;
        Ef1[it] = Ef1_cur;
        Af1[it] = Af1_cur;
        intA11[it] = intA11_cur;
        intA21[it] = intA21_cur;
        int_A0_times_A1[it] = int_A0_times_A1_cur;
      }

      Ef0_old = Ef0_cur;
      Af0_old = Af0_cur;
      intA10_old = intA10_cur;
      intA20_old = intA20_cur;
      Ef1_old = Ef1_cur;
      Af1_old = Af1_cur;
      intA11_old = intA11_cur;
      intA21_old = intA21_cur;
      int_A0_times_A1_old = int_A0_times_A1_cur;
      /// end
    }

    // of_field.close ();
  }

  double action_S0 (double ps0, int itau, int it, int itp) {
    /*
    double tau = itau * dt;
    return (0.5 * ps0 * ps0 + Ip) * tau
      + ps0 * (intA10[it] - intA10[itp])
      + 0.5 * (intA20[it] - intA20[itp]);
    */
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

  double action_dS (double ps0, double ps1, int itau, int it, int itp) {
    double tau = itau * dt;
    // return - ps0 * ps1 * tau + int_A0_times_A1[it] - int_A0_times_A1[itp];


    double phi_t = w0 * (t0 + it * dt - tc); 
    double phi_tau = w0 * itau * dt;
    double phi_sig = w0 * sigma;

    double x = M_SQRT1_2 * phi_t / phi_sig;
    double y = M_SQRT1_2 * phi_sig;
    double eta = M_SQRT1_2 * phi_tau / phi_sig;
/*
    auto z = complex (y, x);
    auto c1 = exp (2.*x*eta - 2.*I*y*eta - eta*eta);
    double dS = - eps * A0 * A1 / w0 * sin (phi1) * sqrt (M_PI) * y * y * exp (-x*x) 
      * real (exp (2.*I*x*y) * ((2.*z-I*eta) * (w_of_z (z) - c1 * w_of_z (z-I*eta)) 
        - I*2./sqrt (M_PI) * (1. - c1)));
    return dS;
*/
    // above approximations, and assume e^{phi_t*phi_tau/(phi_sig^2)} = 1
    // now we have the same Delta S as for a fs laser of plane wave, except for an extra exp (-x*x)
    // Folder: "using_faddeeva_eps_to_1st_order_large_z_limit_2_polynomial_highest_order_exp1"
    double C1 = phi_tau * cos (phi_tau/2) - 2. * sin (phi_tau/2);
    double dS = - A0 * E1 / (w0*w0) * exp (-x*x) * sin (phi1) * C1 * cos (phi_t-phi_tau/2.);

/*
    std::cout << "A0: " << A0 << std::endl;
    std::cout << "E1: " << E1 << std::endl;
    std::cout << "w0: " << w0 << std::endl;
    exit (-1);
*/
    return dS;
  }

  double action_S (double ps0, double ps1, int itau, int it, int itp) {
    return action_S0 (ps0, itau, it, itp) + action_dS (ps0, ps1, itau, it, itp);
  }

  void calc_hhg_dx () {
    const double epsilon = 0.1;
    for (int it = 0; it < nt; it ++) {
      d[it] = 0.;
      // itau loop starts with 1, otherwise tau=0 renders ps singular
      for (int itau = 1; itau <= it; itau ++) {
        double tau = itau * dt;
        int itp = it - itau; // ionization time t'
        // double ps = - (intA1[it] - intA1[itp]) / tau;
        double ps0 = - (intA10[it] - intA10[itp]) / tau;
        double ps1 = - (intA11[it] - intA11[itp]) / tau;
        // double S = action_S (ps, itau, it, itp);
        double S = action_S (ps0, ps1, itau, it, itp);
        d[it] += conj (dipole (ps0 + Af0[it])) * dipole (ps0 + Af0[itp]) * Ef0[itp]
          * exp (- I * S) * pow (M_PI / (epsilon + I * tau / 2.), 1.5);
      }
      d[it] *= I * (-dt);
      if (it % 100 == 0)
        std::cout << "progress: " << 1.0 * it / nt << std::endl;
    }
  }
};

#include "setup.h"

int main (int argc, char* argv[]) {

  set_global_params ();
  E1 = 2e-5;
  global_tag = 4; // distinguish file name
  setup_params ();

  std::string filename_para, filename_data;
  filename_para = "para_4.dat";
  filename_data = "data_4.dat";
  calculate <SFA_HHG_GAUSS> (filename_para, filename_data);

  return 0;
}
