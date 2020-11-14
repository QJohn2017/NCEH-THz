#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <memory>
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

    Pulse_gaussian pulse0 (w0, E0, tc, fwhm);
    Pulse_THz pulse1 (w1, E1, phi1);

    for (int it = 0; it < nt; it ++) {
      time[it] = t0 + it * dt;
      Ef0[it] = pulse0.E (time[it]);
      Af0[it] = pulse0.A (time[it]);
      Ef1[it] = pulse1.E (time[it]);
      Af1[it] = pulse1.A (time[it]);
      if (it == 0) {
        // Af[it] = 0.;          
        intA10[it] = 0.;
        intA20[it] = 0.;
        intA11[it] = 0.;
        intA21[it] = 0.;
      } else {
        // Af[it] = Af[it-1] - Ef[it-1] * dt;
        auto A0 = Af0[it-1], A1 = Af1[it-1];
        intA10[it] = intA10[it-1] + A0 * dt;
        intA20[it] = intA20[it-1] + A0 * A0 * dt;
        intA11[it] = intA11[it-1] + A1 * dt;
        intA21[it] = intA21[it-1] + A1 * A1 * dt;
        int_A0_times_A1[it] = int_A0_times_A1[it-1] + A0 * A1 * dt;
      }
    }
  }

  double action_S0 (double ps0, int itau, int it, int itp) {
    double tau = itau * dt;
    return (0.5 * ps0 * ps0 + Ip) * tau
      + ps0 * (intA10[it] - intA10[itp])
      + 0.5 * (intA20[it] - intA20[itp]);
  }

  double action_dS (double ps0, double ps1, int itau, int it, int itp) {
    double tau = itau * dt;
/*     return - ps0 * ps1 * tau + int_A0_times_A1[it] - int_A0_times_A1[itp]; */

    double phi_t = w0 * (t0 + it * dt - tc);
    double phi_tau = w0 * itau * dt;
    double phi_sig = w0 * sigma;

/*     // final formula
    double dS = - (A0 * A1 / w0) * exp (-phi_t*phi_t/(2.*phi_sig*phi_sig)) * (cos (phi1) - eps * phi_t * sin (phi1))
      * (cos (phi_t) - sin (phi_t) - exp (phi_t*phi_tau/(phi_sig*phi_sig)) * (cos (phi_t-phi_tau) - sin (phi_t-phi_tau)));
    return dS;
 */   


/*     // formula for approximation of f(x, y, eta)
    auto f = [] (double x, double y, double eta) {
      auto z = complex (x, y);
      return (exp (2.*eta*z) - 1.) / (2.*z);
    };

    double x = M_SQRT1_2 * phi_t / phi_sig;
    double y = M_SQRT1_2 * phi_sig;
    double eta = M_SQRT1_2 * phi_tau / phi_sig;

    auto e1 = exp (I * (phi1 + eps * phi_t));

    double dS = A0 * A1 * M_SQRT1_2 / w0 * phi_sig * exp (-phi_t*phi_t/(2.*phi_sig*phi_sig)) 
      * ( - imag (exp (-I*phi_t) * (e1 * f (x, (1-eps)*y, eta) + conj (e1) * f (x, (1+eps)*y, eta)))
          + 4./(eps*phi_tau) * real (exp (-I*phi_t) * f (x, y, eta)) * cos (phi1 + eps * (phi_t-phi_tau/2)) * sin (eps*phi_tau/2) );
    return dS;
 */


    double x = M_SQRT1_2 * phi_t / phi_sig;
    double y = M_SQRT1_2 * phi_sig;
    double eta = M_SQRT1_2 * phi_tau / phi_sig;

    // using Faddeeva function, perfect approximation up to the first order of epsilon

/*     auto z = complex (y, x);
    auto c1 = exp (2.*x*eta - 2.*I*y*eta - eta*eta);
    double dS = - eps * A0 * A1 / w0 * sin (phi1) * sqrt (M_PI) * y * y * exp (-x*x) 
      * real (exp (2.*I*x*y) * ((2.*z-I*eta) * (w_of_z (z) - c1 * w_of_z (z-I*eta)) 
        - I*2./sqrt (M_PI) * (1. - c1)));
    return dS;
 */ 


    // the formula used in the paper, perfect approximation up to the first order of epsilon
/*     auto z = complex (x, -y);
    auto c1 = exp (2. * eta * (z - eta/2.));
    double dS = sqrt (M_PI) * A0 * E1 * sin (phi1) * exp (-x*x) * sigma * sigma
      * imag (exp (I*phi_t) * ((z-eta/2.) * (w_of_z (I*z) - c1 * w_of_z (I*(z-eta))) - 1. / sqrt (M_PI) * (1. - c1)));
    return dS;
 */


/* 
    // approximate the Faddeeva function with i*z / sqrt(pi) / (z^2-0.5), using the *2nd* formula of Table 1
    // in M. R. Zaghloul, ACM Trans. Math. Softw. 44, 22:1 (2017).
    // Folder: "using_faddeeva_eps_to_1st_order_large_z_limit_2"
    auto z = complex (x, -y);
    auto c1 = exp (2. * eta * (z - eta/2.));
    double dS = A0 * E1 * sin (phi1) * exp (-x*x) * sigma * sigma
      * imag (exp (I*phi_t) * ((z-eta/2.) * (z/(z*z+0.5) - c1*(z-eta)/((z-eta)*(z-eta)+0.5)) - (1. - c1)));
    return dS;
 */

/* 
    // approximation neglecting all low order terms in polynomials for prefactors of trigeometric functions
    // Folder: "using_faddeeva_eps_to_1st_order_large_z_limit_2_polynomial_highest_order"     
    double dS = A0 * E1 / (2.*w0*w0) * exp (-x*x) * sin (phi1) * (-phi_tau * cos (phi_t) + 2 * sin (phi_t) - 
      exp (phi_t*phi_tau/(phi_sig*phi_sig)) * (phi_tau * cos (phi_t-phi_tau) + 2 * sin (phi_t-phi_tau)));
    return dS;
 */

    // above approximations, and assume e^{phi_t*phi_tau/(phi_sig^2)} = 1
    // now we have the same Delta S as for a fs laser of plane wave, except for an extra exp (-x*x)
    // Folder: "using_faddeeva_eps_to_1st_order_large_z_limit_2_polynomial_highest_order_exp1"
    double C1 = phi_tau * cos (phi_tau/2) - 2 * sin (phi_tau/2);
    double dS = - A0 * E1 / (w0*w0) * exp (-x*x) * sin (phi1) * C1 * cos (phi_t-phi_tau/2);
    return dS;


/*     auto z = complex (x, -y);
    auto c1 = exp (2. * eta * (z - eta/2.));
    double dS = - A0 * E1 * sin (phi1) * exp (-x*x) * sigma * sigma
      * imag (exp (I*phi_t) * (1. - c1));
    return dS;
 */

/*     double dS = -2.* A0 * A1 * exp (-x*x) * sin (phi1) * eps * (y * cos (2*x*y) + x * sin (2*x*y)) * eta / (w0 * (1 + x*x/(y*y)));
    return dS;
 */
    // double dS = - A0 * A1 * exp (-x*x) / w0 * eps * phi_tau * sin (phi1) * (cos (phi_t) + phi_t * sin (phi_t) / (phi_sig*phi_sig));
    // double dS = - A0 * A1 * exp (-x*x) / w0 * eps * phi_tau * sin (phi1) * cos (phi_t);
    // return dS;

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

template <typename SFA_CALC_MODEL>
void calculate (const std::string& filename_para, const std::string& filename_data) {

  std::ofstream of;
  of.open (filename_para);
  of << Ip << " " << w0 << " " << E0 << std::endl;
  of << t0 << " " << t1 << " " << dt << " " << tc << std::endl;
  of.close ();

  Pulse pulse (w0, E0, tc, fwhm, w1, E1, phi1);
  std::ofstream of_res (filename_data);

  // full calculation
  SFA_CALC_MODEL s;
  // SFA_HHG s;
  s.prepare_pulse (pulse);
  s.calc_hhg_dx ();

  of.open (filename_data);
  for (int it = 0; it < nt; it ++)
    of << s.time[it] << " " << s.Ef[it] << " " << s.Af[it] << " "
       << s.intA1[it] << " "<< real (s.d[it]) << " " << imag (s.d[it])
       << std::endl;
  of.close (); 
}

void scan_over_tc () {
  // scan over tc, reconstructing THz wave form
  double tc_0 = t0 + 12000;
  double tc_1 = t1 - 12000;
  int n_tc = 100;
  double d_tc = (tc_1 - tc_0) / (n_tc - 1);

  std::string filename_para, filename_data;
  std::string dir ("./test-data/scan_test/");

  for (int i_tc = 0; i_tc < n_tc; i_tc ++) {
    tc = tc_0 + i_tc * d_tc;
    std::cout << "current tc: " << tc << std::endl;

    filename_para = dir + "para_" + std::to_string (i_tc) + ".dat";
    filename_data = dir + "data_" + std::to_string (i_tc) + ".dat";
    calculate <SFA_HHG> (filename_para, filename_data);
  }
}

void compare_numerical_and_analytic_solutions () {


  std::string filename_para, filename_data;


  // numerical result
  filename_para = "para.dat";
  filename_data = "data.dat";
  calculate <SFA_HHG> (filename_para, filename_data);


/* 
  // analytically derived result
  filename_para = "para_gauss.dat";
  filename_data = "data_gauss.dat";
  calculate <SFA_HHG_GAUSS> (filename_para, filename_data);
 */
}

int main (int argc, char* argv[]) {
  // scan_over_tc ();
  compare_numerical_and_analytic_solutions ();
  return 0;
}
