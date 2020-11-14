#pragma once
#include <complex>
#include <sstream>
#include "parameters.h"
using namespace parameters;

class SFA_HHG {
 public:
  double *time, *Ef, *Af, *intA1, *intA2;
  complex *d;

  SFA_HHG () {
    std::cout << "nt " << nt  << " @ SFA_HHG allocation" << std::endl;
    time = new double [nt];
    Ef = new double [nt];
    Af = new double [nt];
    intA1 = new double [nt];
    intA2 = new double [nt];
    d = new complex [nt];
  }
  ~SFA_HHG () {
    delete[] time;
    delete[] Ef;
    delete[] Af;
    delete[] intA1;
    delete[] intA2;
    delete[] d;
  }

  void prepare_pulse (Pulse& pulse) {
    double t_cur;
    double Ef_cur, Af_cur, intA1_cur, intA2_cur;  
    double Ef_old, Af_old, intA1_old, intA2_old;

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
    exit (-1);
*/

    long ntp = (tp1_adjusted - tp0_adjusted) / dt + 1;
    std::ostringstream filename;
    filename << "field_" << global_tag << ".dat";
    std::ofstream of_field (filename.str ());

    std::cout << "E1: " << E1 << std::endl;

    for (long itp = 0; itp < ntp; itp ++) {
      /// start
      t_cur = tp0_adjusted + itp * dt;
      Ef_cur = pulse.E (t_cur);
      Af_cur = pulse.A (t_cur);

      if (itp == 0) {
        intA1_cur = 0.;
        intA2_cur = 0.;
      } else {
        intA1_cur = intA1_old + Af_old * dt;
        intA2_cur = intA2_old + Af_old * Af_old * dt;
      }

      of_field << t_cur << " " << Ef_cur << " " << Af_cur << std::endl;

      long it = itp - it_valid_0;
      
      if (it >= 0 && it < nt) {
        time[it] = t_cur;
        Ef[it] = Ef_cur;
        Af[it] = Af_cur;
        intA1[it] = intA1_cur;
        intA2[it] = intA2_cur; 
      }

      Ef_old = Ef_cur;
      Af_old = Af_cur;
      intA1_old = intA1_cur;
      intA2_old = intA2_cur;
      /// end
    }

    of_field.close ();
  }

  virtual double action_S (double ps, int itau, int it, int itp) {
    double tau = itau * dt;
    return (0.5 * ps * ps + Ip) * tau
      + ps * (intA1[it] - intA1[itp])
      + 0.5 * (intA2[it] - intA2[itp]);
  }

  void calc_hhg_dx () {
    const double epsilon = 0.1;
    for (int it = 0; it < nt; it ++) {
      d[it] = 0.;
      // itau loop starts with 1, otherwise tau=0 renders ps singular
      for (int itau = 1; itau <= it; itau ++) {
        double tau = itau * dt;
        int itp = it - itau; // ionization time t'
        double ps = - (intA1[it] - intA1[itp]) / tau;
        double S = action_S (ps, itau, it, itp);
        d[it] += conj (dipole (ps + Af[it])) * dipole (ps + Af[itp]) * Ef[itp]
          * exp (- I * S) * pow (M_PI / (epsilon + I * tau / 2.), 1.5);
      }
      d[it] *= I * (-dt);

      if (it % 1000 == 0) {
        std::cout << "progress (%): " << 1.0 * it / nt << std::endl;
      }

    }
  }

};
