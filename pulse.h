#pragma once

#include <cmath>
#include <complex>
#include "erf.h"

#define complex std::complex<double>

class Pulse_mono {
  double w;   // frequency
  double E0;  // field amplitude
  double A0;
 public:
  Pulse_mono (double w, double E0) :
  w (w), E0 (E0) { A0 = E0 / w; }
  inline double A (double t) {
    return A0 * cos (w * t);
  }
  inline double E (double t) {
    return E0 * sin (w * t);
  }
};

class Pulse_mono_THz {
  double w0;  // frequency
  double E0;  // field amplitude
  double A0;
  double w1;
  double E1;
  double A1;
  double phi1;
 public:
  Pulse_mono_THz (double w0, double E0, double w1, double E1, double phi1) :
  w0 (w0), E0 (E0), w1 (w1), E1 (E1), phi1 (phi1) { 
    A0 = E0 / w0; 
    A1 = E1 / w1;
  }
  inline double A (double t) {
    return A0 * cos (w0 * t) + A1 * cos (w1 * t + phi1);
  }
  inline double E (double t) {
    return E0 * sin (w0 * t) + E1 * sin (w1 * t + phi1);
  }
};

class Pulse_sin2 {
  double w;   // frequency
  double E0;  // field amplitude
  double tc;  // center of pulse
  double nc;  // number of cycles
  double A0;
  double duration;
 public:
  Pulse_sin2 (double w, double E0, double tc, double nc) :
  w (w), E0 (E0), tc (tc), nc (nc) {
    A0 = E0 / w;
    duration = 2. * M_PI * nc / w;
  }
  double A (double t) {
    double result = 0.;
    auto phi = w * (t - tc + duration / 2.);
    if (phi >= 0 && phi < w * duration) {
      auto e = sin (phi / (2.*nc));
      result = A0 * e * e * sin (phi);
    }
    return result;
  }
  double E (double t) {
    double result = 0.;
    auto phi = w * (t - tc + duration / 2.);
    if (phi >= 0 && phi < w * duration) {
      auto phi_m = (1. - 1./nc) * phi;
      auto phi_p = (1. + 1./nc) * phi;
      result = - A0 / (4.*(nc*nc-1)*w)
        * (2.*(1.-nc*nc) * cos (phi)
           + nc * ((nc + 1.) * cos (phi_m) + (nc - 1.) * cos (phi_p)));
    }
    return result;
  }
};

class Pulse_gaussian {
  double w;   // frequency
  double E0;  // field amplitude
  double tc;  // center of pulse
  double fwhm;// FWHM
  double A0;
  double sigma;
 public:
  Pulse_gaussian (double w, double E0, double tc, double fwhm) :
  w (w), E0 (E0), tc (tc), fwhm (fwhm) {
    A0 = E0 / w;
    sigma = fwhm / (2. * sqrt (2. * log (2.)));
  }
  inline double A (double t) {
    auto tt = t - tc;
    double res = A0 * exp (- tt * tt / (2. * sigma * sigma)) * sin (w * tt);
    if (std::isnan (res)) res = 0.;
    return res;
  }
  inline double E (double t) {
    auto tt = t - tc;
    // double res = A0 * exp (-0.5 * sigma * sigma * w * w) * sqrt (0.5*M_PI) * sigma
    //   * imag (cerf (M_SQRT1_2 * complex (tt/sigma, sigma * w)));
    double sigma2 = sigma * sigma;
    double res = A0 / sigma2 * exp (- tt * tt / (2. * sigma2))
      * (tt * sin (w * tt) - w * sigma2 * cos (w * tt));
    if (std::isnan (res)) res = 0.;
    return res;
  }
};

class Pulse_gaussian_and_THz {
  double w0;  // frequency
  double E0;  // field amplitude
  double tc;  // center of pulse
  double fwhm;// FWHM
  double A0;
  double sigma;
  double E1;
  double w1;
  double A1;
  double phi1;
 public:
  Pulse_gaussian_and_THz (double w0, double E0, double tc, double fwhm, double w1, double E1, double phi1) :
    w0 (w0), E0 (E0), tc (tc), fwhm (fwhm), w1 (w1), E1 (E1), phi1 (phi1) {
    A0 = E0 / w0;
    A1 = E1 / w1;
    sigma = fwhm / (2. * sqrt (2. * log (2.)));
  }
  inline double A (double t) {
    auto tt = t - tc;
    auto res = A0 * exp (- tt * tt / (2. * sigma * sigma)) * sin (w0 * tt);
    if (std::isnan (res)) res = 0.;
    auto phi = w1 * t + phi1;
    if (phi >= 0. && phi < 2*M_PI)
      res += A1 * cos (phi) - A1;
    return res;
  }
  inline double E (double t) {
    auto tt = t - tc;
    auto sigma2 = sigma * sigma;
    auto res = A0 / sigma2 * exp (- tt * tt / (2. * sigma2))
      * (tt * sin (w0 * tt) - w0 * sigma2 * cos (w0 * tt));
    if (std::isnan (res)) res = 0.;
    auto phi = w1 * t + phi1;
    if (phi >= 0. && phi < 2*M_PI)
      res += A1 * w1 * sin(phi);
    return res;
  }
};

class Pulse_THz {
  double E1;
  double w1;
  double A1;
  double phi1;
 public:
  Pulse_THz (double w1, double E1, double phi1) :
    w1 (w1), E1 (E1), phi1 (phi1) {
    A1 = E1 / w1;
  }
  inline double A (double t) {
    double res = 0.;
    auto phi = w1 * t + phi1;
    if (phi >= 0. && phi < 2*M_PI)
      res += A1 * cos (phi) - A1;
    return res;
  }
  inline double E (double t) {
    double res = 0.;
    auto phi = w1 * t + phi1;
    if (phi >= 0. && phi < 2*M_PI)
      res += A1 * w1 * sin(phi);
    return res;
  } 
};

class Pulse_gaussian_and_THz_practical {
  double w0;  // frequency
  double E0;  // field amplitude
  double tc;  // center of pulse
  double fwhm;// FWHM
  double A0;
  double sigma;
  double E1;
  double w1;
  double A1;
  double phi1;
  const double a = 50.;
 public:
  Pulse_gaussian_and_THz_practical (double w0, double E0, double tc, double fwhm, double w1, double E1, double phi1) :
    w0 (w0), E0 (E0), tc (tc), fwhm (fwhm), w1 (w1), E1 (E1), phi1 (phi1) {
    A0 = E0 / w0;
    A1 = E1 / w1;
    sigma = fwhm / (2. * sqrt (2. * log (2.)));
  }
  inline double A (double t) {
    auto tt = t - tc;
    auto res = A0 * exp (- tt * tt / (2. * sigma * sigma)) * sin (w0 * tt);
    if (std::isnan (res)) res = 0.;

    double res_t = 0.;
    if (t >= 0) {
      auto phi = w1 * t;
      res_t += - A1 * pow (phi, -10.) / (exp (a/phi) - 1.) * sin (phi + phi1) / 4.67275e-12;
      if (std::isnan (res_t)) res_t = 0.;
    }
    res += res_t;
    return res;
  }
  inline double E (double t) {
    auto tt = t - tc;
    auto sigma2 = sigma * sigma;
    auto res = A0 / sigma2 * exp (- tt * tt / (2. * sigma2))
      * (tt * sin (w0 * tt) - w0 * sigma2 * cos (w0 * tt));
    if (std::isnan (res)) res = 0.;

    double res_t = 0.;
    if (t >= 0) {
      auto phi = w1 * t;
      auto et = exp (a / phi) - 1.;
      res_t += A1 * (cos (phi+phi1) / (et * pow (phi, 9.) * t)) / 4.67275e-12;
      res_t += A1 * (a * exp (a/phi) * sin (phi+phi1) / (et*et * pow (phi, 11.) * t)) / 4.67275e-12;
      res_t += A1 * (-10 * sin (phi+phi1) / (et * pow (phi, 10.) * t)) / 4.67275e-12;
      if (std::isnan (res_t)) res_t = 0.;
    }
    res += res_t;
    return res;
  }
};

class Pulse_THz_practical {
  double E1;
  double w1;
  double A1;
  double phi1;
  const double a = 50.;  
 public:
  Pulse_THz_practical (double w1, double E1, double phi1) :
    w1 (w1), E1 (E1), phi1 (phi1) {
    A1 = E1 / w1;
  }
  inline double A (double t) {
    double res = 0.;
    if (t >= 0) {
      auto phi = w1 * t;
      res += - A1 * pow (phi, -10.) / (exp (a/phi) - 1.) * sin (phi + phi1) / 4.67275e-12;
      if (std::isnan (res)) res = 0.;
    }
    return res;    
  }
  inline double E (double t) {
    double res = 0.;
    if (t >= 0) {
      auto phi = w1 * t;
      auto et = exp (a / phi) - 1.;
      res += A1 * (cos (phi+phi1) / (et * pow (phi, 9.) * t)) / 4.67275e-12;
      res += A1 * (a * exp (a/phi) * sin (phi+phi1) / (et*et * pow (phi, 11.) * t)) / 4.67275e-12;
      res += A1 * (-10 * sin (phi+phi1) / (et * pow (phi, 10.) * t)) / 4.67275e-12;
      if (std::isnan (res)) res = 0.;
    }
    return res;    
  } 
};