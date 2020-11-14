#pragma once

namespace parameters {

#define complex std::complex<double>
const complex I = complex (0., 1.);

// fs laser pulse (mono)
double Ip = 0.5;
double w0 = 0.0353;
double E0 = 0.06;

// THz field
// double w1   = 1e-4;
double w1 = 7.06e-4;
double E1   = 2e-5;
double phi1 = 0.5*M_PI;

double fwhm  = 100 * 41.34;
double tc    = 6000.;
double half_duration = fwhm * 1.478;
double sigma = fwhm / (2. * sqrt (2. * log (2.)));


// used for Gaussian enveloped pulse
int    nt  = int (half_duration * 2.27);
double t0  = tc - half_duration;
double t1  = tc + half_duration;
double tp0 = -120;
double tp1 = 13000;

/*
// used for monochromatic continuous wave
int    nt = 9001;
double t0 = 0;
double t1 = 2.*M_PI/w0 * 50;
double tp0 = 0;
double tp1 = t1;
*/


double dt;
double A0;
double A1;
double Up;

int global_tag = 0;

void setup_params () {
  dt = (t1 - t0) / (nt - 1.);
  A0 = E0 / w0;
  A1 = E1 / w1;
  Up = E0 * E0 / (4. * w0 * w0);
}


inline complex dipole (double p) {
  return I * p;
  // double dnom = (p*p + 2*Ip);
  // return - I * 11.3137 * pow (2.*Ip, 1.25) * p
  //   / (M_PI * dnom * dnom * dnom); // pow(2., 3.5)
  // double alpha = 2.*Ip;
  // return I * pow (1./ (M_PI*alpha), 3./4.) * p / alpha * exp (- p*p / (2.*alpha));
}


}
