#include <complex>
#include <array>
#include <iostream>
#include <fstream>
#include "../parameters.h"

using namespace parameters;

complex bCoef (int n, double t) {
  complex bn;
  switch (n) {
  case 1: // -(2n + 1) = -3
    bn = exp (+2.*I*t) * (2 + t*t - 2 * cos (t) - 2 * t * sin (t)) / (8*t*t);
    break;
  case 0: // -(2n + 1) = -1
    bn = (exp (+2.*I*t) * (2. - t * (t + 2.*I)) - exp (+I*t) * (5. - I*t) + 4. + 2.*I*t -
	  exp (-I*t) * (1. + I * t)) / (8*t*t);
    break;
  case -1:  // -(2n + 1) = 1
    bn = (exp (-2.*I*t) * (2. - t * (t - 2.*I)) - exp (-I*t) * (5. + I*t) + 4. - 2.*I*t -
	  exp (+I*t) * (1. - I * t)) / (8*t*t);
    break;
  case -2:  // -(2n + 1) = 3
    bn = exp (-2.*I*t) * (2 + t*t - 2 * cos (t) - 2 * t * sin (t)) / (8*t*t);
    break;
  default:
    std::cerr << "undefined M in bCoef!!!" << std::endl;
  }
  return bn * A0 * A0 * E0;
}

inline double besselJ (int n, double x) {
  double (*J) (double, double) = std::cyl_bessel_j;
  double nn = (double) n;
  return x > 0 ? J (nn, x) : n % 2 == 0 ? J (nn, -x) : -J (nn, -x);
}

complex integrand_odd (int K, int M, double tau) {
  const double eps = 0.01;
  auto diffuseCoef = pow (M_PI / (eps + 0.5 * I * (tau)), 1.5);
  auto F0 = (Up + Ip) * tau - 2 * Up * (1. - cos (tau)) / tau;
  F0 /= w0;
  auto C0 = sin (tau) - 4 * sin (0.5*tau) * sin (0.5*tau) / tau;
  C0 /= w0;
  complex result = diffuseCoef;
  result *= bCoef (K-M, tau);
  result *= exp (- I * F0);
  result *= pow (I, M);
  result *= besselJ (M, Up * C0);
  result *= exp (I * double(M) * tau);
  return result;
}

int main (int argc, char* argv[]) {

  setup_params (); // calculate Up
  std::cout << "E0: " << E0 << std::endl;
  std::cout << "w0: " << w0 << std::endl;
  std::cout << "Up: " << Up << std::endl;

  std::ofstream outfile ("res.dat");

  const long nt = 4000*2;
  auto dt = 200*2 / (nt - 1.);

  // n = K - M
  std::array<int, 4> n_set = {-2, -1, 0, 1};
  for (int K = 5; K < 50; K ++) {                  // loop over harmonics
    complex s = 0.;
    for (auto n : n_set) {                         // loop over M for the sum; but actually we loop over n = K-M since bM predetermines the range of n, and now M = K-n
      for (long it = 0; it < nt; it ++) {          // loop over tau for the integration
        auto tau = it * dt;
        auto M = K - n;
        auto res = integrand_odd (K, M, tau);
        if (std::isnan (real (res))) res = 0.;
        // if (it == 0 || it == nt - 1) res *= 0.5;
        s += res;
        // outfile << tau << " "<< real (conj (res) * res) << std::endl;
      }
    }
    s *= I * dt;
    outfile << 2*K+1 << " " << abs (s) << std::endl;
  }


  // std::cout << bCoef (-3, 100.) << " " << bCoef (3, 100.) << std::endl;


  outfile.close ();
  // std::cout << integrand (30, 31, 50.) << std::endl;


  // std::cout << std::cyl_bessel_j (30, -5.2) << std::endl;
  return 0;

}
