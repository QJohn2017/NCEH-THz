#include <iostream>
#include <fstream>
#include <complex>
#include "../erf.h"

#define complex std::complex<double>
#define I complex(0., 1.)

inline complex w (complex z) {
    return w_of_z (z);
    return I * z / (sqrt (M_PI) * (z*z - 0.5));
    return I * (z * z - 1.) / (z * sqrt (M_PI) * (z * z - 1.5));

/*     
    auto z2 = z*z;
    auto z4 = z2*z2;
    return (945./(32.*z4*z4*z2) + 105./(16.*z4*z4) + 15./(8.*z4*z2) + 3./(4.*z4) + 1./(2.*z2) + 1.) / (-I * z * sqrt (M_PI)); 
*/
}

complex g (complex z, double e) {
    return w (I * M_SQRT2 * z) - exp (4.*e*(z-e/2.)) * w (I * M_SQRT2 * (z - e));
}

double calc_G (double x, double y, double e) {
    const double sqrt_pi = 2. / M_2_SQRTPI; // sqrt (M_PI);
    complex z = complex (x, -y);
    double im = imag (exp (2.*I*x*y) * g (z*M_SQRT1_2, e*M_SQRT1_2));
    return real (exp (4.*I*x*y) * g (z, e) - g (x, e)) - M_SQRT2 * sqrt_pi / e * im * im;
}

double calc_Phi_of_mono (double phi_t, double phi_tau) {
    auto s = sin (0.5 * phi_tau);
    auto C0 = sin (phi_tau) - 4. * s * s / phi_tau;
    return phi_tau - 2. / phi_tau * (1 - cos (phi_tau)) - C0 * cos (2*phi_t - phi_tau);
}

double calc_Phi_of_gaus (double phi_t, double phi_sig, double phi_tau) {
    return 0.5 * sqrt (M_PI) * phi_sig * calc_G (phi_t/phi_sig*M_SQRT1_2, phi_sig*M_SQRT1_2, phi_tau/phi_sig*M_SQRT1_2);
}


int generate_figure_data () {

    double phi_tau = 4.;
    double phi_sig = 50.;

    const int nt = 4000;
    double t0 = -100.;
    double t1 = 100.;
    double dt = (t1 - t0) / (nt - 1.);

    std::ofstream file;

    // Phi_cw
    file.open ("res_mono.dat");
    for (int it = 0; it < nt; it ++) {            
        double phi_t = t0 + it * dt;
        file << phi_t << " "
            << calc_Phi_of_mono (phi_t, phi_tau)
            << std::endl;
    }
    file.close ();

    // Phi_gauss
    file.open ("res_gaus.dat");
    for (int it = 0; it < nt; it ++) {            
        double phi_t = t0 + it * dt;
        file << phi_t << " "
            << calc_Phi_of_gaus (phi_t, phi_sig, phi_tau)
            << std::endl;
    }
    file.close ();

    // Phi_gauss * exp (-(t/sigma)^2)
    file.open ("res_gaus_env.dat");
    for (int it = 0; it < nt; it ++) {            
        double phi_t = t0 + it * dt;
        file << phi_t << " "
            << calc_Phi_of_gaus (phi_t, phi_sig, phi_tau) * exp (- phi_t*phi_t / (phi_sig*phi_sig))
            << std::endl;
    }
    file.close ();
}


int main (int argc, char* argv[]) {

    // generate data for paper figure
    generate_figure_data ();

/*
    double phi_tau = 4.;
    double phi_sig = 50.;

    std::ofstream file;
    //file.open ("res.dat");
    //file.open ("res_mono.dat");
    //file.open ("res_gaus.dat");
    file.open ("res_gaus_env.dat");

    const int nt = 4000;
    double t0 = -100.;
    double t1 = 100.;
    double dt = (t1 - t0) / (nt - 1.);
*/
/*
    for (int it = 0; it < nt; it ++) {
        double phi_t = t0 + it * dt;
        file << phi_t << " " 
            << calc_G (phi_t/phi_sig*M_SQRT1_2, phi_sig*M_SQRT1_2, phi_tau/phi_sig*M_SQRT1_2) 
            << std::endl;
    }
*/

/*
    // compare the gaussian-enveloped pulse and the monochromatic laser
    auto C0 = [] (double phi_tau) {
        double s = sin (phi_tau / 2.);
        return sin (phi_tau) - 4. / phi_tau * s * s;
    };
    for (int it = 0; it < nt; it ++) {
        double phi_t = t0 + it * dt;
        double term_mono = phi_tau - 2./phi_tau * (1 - cos (phi_tau)) - C0 (phi_tau) * cos (2.*phi_t - phi_tau);
        double term_gaus = 1./M_2_SQRTPI * phi_sig * calc_G (phi_t/phi_sig*M_SQRT1_2, phi_sig*M_SQRT1_2, phi_tau/phi_sig*M_SQRT1_2);
        file << phi_t << " " << term_mono << " " << term_gaus << std::endl;
    }
*/
/*
    file.close ();
*/
    return 0;
}