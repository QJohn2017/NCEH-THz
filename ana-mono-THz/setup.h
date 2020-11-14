void set_global_params () {
  w0 = 0.0353;
  E0 = 0.06;

  w1 = 0.000353;
  E1 = 2e-5;
  phi1 = 0.5*M_PI;

  t0 = 0.;
  t1 = 8899.7 * 4;
  // nt = int ((t1 - t0) / 0.988855) + 1;
  nt = int ((t1 - t0) / 0.8) + 1;
}