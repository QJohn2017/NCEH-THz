void set_global_params () {
  // parameters of THz field
  w1   = 1e-4;
  E1   = 2e-5;
  phi1 = 0.5*M_PI;

  // parameters of fs pulse
  fwhm  = 120 * 41.34;
  tc    = 0.; // 6000
  half_duration = fwhm * 1.478;
  sigma = fwhm / (2. * sqrt (2. * log (2.)));

  // parameters of time grid for integration
  nt  = int (half_duration * 2.27);
  t0  = tc - half_duration;
  t1  = tc + half_duration;
  tp0 = -63000.;
  tp1 =  63000.;
}

template <typename SFA_CALC_MODEL>
void calculate (const std::string& filename_para, const std::string& filename_data) {

  std::ofstream of;
  of.open (filename_para);
  of << Ip << " " << w0 << " " << E0 << std::endl;
  of << t0 << " " << t1 << " " << dt << " " << tc << std::endl;
  of.close ();

  Pulse pulse (w0, E0, tc, fwhm, w1, E1, phi1);
  std::ofstream of_res;

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