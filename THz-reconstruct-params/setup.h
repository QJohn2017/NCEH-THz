void set_global_params () {


  w0 = 0.1139; // 400 nm

  // parameters of THz field
  // w1   = 2e-4;
  w1 = 3e-3; // 20 THz
  E1   = 2e-5;
  phi1 = 0.3*M_PI;

  // parameters of fs pulse
  fwhm  = 40 * 41.34;//120 * 41.34;
  // tc    = 10000.; // 6000
  tc    = 1750.; // for w1 = 3e-3
  half_duration = fwhm * 1.478;
  sigma = fwhm / (2. * sqrt (2. * log (2.)));

  // parameters of time grid for integration
  nt  = int (half_duration * 2.27);
  t0  = tc - half_duration;
  t1  = tc + half_duration;
  // tp0 =  0. - 8000;
  // tp1 =  94247. + 8000.; // for w1 = 2e-4
  tp0 = 0. - 3000.;
  tp1 = 6283. + 3000.; // for w1 = 3e-3
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