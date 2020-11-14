void set_global_params () {
  w0 = 0.0353;
  E0 = 0.06;

  w1 = 0.000353;
  E1 = 2e-5;
  phi1 = 0.5*M_PI;  

  t0 = 0.;
  t1 = 8899.7 * 4;
  nt = int ((t1 - t0) / 0.988855) + 1;
  //nt = int ((t1 - t0) / 0.8) + 1;
  std::cout << "nt " << nt << " @ set_global_params" << std::endl;

  tp0 = t0;
  tp1 = t1;
}

template <typename SFA_CALC_MODEL>
void calculate (const std::string& filename_para, const std::string& filename_data) {

  std::ofstream of;
  of.open (filename_para);
  of << Ip << " " << w0 << " " << E0 << std::endl;
  of << t0 << " " << t1 << " " << dt << std::endl;
  of.close ();

  Pulse pulse (w0, E0, w1, E1, phi1);
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