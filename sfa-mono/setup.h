void set_global_params () {
  w0 = 0.0353;
  t0 = 0.;
  t1 = 8899.7;
  nt = int ((t1 - t0) / 0.988855) + 1;
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

  Pulse pulse (w0, E0);
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