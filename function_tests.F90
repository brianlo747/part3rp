program function_tests
  use mphys_with_ice
  print *, rho_f(0.0_kreal, 0.0_kreal, 0.0_kreal, 1.0_kreal, 0.0_kreal, &
  0.0_kreal, 100000.0_kreal, 298.0_kreal)
  print *, dql_dt__condensation_evaporation(1.5_kreal, 1.2_kreal, &
  0.1_kreal, 0.1_kreal, 298.0_kreal, 100000.0_kreal)
  print *, dqr_dt__condensation_evaporation(0.2_kreal, 0.2_kreal, 1.5_kreal, 298.0_kreal, &
  100000.0_kreal)
  print *, dqh_dt__condensation_evaporation(0.2_kreal, 0.2_kreal, 1.5_kreal, &
  298.0_kreal, 100000.0_kreal)
  print *, dqi_dt__sublimation_deposition(0.2_kreal, 1.5_kreal, 260.0_kreal, &
  100000.0_kreal)
  print *, dqr_dt__autoconversion(0.5_kreal, 0.1_kreal, 1.2_kreal, 298.0_kreal)
  print *, dqh_dt__autoconversion_ice_graupel(0.5_kreal, 0.1_kreal, 1.2_kreal, &
  265.0_kreal)
  print *, dqr_dt__accretion_cloud_rain(0.05_kreal, 1.2_kreal, 0.01_kreal)

  print *, dqh_dt__accretion_ice_graupel(0.01_kreal, 1.2_kreal, 0.01_kreal, &
  250.0_kreal)
  print *, dqh_dt__accretion_cloud_graupel_rain(0.01_kreal, 1.21_kreal, &
  1.20_kreal, 0.01_kreal)
  print *, dqr_dt__accretion_ice_rain_graupel_i(0.01_kreal, 1.2_kreal, &
  1.21_kreal, 0.01_kreal)
  print *, "Fixed?"
  print *, dqr_dt__accretion_ice_rain_graupel_r(0.01_kreal, 0.01_kreal, &
  1.21_kreal, 1.2_kreal, 0.01_kreal)
  print *, dqr_dt__accretion_graupel(0.01_kreal, 1.2_kreal, &
  0.01_kreal, 0.01_kreal, 1.2_kreal, 270.0_kreal, 100000.0_kreal, 0.01_kreal)
  print *, dqr_dt__melting_graupel(0.01_kreal, 1.2_kreal, 0.01_kreal, 0.01_kreal, &
  1.21_kreal, 298.0_kreal, 100000.0_kreal, 0.01_kreal, 0.01_kreal)
  print *, dqi_dt__freezing_graupel(0.01_kreal, 1.2_kreal, 260.0_kreal)
  print *, dqh_dt__freezing_ice(0.1_kreal, 1.6_kreal, 260.0_kreal)
end program function_tests
