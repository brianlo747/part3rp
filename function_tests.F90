program function_tests
  use mphys_with_ice
  print *, "rho_f"
  print *, rho_f(qd=0.95_kreal, qv=0.01_kreal, ql=0.01_kreal, qr=0.01_kreal, &
    qi=0.01_kreal, qh=0.01_kreal, p=101325.0_kreal, temp=298.0_kreal)

  print *, "dqldt_condevap"
  print *, dql_dt__condensation_evaporation(rho=1.22606402885214_kreal, &
    rho_g=1.22615614877405_kreal, qv=0.01_kreal, ql=0.01_kreal, T=298.0_kreal, &
    p=101325.0_kreal)

  print *, "dqrdt_condevap"
  print *, dqr_dt__condensation_evaporation(qv=0.01_kreal, qr=0.01_kreal, &
    rho=1.22606402885214_kreal, T=298.0_kreal, p=101325.0_kreal)

  print *, "dqhdt_condevap"
  print *, dqh_dt__condensation_evaporation(0.2_kreal, 0.2_kreal, 1.5_kreal, &
  298.0_kreal, 100000.0_kreal)

  print *, "dqidt_sublidep"
  print *, dqi_dt__sublimation_deposition(0.2_kreal, 1.5_kreal, 260.0_kreal, &
  100000.0_kreal)

  print *, "dqhdt_sublidep"
  print *, dqh_dt__sublimation_evaporation(0.95_kreal, 0.01_kreal, 0.01_kreal, 300.0_kreal, &
  298.0_kreal, 100000.0_kreal)

  print *, "dqrdt_autoconv"
  print *, dqr_dt__autoconversion(0.5_kreal, 0.1_kreal, 1.2_kreal, 298.0_kreal)

  print *, "dqhdt_autoconv"
  print *, dqh_dt__autoconversion_ice_graupel(0.5_kreal, 0.1_kreal, 1.2_kreal, &
  265.0_kreal)

  print *, "dqrdt_accre_rc"
  print *, dqr_dt__accretion_cloud_rain(ql=0.01_kreal, rho=1.22606402885214_kreal, &
    rho_g=1.22615614877405_kreal, qr=0.01_kreal)

  print *, "dqhdt_accre_hi"
  print *, dqh_dt__accretion_ice_graupel(0.01_kreal, 1.2_kreal, 0.01_kreal, &
  250.0_kreal)

  print *, "dqldt_accre_chr"
  print *, dqh_dt__accretion_cloud_graupel_rain(0.01_kreal, 1.21_kreal, &
  1.20_kreal, 0.01_kreal)

  print *, "dqhdt_accre_hiri"
  print *, dqr_dt__accretion_ice_rain_graupel_i(0.01_kreal, 1.2_kreal, &
  1.21_kreal, 0.01_kreal)

  print *, "dqhdt_accre_hirr"
  print *, dqr_dt__accretion_ice_rain_graupel_r(0.01_kreal, 0.01_kreal, &
  1.21_kreal, 1.2_kreal, 0.01_kreal)

  print *, "dqhdt_accre_hr"
  print *, dqr_dt__accretion_graupel(0.01_kreal, 1.2_kreal, &
  0.01_kreal, 0.01_kreal, 1.2_kreal, 270.0_kreal, 100000.0_kreal, 0.01_kreal)

  print *, "dqrdt_melt_rh"
  print *, dqr_dt__melting_graupel(0.01_kreal, 1.2_kreal, 0.01_kreal, 0.01_kreal, &
  1.21_kreal, 298.0_kreal, 100000.0_kreal, 0.01_kreal, 0.01_kreal)

  print *, "dqldt_melt_ci"
  print *, dqr_dt__melting_ice(0.01_kreal, 1.2_kreal, 0.01_kreal, 0.01_kreal, &
  1.21_kreal, 298.0_kreal, 100000.0_kreal, 0.01_kreal, 0.01_kreal)

  print *, "dqidt_freeze"
  print *, dqi_dt__freezing_ice(0.01_kreal, 1.2_kreal, 260.0_kreal)

  print *, "dqhdt_freeze"
  print *, dqh_dt__freezing_graupel(0.1_kreal, 1.6_kreal, 260.0_kreal)
end program function_tests
