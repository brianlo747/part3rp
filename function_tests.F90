program function_tests
  use mphys_with_ice
  real(kreal) :: rho_f_output
  print *, rho_f(0.1_kreal, 0.1_kreal, 0.1_kreal, 0.1_kreal, 0.1_kreal, &
  0.1_kreal, 100000.0_kreal, 298.0_kreal)
  print *, rho_f_output
end program function_tests
