program integrator_tests
  use integrator_main
  !print *, rho_f(0.0_kreal, 0.0_kreal, 0.0_kreal, 1.0_kreal, 0.0_kreal, &
  !0.0_kreal, 100000.0_kreal, 298.0_kreal)
  real(kreal), dimension (1 : 7) :: y
  real(kreal) :: t0 = 0.0_kreal
  real(kreal) :: t_end=10.0_kreal
  y(1) = 298.0_kreal
  y(2) = 101325.0_kreal
  y(3) = 0.0_kreal
  y(4) = 0.0_kreal
  y(5) = 0.01_kreal
  y(6) = 0.0_kreal
  y(7) = 0.0_kreal
  call integrate(y=y, t0=t0, t_end=t_end)
end program integrator_tests
