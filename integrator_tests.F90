program integrator_tests
  use integrator_main
  !print *, rho_f(0.0_kreal, 0.0_kreal, 0.0_kreal, 1.0_kreal, 0.0_kreal, &
  !0.0_kreal, 100000.0_kreal, 298.0_kreal)
  real(8), dimension (1 : 7) :: y
  real(8) :: t0 = 0.0
  real(8) :: t_end=1000.0
  y(1) = 298.0
  y(2) = 101325.0
  y(3) = 0.01
  y(4) = 0.01
  y(5) = 0.01
  y(6) = 0.0
  y(7) = 0.0
  call integrate(y=y, t0=t0, t_end=t_end)
end program integrator_tests
