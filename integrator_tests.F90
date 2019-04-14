program integrator_tests
  use integrator_main
  !print *, rho_f(0.0_kreal, 0.0_kreal, 0.0_kreal, 1.0_kreal, 0.0_kreal, &
  !0.0_kreal, 100000.0_kreal, 298.0_kreal)
  real(kreal), dimension (1 : 7) :: y
  real(kreal) :: t_0 = 0.0_kreal
  real(kreal) :: t_end=1000.0_kreal
  real(kreal) :: t_step=10.0_kreal
  real(kreal) :: t_curr = 0.0_kreal

  y(1) = 298.0_kreal
  y(2) = 101325.0_kreal
  y(3) = 0.00_kreal
  y(4) = 0.00_kreal
  y(5) = 0.02_kreal
  y(6) = 0.00_kreal
  y(7) = 0.00_kreal

  do
    print *, t_curr
    print *, y
    if (t_curr >= t_end) exit
    call integrate(y=y, t0=t_curr, t_end=t_curr+t_step)
    y(1) = y(1) - 1.96_kreal
    t_curr = t_curr + t_step
  end do

end program integrator_tests
