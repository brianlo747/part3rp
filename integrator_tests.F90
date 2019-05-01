program integrator_tests
  use integrator_main
  use mphys_with_ice, only: rho_f
  !print *, rho_f(0.0_kreal, 0.0_kreal, 0.0_kreal, 1.0_kreal, 0.0_kreal, &
  !0.0_kreal, 100000.0_kreal, 298.0_kreal)
  real(kreal), dimension (1 : 7) :: y
  real(kreal) :: t_0 = 0.0_kreal
  real(kreal) :: t_end=7000.0_kreal
  real(kreal) :: t_step=10.0_kreal
  real(kreal) :: t_curr = 0.0_kreal
  real(kreal) :: qd, rho_curr

  y(1) = 298.0_kreal
  y(2) = 101325.0_kreal
  y(3) = 0.00_kreal
  y(4) = 0.00_kreal
  y(5) = 0.02_kreal
  y(6) = 0.00_kreal
  y(7) = 0.00_kreal

  open(unit=1,file='output/test.txt',status='unknown')

  do
    print *, "t =>", t_curr, sum(y(3:))
    print *, y
    write(1,*) &
      t_curr,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',',y(6),',',y(7)

    if (t_curr >= t_end) exit

    call integrate(y=y, t0=t_curr, t_end=t_curr+t_step)
    y(1) = y(1) - 0.196_kreal
    qd = 1.0_kreal - y(3) - y(4) - y(5) - y(6) - y(7)
    rho_curr = rho_f(qd=qd, qv=y(5), ql=y(3), qr=y(4), qi=y(6), qh=y(7), p=y(2), temp=y(1))
    y(2) = y(2) - (20.0_kreal*9.80665_kreal*rho_curr)
    t_curr = t_curr + t_step
  end do

end program integrator_tests
