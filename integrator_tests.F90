program integrator_tests
  use integrator_main
  use microphysics_constants, only: rho_l => rho_w, g_earth, pi, rho_h
  use integrator_config
  use mphys_with_ice, only: rho_f
  use integrator_tests_precip, only: rain_after_precip, graupel_after_precip

  real(kreal), dimension (1 : 7) :: y
  real(kreal) :: t_0 = 0.0_kreal
  !real(kreal) :: t_end
  real(kreal) :: t_step = 2.0_kreal / updraft_vel
  real(kreal) :: t_curr = 0.0_kreal
  real(kreal) :: qd, rho_curr, rho_g, height

  y(1) = 298.0_kreal      ! Temperature
  y(2) = 101325.0_kreal   ! Pressure
  y(3) = 0.00_kreal       ! Cloud water
  y(4) = 0.00_kreal       ! Rain
  y(5) = 0.02_kreal       ! Water vapour
  y(6) = 0.00_kreal       ! Ice
  y(7) = 0.00_kreal       ! Graupel

  qv_sat = 0.00_kreal     ! For testing saturation vapour concentration

  open(unit=1,file='output/05.txt',status='unknown') !TODO: Change the output name!

  do
    height = t_curr*updraft_vel
    print *, "Height =>", t_curr*updraft_vel, sum(y(3:))
    print *, y
    qd = 1.0_kreal - y(3) - y(4) - y(5) - y(6) - y(7)
    rho_curr = rho_f(qd=qd, qv=y(5), ql=y(3), qr=y(4), qi=y(6), qh=y(7), p=y(2), temp=y(1))
    write(1,*) &
      height,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',',y(6),',',y(7),',',rho_curr

    if (height >= 15000.0_kreal .or. y(1) < 233.15_kreal) exit

    call integrate(y=y, t0=t_curr, t_end=t_curr+t_step)

    ! Precipitation
    rho_curr = rho_f(qd=qd, qv=y(5), ql=y(3), qr=y(4), qi=y(6), qh=y(7), p=y(2), temp=y(1))
    rho_g = rho_f(qd=qd, qv=y(5), ql=0.0_kreal, qr=0.0_kreal, qi= 0.0_kreal, qh=0.0_kreal, p=y(2), temp=y(1))

    if (rain_precip) then
      y(4) = rain_after_precip(rho_g=rho_g, rho=rho_curr, temp=y(1), p=y(2), qr=y(4), del_t=t_step, del_z=cloud_thickness)
    endif
    if (graupel_precip) then
      y(7) = graupel_after_precip(rho_g=rho_g, rho=rho_curr, temp=y(1), p=y(2), qh=y(7), del_t=t_step, del_z=cloud_thickness)
    endif

    y(1) = y(1) - lapse_rate*updraft_vel*t_step
    qd = 1.0_kreal - y(3) - y(4) - y(5) - y(6) - y(7)
    rho_curr = rho_f(qd=qd, qv=y(5), ql=y(3), qr=y(4), qi=y(6), qh=y(7), p=y(2), temp=y(1))
    y(2) = y(2) - (updraft_vel*t_step*g_earth*rho_curr)

    ! Timestep
    t_curr = t_curr + t_step
  end do

end program integrator_tests
