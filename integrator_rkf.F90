module integrator_rkf
  use microphysics_constants, only: kreal
  use microphysics_constants, only: abs_tol => integration_abs_tol, rel_tol => integration_rel_tol
  use integrator_species_rate, only: dydt => dydt_mphys_with_ice

  implicit none

  public integrate_with_message

  logical, parameter :: debug = .false.
  real(8), parameter :: dt_min = 1.0e-10
  integer, parameter :: max_steps = 10

contains
  !> Outer loop on integrator that guarantees that we reach `t_end` and
  !> keeps checking `msg` to see if we have an error

  !> Calculate density of mixture from equation of state
  pure function density_eos(y) result(rho)
    !use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain
    !use microphysics_register, only: idx_temp, idx_pressure
    use microphysics_constants, only: R_v, R_d, rho_h, rho_i, rho_l => rho_w

    real(8), intent(in), dimension(:) :: y

    real(8) :: temp, p, qd, qv, ql, qr, qi, qh
    real(8) :: rho

    temp = y(1)
    p = y(2)

    ql = y(3)
    !if (idx_rain /= 0) then
    qr = y(4)
    !else
    !  qr = 0.
    !endif
    qv = y(5)
    qi = y(6)
    qh = y(7)
    qd = 1.0 - ql - qr - qv - qi - qh

    rho = 1.0/( (qd*R_d + qv*R_v)*temp/p + (ql+qr)/rho_l + qi/rho_i + qh/rho_h )
  end function density_eos

  !> Compute pressure from equation of state.
  ! Because state representation doesn't contain density (it is diagnosed)
  ! we need to supply a density so that pressure can be calculated from the
  ! other state variables
  pure function pressure_eos(y, rho) result(p)
    !use microphysics_register, only: idx_temp
    !use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain
    use microphysics_constants, only: R_v, R_d, rho_h, rho_i, rho_l => rho_w

    real(8), intent(in) :: rho
    real(8), intent(in), dimension(:) :: y

    real(8) :: temp, qd, qv, ql, qr, qi, qh
    real(8) :: p

    temp = y(1)

    ql = y(3)
    !if (idx_rain /= 0) then
    qr = y(4)
    !else
    !   qr = 0.
    !endif
    qv = y(5)
    qi = y(6)
    qh = y(7)
    qd = 1.0 - ql - qr - qv - qi - qh

    p = temp*(qd*R_d + qv*R_v)/(1./rho - (ql + qr)/rho_l - qi/rho_i - qh/rho_h)
  end function pressure_eos

  function fix_y_isometric(y, dy) result(y_new)
    !use microphysics_register, only: idx_pressure

    real(8), intent(in), dimension(:) :: y, dy
    real(8) :: y_new(size(y)), rho_old

    rho_old = density_eos(y)
    y_new = y + dy

    ! need to update pressure, since we assume constant volume density is
    ! unchanged, use old density to calculate update pressure
    y_new(2) = pressure_eos(y_new, rho_old)
  end function fix_y_isometric

  subroutine mass_scale(y, msg)

    real(8), intent(inout), dimension(:) :: y
    character(len=100), intent(out) :: msg
    real(8) :: q_negative_sum, q_positive_sum

    q_negative_sum = sum(y(3:), mask = y(3:) < 0.0_kreal)
    q_positive_sum = sum(y(3:), mask = y(3:) >= 0.0_kreal)

    if (q_negative_sum < -1.0_kreal) then
      msg = "failed to scale mass"
      return
    endif

    where (y(3:) < 0.0_kreal)
      y(3:) = 0.0_kreal
    elsewhere
      y(3:) = y(3:) + y(3:) * (-1.0_kreal*q_negative_sum) / q_positive_sum
    endwhere
  end subroutine

  subroutine integrate_with_message(y, t0, t_end, msg, m_total)
    real(8), intent(inout), dimension(:) :: y
    real(8), intent(in) :: t_end, t0
    real(8) :: dt_s, t  ! sub-step timestep
    integer, intent(out) :: m_total

    character(len=100), intent(out) :: msg
    integer :: m
    !f2py raise_python_exception msg

    msg = " "

    ! make initial sub-cycling full timestep, in case we can do that
    dt_s = t_end - t0
    t = t0

    do while (t < t_end)
      !print *, "Step start", t, dt_s
      if (t + dt_s > t_end) then
        dt_s = t_end - t
      endif

      ! evolve y, t and dt_s using the integrator
      ! dt_s will be the time-step that was actually used, so that we can
      ! use that for the next integration step

      if (debug) then
        print *, "=> t", t
      endif

      m = 0
      call rkf34_original(y, t, dt_s, msg, m)
      m_total = m_total + m

      if (msg(1:1) /= " ") then
        exit
      endif

    enddo
  end subroutine integrate_with_message

  recursive subroutine rkf34_original(y, t, dt, msg, m)
    character(len=100), intent(out) :: msg
    real(8), intent(inout) :: dt
    real(8), intent(inout), dimension(:) :: y
    real(8), intent(inout) :: t
    integer, intent(in) :: m

    real(8) :: max_abs_err
    real(8) :: dt_min__posdef
    real(8) :: max_rel_err, err

    !f2py raise_python_exception msg

    ! Coefficients for RKFelhberg4(5)
     ! real(8), parameter :: &
     ! a2=1./4.,   b21=1./4., &
     ! a3=3./8.,  b31=3./32.,   b32= 9./32., &
     ! a4=12./13., b41=1932./2197., b42=-7200./2197., b43=7296./2197., &
     ! a5=1.0,     b51=439./216.,   b52= -8.0,        b53=3680./513.,   b54=-845./4104., &
     ! a6=1./2.,  b61=-8./27.,  b62=2.0, b63=-3544./2565., b64=1859./4104., b65=-11./40.
     ! real(8) :: &
     ! c1_1=25./216.,    c2_1=0.0,  c3_1=1408./2565.,  c4_1=2197./4104., c5_1=-1./5., c6_1=0.0, &
     ! c1_2=16./135.,  c2_2=0.0,  c3_2=6656./12825.,   c4_2=28561./56430., c5_2=-9./50., c6_2=2./55.

    ! Coefficients for RKFelhberg3(4)
    real(8), parameter :: &
    a2=2./7.,   b21=2./7., &
    a3=7./15.,  b31=77./900.,   b32= 343./900., &
    a4=35./38., b41=805./1444., b42=-77175./54872., b43=97125./54872., &
    a5=1.0,     b51=79./490.,   b52= 0.0,           b53=2175./3626.,   b54=2166./9065.
    real(8) :: &
    c1_1=79./490.,    c2_1=0.0,  c3_1=2175./3626.,  c4_1=2166./9065., c5_1=0.0, &
    c1_2=229./1470.,  c2_2=0.0,  c3_2=1125./1813.,   c4_2=13718./81585., c5_2=1./18.

    ! Coefficients for RK3(4) method
    ! real(8), parameter :: &
    ! !a1=0.0, &
    ! a2=0.5, b21=0.5, &
    ! a3=0.5, b31=0.0,        b32= 0.5, &
    ! a4=1.0, b41=0.0,        b42= 0.0,           b43=1.0, &
    ! a5=1.0, b51=1./6.,      b52= 1./3.,         b53=1./3.,   b54= 1./6.
    ! real(8), parameter :: &
    ! c1_1=0.0, c2_1=0.0,  c3_1=0.0,  c4_1=0.0,  c5_1=1.,&
    ! c1_2=1./6., c2_2=1./3.,    c3_2=1./3.,     c4_2=0.,   c5_2=1./6.

    real(8), dimension(size(y)) :: k1, k2, k3, k4, k5, k6
    real(8), dimension(size(y)) :: abs_err, y_n1, y_n2, rel_err
    real(8), dimension(size(y)) :: dydt0, max_total_err
    real(8) :: s

    logical :: done = .false.

    if (debug) then
      print *, ""
      print *, " -- substep -- m dt", m, dt
      print *, 'y', y
    endif

    done = .false.

    ! Calculate rate of change based on initial t and y
    dydt0 = dydt(t, y)

    ! TODO: use a vector for abs error here, so that we can have a
    ! different abs error for say specific concentration and temperature

    ! State components that are already zero do not contribute to the
    ! relative error calculation

    ! Calculate a timestep that would lead to a species reaching 0
    ! This is to ensure positive definite species values when integrating
    dt_min__posdef = minval(abs(y/dydt0), y /= 0.0_kreal)

    ! Check if positive definite timestep is not too small
    if (dt_min__posdef > dt_min) then

      ! Check if current timestep is too large for positive definite condition
      if (dt > dt_min__posdef) then

        ! Timestep too large for positive definite condition
        if (debug) then
          print *, "Adjusting integration step down, too big to be pos def"
          print *, "dt dt_min__posdef, m", dt, dt_min__posdef, m
          print *, ""
        endif

        ! Scaling includes half for safety
        s = 0.5*dt_min__posdef/dt
        dt = s*dt

        if (debug) then
          print *, "Initial timestep risks solution becoming negative, scaling timestep"
          print *, "by s=", s
          print *, "dt_new", dt
        endif
      endif
    else

      ! Timestep too small, smaller than the universally defined limit
      if (debug) then
        print *, "Timestep required for pos def is very small so we just take a single forward"
        print *, "Euler step before runge-kutta integration"
        print *, "dt dt_min, dt_min__posdef", dt, dt_min, dt_min__posdef
        print *, "y=", y
        print *, "dt__pd=", y/dydt0
        print *, "dy=", dt_min__posdef*dydt0
        print *, ""
      endif

      ! Euler step
      y = y + dt_min__posdef*dydt0
      t = t + dt_min__posdef

      if (debug) then
        print *, "y_new=", y
      endif

      ! TODO: There's got to be a better way than this, we want to make
      ! sure we get exactly to zero
      where (y(:) < tiny(y(1)))
        y(:) = 0.0_kreal
      endwhere
    endif

    if (.not. done) then
      k1 = 0.0
      k2 = 0.0
      k3 = 0.0
      k4 = 0.0
      k5 = 0.0

      k1 = dt*dydt(t,       y)
      k2 = dt*dydt(t+a2*dt, fix_y_isometric(y, b21*k1))
      k3 = dt*dydt(t+a3*dt, fix_y_isometric(y, b31*k1 + b32*k2))
      k4 = dt*dydt(t+a4*dt, fix_y_isometric(y, b41*k1 + b42*k2 + b43*k3))
      k5 = dt*dydt(t+a5*dt, fix_y_isometric(y, b51*k1 + b52*k2 + b53*k3 + b54*k4))
      !k6 = dt*dydt(t+a6*dt, fix_y_isometric(y, b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5))

      y_n1 = fix_y_isometric(y, c1_1*k1 + c2_1*k2 + c3_1*k3 + c4_1*k4 + c5_1*k5 )
      y_n2 = fix_y_isometric(y, c1_2*k1 + c2_2*k2 + c3_2*k3 + c4_2*k4 + c5_2*k5 )

      abs_err = abs(y_n1 - y_n2)
      !abs_err = abs(1./6.*(k4 - k5))

      ! TODO: make abs_tol and rel_tol vectors
      max_total_err = (abs_tol + rel_tol*abs(y_n2)) * dt
      s = 0.84*(minval(max_total_err/abs_err, abs_err > 0.0))**0.25

      ! Check to see if any solution is negative
      if (any(y_n2 < 0.0)) then
        !msg = "Solution became negative"
        if (debug) then
          print *, "=> solution became negative"
          print *, "dt s", dt, s
          print *, "y", y
          print *, "y_n2", y_n2
          print *, "max_tot_err", max_total_err
          print *, "abs_err", abs_err
        endif

        ! Check if s is large
        ! Large s suggests the current or a bigger timestep can be taken
        if (s > 1.0) then
          msg = "s was huge"
          done = .true.
        endif

        ! Finally, mass scale to deal with negative q species
        call mass_scale(y_n2, msg)
      endif

      if (debug) then
        print *, "max_tot_err=", max_total_err
        print *, "abs_err    =", abs_err
        !print *, "s_all=", max_total_err/abs_err
        !print *, "s=", s
        !print *, "dx=", dx
        !print *, abs_err < max_total_err

        !s = 0.84*(rel_tol*dx/max_abs_err)**0.25
        print *, ":: scaling by s", s
      endif

      ! Check if any solution is not a number...
      if (any(isnan(y_n1)) .or. any(isnan(y_n2))) then
        if (debug) then
          print *, "Solution became nan"
          print *, "s dt", s, dt
          print *, "y=", y
        endif
        if (s > 1.0) then
          s = 0.1
        endif
        !msg = "solution became nan"
        !done = .true.
      ! Success case
      else
        if (all(abs_err < max_total_err)) then
          ! Clear the error message, in case it was set earlier
          msg = " "
          done = .true.
          y = y_n2
          t = t + dt
        endif
      endif

      dt = dt*s

      ! Recursive call if integration isn't complete yet
      if (.not. done) then
        if (m > max_steps) then
          msg = "Didn't converge"
          print *, "last s=", s
          y = y_n2
        else
          !if (s >= 1.0) then
          !print *, "s=", s
          !print *, "Warning: Incorrect scaling, timestep is growing."
          !!s = 0.1
          !endif
          call rkf34_original(y, t, dt, msg, m+1)
        endif
      endif
    endif

  end subroutine rkf34_original
end module integrator_rkf
