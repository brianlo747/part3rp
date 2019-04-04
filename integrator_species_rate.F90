module integrator_species_rate
  use mphys_with_ice

  implicit none

  contains
    pure function dydt_mphys_with_ice(t, y, c_m) result(dydt)
       !use microphysics_register, only: idx_temp, idx_pressure
       !use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain, idx_cice, idx_graupel
       use microphysics_constants, only: T0, L_v => L_cond
       use mphys_with_ice

       real(kreal), dimension(:), intent(in) :: y
       real(kreal), intent(in) :: t, c_m
       real(kreal) :: dydt(size(y))

       real(kreal) :: ql, qv, qg, qr, qd, qi, qh
       real(kreal) :: rho, rho_g
       real(kreal) :: dqldt_condevap, dqrdt_condevap, dqhdt_condevap
       real(kreal) :: dqidt_sublidep, dqhdt_sublidep
       real(kreal) :: dqrdt_autoconv, dqhdt_autoconv
       real(kreal) :: dqrdt_accre_rc, dqhdt_accre_hi, dqldt_accre_chr
       real(kreal) :: dqhdt_accre_hiri, dqhdt_accre_hirr, dqhdt_accre_hr
       real(kreal) :: dqrdt_melt_rh, dqldt_melt_lc
       real(kreal) :: dqidt_freeze , dqhdt_freeze
       real(kreal) :: temp, p

       ! OBS: it's important to make sure that return variable is initiated to
       ! zero
       dydt = 0.0_kreal

       temp = y(1)
       p = y(2)

       ! pick out specific concentrations from state vectors
       ql = y(3)           !idx_cwater
       qr = y(4)             !idx_rain
       qv = y(5)     !idx_water_vapour
       qi = y(6)             !idx_cice
       qh = y(7)          !idx_graupel
       qd = 1.0_kreal - ql - qr - qv - qi - qh
       qg = qv + qd

       ! compute gas and mixture density using equation of state
       rho = rho_f(qd, qv, ql, qr, qi, qh, p, T)
       rho_g = rho_f(qd, qv, 0.0_kreal, 0.0_kreal, 0.0_kreal, 0.0_kreal, p, T)

       ! compute time derivatives for each process TODO: Functions for each function

       dqldt_condevap = dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=temp, p=p)
       dqrdt_condevap = dqr_dt__condensation_evaporation(qv=qv, qr=qr, rho=rho, T=temp, p=p)
       dqhdt_condevap = dqh_dt__condensation_evaporation(qv=qv, qh=qh, rho=rho, T=temp, p=p)
       dqidt_sublidep = dqi_dt__sublimation_deposition(qi=qi, rho=rho, T=temp, p=p)
       dqhdt_sublidep = dqh_dt__sublimation_evaporation(qg=qg, qv=qv, qh=qh, rho=rho, T=temp, p=p)
       dqrdt_autoconv = dqr_dt__autoconversion(ql=ql, qg=qg, rho_g=rho_g, T=temp)
       dqhdt_autoconv = dqh_dt__autoconversion_ice_graupel(qi=qi, qg=qg, rho_g=rho_g, T=temp)
       dqrdt_accre_rc = dqr_dt__accretion_cloud_rain(ql=ql, rho_g=rho_g, qr=qr)
       dqhdt_accre_hi = dqh_dt__accretion_ice_graupel(qi=qi, rho_g=rho_g, qh=qh, T=temp)
       dqldt_accre_chr= dqh_dt__accretion_cloud_graupel_rain(ql, rho_g, rho, qh)
       dqhdt_accre_hiri= dqr_dt__accretion_ice_rain_graupel_i(qi, rho_g, rho, qr)
       dqhdt_accre_hirr= dqr_dt__accretion_ice_rain_graupel_r(qi, ql, rho, rho_g, qr)
       dqhdt_accre_hr = dqr_dt__accretion_graupel(qg=qg, rho_g=rho_g, qv=qv, qh=qh, rho=rho, T=temp, p=p, qr=qr)
       dqrdt_melt_rh  = dqr_dt__melting_graupel(qg=qg, rho_g=rho_g, qv=qv, qh=qh, rho=rho, T=temp, p=p, qr=qr, ql=ql)
       dqldt_melt_lc  = dqr_dt__melting_ice(qg=qg, rho_g=rho_g, qv=qv, qh=qh, rho=rho, T=temp, p=p, qr=qr, ql=ql)
       dqidt_freeze   = dqi_dt__freezing_graupel(qh=qh, rho=rho, T=temp)
       dqhdt_freeze   = dqh_dt__freezing_ice(ql=ql, rho=rho, T=temp)



       ! combine to create time derivatives for species TODO: Include many more terms!

       if (temp > T0) then
         dydt(3) =  dqldt_condevap - dqrdt_autoconv - dqrdt_accre_rc - dqidt_freeze
         dydt(4) =                   dqrdt_autoconv + dqrdt_accre_rc + dqrdt_condevap
         dydt(5) = -dqldt_condevap                                - dqrdt_condevap
         dydt(6) = 0.0
         dydt(7) = 0.0

         dydt(1) = L_v/c_m*dqldt_condevap
       else
         dydt(3) =  dqldt_condevap - dqrdt_autoconv - dqrdt_accre_rc - dqidt_freeze
         dydt(4) =                   dqrdt_autoconv + dqrdt_accre_rc + dqrdt_condevap
         dydt(5) = -dqldt_condevap                                - dqrdt_condevap
         dydt(6) = 0.0
         dydt(7) = 0.0

         dydt(1) = L_v/c_m*dqldt_condevap
       endif

    end function

end module integrator_species_rate
