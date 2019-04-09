module integrator_helpers
   !use microphysics_register, only: dydt_mphys => dydt_registered
   use integrator_rkf, only: integrate_with_message

   implicit none

   public integrate_isometric

   contains
      ! subroutine integrate_isometric(y, t0, t_end, msg_out, m_total, n)
      !    integer, intent(in) :: n
      !    real(8), intent(inout), dimension(n) :: y
      !    real(8), intent(in) :: t_end, t0
      !    integer, intent(inout) :: m_total
      !    character(len=100), intent(inout), optional :: msg_out
      !
      !    call integrate_with_message(dydt_isometric, y, increment_state_isometric, t0, t_end, msg_out, m_total)
      ! end subroutine integrate_isometric

      ! pure function dydt_isometric(t, y) result(dydt)
      !    use microphysics_common, only: cv_mixture
      !    use integrator_species_rate, only: dydt_mphys_with_ice
      !
      !    real(8), intent(in) :: t
      !    real(8), intent(in), dimension(:) :: y
      !    real(8) :: c_m, dydt(size(y))
      !
      !    c_m = cv_mixture(y)
      !
      !    dydt = dydt_mphys_with_ice(t, y, c_m)
      ! end function dydt_isometric

      !> Return a new state changed by increment `dy` keeping constant volume,
      !method is intended to guarantee self-consistency of new state `y_new`
      ! TODO: Implement mass-scaling
      
end module
