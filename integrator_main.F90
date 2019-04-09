module integrator_main
  use microphysics_constants, only: kreal
  use microphysics_constants, only: abs_tol => integration_abs_tol, rel_tol => integration_rel_tol
  !use microphysics_register, only: n_variables
  use integrator_rkf, only: integrate_with_message
  !#ifdef MPI
  !use mpi !TODO: What is mpi
  !#endif

  implicit none

  ! interface
  !   subroutine constraint_based_integrator(y, t0, t_end, msg_out, m_total, n)
  !     integer, intent(in) :: n
  !     real(8), intent(inout), dimension(n) :: y
  !     real(8), intent(in) :: t_end, t0
  !     integer, intent(inout) :: m_total
  !     character(len=100), intent(inout), optional :: msg_out
  !   end subroutine
  ! end interface

  !> Will be assigned to one of the "integration helpers" which gaurantee
  !> either isometric or isobaric integration
  !procedure(constraint_based_integrator), pointer :: integrate_with_constraint => null()

contains
  !> Public subroutine that will be called by ATHAM/python-wrapper etc.
  subroutine integrate(y, t0, t_end, msg_out)
    integer, parameter :: n_variables = 7
    real(8), intent(inout), dimension(n_variables) :: y
    real(8), intent(in) :: t_end, t0
    character(len=100), intent(inout), optional :: msg_out
    character(len=100) :: msg
    real(8) :: y0(size(y))
    integer :: mpi_rank, ierror
    integer :: m_total

    ! Copy the initial state for future reference
    y0(:) = y

    msg = " "
    m_total = 0

    call integrate_with_message(y, t0, t_end, msg, m_total)

    if (present(msg_out)) then
      !TODO: when calling from ATHAM this "optional" value is set although
      !it it shouldn't be. Can't use it right now
      !msg_out = msg
    else
      if (msg(1:1) /= " ") then
        !#ifdef MPI
        !call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierror)
        !#endif
        print *, "==============================="
        print *, "integration failed"
        print *, mpi_rank, ":", y0
        print *, msg
        call exit(-2)
      endif
    endif
  end subroutine integrate
end module integrator_main
