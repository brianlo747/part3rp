!-*- F90 -*- so emacs thinks this is an f90 file
module precision
!----------------------------------------------------------------------!
!                                                                      !
! Purpose:                                                             !
!       Define the precision for real and integer numbers              !
!       used throughout the model.                                     !
!                                                                      !
! Author:                                                              !
!       Michael Herzog (Michael.Herzog@noaa.gov)                       !
!                                                                      !
!       kreal, kint are the default kind values for real and integer   !
!----------------------------------------------------------------------!
  implicit none

  private ! this private is needed so mpif.h doesn't clash with
          ! later include "mpif.h"'s

#ifdef MPI
  include 'mpif.h'
  integer, public, parameter :: my_real=MPI_REAL8
#endif

  integer, public, parameter :: real2 = selected_real_kind(3),            &
                                real4 = selected_real_kind(6),            &
                                real8 = selected_real_kind(12),           &
                                kreal = real8
  integer, public, parameter :: int4 = selected_int_kind(7),              &
                                int8 = selected_int_kind(13),             &
                                kint = int4

end module precision
