module integrator_config

  use microphysics_constants, only: kreal

  implicit none

  ! ! RUN 1 (Default w/ precip)
  ! ! CONFIG
  ! logical, parameter :: rain_precip = .true.
  ! logical, parameter :: graupel_precip = .true.
  ! logical, parameter :: ice_graupel_processes = .true.
  !
  ! ! UPDRAFT CONFIG
  ! real(kreal), parameter :: updraft_vel = 2.0_kreal
  ! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
  !
  ! ! PRECIP CONFIG
  ! real(kreal), parameter :: cloud_thickness = 10.0_kreal
  !
  ! ! CLOUD DROPLET NUMBER CONFIG
  ! real(kreal), parameter :: N_cloud = 2.0e8_kreal

  ! ! RUN 2 (Updraft vel 0.5)
  ! ! CONFIG
  ! logical, parameter :: rain_precip = .true.
  ! logical, parameter :: graupel_precip = .true.
  ! logical, parameter :: ice_graupel_processes = .true.
  !
  ! ! UPDRAFT CONFIG
  ! real(kreal), parameter :: updraft_vel = 0.5_kreal
  ! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
  !
  ! ! PRECIP CONFIG
  ! real(kreal), parameter :: cloud_thickness = 10.0_kreal
  !
  ! ! CLOUD DROPLET NUMBER CONFIG
  ! real(kreal), parameter :: N_cloud = 2.0e8_kreal

  ! ! RUN 3 (Updraft vel 10)
  ! ! CONFIG
  ! logical, parameter :: rain_precip = .true.
  ! logical, parameter :: graupel_precip = .true.
  ! logical, parameter :: ice_graupel_processes = .true.
  !
  ! ! UPDRAFT CONFIG
  ! real(kreal), parameter :: updraft_vel = 10.0_kreal
  ! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
  !
  ! ! PRECIP CONFIG
  ! real(kreal), parameter :: cloud_thickness = 10.0_kreal
  !
  ! ! CLOUD DROPLET NUMBER CONFIG
  ! real(kreal), parameter :: N_cloud = 2.0e8_kreal

  ! ! RUN 4 (CDNC 1.0e9)
  ! ! CONFIG
  ! logical, parameter :: rain_precip = .true.
  ! logical, parameter :: graupel_precip = .true.
  ! logical, parameter :: ice_graupel_processes = .true.
  !
  ! ! UPDRAFT CONFIG
  ! real(kreal), parameter :: updraft_vel = 2.0_kreal
  ! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
  !
  ! ! PRECIP CONFIG
  ! real(kreal), parameter :: cloud_thickness = 10.0_kreal
  !
  ! ! CLOUD DROPLET NUMBER CONFIG
  ! real(kreal), parameter :: N_cloud = 1.0e9_kreal

  ! RUN 5 (CDNC 1.0e11)
  ! CONFIG
  logical, parameter :: rain_precip = .true.
  logical, parameter :: graupel_precip = .true.
  logical, parameter :: ice_graupel_processes = .true.

  ! UPDRAFT CONFIG
  real(kreal), parameter :: updraft_vel = 2.0_kreal
  real(kreal), parameter :: lapse_rate = 9.8e-3_kreal

  ! PRECIP CONFIG
  real(kreal), parameter :: cloud_thickness = 10.0_kreal

  ! CLOUD DROPLET NUMBER CONFIG
  real(kreal), parameter :: N_cloud = 1.0e11_kreal

  ! ! RUN 6 (Default w/o ice)
  ! ! CONFIG
  ! logical, parameter :: rain_precip = .true.
  ! logical, parameter :: graupel_precip = .true.
  ! logical, parameter :: ice_graupel_processes = .false.
  !
  ! ! UPDRAFT CONFIG
  ! real(kreal), parameter :: updraft_vel = 2.0_kreal
  ! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
  !
  ! ! PRECIP CONFIG
  ! real(kreal), parameter :: cloud_thickness = 10.0_kreal
  !
  ! ! CLOUD DROPLET NUMBER CONFIG
  ! real(kreal), parameter :: N_cloud = 2.0e8_kreal


!!!!! NO PRECIP!


! ! RUN 7 (Default)
! ! CONFIG
! logical, parameter :: rain_precip = .false.
! logical, parameter :: graupel_precip = .false.
! logical, parameter :: ice_graupel_processes = .true.
!
! ! UPDRAFT CONFIG
! real(kreal), parameter :: updraft_vel = 2.0_kreal
! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
!
! ! PRECIP CONFIG
! real(kreal), parameter :: cloud_thickness = 10.0_kreal
!
! ! CLOUD DROPLET NUMBER CONFIG
! real(kreal), parameter :: N_cloud = 2.0e8_kreal

! ! RUN 8 (Updraft vel 0.5)
! ! CONFIG
! logical, parameter :: rain_precip = .false.
! logical, parameter :: graupel_precip = .false.
! logical, parameter :: ice_graupel_processes = .true.
!
! ! UPDRAFT CONFIG
! real(kreal), parameter :: updraft_vel = 0.5_kreal
! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
!
! ! PRECIP CONFIG
! real(kreal), parameter :: cloud_thickness = 10.0_kreal
!
! ! CLOUD DROPLET NUMBER CONFIG
! real(kreal), parameter :: N_cloud = 2.0e8_kreal

! ! RUN 9 (Updraft vel 10)
! ! CONFIG
! logical, parameter :: rain_precip = .false.
! logical, parameter :: graupel_precip = .false.
! logical, parameter :: ice_graupel_processes = .true.
!
! ! UPDRAFT CONFIG
! real(kreal), parameter :: updraft_vel = 10.0_kreal
! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
!
! ! PRECIP CONFIG
! real(kreal), parameter :: cloud_thickness = 10.0_kreal
!
! ! CLOUD DROPLET NUMBER CONFIG
! real(kreal), parameter :: N_cloud = 2.0e8_kreal

! ! RUN 10 (CDNC 1.0e9)
! ! CONFIG
! logical, parameter :: rain_precip = .false.
! logical, parameter :: graupel_precip = .false.
! logical, parameter :: ice_graupel_processes = .true.
!
! ! UPDRAFT CONFIG
! real(kreal), parameter :: updraft_vel = 2.0_kreal
! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
!
! ! PRECIP CONFIG
! real(kreal), parameter :: cloud_thickness = 10.0_kreal
!
! ! CLOUD DROPLET NUMBER CONFIG
! real(kreal), parameter :: N_cloud = 1.0e9_kreal

! ! RUN 11 (CDNC 1.0e11)
! ! CONFIG
! logical, parameter :: rain_precip = .false.
! logical, parameter :: graupel_precip = .false.
! logical, parameter :: ice_graupel_processes = .true.
!
! ! UPDRAFT CONFIG
! real(kreal), parameter :: updraft_vel = 2.0_kreal
! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
!
! ! PRECIP CONFIG
! real(kreal), parameter :: cloud_thickness = 10.0_kreal
!
! ! CLOUD DROPLET NUMBER CONFIG
! real(kreal), parameter :: N_cloud = 1.0e11_kreal

! ! RUN 12 (Default w/o ice)
! ! CONFIG
! logical, parameter :: rain_precip = .false.
! logical, parameter :: graupel_precip = .false.
! logical, parameter :: ice_graupel_processes = .false.
!
! ! UPDRAFT CONFIG
! real(kreal), parameter :: updraft_vel = 2.0_kreal
! real(kreal), parameter :: lapse_rate = 9.8e-3_kreal
!
! ! PRECIP CONFIG
! real(kreal), parameter :: cloud_thickness = 10.0_kreal
!
! ! CLOUD DROPLET NUMBER CONFIG
! real(kreal), parameter :: N_cloud = 2.0e8_kreal
end module integrator_config
