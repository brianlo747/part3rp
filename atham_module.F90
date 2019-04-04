!-*- F90 -*- so emacs thinks this is an f90 file
module atham_module
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    Dec 2004                                                  !
  ! version: 0.1                                                       !
  !                                                                    !
  ! definition of model data                                           !
  ! collection utility routines                                        !
  !--------------------------------------------------------------------!
  use precision
  use phys_constants, only: r0

  implicit none
  private
  public :: volcano_setup, coignimbrite_setup, procsconfig_setup,       &
            nuclear_setup, convection_setup,                            &
			config_surf_from_file
  public :: rh_init_water
  public :: cylindric_geometry, cyclic_boundary, enforce_sym_boundary
  public :: spmax, spstep, spfrac, sponge_in_vertical, nudge_sponge
  public :: netcdf, grid2geo, timezone,                                 &
            sim_d_YYYY, sim_d_MM, sim_d_DD, sim_t_hh, sim_t_mm, sim_t_ss, &
            current_year, current_month, current_jday, current_lst
  public :: nrep, periodt, dtmin, dtmax, cfllim, cpumax, spinup, cfl_diflim
  public :: ntx, nty, nx, ny, nz, npx, npy, nxl, nxr, nyv, nyh,         &
            icenter, jcenter, kcenter,nxtrans, nytrans, nztrans,        &
            nxtrans_boundary, nytrans_boundary, nztrans_boundary,       &
            xstart, xtotal, dxzoom, ystart, ytotal, dyzoom,             &
            zstart, ztotal, dzzoom, xg, yg, iindex, jindex,             &
            x, y, z, xv, yv, zv, scaler, scalel, scalerv, scalelv,      &
            iflgs, ifeld, landseamask, read_zgrid
  public :: deglat, deglon, corsin, corcos, corsin_cyl
  public :: include_boundary, no_boundary
  public :: mpi, mpi_initialized, nprocs, myrank,                       &
            iprocxr, iprocxl, iprocyh, iprocyv,                         &
            my_cart, mycoor, my_planexz, my_planeyz,                    &
            my_rowx, my_rowy, myrank_rowx, myrank_rowy,                 &
            scatter_displs_x, scatter_displs_y, scatter_x, scatter_y,   &
             gather_displs_x,  gather_displs_y,  gather_x,  gather_y
  public :: ntgas, ntrac, ntpas
  public :: unew, vnew, wnew, pnew, tetnew, tetflx, wtet, pflx,         &
            turbhor, turbver, turblen, horflx, verflx, scaflx,          &
            turbheight, turbrough,                                      &
            tgasnew, tgasflx, tracnew, tracflx, wtrac,                  &
            tpasnew, tpasflx, wtpas, wtgas, tgas, trac, tpas,           &
            tracdep, tpasdep, p0, temp, tss, relh, ronull,              &
            ugeos, vgeos, viskos, freewy, ozone, spech,                 &
            density, tempnew, vratnew, rgasnew, cpgasnew, cptot,         &
            wflux, wgasflux,                                            &
            cptgas, cvtgas, cptrac, rhotrac, radtrac, gas_tpas,         &
            c0, alpha_tet, alpha_tgas, alpha_trac, alpha_tpas,          &
            arith_mean, restu, restv, restw
  public :: restart_file, movie_file, picture_file,                     &
            movie_desfile, picture_desfile, mask_file, mask_desfile,    &
            movie_datfile, picture_datfile, mask_datfile
  public :: var_tgas_picture, var_tgas_movie, des_tgas_picture, des_tgas_movie, &
            var_trac_picture, var_trac_movie, des_trac_picture, des_trac_movie, &
            var_tpas_picture, var_tpas_movie, des_tpas_picture, des_tpas_movie
  public :: ntgas_movie, ntrac_movie, nwtrac_movie, ntpas_movie,        &
            ntgas_picture, ntrac_picture, nwtrac_picture, ntpas_picture
  public :: itgas_movie,itrac_movie,iwtrac_movie,itpas_movie,           &
            itgas_picture,itrac_picture,iwtrac_picture,itpas_picture
  public :: des_unew, des_vnew, des_wnew, des_pnew, des_tetnew,         &
            des_tempnew, des_density, des_turbhor, des_turbver, des_turblen
  public :: longname_gas, longname_trac, longname_pas, longname_wtrac,  &
            varname_gas, varname_trac, varname_pas, varname_wtrac,      &
            dyn_movie

  public :: tracer_allocate, atham_stop, atham_filenames, makename

  !--------------------------------------------------------------------!
  ! configuratione control flags:                                      !
  !                                                                    !
  ! volcano        volcanic eruption plume                             !
  ! procsconfig    process configuration system                        !
  !--------------------------------------------------------------------!
  logical, save :: volcano_setup      =.false.,                         &
                   coignimbrite_setup =.false.,                         &
                   convection_setup   =.false.,                         &
                   nuclear_setup      =.false.,                         &
                   procsconfig_setup  =.false.,                         &
                   config_surf_from_file = .false.
  !--------------------------------------------------------------------!
  ! initialization for water vapor                                     !
  ! rh_init_water  treat initial relative humidity in INPUT_profile    !
  !                always as RH over water if true                     !
  !                as RH over water or ice if false                    !
  !--------------------------------------------------------------------!
  logical, save :: rh_init_water=.true.
  !--------------------------------------------------------------------!
  ! run control flags:                                                 !
  !                                                                    !
  ! cylindric_geometry    2d version in cylindrical geometry           !
  ! cyclic_boundary       cyclic lateral boundaries                    !
  ! enforce_sym_boundary  enforce identical spatial resolution         !
  !                       near lateral boundaries                      !
  !--------------------------------------------------------------------!
  logical, save :: cylindric_geometry   =.false.,                      &
                   cyclic_boundary      =.false.,                      &
                   enforce_sym_boundary =.false.
  !--------------------------------------------------------------------!
  ! sponge layer                                                       !
  !--------------------------------------------------------------------!
  real(kreal), save :: spmax, spstep, spfrac
  logical, save :: sponge_in_vertical=.false.
  logical, save :: nudge_sponge
  !--------------------------------------------------------------------!
  ! output control flags:                                              !
  !                                                                    !
  ! netcdf                netcdf output                                !
  !                       (single precision binary if false)           !
  ! grid2geo              linear approximation of xyz grid in          !
  !                       geographical (lon/lat/lev) frame             !
  ! deglat/lon_start      lat/lon coordinates of grid origin (grid2geo)!
  ! sim_d_*               simulation onset date YYYY-MM-DD   (grid2geo)!
  ! sim_t_*               simulation onset time hh:mm:ss     (grid2geo)!
  ! current_year/month    current year and month                       !
  ! current_jday          current Julian day (assuming 365 days/year   !
  ! current_lst           current Local Solar Time                     !
  !--------------------------------------------------------------------!
  logical, save :: netcdf=.false.
  logical, save :: grid2geo=.true.
  integer(kint), save :: sim_d_YYYY
  integer(kint), save :: sim_d_MM
  integer(kint), save :: sim_d_DD
  integer(kint), save :: sim_t_hh
  integer(kint), save :: sim_t_mm
  integer(kint), save :: sim_t_ss
  real(kreal), save   :: timezone
  integer(kint), save :: current_year
  integer(kint), save :: current_month
  integer(kint), save :: current_jday
  real(kreal), save   :: current_lst
  !--------------------------------------------------------------------!
  ! run time control variables:                                        !
  !                                                                    !
  ! nrep                    number of periods for i/o                  !
  ! periodt      [sec]      length of period between i/o               !
  ! dtmin        [sec]      smallest allowed time step                 !
  ! dtmax        [sec]      largest allowed time step                  !
  ! cfllim       [1]        CFL limit for advection                    !
  ! cfl_diflim   [1]        CFL limit for vert diffusion coef          !
  ! cpumax       [sec]      maximal allowed cpu-time per run           !
  ! spinup       [sec]      time for spinup                            !
  !--------------------------------------------------------------------!
  integer(kint), save :: nrep
  real(kreal), save   :: periodt,dtmin,dtmax,cfllim,cpumax,spinup=r0,   &
                         cfl_diflim
  !--------------------------------------------------------------------!
  ! model geometry:                                                    !
  !                                                                    !
  ! ntx,nty,nz         total number of grid points in each direction   !
  ! nx,ny,nz           number of grid points per processor             !
  !                    nty/ny=4 in 2d versions                         !
  ! nxl,nxr            index range in x-direction                      !
  ! nyv,nyh            index range in y-direction                      !
  !                                                                    !
  ! npx,npy      number of processors in x/y direction (MPI version)   !
  !                                                                    !
  ! xstart       [m]   start of x-domain                               !
  ! xtotal       [m]   total x-domain                                  !
  ! dxzoom       [m]   spatial resolution in center of zoom in x       !
  ! icenter      [1]   location of zoom in x in grid space             !
  ! nxtrans      [1]   width of zoom in x in grid space                !
  ! nxtrans_boundary [1]   size of transition zone at lateral boundary !
  !                                                                    !
  ! ystart       [m]   start of y-domain                               !
  ! ytotal       [m]   total y-domain                                  !
  ! dyzoom       [m]   spatial resolution in center of zoom in y       !
  ! jcenter      [1]   location of zoom in y in grid space             !
  ! nytrans      [1]   width of zoom in y in grid space                !
  ! nytrans_boundary [1]   size of transition zone at lateral boundary !
  !                                                                    !
  ! zstart       [m]   start of z-domain                               !
  ! ztotal       [m]   total x-domain                                  !
  ! dzzoom       [m]   spatial resolution in center of zoom in z       !
  ! kcenter      [1]   location of zoom in z in grid space             !
  ! nztrans      [1]   width of zoom in z in grid space                !
  ! nztrans_boundary [1]   size of transition zone at top boundary     !
  !                                                                    !
  ! xg, yg       [m]   global coordinate for scalars                   !
  ! x, y ,z      [m]   local coordinate for scalars                    !
  ! xv, yv, zv   [m]   local coordinate for vectors                    !
  ! scaler/v, scalel/v [1]  scaling factors for cylindric coordinate   !
  ! iflgs, ifeld [1]   mask fields for topography                      !
  ! landseamask  [1]   land (=1) sea (=0) mask                         !
  !                                                                    !
  ! read_zgrid   [ ]   switch to prescribe the z-Grid                  !
  !--------------------------------------------------------------------!
  integer(kint), save :: ntx, nty, nx, ny, nz, npx, npy
  integer(kint), save :: nxl, nxr, nyv, nyh
  integer(kint), save :: icenter, jcenter, kcenter,nxtrans, nytrans, nztrans, &
                         nxtrans_boundary, nytrans_boundary, nztrans_boundary
  integer(kint), save :: iindex, jindex
  real(kreal), save :: xstart, xtotal, dxzoom,                          &
                       ystart, ytotal, dyzoom,                          &
                       zstart, ztotal, dzzoom
  real(kreal), allocatable, save, dimension(:) :: xg, yg
  real(kreal), allocatable, save, dimension(:) ::                       &
       x, y, z, xv, yv, zv, scaler, scalel, scalerv, scalelv
  integer(kint), allocatable, save, dimension(:,:,:) :: iflgs
  integer(kint), allocatable, save, dimension(:,:)   :: ifeld, landseamask
  logical, save :: read_zgrid
  !--------------------------------------------------------------------!
  ! deglat   [degree]   latitudinal location of the grid               !
  ! deglon   [degree]   longitude   location of the grid               !
  ! corsin   [1]        sine term for coriolis force                   !
  ! corcos   [1]        cosine term for coriolis force                 !
  !--------------------------------------------------------------------!
  real(kreal), save :: deglat, deglon
  real(kreal), allocatable, save, dimension(:,:) :: corsin, corcos, corsin_cyl
  !--------------------------------------------------------------------!
  ! flags for boundary conditions                                      !
  !--------------------------------------------------------------------!
  logical, save :: include_boundary=.true.,                             &
                   no_boundary=.false.
  !--------------------------------------------------------------------!
  ! MPI related quantities                                             !
  !                                                                    !
  ! mpi                   logical that flags use of MPI                !
  ! nprocs                number of processors                         !
  ! myrank, myrank_rowx/y rank / processor number                      !
  ! my_cart, my_rowx/y    communicator handle                          !
  ! iprocxr/xl/yh/yv      processor number of neighbor in x/y direction!
  ! mycoor                coordinate in 2d domain decomposition        !
  ! my_planexz/yz         derived data types for communication         !
  ! scatter_displs_x/y    displacements for scatterv operation         !
  ! scatter_x/y           send counts for scatterv operation           !
  ! gather_displs_x/y     displacements for gatherv operation          !
  ! gather_x/y            send counts for gatherv operation            !
  !--------------------------------------------------------------------!
  logical, save :: mpi, mpi_initialized=.false.
  integer(kint), save :: nprocs, myrank, myrank_rowx, myrank_rowy,      &
       my_cart, iprocxr, iprocxl, iprocyh, iprocyv,                     &
       mycoor(2), my_rowx, my_rowy, my_planexz, my_planeyz
  integer(kint), allocatable, save, dimension(:) ::                     &
       scatter_displs_x, scatter_displs_y, scatter_x, scatter_y,        &
        gather_displs_x,  gather_displs_y,  gather_x,  gather_y
  !--------------------------------------------------------------------!
  ! number of tracers                                                  !
  ! ntgas        [1]   number of gaseous tracers                       !
  ! ntrac        [1]   number of incompressible (active tracers)       !
  ! ntpas        [1]   number of passive tracers                       !
  !--------------------------------------------------------------------!
  integer(kint), save :: ntgas, ntrac, ntpas
  !--------------------------------------------------------------------!
  ! model data                                                         !
  ! u/v/wnew      [m/s]    components of wind vector (c-grid)          !
  ! pnew          [Pa]     pressure anomaly                            !
  ! tetnew        [K]      potential temperature                       !
  ! turbhor/ver   [m2/s2]  vertical/horizontal turbulent kinetic energy!
  ! turblen       [m]      turbulent length scale                      !
  ! tgasnew       [kg/kg]  gaseous tracer                              !
  ! tracnew       [kg/kg]  incompressible tracer                       !
  ! tpasnew       [#/kg]   number concentration (passive tracer)       !
  !                                                                    !
  ! density       [kg/m3]  density                                     !
  ! tempnew       [K]      in situ temperature                         !
  ! vratnew       [1]      ratio if specific heats cp/cv               !
  ! rgasnew       [J/kg/K] gas constant cp-cv                          !
  ! cptot         [J/kg/K] specific heat at constant pressure          !
  !                                                                    !
  ! wflux         [m/s]    velocity for surface forcing                !
  ! wgasflux      [m/s]    velocity anomaly for gas for surface forcing!
  !--------------------------------------------------------------------!
  real(kreal), allocatable, save, dimension(:,:,:) ::                   &
       unew, vnew, wnew, pnew, tetnew, tetflx, wtet, pflx,              &
       restu, restv, restw
  real(kreal),allocatable, save, dimension(:,:,:) ::                    &
       turbhor, turbver, turblen, horflx, verflx, scaflx
  real(kreal),allocatable, save, dimension(:,:) ::                      &
       turbheight, turbrough
  real(kreal), allocatable, target, save, dimension(:,:,:,:) ::         &
       tgasnew, tgasflx,                                                &
       tracnew, tracflx, wtrac,                                         &
       tpasnew, tpasflx, wtpas
  real(kreal), allocatable, target, save, dimension(:,:,:) ::           &
       wtgas
  real(kreal), allocatable, target, save, dimension(:,:) ::             &
       tgas, trac, tpas
  real(kreal), allocatable, target, save, dimension(:,:,:) ::           &
       tracdep, tpasdep
  real(kreal), allocatable, save, dimension(:) ::                       &
       p0, temp, tss, relh, ronull, ugeos, vgeos, viskos, freewy, ozone, spech
  real(kreal), allocatable, save, dimension(:,:,:) ::                   &
       density,tempnew,vratnew,rgasnew,cpgasnew,cptot
  real(kreal), allocatable, save, dimension(:,:) ::                     &
       wflux, wgasflux

  real(kreal), save :: c0=0.32_kreal, alpha_tet=2._kreal
  real(kreal), target, save, dimension(99) ::                           &
       cptgas, cvtgas, cptrac, rhotrac, radtrac,                        &
       alpha_tgas, alpha_trac, alpha_tpas
  logical, save, dimension(99) ::                                       &
       gas_tpas=.true.
  logical, save :: arith_mean=.true.
  !--------------------------------------------------------------------!
  ! output                                                             !
  !--------------------------------------------------------------------!
  character(len=25), save :: restart_file, movie_file, picture_file,    &
       movie_desfile, picture_desfile, mask_file, mask_desfile,         &
       movie_datfile, picture_datfile, mask_datfile
!!$  character(len=25), save ::                                            &
!!$       des_unew    ='x_wind                    ',                       &
!!$       des_vnew    ='y_wind                    ',                       &
!!$       des_wnew    ='upward_air_velocity       ',                       &
!!$       des_pnew    ='air_pressure_anomaly      ',                       &
!!$       des_tetnew  ='air__potential_temperature',                       &
!!$       des_tempnew ='air_temperature           ',                       &
!!$       des_density ='air_density               ',                       &
!!$       des_turbhor ='hor_turb_energy           ',                       &
!!$       des_turbver ='ver_turb_energy           ',                       &
!!$       des_turblen ='turb_length_scale         '
  character(len=25), save ::                                            &
       des_unew    ='u_wind                   ',                        &
       des_vnew    ='v_wind                   ',                        &
       des_wnew    ='w_wind                   ',                        &
       des_pnew    ='pressure_anomaly         ',                        &
       des_tetnew  ='potential_temperature    ',                        &
       des_tempnew ='in_situ_temperature      ',                        &
       des_density ='bulk_density             ',                        &
       des_turbhor ='hor_turb_kinetic_energy  ',                        &
       des_turbver ='ver_turb_kinetic_energy  ',                        &
       des_turblen ='turb_length_scale        '

  character(len=25), save, dimension(99) ::                             &
       var_tgas_picture,var_tgas_movie,des_tgas_picture,des_tgas_movie, &
       var_trac_picture,var_trac_movie,des_trac_picture,des_trac_movie, &
       var_tpas_picture,var_tpas_movie,des_tpas_picture,des_tpas_movie

  integer(kint), save ::                                                &
       ntgas_movie, ntrac_movie, nwtrac_movie, ntpas_movie,             &
       ntgas_picture, ntrac_picture, nwtrac_picture, ntpas_picture
  integer(kint), save, dimension(99) ::                                 &
       itgas_movie, itrac_movie, iwtrac_movie, itpas_movie,             &
       itgas_picture, itrac_picture, iwtrac_picture, itpas_picture

  logical , save :: dyn_movie=.false.

  !--------------------------------------------------------------------!
  ! description names for tracer output (netcdf)                       !
  !--------------------------------------------------------------------!
  character(len=50), dimension(99), save :: longname_gas, longname_trac, &
       longname_pas, longname_wtrac
  character(len=10), dimension(99), save :: varname_gas, varname_trac, &
       varname_pas, varname_wtrac

contains
  !=====================================================================
  subroutine tracer_allocate()
    !------------------------------------------------------------------!
    ! allocate and initialize tracer data for process modules          !
    !------------------------------------------------------------------!

    allocate(tgasnew(nx,ny,nz,ntgas),tgasflx(nx,ny,nz,ntgas),           &
             wtgas(nx,ny,nz),tgas(nz,ntgas))
    allocate(tracnew(nx,ny,nz,ntrac),tracflx(nx,ny,nz,ntrac),           &
             wtrac(nx,ny,nz,ntrac),trac(nz,ntrac),tracdep(nx,ny,ntrac))
    allocate(tpasnew(nx,ny,nz,ntpas),tpasflx(nx,ny,nz,ntpas),           &
             wtpas(nx,ny,nz,ntpas),tpas(nz,ntpas),tpasdep(nx,ny,ntpas))

    tgasnew(:,:,:,:)=r0
    tgasflx(:,:,:,:)=r0
    wtgas(:,:,:)=r0
    tgas(:,:)=r0

    tracnew(:,:,:,:)=r0
    tracflx(:,:,:,:)=r0
    wtrac(:,:,:,:)=r0
    trac(:,:)=r0
    tracdep(:,:,:)=r0

    tpasnew(:,:,:,:)=r0
    tpasflx(:,:,:,:)=r0
    wtpas(:,:,:,:)=r0
    tpas(:,:)=r0
    tracdep(:,:,:)=r0

  end subroutine tracer_allocate
  !====================================================================!
  subroutine atham_stop(message)
    !------------------------------------------------------------------!
    ! print message, stop prgram execution                             !
    !------------------------------------------------------------------!
    character(len=*), intent(in) :: message
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
#ifdef MPI
    include 'mpif.h'
    integer(kint) :: ierror
#endif
    print 100, trim(message)
#ifdef MPI
    call mpi_abort(mpi_comm_world,ierror)
#endif
    stop

100 format(/,a50,/)

  end subroutine atham_stop
  !=====================================================================
  subroutine atham_filenames()
    !------------------------------------------------------------------!
    ! build filename for restart and output                            !
    !------------------------------------------------------------------!
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    character(len=25), save ::                                          &
         restart_prefix  ='output/ATDAT             ',                  &
         gradsmov_prefix ='gradsmov                 ',                  &
         netcdfmov_prefix='netCDF_MOV               ',                  &
         gradspic_prefix ='gradspic                 ',                  &
         netcdfpic_prefix='netCDF_PIC               ',                  &
         gradsmask_prefix='mask                     '
    integer(kint) :: nchar
    !------------------------------------------------------------------!
    ! if mpi use coordinates in file names                             !
    !------------------------------------------------------------------!
    if (mpi) then
       call makename(restart_prefix,mycoor,restart_file)
       if (netcdf) then
          movie_file=netcdfmov_prefix
          picture_file=netcdfpic_prefix
!!$          call makename(netcdfmov_prefix,mycoor,movie_file)
!!$          call makename(netcdfpic_prefix,mycoor,picture_file)
       else
          call makename(gradsmov_prefix,mycoor,movie_file)
          call makename(gradspic_prefix,mycoor,picture_file)
          call makename(gradsmask_prefix,mycoor,mask_file)
       endif
    else
       restart_file=restart_prefix
       if (netcdf) then
          movie_file=netcdfmov_prefix
          picture_file=netcdfpic_prefix
       else
          movie_file=gradsmov_prefix
          picture_file=gradspic_prefix
          mask_file=gradsmask_prefix
       endif
    endif
    !------------------------------------------------------------------!
    ! add .dat to movie and picture file                               !
    !------------------------------------------------------------------!
    movie_file=adjustl(movie_file)
    nchar=len_trim(movie_file)
    if (netcdf) then
       movie_desfile='output/'//movie_file(1:nchar)//'.des'
       movie_datfile=movie_file(1:nchar)//'.nc'
       movie_file='output/'//movie_file(1:nchar)//'.nc'
    else
       movie_desfile='output/'//movie_file(1:nchar)//'.des'
       movie_datfile=movie_file(1:nchar)//'.dat'
       movie_file='output/'//movie_file(1:nchar)//'.dat'
    endif
    picture_file=adjustl(picture_file)
    nchar=len_trim(picture_file)
    if (netcdf) then
       picture_desfile='output/'//picture_file(1:nchar)//'.des'
       picture_datfile=picture_file(1:nchar)//'.nc'
       picture_file='output/'//picture_file(1:nchar)//'.nc'
    else
       picture_desfile='output/'//picture_file(1:nchar)//'.des'
       picture_datfile=picture_file(1:nchar)//'.dat'
       picture_file='output/'//picture_file(1:nchar)//'.dat'

       nchar=len_trim(mask_file)
       mask_desfile='output/'//mask_file(1:nchar)//'.des'
       mask_datfile=mask_file(1:nchar)//'.dat'
       mask_file='output/'//mask_file(1:nchar)//'.dat'
    endif
  end subroutine atham_filenames
  !====================================================================!
  subroutine makename(prefix,icoor,name)
    !------------------------------------------------------------------!
    ! create character name=prefix_(icoor(1)+1)_(icoor(2)+1)           !
    !------------------------------------------------------------------!
    character(len=*), intent(in) :: prefix
    integer(kint), dimension(2), intent(in) :: icoor
    character(len=*), intent(out) :: name
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    character(len=25) numx,numy
    integer(kint) :: lenx,leny,lenf

    write(numx,'(i5)') icoor(1)+1
    numx=adjustl(numx)
    lenx=len_trim(numx)

    write(numy,'(i5)') icoor(2)+1
    numy=adjustl(numy)
    leny=len_trim(numy)

    name=adjustl(prefix)
    lenf=len_trim(name)

    name=name(1:lenf)//'_'//numx(1:lenx)//'_'//numy(1:leny)

  end subroutine makename
  !====================================================================!
end module atham_module
