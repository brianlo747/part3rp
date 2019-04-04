!-*- F90 -*- so emacs thinks this is an f90 file
module process_data
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  !          Gunnar Luderer                                            !
  !          Alex Hoffmann                                             !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    July 2005 / April 2009                                    !
  ! version: 0.2                                                       !
  !                                                                    !
  ! collection pointers and data for process modules                   !
  !--------------------------------------------------------------------!
  use precision, only: kreal, kint
  !--------------------------------------------------------------------!
  ! tracer number and setup flags                                      !
  !--------------------------------------------------------------------!
  integer(kint), save ::                                                &
       ntgas_forcing,ntgas_micro,ntgas_chem,                            &
       ntrac_forcing,ntrac_micro,ntrac_chem,                            &
       nwtrac_forcing,nwtrac_micro,nwtrac_chem,                         &
       ntpas_forcing,ntpas_micro,ntpas_chem
  integer(kint), save ::                                                &
       ntgas_forc_movie,ntgas_forc_picture,                             &
       ntrac_forc_movie,ntrac_forc_picture,                             &
       nwtrac_forc_movie,nwtrac_forc_picture,                           &
       ntpas_forc_movie,ntpas_forc_picture            
  integer(kint), save, dimension(99) ::                                 &
       itgas_forc_movie,itgas_forc_picture,                             &
       itrac_forc_movie,itrac_forc_picture,                             &
       iwtrac_forc_movie,iwtrac_forc_picture,                           &
       itpas_forc_movie,itpas_forc_picture          
  integer(kint), save ::                                                &
       ntgas_mic_movie,ntgas_mic_picture,                               &
       ntrac_mic_movie,ntrac_mic_picture,                               &
       nwtrac_mic_movie,nwtrac_mic_picture,                             &
       ntpas_mic_movie,ntpas_mic_picture          
  integer(kint), dimension(99) ::                                       &
       itgas_mic_movie,itgas_mic_picture,                               &
       itrac_mic_movie,itrac_mic_picture,                               &
       iwtrac_mic_movie,iwtrac_mic_picture,                             &
       itpas_mic_movie,itpas_mic_picture

  logical, save ::                                                      &
       fireforcing_setup,volcforcing_setup,warmbubble_setup,            &
       kessler_setup,twomicrophys_setup,sbm_setup,chem_setup,           &
       scav_setup		   

  !--------------------------------------------------------------------!
  ! arrays, fluxes and data for bulk cloud microphysics                !
  ! kessler_microphysics.F90                                           !
  !--------------------------------------------------------------------!
  real(kreal), pointer, save, dimension(:,:,:) ::                       &
       wetnew,watcnew,watpnew,icenew,granew,                            &
       wetflx,watcflx,watpflx,iceflx,graflx,                            &
       wwatc,wwatp,wice,wgra,                                           &
       wawatc,wawatp,waice,wagra

  real(kreal), pointer, save :: radwatc,radice,rhowat,rhoice,cpwat,cpice

  !--------------------------------------------------------------------!
  ! arrays, fluxes and data for particles  (forcing)                   !
  !--------------------------------------------------------------------!
  real(kreal),pointer, save, dimension(:,:,:) ::                        &
       ashnew,stonew,intaer,ashflx,stoflx,intaerflx,                    &
       wint

  !--------------------------------------------------------------------!
  ! additional arrays, fluxes and data for twomoment microphysics      !
  !--------------------------------------------------------------------!
  real(kreal),pointer, save, dimension(:,:,:) ::                        &
       xnwatc,xnwatcflx,xnwatp,xnwatpflx,xnice,xniceflx,xngra,xngraflx, &
       wnwatc, wnwatp, wnice, wngra
  real(kreal),pointer, save, dimension(:,:,:) ::                        &
       awatc,awatp,aice,agra,awatcflx,awatpflx,aiceflx,agraflx
  real(kreal), pointer, save :: radwatp, radgra, radash, radsto,        &
       radashc, radstoc,                                                &
       rhowatp, rhogra, rhoash, rhosto, rhoashc, rhostoc,               &
       cpwatp, cpgra, cpash, cpsto, cpashc, cpstoc,                     &
       cpint, radint, rhoint
  !--------------------------------------------------------------------!
  ! additional passive tracer arrays and fluxes                        !
  !--------------------------------------------------------------------!
  real(kreal),pointer, save, dimension(:,:,:) ::                        &
       frednew, fredflx, wfred
  
  !--------------------------------------------------------------------!
  ! arrays and fluxes for chemistry                                    !
  !--------------------------------------------------------------------!
  real(kreal),allocatable, save, dimension(:,:,:):: o3new, o3flx, conew,&
       & coflx, co2enew, co2eflx, co2bnew, co2bflx, coenew, coeflx,&
       & cobnew, cobflx,&
       & hchoenew, hchoeflx, hchobnew, hchobflx, noenew, noeflx,&
       & nobnew, nobflx, no2enew, no2eflx, no2bnew, no2bflx,&
       & hno3enew, hno3eflx, hno3bnew, hno3bflx,&
       & h2o2new, h2o2flx, ho2new, ho2flx, o1dnew, o1dflx, &
       & ohnew, ohflx, ch4enew, ch4eflx, ch4bnew, ch4bflx,&
       & ch3o2new, ch3o2flx, ch3oohnew, ch3oohflx, ch3ohenew, &
       & ch3oheflx, ch3ohbnew, ch3ohbflx, ch3choenew, ch3choeflx,&
       & ch3chobnew, ch3chobflx, panenew, paneflx, panbnew, panbflx,&
       & ch3co3new, ch3co3flx, ch3co3hnew, ch3co3hflx, hcoohenew,&
       & hcooheflx, hcoohbnew, hcoohbflx, ch3coohenew, ch3cooheflx,&
       & ch3coohbnew, ch3coohbflx, c2h6enew, c2h6eflx, c2h6bnew,&
       & c2h6bflx, c2h4enew, c2h4eflx, c2h4bnew, c2h4bflx,&
       & c2h4o2new, c2h4o2flx, c2h4onew, c2h4oflx, c2h4oohnew,&
       & c2h4oohflx, c3h6enew, c3h6eflx, c3h6bnew, c3h6bflx, &
       & acetenew, aceteflx, acetbnew, acetbflx, c2h6o2new, c2h6o2flx,&
       & c2h6oohnew, c2h6oohflx, c3h6o2new, c3h6o2flx, c3h6oohnew,&
       & c3h6oohflx, aceto2new, aceto2flx, acetpnew, acetpflx,&
       & no3enew, no3eflx, no3bnew, no3bflx, n2o5enew, n2o5eflx, &
       & n2o5bnew, n2o5bflx, hno4enew, hno4eflx, hno4bnew, hno4bflx, &
       & honoenew, honoeflx, honobnew, honobflx, ohch2choenew, &
       & ohch2choeflx, ohch2chobnew, ohch2chobflx, ohch2co3new, &
       & ohch2co3flx, ch3o2no2enew, ch3o2no2eflx, ch3o2no2bnew, ch3o2no2bflx

  !--------------------------------------------------------------------!
  ! arrays, fluxes and data for surface model & related interfaces     !
  ! surfaceflux.F90 (and related files)                                !
  ! variables with _ATHAM extension are MODEL-ATHAM coupling variables !
  ! variables with _OUT extension are pure diagnostics                 !
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! general:                                                           !
  !                                                                    !
  ! surface_model          flag is .TRUE. if surface model to be used  !
  ! use_dynamic_osa        flag is .TRUE. if ocean surface albedo is   !
  !                        calculated dynamically (surface_model=.T.)  !
  !                                                                    !
  ! ground properties:                                                 !
  !                                                                    !
  ! es_ATHAM     [/]   ground surface emissivity (avg over gridbox)    !
  ! al_ATHAM     [/]   ground surface vis albedo (avg over gridbox)    !
  ! Ts_ATHAM     [K]   ground surface skin temperature                 !
  ! urbanvegmask [/]   urban-vegetation mask, where landseamask=1:     !
  !                    veg=0, urbld=1, urbhd=2                         !
  ! fh_OUT       [W/m2]ground sensible heat flux                       !
  ! ET_OUT       [mm/s]ground evapo-transpiration (prop to latent heat)!
  !                                                                    !
  ! parameters for dynamic_osa                                         !
  !                                                                    !
  ! sea_chlorophyll [mg/m3] ocean chlorophyll content                  !
  ! tau_aeros       [/]     aerosol optical thickness, for estimating  !
  !                         diffuse radiation component                !
  !--------------------------------------------------------------------!
  logical, save :: surface_model=.false., mean_flux_out=.false.
  logical, save :: use_dynamic_osa

  real(kreal),allocatable, save, dimension(:,:) ::                      &
       & es_ATHAM, al_ATHAM, Ts_ATHAM, urbanvegmask,                    &
       & fh_OUT, ET_OUT, fh_mean_OUT, ET_mean_OUT
  real(kreal), save :: sea_chlorophyll, tau_aeros

end module process_data
