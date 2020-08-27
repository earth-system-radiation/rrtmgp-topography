
subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "test_flux_compute stopping"
    stop
  end if

end subroutine stop_on_err
!-----------------------------
program flux_compute
  use mo_rte_kind,        only: wp
  use mo_gas_optics_rrtmgp, &
                        only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,  &
                        only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props_2str, ty_optical_props_1scl
  use mo_fluxes_byband, only: ty_fluxes_byband
  ! ---- RRTMPG driver
  use mo_rrtmgp_clr_all_sky,  &
                        only: rte_lw, rte_sw
  use mo_heating_rates, only: compute_heating_rate

  ! ---- I/O for test format files.
  use mo_garand_atmos_io,   only: read_atmos
  use mo_test_files_io,     only: is_lw, is_sw, &
                                  read_lw_bc, read_sw_bc, read_lw_rt,  &
                                  write_fluxes, write_dir_fluxes, write_heating_rates, &
                                  write_spectral_disc
  use mo_load_coefficients, only: load_and_init

  implicit none
  real, parameter :: pi = acos(-1._wp)
  ! ----------------------------------------------------------------------------------

  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay
  real(wp), dimension(:,:),   allocatable :: p_lev, t_lev
  real(wp), dimension(:,:),   allocatable :: col_dry
  ! Longwave only
  real(wp), dimension(:),     allocatable :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band
  ! Shortwave only
  real(wp), dimension(:),     allocatable :: sza, tsi, mu0
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  real(wp)                                :: tsi_scaling = -999._wp


  real(wp), dimension(:,:  ), target, &
                              allocatable ::     flux_up,      flux_dn, &
                                                 flux_net,     flux_dir, &
                                             clr_flux_up,  clr_flux_dn
  real(wp), dimension(:,:,:), target, &
                              allocatable :: bnd_flux_up,  bnd_flux_dn, &
                                             bnd_flux_net, bnd_flux_dir


  real(wp), dimension(:,:  ), target, &
                  allocatable ::     plux_up,      plux_dn, &
                                    plux_net,     plux_dir
  real(wp), dimension(:,:,:), target, &
                  allocatable :: bnd_plux_up,  bnd_plux_dn, &
                                  bnd_plux_net, bnd_plux_dir
                                                                  
  real(wp), dimension(:,:),   allocatable :: heating_rate
  real(wp), dimension(:,:,:), allocatable :: bnd_heating_rate

  ! topography related
  real(wp), dimension(:,:),   allocatable :: Vlw  ! LW factors
  real(wp), dimension(:),     allocatable :: Vsw  ! SW factor 

  !
  ! Derived types for interacting with RRTMGP
  !
  type(ty_gas_optics_rrtmgp) :: k_dist

  type(ty_optical_props_2str)       :: clouds,    clouds_subset
  ! type(ty_optical_props_1scl)       :: clouds,    clouds_subset
  type(ty_gas_concs)                :: gas_concs, gas_concs_subset
  type(ty_fluxes_byband)            :: allsky_fluxes, clrsky_fluxes

  integer :: ncol, nlay, nbnd, ngpt, nlev
  integer :: b, nBlocks, colS, colE, nSubcols, nang
  integer :: bottom_level
  integer, parameter :: blockSize = 4
  logical, parameter :: effectTopInclude = .true.

  character(len=256) :: k_dist_file, input_file
  ! ----------------------------------------------------------------------------------
  !
  ! k-distribution file and input-output files must be paired: LW or SW
  !
  k_dist_file = 'coefficients.nc'
  input_file = "rrtmgp-flux-inputs-outputs.nc"

  !
  ! Read temperature, pressure, gas concentrations, then variables specific
  !  to LW or SW problems. Arrays are allocated as they are read
  !
  call read_atmos(input_file,                         &
                   p_lay, t_lay, p_lev, t_lev,        &
                   gas_concs, col_dry)
  if(is_lw(input_file)) then
    call read_lw_bc(input_file, t_sfc, emis_sfc)
    ! Number of quadrature angles
    call read_lw_rt(input_file, nang)
  else
    call read_sw_bc(input_file, sza, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    allocate(mu0(size(sza)))
    mu0 = cos(sza * pi/180.)
  end if

  call load_and_init(k_dist, k_dist_file, gas_concs)
  if(k_dist%source_is_internal() .neqv. is_lw(input_file)) &
    call stop_on_err("flux_compute: gas optics and input-output file disagree about SW/LW")

  !
  ! Problem sizes; allocation of output arrays for full problem
  !
  ncol = size(p_lay, 1)
  nlay = size(p_lay, 2)
  nlev = nlay + 1
  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()

  if (effectTopInclude) then
    if(is_lw(input_file)) then
      ! topography correction
      call read_lw_topography()

      do b=1, ncol
        print '(I3, 1p10E12.4)', b, Vlw(b, :)
      enddo
    else
      call read_sw_topography ()
      call adjustDirAlbedo(ncol, nbnd, sfc_alb_dir, Vsw)
    end if
  end if

  ! Need to ask for some clear-sky fluxes or will get an error.
  allocate(clr_flux_up (ncol,nlev     ), clr_flux_dn(ncol,nlev     ))
  allocate(    flux_up (ncol,nlev     ),     flux_dn(ncol,nlev     ), &
               flux_net(ncol,nlev     ))
  allocate(bnd_flux_up (ncol,nlev,nbnd), bnd_flux_dn(ncol,nlev,nbnd), &
           bnd_flux_net(ncol,nlev,nbnd))
  allocate(heating_rate(ncol,nlay), bnd_heating_rate(ncol,nlay,nbnd))

  allocate(    plux_up (ncol,nlev     ),     plux_dn(ncol,nlev     ), &
               plux_net(ncol,nlev     ))
  allocate(bnd_plux_up (ncol,nlev,nbnd), bnd_plux_dn(ncol,nlev,nbnd), &
           bnd_plux_net(ncol,nlev,nbnd))

  if(is_sw(input_file)) allocate(flux_dir(ncol,nlev), bnd_flux_dir(ncol,nlev,nbnd))

  call stop_on_err(clouds%init(k_dist%get_band_lims_wavenumber(), &
                               k_dist%get_band_lims_gpoint()))

! 1scl
  ! call stop_on_err(clouds%alloc_1scl(ncol,nlay))

! 2str
  call stop_on_err(clouds%alloc_2str(ncol,nlay))
    clouds%ssa(:,:,:) = 1._wp
    clouds%g  (:,:,:) = 0._wp
                         
  ! clouds
  clouds%tau(:,:,:) = 0._wp

  if (p_lay(1, 1) < p_lay(1, nlay)) then
    bottom_level = nlev
  else
    bottom_level = 1
  endif

  !
  ! Loop over subsets of the problem
  !
  nBlocks = ncol/blockSize ! Integer division
  ! print *, "Doing ", nBlocks, "blocks of size ", blockSize
  do b = 1, nBlocks
    colS = (b-1) * blockSize + 1
    colE =  b    * blockSize
    nSubcols = blockSize
    clrsky_fluxes%flux_up      => clr_flux_up(colS:colE,:)
    clrsky_fluxes%flux_dn      => clr_flux_dn(colS:colE,:)
    allsky_fluxes%flux_up      => flux_up(colS:colE,:)
    allsky_fluxes%flux_dn      => flux_dn(colS:colE,:)
    allsky_fluxes%flux_net     => flux_net(colS:colE,:)
    allsky_fluxes%bnd_flux_up  => bnd_flux_up(colS:colE,:,:)
    allsky_fluxes%bnd_flux_dn  => bnd_flux_dn(colS:colE,:,:)
    allsky_fluxes%bnd_flux_net => bnd_flux_net(colS:colE,:,:)
    if(is_sw(input_file)) then
      allsky_fluxes%flux_dn_dir     => flux_dir(colS:colE,:)
      allsky_fluxes%bnd_flux_dn_dir => bnd_flux_dir(colS:colE,:,:)
    end if
    call stop_on_err(gas_concs%get_subset(colS, nSubcols, gas_concs_subset))
    call stop_on_err(   clouds%get_subset(colS, nSubcols, clouds_subset))

    if(is_sw(input_file)) then
      if(allocated(col_dry)) then
        if(tsi_scaling > 0.0_wp) then
          call stop_on_err(rte_sw(k_dist, gas_concs_subset,             &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:),                     &
                                   mu0(colS:colE),                         &
                                   sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
                                   clouds_subset,                          &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   col_dry = col_dry(colS:colE,:),         &
                                   tsi_scaling = tsi_scaling ))
        else
          call stop_on_err(rte_sw(k_dist, gas_concs_subset,             &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:),                     &
                                   mu0(colS:colE),                         &
                                   sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
                                   clouds_subset,                          &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   col_dry = col_dry(colS:colE,:) ))
        endif
      else
        if(tsi_scaling > 0.0_wp) then
          call stop_on_err(rte_sw(k_dist, gas_concs_subset,               &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:),                     &
                                   mu0(colS:colE),                         &
                                   sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
                                   clouds_subset,                          &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   tsi_scaling = tsi_scaling ))
        else
          call stop_on_err(rte_sw(k_dist, gas_concs_subset,               &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:),                     &
                                   mu0(colS:colE),                         &
                                   sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
                                   clouds_subset,                          &
                                   allsky_fluxes, clrsky_fluxes ))
        end if
      end if
    else
      if(allocated(col_dry)) then
        call stop_on_err(rte_lw(k_dist, gas_concs_subset,               &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:), t_sfc(colS:colE  ), &
                                   emis_sfc(:,colS:colE), clouds_subset,   &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   col_dry = col_dry(colS:colE,:),         &
                                   t_lev   = t_lev  (colS:colE,:), n_gauss_angles = nang))
      else
        call stop_on_err(rte_lw(k_dist, gas_concs_subset,               &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:), t_sfc(colS:colE  ), &
                                   emis_sfc(:,colS:colE), clouds_subset,   &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   t_lev = t_lev(colS:colE,:), n_gauss_angles = nang))
      end if
      if (effectTopInclude) call adjustFluxTopography_lw(nSubcols,bottom_level,  emis_sfc(:,colS:colE),  Vlw(colS:colE, :), allsky_fluxes)
    end if
  end do
  if(mod(ncol, blockSize) /= 0) then
    colS = ncol/blockSize * blockSize + 1  ! Integer arithmetic
    colE = ncol
    nSubcols = colE-colS+1
    ! print *, "Doing ", nSubcols, "extra columns"
    call stop_on_err(   clouds%get_subset(colS, nSubcols, clouds_subset))
    clrsky_fluxes%flux_up      => clr_flux_up(colS:colE,:)
    clrsky_fluxes%flux_dn      => clr_flux_dn(colS:colE,:)
    allsky_fluxes%flux_up      => flux_up(colS:colE,:)
    allsky_fluxes%flux_dn      => flux_dn(colS:colE,:)
    allsky_fluxes%flux_net     => flux_net(colS:colE,:)
    allsky_fluxes%bnd_flux_up  => bnd_flux_up(colS:colE,:,:)
    allsky_fluxes%bnd_flux_dn  => bnd_flux_dn(colS:colE,:,:)
    allsky_fluxes%bnd_flux_net => bnd_flux_net(colS:colE,:,:)
    if(is_sw(input_file)) then
      allsky_fluxes%flux_dn_dir     => flux_dir(colS:colE,:)
      allsky_fluxes%bnd_flux_dn_dir => bnd_flux_dir(colS:colE,:,:)
    end if
    call stop_on_err(gas_concs%get_subset(colS, nSubcols, gas_concs_subset))
    call stop_on_err(clouds%get_subset(colS, nSubcols, clouds_subset))

    if(is_sw(input_file)) then
      if(allocated(col_dry)) then
        if(tsi_scaling > 0.0_wp) then
          call stop_on_err(rte_sw(k_dist, gas_concs_subset,             &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:),                     &
                                   mu0(colS:colE),                         &
                                   sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
                                   clouds_subset,                          &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   col_dry = col_dry(colS:colE,:),         &
                                   tsi_scaling = tsi_scaling ))
        else
          call stop_on_err(rte_sw(k_dist, gas_concs_subset,             &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:),                     &
                                   mu0(colS:colE),                         &
                                   sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
                                   clouds_subset,                          &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   col_dry = col_dry(colS:colE,:) ))
        endif
      else
        if(tsi_scaling > 0.0_wp) then
          call stop_on_err(rte_sw(k_dist, gas_concs_subset,             &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:),                     &
                                   mu0(colS:colE),                         &
                                   sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
                                   clouds_subset,                          &
                                   allsky_fluxes, clrsky_fluxes ))
        end if
      end if
    else
      if(allocated(col_dry)) then
        call stop_on_err(rte_lw(k_dist, gas_concs_subset,               &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:), t_sfc(colS:colE  ), &
                                   emis_sfc(:,colS:colE), clouds_subset,   &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   col_dry = col_dry(colS:colE,:),         &
                                   t_lev   = t_lev  (colS:colE,:), n_gauss_angles = nang))
      else
        call stop_on_err(rte_lw(k_dist, gas_concs_subset,               &
                                   p_lay(colS:colE,:), t_lay(colS:colE,:), &
                                   p_lev(colS:colE,:), t_sfc(colS:colE  ), &
                                   emis_sfc(:,colS:colE), clouds_subset,   &
                                   allsky_fluxes, clrsky_fluxes,           &
                                   t_lev = t_lev(colS:colE,:), n_gauss_angles = nang))
      end if
      if (effectTopInclude) call adjustFluxTopography_lw(nSubcols,bottom_level,  emis_sfc(:,colS:colE),  Vlw(colS:colE, :), allsky_fluxes)
    end if
  end if
  !
  ! Heating rates
  !
  call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))
  do b = 1, nbnd
      call stop_on_err(compute_heating_rate(bnd_flux_up(:,:,b), bnd_flux_dn(:,:,b), p_lev, bnd_heating_rate(:,:,b)))
  end do

  !
  ! ... and write everything out
  !
  call write_spectral_disc(input_file, clouds)
  call write_fluxes(input_file, flux_up, flux_dn, flux_net, bnd_flux_up, bnd_flux_dn, bnd_flux_net)
  call write_heating_rates(input_file, heating_rate, bnd_heating_rate)
  if(k_dist%source_is_external()) &
    call write_dir_fluxes(input_file, flux_dir, bnd_flux_dir)

    if (allocated(Vlw)) deallocate(Vlw)
    if (allocated(Vsw)) deallocate(Vsw)
    contains

    subroutine read_lw_topography()
      use mo_simple_netcdf, only: read_field, read_char_vec, read_logical_vec, var_exists, get_dim_size
      use netcdf
      ! local
      integer ncid
      integer ncol
      character(len=*), parameter :: fname='ref/rrtmgp-lw-topography-sample.nc'
      print  * , 'reading ', fname
      if(nf90_open(fname, NF90_NOWRITE, ncid) /= NF90_NOERR) &
          call stop_on_err("read_lw_topography(): can't open file " // trim(fname))

      ncol = get_dim_size(ncid,'col')
      allocate ( Vlw(ncol, 3) ) 
      Vlw(:,1) = read_field(ncid, 'V', ncol)
      Vlw(:,2) = read_field(ncid, 'Vav', ncol)
      Vlw(:,3) = read_field(ncid, 'InvCosA', ncol)

      ncid = nf90_close(ncid)

    end subroutine read_lw_topography

    subroutine adjustDirAlbedo(ncol, nband, albedo, factor)
      integer,                      intent(in   ) :: ncol, nband
      real(wp), dimension(nband, ncol),   intent(inout) :: albedo          
      real(wp), dimension(ncol),          intent(in   ) :: factor          

      integer :: icol
      do icol=1, ncol
        albedo(:, icol) = albedo(:, icol)*factor(icol)
      enddo
    end subroutine adjustDirAlbedo

    subroutine read_sw_topography()
      use mo_simple_netcdf, only: read_field, read_char_vec, read_logical_vec, var_exists, get_dim_size
      use netcdf
      ! local
      integer ncid
      integer ncol, nzen, nazm, icol
      character(len=*), parameter :: fname='ref/rrtmgp-sw-topography-sample.nc'

      real(wp), dimension(:),      allocatable :: solarAzm
      real(wp), dimension(:),      allocatable :: solarZen
      real(wp), dimension(:),      allocatable :: solarAzmGrid 
      real(wp), dimension(:),      allocatable :: solarZenGrid 
      real(wp), dimension(:,:,:),  allocatable :: oxAvMask 

      ! reading solar zenith and azimuth fro a computational cells
      if(nf90_open(trim(input_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
           call stop_on_err("read_lw_topography(): can't open file " // trim(input_file))
      ncol = get_dim_size(ncid,'col')
      allocate ( solarAzm(ncol), solarZen(ncol) ) 
      solarAzm= read_field(ncid, 'solar_zenith_angle', ncol)
      solarZen= read_field(ncid, 'solar_azimuth_angle', ncol)
      ncid = nf90_close(ncid)

      ! reading topography parameter tables 
      if(nf90_open(fname, NF90_NOWRITE, ncid) /= NF90_NOERR) &
           call stop_on_err("read_lw_topography(): can't open file " // trim(fname))
      nzen = get_dim_size(ncid,'zen')
      nazm = get_dim_size(ncid,'azm')

      allocate ( solarAzmGrid(nazm), solarZenGrid(nzen) ) 
      allocate ( Vsw(ncol) ) 

      allocate ( oxAvMask(nzen, nazm, ncol)) 
      oxAvMask = read_field(ncid, 'oxAvMask', nzen, nazm, ncol)

      solarAzmGrid = read_field(ncid, 'solAzm', nazm)
      solarZenGrid = read_field(ncid, 'solZen', nzen)

      ncid = nf90_close(ncid)
      ! reading topography parameter iterpolation to local solar azimuth and zenith 

      do icol=1, ncol
        Vsw(icol) = interpolate2D(nzen, nazm, solarZen(icol), solarAzm(icol), &
                                  solarZenGrid, solarAzmGrid, oxAvMask(:,:,icol))
      enddo

      deallocate ( oxAvMask ) 
      deallocate ( solarAzm, solarZen) 
      deallocate ( solarAzmGrid, solarZenGrid) 
    end subroutine read_sw_topography

!--------------------------------------------------------------------------------------------------------------------
  !
  ! linear interpolation in 2D array
  !
    function  interpolate2D(nx, ny, x, y, xGrid, yGrid, dat) result(factor)
      integer,                      intent(in   ) :: nx, ny
      real(wp),                     intent(in   ) :: x, y
      real(wp), dimension(nx),      intent(in   ) :: xGrid          
      real(wp), dimension(ny),      intent(in   ) :: yGrid          
      real(wp), dimension(nx, ny),  intent(in   ) :: dat     

      real(wp)                                    :: factor

      ! -------------
      integer :: tx, ty
      integer :: ix, iy
      do ix=2, nx-1
        if (xGrid(ix) > x) exit
      enddo
      tx = (xGrid(ix) - x)/(xGrid(ix)-xGrid(ix-1))
      do iy=2, ny-1
        if (yGrid(iy) > y) exit
      enddo
      ty = (yGrid(iy) - y)/(yGrid(iy)-yGrid(iy-1))
      factor = ( dat(ix-1, iy-1)*tx+dat(ix, iy-1)*(1_wp-tx) )*ty + ( dat(ix-1, iy)*tx+dat(ix, iy)*(1_wp-tx) )*(1_wp-ty)
    end function  interpolate2D

 !--------------------------------------------------------------------------------------------------------------------
    subroutine adjustFluxTopography_lw(ncol,bottom_level, &
      sfc_emis,   V, fluxes)
use mo_fluxes_broadband_kernels ,  only : sum_broadband

integer,                      intent(in   ) :: ncol           ! column number
integer,                      intent(in   ) :: bottom_level          ! bottom level index
                         ! (if not, ordering is bottom-to-top)
real(wp), dimension(:,:),     intent(in   ) :: sfc_emis       ! emissivity at surface [] (nband, ncol)
real(wp), dimension(:,:),     intent(in   ) :: V              ! V(:,1) - view factor at surface [] <V>  (ncol)
                                 !V(:,2) -  weighted view factor at surface <V/cosi> [] (ncol)
                                 !V(:,3) -  inverse surface slope cosine  <1/cosi> [] (ncol)
class(ty_fluxes_byband),      intent(inout) :: fluxes
! ------

integer                                     :: nband, ngpt
real(wp)                                    :: factor
integer                                     :: icol, ibnd
nband = size(sfc_emis, DIM=1)

! band fluxes
!$acc parallel loop collapse(2) copyin(fluxes, sfc_emis, V) copyout(fluxes)
do icol=1, ncol
do ibnd = 1, nband
factor = 1_wp - (1_wp - sfc_emis(ibnd,icol))*(1. - V(icol,1))
fluxes%bnd_flux_net(icol, bottom_level, ibnd)=(fluxes%bnd_flux_dn(icol, bottom_level, ibnd) &
       - fluxes%bnd_flux_up(icol, bottom_level, ibnd))/factor
enddo
enddo

!$acc parallel loop collapse(2) copyin(fluxes, sfc_emis, V) copyout(fluxes)
do ibnd = 1, nband
do icol=1, ncol
fluxes%bnd_flux_up(icol, bottom_level, ibnd) = V(icol, 3)*( &
fluxes%bnd_flux_dn(icol, bottom_level, ibnd) - fluxes%bnd_flux_net(icol, bottom_level,ibnd))

fluxes%bnd_flux_dn(icol, bottom_level, ibnd) = V(icol, 3)*fluxes%bnd_flux_dn(icol, bottom_level, ibnd) &
-(V(icol, 3) - V(icol, 2))*fluxes%bnd_flux_net(icol, bottom_level, ibnd)

fluxes%bnd_flux_net(icol, bottom_level, ibnd) = V(icol, 2)*fluxes%bnd_flux_net(icol, bottom_level, ibnd)

enddo
enddo

! total fluxes
call sum_broadband(ncol, 1, nband, fluxes%bnd_flux_up (:, bottom_level, :), fluxes%flux_up (:, bottom_level) )
call sum_broadband(ncol, 1, nband, fluxes%bnd_flux_dn (:, bottom_level, :), fluxes%flux_dn (:, bottom_level) )
call sum_broadband(ncol, 1, nband, fluxes%bnd_flux_net(:, bottom_level, :), fluxes%flux_net(:, bottom_level) )

end subroutine adjustFluxTopography_lw
! --------------------------------      
end program flux_compute
