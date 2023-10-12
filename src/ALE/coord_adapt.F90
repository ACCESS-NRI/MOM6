! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Regrid columns for the adaptive coordinate
module coord_adapt

use MOM_coms,          only : reproducing_sum
use MOM_EOS,           only : calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use filter_utils,      only : filter_CS, filtered_grid_motion
use coord_zlike,       only : init_coord_zlike, end_coord_zlike, zlike_CS, set_zlike_params, build_zstar_column

implicit none ; private

#include <MOM_memory.h>

!> AG coordinate diagnostic control structure
type, public :: adapt_diag_CS
  !> Along-coordinate i-gradient of density (used for density term)
  real, dimension(:,:,:), allocatable :: slope_u
  !> Along-coordinate j-gradient of density (used for density term)
  real, dimension(:,:,:), allocatable :: slope_v

  !> Denominator used for calculating density displacement, i-direction
  real, dimension(:,:,:), allocatable :: denom_u
  !> Denominator used for calculating density displacement, j-direction
  real, dimension(:,:,:), allocatable :: denom_v

  !> Physical-space slope of interface along i-direction (used for density weighting)
  real, dimension(:,:,:), allocatable :: phys_u
  !> Physical-space slope of interface along j-direction (used for density weighting)
  real, dimension(:,:,:), allocatable :: phys_v

  !> Coordinate-space slope of interface along i-direction (used for density weighting)
  real, dimension(:,:,:), allocatable :: coord_u
  !> Coordinate-space slope of interface along j-direction (used for density weighting)
  real, dimension(:,:,:), allocatable :: coord_v

  !> Amount of limiting applied to density (before weighting)
  real, dimension(:,:,:), allocatable :: limiting_density
  !> Amount of limiting applied to smoothing (before weighting)
  real, dimension(:,:,:), allocatable :: limiting_smoothing

  !> The adjustment provided by the convective adjustment term
  real, dimension(:,:,:), allocatable :: w_adjust

  !> Interface displacement due to density term
  real, dimension(:,:,:), allocatable :: disp_density
  !> Interface displacement due to smoothing term
  real, dimension(:,:,:), allocatable :: disp_smoothing
  !> Interface displacement due to unlimited smoothing term
  real, dimension(:,:,:), allocatable :: disp_unlimited

  !>@{ Diagnostics IDs
  integer :: id_slope_u, id_slope_v
  integer :: id_denom_u, id_denom_v
  integer :: id_phys_u, id_phys_v
  integer :: id_coord_u, id_coord_v
  integer :: id_limiting_density, id_limiting_smoothing
  integer :: id_w_adjust
  integer :: id_disp_density, id_disp_smoothing, id_disp_unlimited
  !>@}
end type adapt_diag_CS

!> Control structure for adaptive coordinates (coord_adapt).
type, public :: adapt_CS ; private

  !> Number of layers/levels
  integer :: nk

  !> Nominal near-surface resolution [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: coordinate_resolution

  !> If positive, a manual coefficient for the density adaptivity term.
  !! If negative, either density or pressure adaptivity are chosen,
  !! depending on the local coordinate slope, with a minimum of min_smooth
  !! going toward the pressure term.
  real :: alpha_rho

  !> The complement of alpha_rho: a positive value is a manually-specified
  !! coefficient; a negative value is automatically-determined, with a
  !! value of at least min_smooth.
  real :: alpha_p

  !> Minimum weighting of the pressure adaptivity (smoothing) term, used
  !! when alpha_rho and alpha_p are negative.
  real :: min_smooth

  !> The timescale over which to apply the diffusive adaptivity terms. [T ~> s]
  real :: adaptivity_timescale

  !> The timescale over which to restore towards the calculated
  !! or pre-defined target coordinate. [T ~> s]
  real :: restoring_timescale

  !> Interface slope cutoff for defining stratified/unstratified regions.
  real :: slope_cutoff

  !> If true, use the uniform mean of thicknesses where required.
  !! Otherwise, use the "upstream" thickness in the direction of
  !! interface movement due to adaptivity.
  logical :: use_mean_h

  !> If true, the on-interface density gradient is calculated in the layers
  !! above and below. They must agree on sign to prevent a null mode, and the
  !! minimum is chosen, to prefer smoothing.
  !! Otherwise, the gradient is calculated directly on the interface.
  logical :: use_twin_gradient

  !> If true, calculate the slope in physical space (taking into account the
  !! vertical distance between adjacent points). Otherwise, the slope is only
  !! calculated along the interface.
  logical :: use_physical_slope

  !> If true, restore towards the dynamically-determined mean position of
  !! a given interface. Otherwise, use the specified coordinate locations.
  logical :: do_restore_mean

  !> The non-dimensional scale for the adjustment performed for diagonal
  !! convective instabilities.
  real :: adjustment_scale

  !> Used if do_restore_mean is .false.: delegate to a zlike coordinate
  !! for the restoring term target.
  type(zlike_CS), pointer :: zlike_CS => null()

  !> Used for outputting diagnostics from within the regridding routine.
  type(adapt_diag_CS), pointer :: diag_CS => null()
end type adapt_CS

public init_coord_adapt, set_adapt_params, build_adapt_grid, end_coord_adapt
public associate_adapt_diag, get_adapt_diag_CS

contains

!> Initialise an adapt_CS with parameters
subroutine init_coord_adapt(CS, nk, coordinate_resolution)
  type(adapt_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk !< Number of layers in the grid
  real, dimension(:), intent(in) :: coordinate_resolution !< Nominal near-surface resolution [m] or
                                       !! other units specified with m_to_H

  if (associated(CS)) call MOM_error(FATAL, "init_coord_adapt: CS already associated")
  allocate(CS)
  allocate(CS%coordinate_resolution(nk))

  CS%nk = nk
  CS%coordinate_resolution(:) = coordinate_resolution(:)

  CS%alpha_rho = -1.0
  CS%alpha_p   = -1.0

  CS%use_mean_h = .false.
  CS%use_twin_gradient = .true.
  CS%use_physical_slope = .true.
  CS%do_restore_mean = .false.

  call init_coord_zlike(CS%zlike_CS, nk, coordinate_resolution)

end subroutine init_coord_adapt

!> Clean up the coordinate control structure
subroutine end_coord_adapt(CS)
  type(adapt_CS), pointer :: CS  !< The control structure for this module

  ! nothing to do
  if (.not. associated(CS)) return

  call end_coord_zlike(CS%zlike_CS)

  if (associated(CS%diag_CS)) deallocate(CS%diag_CS)

  deallocate(CS%coordinate_resolution)
  deallocate(CS)
end subroutine end_coord_adapt

!> This subtroutine can be used to set the parameters for coord_adapt module
subroutine set_adapt_params(CS, alpha_rho, alpha_p, adaptivity_timescale, use_mean_h, &
     use_twin_gradient, slope_cutoff, min_smooth, use_physical_slope, restoring_timescale, do_restore_mean, &
     adjustment_scale)

  type(adapt_CS),    pointer    :: CS  !< The control structure for this module
  real,    optional, intent(in) :: alpha_rho !< Density adaptivity coefficient
  real,    optional, intent(in) :: alpha_p !< Pressure adaptivity coefficient
  real,    optional, intent(in) :: adaptivity_timescale !< Adaptivity timescale
  logical, optional, intent(in) :: use_mean_h !< Use uniform or "upstream" mean thickness?
  logical, optional, intent(in) :: use_twin_gradient !< Calculate interface density gradient layers above and below
  real,    optional, intent(in) :: slope_cutoff !< Stratified/unstratified cutoff
  real,    optional, intent(in) :: min_smooth !< Minimum pressure adaptivity contribution
  logical, optional, intent(in) :: use_physical_slope !< Use physical or along-interface slope
  real,    optional, intent(in) :: restoring_timescale !< Timescale for restoring term
  logical, optional, intent(in) :: do_restore_mean !< Restore to the mean height?
  real,    optional, intent(in) :: adjustment_scale !< Hydrostatic adjustment scale

  if (.not. associated(CS)) call MOM_error(FATAL, "set_adapt_params: CS not associated")

  if (present(alpha_rho))            CS%alpha_rho = alpha_rho
  if (present(alpha_p))              CS%alpha_p = alpha_p
  if (present(adaptivity_timescale)) CS%adaptivity_timescale = adaptivity_timescale
  if (present(use_mean_h))           CS%use_mean_h = use_mean_h
  if (present(use_twin_gradient))    CS%use_twin_gradient = use_twin_gradient
  if (present(slope_cutoff))         CS%slope_cutoff = slope_cutoff
  if (present(min_smooth))           CS%min_smooth = min_smooth
  if (present(use_physical_slope))   CS%use_physical_slope = use_physical_slope
  if (present(restoring_timescale))  CS%restoring_timescale = restoring_timescale
  if (present(do_restore_mean))      CS%do_restore_mean = do_restore_mean
  if (present(adjustment_scale))     CS%adjustment_scale = adjustment_scale
end subroutine set_adapt_params

!> Associate a diagnostic control structure with an existing
!! AG control structure -- used to get around the circular
!! dependency of diagnostics depending on coordinates.
subroutine associate_adapt_diag(CS, diag_CS)
  type(adapt_CS), pointer :: CS
  type(adapt_diag_CS), target :: diag_CS

  if (associated(CS%diag_CS)) deallocate(CS%diag_CS)
  CS%diag_CS => diag_CS
end subroutine associate_adapt_diag

!> Return the associated diagnostic control structure for an
!! AG control structure
function get_adapt_diag_CS(CS)
  type(adapt_CS), intent(in) :: CS
  type(adapt_diag_CS), pointer :: get_adapt_diag_CS

  get_adapt_diag_CS => CS%diag_CS
end function get_adapt_diag_CS

!> Calculate the along-coordinate density derivatives
!! and the physical analogue thereof. The derivatives can
!! be calculated in the i- or j-direction, depending on the
!! value of di/dj.
subroutine calc_derivs(G, GV, CS, US, h, z_int, tv, i, j, k, &
     di, dj, dk_sig_int, alpha, beta, Idx, mask, hd_sig, hd_sig_phys)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  type(adapt_CS), intent(in) :: CS
  type(unit_scale_type), intent(in) :: US
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)), intent(in) :: h
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1), intent(in) :: z_int
  type(thermo_var_ptrs), intent(in) :: tv
  integer, intent(in) :: i, j, k, di, dj
  real, dimension(SZI_(G), SZJ_(G)), intent(in) :: dk_sig_int
  real, intent(in) :: alpha, beta, Idx, mask
  real, intent(out) :: hd_sig, hd_sig_phys

  real :: H_to_L
  real :: d_sig_up, d_sig_dn, d_sig, dk_sig, h_interp

  H_to_L = GV%H_to_Z * US%Z_to_L

  if (CS%use_twin_gradient) then
    d_sig_up = alpha * (tv%t(i+di,j+dj,k-1) - tv%t(i,j,k-1)) &
         + beta * (tv%s(i+di,j+dj,k-1) - tv%s(i,j,k-1))
    d_sig_dn = alpha * (tv%t(i+di,j+dj,k) - tv%t(i,j,k)) &
         + beta * (tv%s(i+di,j+dj,k) - tv%s(i,j,k))

    if (d_sig_up * d_sig_dn <= 0.) then
      d_sig = 0.
    else
      d_sig = sign(min(abs(d_sig_up), abs(d_sig_dn)), d_sig_up)
    end if
  end if

  dk_sig = 0.5 * (dk_sig_int(i,j) + dk_sig_int(i+di,j+dj))

  if (d_sig * dk_sig < 0.) then
    h_interp = 0.5 * (h(i,j,k-1) + h(i+di,j+dj,k))
  else
    h_interp = 0.5 * (h(i,j,k) + h(i+di,j+dj,k-1))
  end if

  if (CS%use_mean_h) &
       h_interp = 0.25 * ((h(i,j,k-1) + h(i+di,j,k)) + (h(i,j,k) + h(i+di,j,k-1)))

  hd_sig = h_interp * d_sig * Idx * H_to_L * mask
  hd_sig_phys = hd_sig - Idx * dk_sig * (z_int(i+di,j+dj,K) - z_int(i,j,K)) * H_to_L * mask
end subroutine calc_derivs

subroutine build_adapt_grid(G, GV, US, h, tv, dzInterface, CS, fCS, min_thickness, dt)
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US   !< The dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h    !< Layer thicknesses, in H (usually m or kg m-2)
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: dzInterface !< The interface changes
  type(adapt_CS),                              intent(in)    :: CS  !< Regridding control structure
  type(filter_CS),                             intent(in)    :: fCS !< Filtering control structure
  real,                                        intent(in)    :: min_thickness !< ALE layer minimum thickness
  real,                              optional, intent(in)    :: dt !< The intended timestep over which this
                                                                   !! regridding operation applies

  ! local variables
  integer :: i, j, k, k2, kt, nz ! indices and dimension lengths

  ! interface heights
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: z_int, z_new, h_int
  ! drho/dt and drho/ds on interfaces
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: alpha_int
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: beta_int
  ! vertical gradient in sigma
  real, dimension(SZI_(G),SZJ_(G)) :: dk_sig_int
  ! final change in interface height
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: dz_a, dz_p, dz_r

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_upd
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: w

  ! interface position after adaptivity, mean interface position across basin
  real, dimension(SZK_(GV)+1) :: z_mean, h_col, z_col, z_upd, dz_col

  ! numerator of density term and upstreamed h
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: hdi_sig, hdi_sig_phys
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: hdj_sig, hdj_sig_phys
  ! temporary alpha/beta interpolated to velocity points
  real :: alpha, beta
  ! some temporary quantities
  real :: eps, weight, weight2, h_interp, i_denom, j_denom
  ! numerator (and intermediates) of density term before multiplication by h
  real :: di_sig, di_sig_up, di_sig_dn
  real :: dj_sig, dj_sig_up, dj_sig_dn
  ! difference quantities interpolated to other locations
  real :: hdi_sig_u, hdj_sig_u, hdi_sig_v, hdj_sig_v, dk_sig_u, dk_sig_v
  real :: ts_ratio, slope, phys_slope
  real :: global_z_sum, global_h_sum
  real :: dz_p_unlim
  real :: tmp, dir, CFL
  real :: dsig_horiz, dsig_vert_up, dsig_vert_down
  real :: H_to_L, L_to_H

  logical :: do_diag

  character(len=11) :: fname

  ! we could probably assume some limit without a specified timestep
  if (.not. present(dt)) then
    dzInterface(:,:,:) = 0.0
    return
  end if

  eps = 1. ; eps = epsilon(eps)
  nz = GV%ke

  L_to_H = US%L_to_Z * GV%Z_to_H

  call set_zlike_params(CS%zlike_CS, min_thickness=min_thickness)

  ! zero out diagnostic arrays
  do_diag = .true.
  if (.not. associated(CS%diag_CS)) then
    call MOM_error(WARNING, 'build_adapt_grid expected diag_CS associated')
    do_diag = .false.
  end if

  if (do_diag) then
    if (allocated(CS%diag_CS%phys_u)) CS%diag_CS%phys_u(:,:,:) = 0.0
    if (allocated(CS%diag_CS%phys_v)) CS%diag_CS%phys_v(:,:,:) = 0.0
    if (allocated(CS%diag_CS%slope_u)) CS%diag_CS%slope_u(:,:,:) = 0.0
    if (allocated(CS%diag_CS%slope_v)) CS%diag_CS%slope_v(:,:,:) = 0.0
    if (allocated(CS%diag_CS%denom_u)) CS%diag_CS%denom_u(:,:,:) = 0.0
    if (allocated(CS%diag_CS%denom_v)) CS%diag_CS%denom_v(:,:,:) = 0.0
    if (allocated(CS%diag_CS%coord_u)) CS%diag_CS%coord_u(:,:,:) = 0.0
    if (allocated(CS%diag_CS%coord_v)) CS%diag_CS%coord_v(:,:,:) = 0.0
    if (allocated(CS%diag_CS%limiting_smoothing)) CS%diag_CS%limiting_smoothing(:,:,:) = 0.0
    if (allocated(CS%diag_CS%limiting_density)) CS%diag_CS%limiting_density(:,:,:) = 0.0
    if (allocated(CS%diag_CS%w_adjust)) CS%diag_CS%w_adjust(:,:,:) = 0.0
    if (allocated(CS%diag_CS%disp_density)) CS%diag_CS%disp_density(:,:,:) = 0.0
    if (allocated(CS%diag_CS%disp_smoothing)) CS%diag_CS%disp_smoothing(:,:,:) = 0.0
    if (allocated(CS%diag_CS%disp_unlimited)) CS%diag_CS%disp_unlimited(:,:,:) = 0.0
  end if

  ! sum from free surface downward
  z_int(:,:,1) = sum(h, 3) - G%bathyT(:,:) * GV%Z_to_H ! free-surface
  do K = 1,nz
    z_int(:,:,K+1) = z_int(:,:,K) - h(:,:,k)
  enddo

  if (CS%do_restore_mean) then
    ! calculate geometric mean of thicknesses on interfaces
    ! we only need to do this in our own domain because this
    ! is a global sum
    z_new(:,:,:) = 0. ;  h_int(:,:,:) = 0.
    do j = G%jsc,G%jec
      do i = G%isc,G%iec
        h_int(i,j,2:nz) = (h(i,j,2:nz) * h(i,j,1:nz-1)) / &
             (h(i,j,2:nz) + h(i,j,1:nz-1) + GV%H_subroundoff)
        ! we don't really want to volume-weight this, we just want to discount vanished layers
        ! this way, we won't bias towards thick layers
        h_int(i,j,2:nz) = max(GV%H_to_m * h_int(i,j,2:nz), 1.0)
        h_int(i,j,2:nz) = h_int(i,j,2:nz) * (G%areaT(i,j) * G%mask2dT(i,j))
        ! weight height by thickness
        z_new(i,j,2:nz) = z_int(i,j,2:nz) * h_int(i,j,2:nz)
      enddo
    enddo
    global_z_sum = reproducing_sum(z_new, G%isc, G%iec, G%jsc, G%jec, sums=z_mean)
    global_h_sum = reproducing_sum(h_int, G%isc, G%iec, G%jsc, G%jec, sums=h_col)
    z_mean(2:nz) = z_mean(2:nz) / h_col(2:nz)

    do K = 2,nz-1
      if (z_mean(K) < z_mean(K+1)) then
        print *, z_mean
        call MOM_error(FATAL, 'tangled z_mean')
      endif
    enddo
  else
    ! we'll restore to the predefined coordinate resolution
    z_mean(1) = 0.
    do K = 2,nz
      z_mean(K) = z_mean(K-1) - CS%coordinate_resolution(k-1) * GV%Z_to_H
    end do
  end if

  ! the top and bottom interfaces don't move
  dz_a(:,:,1) = 0. ; dz_a(:,:,nz+1) = 0.
  dz_p(:,:,1) = 0. ; dz_p(:,:,nz+1) = 0.
  dz_r(:,:,1) = 0. ; dz_r(:,:,nz+1) = 0.
  w(:,:,1) = 0. ; w(:,:,nz+1) = 0.

  h_upd(:,:,:) = 0.

  ! nondimensionalise the adaptivity timescale wrt. the ALE timescale
  ! to get a scaling for diffusive adaptivity
  ts_ratio = dt / CS%adaptivity_timescale
  ts_ratio = min(ts_ratio, 1.0)

  !$omp parallel default(none) &
  !$omp          shared(tv, GV, G, CS, US, z_int, h, alpha_int, beta_int) &
  !$omp          shared(hdi_sig, hdj_sig, hdi_sig_phys, hdj_sig_phys) &
  !$omp          shared(L_to_H, ts_ratio, dz_a, dz_p, do_diag, eps, nz) &
  !$omp          private(i, j, k, dk_sig_int, alpha, beta) &
  !$omp          private(hdi_sig_u, hdj_sig_u, dk_sig_u, hdi_sig_v, hdj_sig_v, dk_sig_v, i_denom, j_denom, dz_p_unlim, slope, phys_slope, weight, weight2)
  block
    ! for some reason we get a segfault if these are brought in as private to the
    ! parallel block, so instead we allocate them locally (they'll be deallocated at the
    ! end of the block anyway, but annoying to have to use heap space)
    real, allocatable, dimension(:,:) :: t_int, s_int
    real, allocatable, dimension(:,:) :: dz_s_i, dz_s_j, dz_p_i, dz_p_j, dz_i, dz_j

    allocate(t_int(SZI_(G),SZJ_(G)), s_int(SZI_(G),SZJ_(G)))
    allocate(dz_s_i(SZIB_(G),SZJ_(G)), dz_s_j(SZI_(G),SZJB_(G)))
    allocate(dz_p_i(SZIB_(G),SZJ_(G)), dz_p_j(SZI_(G),SZJB_(G)))
    allocate(dz_i(SZIB_(G),SZJ_(G)), dz_j(SZI_(G),SZJB_(G)))

  !$omp do
  do K = 2,nz
    dz_s_i(:,:) = 0. ; dz_s_j(:,:) = 0.
    dz_p_i(:,:) = 0. ; dz_p_j(:,:) = 0.
    dz_i(:,:) = 0. ; dz_j(:,:) = 0.

    do j = G%jsc-2,G%jec+2
      do i = G%isc-2,G%iec+2
        t_int(i,j) = ( &
             tv%t(i,j,k-1) * (h(i,j,k) + GV%H_subroundoff) + &
             tv%t(i,j,k) * (h(i,j,k-1) + GV%H_subroundoff)) / &
             (h(i,j,k-1) + h(i,j,k) + 2*GV%H_subroundoff)
        s_int(i,j) = ( &
             tv%s(i,j,k-1) * (h(i,j,k) + GV%H_subroundoff) + &
             tv%s(i,j,k) * (h(i,j,k-1) + GV%H_subroundoff)) / &
             (h(i,j,k-1) + h(i,j,k) + 2*GV%H_subroundoff)
      enddo

      call calculate_density_derivs(t_int(:,j), s_int(:,j), -z_int(:,j,K) * GV%H_to_Pa, &
           alpha_int(:,j,K), beta_int(:,j,K), G%isc-2, G%iec+2 - (G%isc-2) + 1, tv%eqn_of_state)

      do i = G%isc-2,G%iec+2
        dk_sig_int(i,j) = alpha_int(i,j,K) * (tv%t(i,j,k) - tv%t(i,j,k-1)) + &
             beta_int(i,j,K) * (tv%s(i,j,k) - tv%s(i,j,k-1))
     enddo
    enddo

    ! calculate horizontal derivatives on i-points
    ! reduce I-halo 2 -> 1
    do j = G%jsc-2,G%jec+2
      do I = G%IscB-1,G%IecB+1
        alpha = 0.5 * (alpha_int(i,j,K) + alpha_int(i+1,j,K))
        beta = 0.5 * (beta_int(i,j,K) + beta_int(i+1,j,K))

        call calc_derivs(G, GV, CS, US, h, z_int, tv, I, j, k, 1, 0, dk_sig_int, alpha, beta, G%IdxCu(I,j), &
             G%mask2dCu(I,j), hdi_sig(I,j,K), hdi_sig_phys(I,j,K))
      enddo
    enddo

    ! calculate horizontal derivatives on j-points
    ! reduce J-halo 2 -> 1
    do J = G%JscB-1,G%JecB+1
      do i = G%isc-2,G%iec+2
        alpha = 0.5 * (alpha_int(i,j,K) + alpha_int(i,j+1,K))
        beta = 0.5 * (beta_int(i,j,K) + beta_int(i,j+1,K))

        call calc_derivs(G, GV, CS, US, h, z_int, tv, i, J, k, 0, 1, dk_sig_int, alpha, beta, G%IdyCv(i,J), &
             G%mask2dCv(i,J), hdj_sig(i,J,K), hdj_sig_phys(i,J,K))
      enddo
    enddo

    ! u-points
    do j = G%jsc-1,G%jec+1
      do I = G%IscB-1,G%IecB+1
        if (G%mask2dCu(I,j) < 0.5) then
          dz_i(I,j) = 0.
          dz_s_i(I,j) = 0.
          dz_p_i(I,j) = 0.
          cycle
        endif

        ! interpolate terms in the denominator onto the u-point
        hdi_sig_u = hdi_sig(I,j,K)**2
        hdj_sig_u = 0.25 * ((hdj_sig(i,J,K)**2 + hdj_sig(i+1,J-1,K)**2) + &
             (hdj_sig(i+1,J,K)**2 + hdj_sig(i,J-1,K)**2))
        dk_sig_u = 0.5 * (dk_sig_int(i,j)**2 + dk_sig_int(i+1,j)**2)

        i_denom = hdi_sig_u + hdj_sig_u + dk_sig_u
        if (abs(i_denom) < eps) then
          ! if gradients in all directions are exactly zero, we don't want any flux
          dz_s_i(I,j) = 0.
        else
          dz_s_i(I,j) = hdi_sig(I,j,K) / sign(sqrt(i_denom), dk_sig_u)
        end if

        if (do_diag) then
          ! DIAG: slope_u
          if (allocated(CS%diag_CS%slope_u)) CS%diag_CS%slope_u(I,j,K) = dz_s_i(I,j)
          ! DIAG: denom_u
          if (allocated(CS%diag_CS%denom_u)) CS%diag_CS%denom_u(I,j,K) = sqrt(i_denom)
        end if

        ! to convert from the density gradient to the flux, flip the sign and multiply by
        ! kappa*dt
        dz_s_i(I,j) = -dz_s_i(I,j) * G%dxCu(I,j)**2 * ts_ratio * L_to_H**2

        dz_p_unlim = dz_s_i(I,j)

        ! limit slope based on adjacent layers
        ! dz_s_i has opposite sign to hdi_sig
        if (dz_s_i(I,j) < 0.) then
          ! hdi_sig positive -- left down, right up
          dz_s_i(I,j) = max(dz_s_i(I,j), -0.125 * min( &
               h(i,j,k) * G%areaT(i,j), &
               h(i+1,j,k-1) * G%areaT(i+1,j)) * G%IdyCu(I,j) * L_to_H)
        else
          ! hdi_sig negative -- left up, right down
          dz_s_i(I,j) = min(dz_s_i(I,j), 0.125 * min( &
               h(i,j,k-1) * G%areaT(i,j), &
               h(i+1,j,k) * G%areaT(i+1,j)) * G%IdyCu(I,j) * L_to_H)
        end if

        if (do_diag) then
          ! DIAG: limiting_density
          ! difference between the unlimited slope flux and the limited, across the participating adjacent cells
          if (allocated(CS%diag_CS%limiting_density)) then
            CS%diag_CS%limiting_density(i,j,K) = CS%diag_CS%limiting_density(i,j,K) + (dz_s_i(I,j) - dz_p_unlim)
            CS%diag_CS%limiting_density(i+1,j,K) = CS%diag_CS%limiting_density(i+1,j,K) + (dz_s_i(I,j) - dz_p_unlim)
          end if
        end if

        ! we also calculate the difference in pressure (interface position)
        dz_p_i(I,j) = (z_int(i+1,j,K) - z_int(i,j,K)) * G%dxCu(I,j) * ts_ratio * L_to_H
        dz_p_unlim = dz_p_i(I,j)
        ! dz_p_i positive => left is further down than right
        ! => move left up, right down

        if (dz_p_i(I,j) < 0.) then
          ! dz_p_i negative -- right up, left down
          dz_p_i(I,j) = max(dz_p_i(I,j), -0.125 * min( &
               h(i,j,k) * G%areaT(i,j), &
               h(i+1,j,k-1) * G%areaT(i+1,j)) * G%IdyCu(I,j) * L_to_H)
        else
          ! dz_p_i positive -- left up, right down
          dz_p_i(I,j) = min(dz_p_i(I,j), 0.125 * min( &
               h(i,j,k-1) * G%areaT(i,j), &
               h(i+1,j,k) * G%areaT(i+1,j)) * G%IdyCu(I,j) * L_to_H)
        end if

        if (do_diag) then
          ! DIAG: limiting_smoothing
          ! similar to limiting_density, but applied on the pressure (smoothing) term
          if (allocated(CS%diag_CS%limiting_smoothing)) then
            CS%diag_CS%limiting_smoothing(i,j,K) = CS%diag_CS%limiting_smoothing(i,j,K) + (dz_p_i(I,j) - dz_p_unlim)
            CS%diag_CS%limiting_smoothing(i+1,j,K) = CS%diag_CS%limiting_smoothing(i+1,j,K) + (dz_p_i(I,j) - dz_p_unlim)
          end if
        end if

        ! calculate and diagnose along-coordinate slope
        if (abs(i_denom) < eps) then
          slope = 1.0
        else
          slope = (hdi_sig_u + hdj_sig_u) / i_denom
        endif

        ! calculate physical slope
        hdi_sig_u = hdi_sig_phys(I,j,K)**2
        hdj_sig_u = 0.25 * ((hdj_sig_phys(i,J,K)**2 + hdj_sig_phys(i+1,J-1,K)**2) + &
             (hdj_sig_phys(i+1,J,K)**2 + hdj_sig_phys(i,J-1,K)**2))
        i_denom = hdi_sig_u + hdj_sig_u + dk_sig_u

        if (abs(i_denom) < eps) then
          ! unstratified limit
          phys_slope = 1.0
        else
          phys_slope = (hdi_sig_u + hdj_sig_u) / i_denom
        endif

        if (do_diag) then
          ! DIAG: coord_u
          if (allocated(CS%diag_CS%coord_u)) CS%diag_CS%coord_u(I,j,K) = slope
          ! DIAG: phys_u
          if (allocated(CS%diag_CS%phys_u)) CS%diag_CS%phys_u(I,j,K) = phys_slope
        end if

        ! use physical slope or not?
        if (CS%use_physical_slope) slope = phys_slope

        ! calculate weighting between density and pressure terms
        ! by a cutoff value on the local normalised stratification
        if (slope <= CS%slope_cutoff**2 .and. k > 2) then
          weight = 1.0 - CS%min_smooth ; weight2 = 0.
        else
          weight = 0.0 ; weight2 = 1.0 - CS%min_smooth
        endif

        ! override weights if required
        if (CS%alpha_rho >= 0.) then
          weight = CS%alpha_rho

          if (CS%alpha_p < 0.) then
            weight2 = 1.0 - CS%alpha_rho
          else
            weight2 = CS%alpha_p
          endif
        else if (CS%alpha_p >= 0.) then
          weight2 = CS%alpha_p
          weight = 1.0 - CS%alpha_p
        endif

        dz_s_i(I,j) = dz_s_i(I,j) * weight
        dz_p_i(I,j) = dz_p_i(I,j) * weight2

        ! combining density and pressure fluxes
        ! and re-apply limiter -- with a full cut-off this isn't necessary
        dz_i(I,j) = dz_s_i(I,j) + dz_p_i(I,j)
        if (dz_i(I,j) < 0.) then
          ! hdi_sig positive -- left down, right up
          dz_i(I,j) = max(dz_i(I,j), -0.125 * min( &
               h(i,j,k) * G%areaT(i,j), &
               h(i+1,j,k-1) * G%areaT(i+1,j)) * G%IdyCu(I,j) * L_to_H)
        else
          ! hdi_sig negative -- left up, right down
          dz_i(I,j) = min(dz_i(I,j), 0.125 * min( &
               h(i,j,k-1) * G%areaT(i,j), &
               h(i+1,j,k) * G%areaT(i+1,j)) * G%IdyCu(I,j) * L_to_H)
        end if
      end do
    end do

    ! v-points
    do J = G%JscB-1,G%JecB+1
      do i = G%isc-1,G%iec+1
        if (G%mask2dCv(i,J) < 0.5) then
          dz_j(i,J) = 0.
          dz_s_j(i,J) = 0.
          dz_p_j(i,J) = 0.
          cycle
        endif

        hdj_sig_v = hdj_sig(i,J,K)**2
        hdi_sig_v = 0.25 * ((hdi_sig(I,j,K)**2 + hdi_sig(I-1,j+1,K)**2) + &
             (hdi_sig(I,j+1,K)**2 + hdi_sig(I-1,j,K)**2))
        dk_sig_v = 0.5 * (dk_sig_int(i,j)**2 + dk_sig_int(i,j+1)**2)

        j_denom = hdj_sig_v + hdi_sig_v + dk_sig_v
        if (abs(j_denom) < eps) then
          dz_s_j(i,J) = 0.
        else
          dz_s_j(i,J) = hdj_sig(i,J,K) / sign(sqrt(j_denom), dk_sig_v)
        end if

        if (do_diag) then
          ! DIAG: slope_v
          if (allocated(CS%diag_CS%slope_v)) CS%diag_CS%slope_v(i,J,K) = dz_s_j(i,J)
          ! DIAG: denom_v
          if (allocated(CS%diag_CS%denom_v)) CS%diag_CS%denom_v(i,J,K) = sqrt(j_denom)
        end if

        ! dz_s_j beforehand is unitless (ratio of densities)
        dz_s_j(i,J) = -dz_s_j(i,J) * G%dyCv(i,J)**2 * ts_ratio * L_to_H**2
        ! dz_s_j is now [m2]

        dz_p_unlim = dz_s_j(i,J)

        ! density limiter
        ! dz_s_j [m2]
        if (dz_s_j(i,J) < 0.) then
          ! hdj_sig positive -- left down, right up
          dz_s_j(i,J) = max(dz_s_j(i,J), -0.125 * min( &
               h(i,j,k) * G%areaT(i,j), &
               h(i,j+1,k-1) * G%areaT(i,j+1)) * G%IdxCv(i,J) * L_to_H)
        else
          ! hdj_sig negative -- left up, right down
          dz_s_j(i,J) = min(dz_s_j(i,J), 0.125 * min( &
               h(i,j,k-1) * G%areaT(i,j), &
               h(i,j+1,k) * G%areaT(i,j+1)) * G%IdxCv(i,J) * L_to_H)
        end if

        if (do_diag) then
          ! DIAG: limiting_density
          ! see u-point loop for explanation
          if (allocated(CS%diag_CS%limiting_density)) then
            CS%diag_CS%limiting_density(i,j,K) = CS%diag_CS%limiting_density(i,j,K) + (dz_s_j(i,J) - dz_p_unlim)
            CS%diag_CS%limiting_density(i,j+1,K) = CS%diag_CS%limiting_density(i,j+1,K) + (dz_s_j(i,J) - dz_p_unlim)
          end if
        end if

        dz_p_j(i,J) = (z_int(i,j+1,K) - z_int(i,j,K)) * G%dyCv(i,J) * ts_ratio * L_to_H
        dz_p_unlim = dz_p_j(i,J)

        if (dz_p_j(i,J) < 0.) then
          dz_p_j(i,J) = max(dz_p_j(i,J), -0.125 * min( &
               h(i,j,k) * G%areaT(i,j), &
               h(i,j+1,k-1) * G%areaT(i,j+1)) * G%IdxCv(i,J) * L_to_H)
        else
          dz_p_j(i,J) = min(dz_p_j(i,J), 0.125 * min( &
               h(i,j,k-1) * G%areaT(i,j), &
               h(i,j+1,k) * G%areaT(i,j+1)) * G%IdxCv(i,J) * L_to_H)
        end if

        if (do_diag) then
          ! DIAG: limiting_smoothing
          if (allocated(CS%diag_CS%limiting_smoothing)) then
            CS%diag_CS%limiting_smoothing(i,j,K) = CS%diag_CS%limiting_smoothing(i,j,K) + (dz_p_j(i,J) - dz_p_unlim)
            CS%diag_CS%limiting_smoothing(i,j+1,K) = CS%diag_CS%limiting_smoothing(i,j+1,K) + (dz_p_j(i,J) - dz_p_unlim)
          end if
        end if

        ! diagnose along-coordinate slope
        if (abs(j_denom) < eps) then
          slope = 1.0
        else
          slope = (hdi_sig_v + hdj_sig_v) / j_denom
        endif

        hdj_sig_v = hdj_sig_phys(i,J,K)**2
        hdi_sig_v = 0.25 * ((hdi_sig_phys(I,j,K)**2 + hdi_sig_phys(I-1,j+1,K)**2) + &
             (hdi_sig_phys(I,j+1,K)**2 + hdi_sig_phys(I-1,j,K)**2))
        j_denom = hdi_sig_v + hdj_sig_v + dk_sig_v

        if (abs(j_denom) < eps) then
          phys_slope = 1.0
        else
          phys_slope = (hdi_sig_v + hdj_sig_v) / j_denom
        endif

        if (do_diag) then
          ! DIAG: coord_v
          if (allocated(CS%diag_CS%coord_v)) CS%diag_CS%coord_v(i,J,K) = slope
          ! DIAG: phys_v
          if (allocated(CS%diag_CS%phys_v)) CS%diag_CS%phys_v(i,J,K) = phys_slope
        end if

        if (CS%use_physical_slope) slope = phys_slope

        if (slope <= CS%slope_cutoff**2 .and. k > 2) then
          weight = 1.0 - CS%min_smooth ; weight2 = 0.
        else
          weight = 0.0 ; weight2 = 1.0 - CS%min_smooth
        endif

        ! override weights if required
        if (CS%alpha_rho >= 0.) then
          weight = CS%alpha_rho

          if (CS%alpha_p < 0.) then
            weight2 = 1.0 - CS%alpha_rho
          else
            weight2 = CS%alpha_p
          endif
        else if (CS%alpha_p >= 0.) then
          weight2 = CS%alpha_p
          weight = 1.0 - CS%alpha_p
        endif

        dz_s_j(i,J) = dz_s_j(i,J) * weight
        dz_p_j(i,J) = dz_p_j(i,J) * weight2

        dz_j(i,J) = dz_s_j(i,J) + dz_p_j(i,J)
        if (dz_j(i,J) < 0.) then
          ! hdj_sig positive -- left down, right up
          dz_j(i,J) = max(dz_j(i,J), -0.125 * min( &
               h(i,j,k) * G%areaT(i,j), &
               h(i,j+1,k-1) * G%areaT(i,j+1)) * G%IdxCv(i,J) * L_to_H)
        else
          ! hdj_sig negative -- left up, right down
          dz_j(i,J) = min(dz_j(i,J), 0.125 * min( &
               h(i,j,k-1) * G%areaT(i,j), &
               h(i,j+1,k) * G%areaT(i,j+1)) * G%IdxCv(i,J) * L_to_H)
        end if
      end do
    end do

    do j = G%jsc-1,G%jec+1
      do i = G%isc-1,G%iec+1
        ! prior to this point, dz_a and dz_p should be limited such that they
        ! can't cause any tangling. however, they may still lead to some grid-scale
        ! checkerboarding, so we reduce by another factor of 2
        dz_a(i,j,K) = 0.25 * G%IareaT(i,j) / L_to_H &
             * ((G%dyCu(I,j) * dz_i(I,j) - G%dyCu(I-1,j) * dz_i(I-1,j)) &
              + (G%dxCv(i,J) * dz_j(i,J) - G%dxCv(i,J-1) * dz_j(i,J-1)))

        ! apply the change in interface position due to this flux immediately
        z_int(i,j,K) = z_int(i,j,K) + dz_a(i,j,K)
      end do
    end do

    if (do_diag) then
      ! DIAG: disp_density
      if (allocated(CS%diag_CS%disp_density)) then
        do j = G%jsc-1,G%jec+1
          do i = G%isc-1,G%iec+1
            CS%diag_CS%disp_density(i,j,K) = 0.25 * G%IareaT(i,j) / L_to_H &
                 * ((G%dyCu(I,j) * dz_s_i(I,j) - G%dyCu(I-1,j) * dz_s_i(I-1,j)) &
                 +  (G%dxCv(i,J) * dz_s_j(i,J) - G%dxCv(i,J-1) * dz_s_j(i,J-1)))
          end do
        end do
      end if
      ! DIAG: disp_smoothing
      if (allocated(CS%diag_CS%disp_smoothing)) then
        do j = G%jsc-1,G%jec+1
          do i = G%isc-1,G%iec+1
            CS%diag_CS%disp_smoothing(i,j,K) = 0.25 * G%IareaT(i,j) / L_to_H &
                 * ((G%dyCu(I,j) * dz_p_i(I,j) - G%dyCu(I-1,j) * dz_p_i(I-1,j)) &
                 +  (G%dxCv(i,J) * dz_p_j(i,J) - G%dxCv(i,J-1) * dz_p_j(i,J-1)))
          end do
        end do
      end if
    end if

    ! calculate the z-smoothing fluxes and apply in a second step
    ! this lets us use a "barotropic" limiter, which should be much less
    ! restrictive than the layer-based one
    do j = G%jsc-1,G%jec+1
      do I = G%IscB-1,G%IecB+1
        if (G%mask2dCu(I,j) < 0.5) then
          dz_p_i(I,j) = 0.
          cycle
        endif

        dz_p_i(I,j) = (z_int(i+1,j,K) - z_int(i,j,K)) * G%dxCu(I,j) * ts_ratio * L_to_H
        ! dz_p_i positive => left is further down than right
        ! => move left up, right down

        ! XXX this becomes a barotropic limiter
        if (dz_p_i(I,j) < 0.) then
          ! dz_p_i negative -- right up, left down
          dz_p_i(I,j) = max(dz_p_i(I,j), -min( &
               (z_int(i,j,K) - z_int(i,j,nz+1)) * G%areaT(i,j), &
               (z_int(i+1,j,1) - z_int(i+1,j,K)) * G%areaT(i+1,j)) * G%IdyCu(I,j) * L_to_H)
        else
          ! dz_p_i positive -- left up, right down
          dz_p_i(I,j) = min(dz_p_i(I,j), min( &
               (z_int(i,j,1) - z_int(i,j,K)) * G%areaT(i,j), &
               (z_int(i+1,j,K) - z_int(i+1,j,nz+1)) * G%areaT(i+1,j)) * G%IdyCu(I,j) * L_to_H)
        end if
        dz_p_i(I,j) = dz_p_i(I,j) * CS%min_smooth
      end do
    end do

    do J = G%JscB-1,G%JecB+1
      do i = G%isc-1,G%iec+1
        if (G%mask2dCv(i,J) < 0.5) then
          dz_p_j(i,J) = 0.
          cycle
        endif

        dz_p_j(i,J) = (z_int(i,j+1,K) - z_int(i,j,K)) * G%dyCv(i,J) * ts_ratio * L_to_H

        if (dz_p_j(i,J) < 0.) then
          dz_p_j(i,J) = max(dz_p_j(i,J), -min( &
               (z_int(i,j,K) - z_int(i,j,nz+1)) * G%areaT(i,j), &
               (z_int(i,j+1,1) - z_int(i,j+1,K)) * G%areaT(i,j+1)) * G%IdxCv(i,J) * L_to_H)
        else
          dz_p_j(i,J) = min(dz_p_j(i,J), min( &
               (z_int(i,j,1) - z_int(i,j,K)) * G%areaT(i,j), &
               (z_int(i,j+1,K) - z_int(i,j+1,nz+1)) * G%areaT(i,j+1)) * G%IdxCv(i,J) * L_to_H)
        end if
        dz_p_j(i,J) = dz_p_j(i,J) * CS%min_smooth
      end do
    end do

    ! calculate flux due to barotropically-limited smoothing term
    do j = G%jsc-1,G%jec+1
      do i = G%isc-1,G%iec+1
        dz_p(i,j,K) = 0.5 * 0.25 * G%IareaT(i,j) / L_to_H &
             * ((G%dyCu(I,j) * dz_p_i(I,j) - G%dyCu(I-1,j) * dz_p_i(I-1,j)) &
              + (G%dxCv(i,J) * dz_p_j(i,J) - G%dxCv(i,J-1) * dz_p_j(i,J-1)))
      end do
    end do
  end do
  !$omp end do
  end block
  !$omp end parallel

  if (do_diag) then
    ! DIAG: disp_unlimited
    if (allocated(CS%diag_CS%disp_unlimited)) &
         CS%diag_CS%disp_unlimited(:,:,:) = dz_p(:,:,:)
  end if

  ts_ratio = dt / CS%restoring_timescale
  !$omp parallel do private(z_upd, z_col, i, j, k)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1
      dzInterface(i,j,:) = 0.
      ! for land points, leave interfaecs undisturbed (possibly doesn't matter)
      if (G%mask2dT(i,j) < 0.5) cycle

      ! calculate change in interface position due to restoring term
      ! z_int has already been updated by layer-limited fluxes
      ! add the barotropically limited flux too
      z_upd(:) = z_int(i,j,:) + dz_p(i,j,:)

      if (fCS%depth_of_time_filter_shallow > 0.) then
        ! build a z-star column
        call build_zstar_column(CS%zlike_CS, G%bathyT(i,j) * GV%Z_to_H, sum(h(i,j,:)), z_mean, zScale=GV%Z_to_H)

        ! filtered_grid_motion will fail if z_upd and z_mean are tangled with each other
        ! this basically means that every pair (z_upd(K),z_mean(K)) should be adjacent in a sorted list
        ! we can't (shouldn't?) change z_upd, so we can only tweak z_mean to ensure this condition is met
        ! restore with depth-dependent profile
        z_col(:) = z_mean(:)

        call filtered_grid_motion(fCS, CS%nk, nz, z_upd, z_col, dz_col)
        ! dz_col is the additional displacement on top of the interface displacement we already had
        dzInterface(i,j,2:nz) = dz_a(i,j,2:nz) + dz_p(i,j,2:nz) + dz_col(2:nz)
      else
        do K = 2,nz
          dz_r(i,j,K) = ts_ratio * (max(min(z_mean(K), z_upd(1)), z_upd(nz+1)) - z_upd(K)) &
               / (1.0 + ts_ratio)

          ! using filtered_grid_motion to obtain our dzInterface leads to a loss of precision:
          ! we effectively add the depth of the ocean and immediately subtract it out, losing
          ! about 4-5 orders of magnitude!
          ! instead, we just apply the calculated value directly
          ! combine both the layer-limited and barotropically-limited fluxes
          dzInterface(i,j,K) = dz_a(i,j,K) + dz_p(i,j,K)

          if (CS%restoring_timescale > 0.) &
               dzInterface(i,j,K) = dzInterface(i,j,K) + dz_r(i,j,K)
        enddo
      endif

      ! update h from previous steps in preparation for adjustment
      do k = 1,nz
        h_upd(i,j,k) = h(i,j,k) + (dzInterface(i,j,K) - dzInterface(i,j,K+1))
      enddo
    enddo
  enddo
end subroutine build_adapt_grid

end module coord_adapt
