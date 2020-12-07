!> Contains utility functions for filtering ALE grids
module filter_utils

use MOM_error_handler, only : MOM_error, FATAL

implicit none

type, public :: filter_CS
  !> Weight given to old coordinate when blending between new and old grids [nondim]
  !! Used only below depth_of_time_filter_shallow, with a cubic variation
  !! from zero to full effect between depth_of_time_filter_shallow and
  !! depth_of_time_filter_deep.
  real :: old_grid_weight = 0.

  !> Depth above which no time-filtering of grid is applied [H ~> m or kg m-2]
  real :: depth_of_time_filter_shallow = 0.

  !> Depth below which time-filtering of grid is applied at full effect [H ~> m or kg m-2]
  real :: depth_of_time_filter_deep = 0.
end type filter_CS

public filtered_grid_motion

contains

!> Returns the change in interface position motion after filtering and
!! assuming the top and bottom interfaces do not move.  The filtering is
!! a function of depth, and is applied as the integrated average filtering
!! over the trajectory of the interface.  By design, this code can not give
!! tangled interfaces provided that z_old and z_new are not already tangled.
subroutine filtered_grid_motion( CS, nkt, nk, z_old, z_new, dz_g )
  type(filter_CS),          intent(in)    :: CS !< Regridding control structure
  integer,                  intent(in)    :: nkt !< Number of cells in target grid
  integer,                  intent(in)    :: nk !< Number of cells in source grid
  real, dimension(nk+1),    intent(in)    :: z_old !< Old grid position [H ~> m or kg m-2]
  real, dimension(nkt+1), intent(in)    :: z_new !< New grid position before filtering [H ~> m or kg m-2]
  real, dimension(nkt+1), intent(inout) :: dz_g !< Change in interface positions
                                                !! including the effects of filtering [H ~> m or kg m-2]
  ! Local variables
  real :: sgn     ! The sign convention for downward [nondim].
  real :: dz_tgt  ! The target grid movement of the unfiltered grid [H ~> m or kg m-2]
  real :: zr1     ! The old grid position of an interface relative to the surface [H ~> m or kg m-2]
  real :: z_old_k ! The corrected position of the old grid [H ~> m or kg m-2]
  real :: Aq      ! A temporary variable related to the grid weights [nondim]
  real :: Bq      ! A temporary variable used in the linear term in the quadratic expression for the
                  ! filtered grid movement [H ~> m or kg m-2]
  real :: z0, dz0 ! Together these give the position of an interface relative to a reference hieght
                  ! that may be adjusted for numerical accuracy in a solver [H ~> m or kg m-2]
  real :: F0      ! An estimated grid movement  [H ~> m or kg m-2]
  real :: zs      ! The depth at which the shallow filtering timescale applies [H ~> m or kg m-2]
  real :: zd      ! The depth at which the deep filtering timescale applies [H ~> m or kg m-2]
  real :: dzwt    ! The depth range over which the transition in the filtering timescale occurs [H ~> m or kg m-2]
  real :: Idzwt   ! The Adcroft reciprocal of dzwt [H-1 ~> m-1 or m2 kg-1]
  real :: wtd     ! The weight given to the new grid when time filtering [nondim]
  real :: Iwtd    ! The inverse of wtd [nondim]
  real :: Int_zs  ! A depth integral of the weights in [H ~> m or kg m-2]
  real :: Int_zd  ! A depth integral of the weights in [H ~> m or kg m-2]
  real :: dInt_zs_zd ! The depth integral of the weights between the deep and shallow depths in [H ~> m or kg m-2]
! For debugging:
  real, dimension(nk+1) :: z_act ! The final grid positions after the filtered movement [H ~> m or kg m-2]
!  real, dimension(nk+1) :: ddz_g_s, ddz_g_d
  logical :: debug = .false.
  integer :: k

  if ((z_old(nk+1) - z_old(1)) * (z_new(nkt+1) - z_new(1)) < 0.0) then
    call MOM_error(FATAL, "filtered_grid_motion: z_old and z_new use different sign conventions.")
  elseif ((z_old(nk+1) - z_old(1)) * (z_new(nkt+1) - z_new(1)) == 0.0) then
    ! This is a massless column, so do nothing and return.
    do k=1,nkt+1 ; dz_g(k) = 0.0 ; enddo ; return
  elseif ((z_old(nk+1) - z_old(1)) + (z_new(nkt+1) - z_new(1)) > 0.0) then
    sgn = 1.0
  else
    sgn = -1.0
  endif

  if (debug) then
    do k=2,nkt+1
      if (sgn*(z_new(k)-z_new(k-1)) < -5e-16*(abs(z_new(k))+abs(z_new(k-1))) ) &
        call MOM_error(FATAL, "filtered_grid_motion: z_new is tangled.")
    enddo
    do k=2,nk+1
      if (sgn*(z_old(k)-z_old(k-1)) < -5e-16*(abs(z_old(k))+abs(z_old(k-1))) ) &
        call MOM_error(FATAL, "filtered_grid_motion: z_old is tangled.")
    enddo
    ! ddz_g_s(:) = 0.0 ; ddz_g_d(:) = 0.0
  endif

  zs = CS%depth_of_time_filter_shallow
  zd = CS%depth_of_time_filter_deep
  wtd = 1.0 - CS%old_grid_weight
  Iwtd = 1.0 / wtd

  dzwt = (zd - zs)
  Idzwt = 0.0 ; if (abs(zd - zs) > 0.0) Idzwt = 1.0 / (zd - zs)
  dInt_zs_zd = 0.5*(1.0 + Iwtd) * (zd - zs)
  Aq = 0.5*(Iwtd - 1.0)

  dz_g(1) = 0.0
  z_old_k = z_old(1)
  do k = 2,nkt+1
    if (k<=nk+1) z_old_k = z_old(k) ! This allows for virtual z_old interface at bottom of the model
    ! zr1 is positive and increases with depth, and dz_tgt is positive downward.
    dz_tgt = sgn*(z_new(k) - z_old_k)
    zr1 = sgn*(z_old_k - z_old(1))

    !   First, handle the two simple and common cases that do not pass through
    ! the adjustment rate transition zone.
    if ((zr1 > zd) .and. (zr1 + wtd * dz_tgt > zd)) then
      dz_g(k) = sgn * wtd * dz_tgt
    elseif ((zr1 < zs) .and. (zr1 + dz_tgt < zs)) then
      dz_g(k) = sgn * dz_tgt
    else
      ! Find the new value by inverting the equation
      !   integral(0 to dz_new) Iwt(z) dz = dz_tgt
      ! This is trivial where Iwt is a constant, and agrees with the two limits above.

      ! Take test values at the transition points to figure out which segment
      ! the new value will be found in.
      if (zr1 >= zd) then
        Int_zd = Iwtd*(zd - zr1)
        Int_zs = Int_zd - dInt_zs_zd
      elseif (zr1 <= zs) then
        Int_zs = (zs - zr1)
        Int_zd = dInt_zs_zd + (zs - zr1)
      else
!        Int_zd = (zd - zr1) * (Iwtd + 0.5*(1.0 - Iwtd) * (zd - zr1) / (zd - zs))
        Int_zd = (zd - zr1) * (Iwtd*(0.5*(zd+zr1) - zs) + 0.5*(zd - zr1)) * Idzwt
        Int_zs = (zs - zr1) * (0.5*Iwtd * ((zr1 - zs)) + (zd - 0.5*(zr1+zs))) * Idzwt
        ! It has been verified that  Int_zs = Int_zd - dInt_zs_zd to within roundoff.
      endif

      if (dz_tgt >= Int_zd) then ! The new location is in the deep, slow region.
        dz_g(k) = sgn * ((zd-zr1) + wtd*(dz_tgt - Int_zd))
      elseif (dz_tgt <= Int_zs) then ! The new location is in the shallow region.
        dz_g(k) = sgn * ((zs-zr1) + (dz_tgt - Int_zs))
      else  ! We need to solve a quadratic equation for z_new.
        ! For accuracy, do the integral from the starting depth or the nearest
        ! edge of the transition region.  The results with each choice are
        ! mathematically equivalent, but differ in roundoff, and this choice
        ! should minimize the likelihood of inadvertently overlapping interfaces.
        if (zr1 <= zs) then ; dz0 = zs-zr1 ; z0 = zs ; F0 = dz_tgt - Int_zs
        elseif (zr1 >= zd) then ; dz0 = zd-zr1 ; z0 = zd ; F0 = dz_tgt - Int_zd
        else ; dz0 = 0.0 ; z0 = zr1 ; F0 = dz_tgt ; endif

        Bq = (dzwt + 2.0*Aq*(z0-zs))
        ! Solve the quadratic: Aq*(zn-z0)**2 + Bq*(zn-z0) - F0*dzwt = 0
        ! Note that b>=0, and the two terms in the standard form cancel for the right root.
        dz_g(k) = sgn * (dz0 + 2.0*F0*dzwt / (Bq + sqrt(Bq**2 + 4.0*Aq*F0*dzwt) ))

!       if (debug) then
!         dz0 = zs-zr1 ; z0 = zs ; F0 = dz_tgt - Int_zs ; Bq = (dzwt + 2.0*Aq*(z0-zs))
!         ddz_g_s(k) = sgn * (dz0 + 2.0*F0*dzwt / (Bq + sqrt(Bq**2 + 4.0*Aq*F0*dzwt) )) - dz_g(k)
!         dz0 = zd-zr1 ; z0 = zd ; F0 = dz_tgt - Int_zd ; Bq = (dzwt + 2.0*Aq*(z0-zs))
!         ddz_g_d(k) = sgn * (dz0 + 2.0*F0*dzwt / (Bq + sqrt(Bq**2 + 4.0*Aq*F0*dzwt) )) - dz_g(k)
!
!         if (abs(ddz_g_s(k)) > 1e-12*(abs(dz_g(k)) + abs(dz_g(k)+ddz_g_s(k)))) &
!           call MOM_error(WARNING, "filtered_grid_motion: Expect z_output to be tangled (sc).")
!         if (abs(ddz_g_d(k) - ddz_g_s(k)) > 1e-12*(abs(dz_g(k)+ddz_g_d(k)) + abs(dz_g(k)+ddz_g_s(k)))) &
!           call MOM_error(WARNING, "filtered_grid_motion: Expect z_output to be tangled.")
!       endif
      endif

    endif
  enddo
 !dz_g(nkt+1) = 0.0

  if (debug) then
    z_old_k = z_old(1)
    do k=1,nkt+1
      if (k<=nk+1) z_old_k = z_old(k) ! This allows for virtual z_old interface at bottom of the model
      z_act(k) = z_old_k + dz_g(k)
    enddo
    do k=2,nkt+1
      if (sgn*((z_act(k))-z_act(k-1)) < -1e-15*(abs(z_act(k))+abs(z_act(k-1))) ) &
        call MOM_error(FATAL, "filtered_grid_motion: z_output is tangled.")
    enddo
  endif

end subroutine filtered_grid_motion

end module filter_utils
