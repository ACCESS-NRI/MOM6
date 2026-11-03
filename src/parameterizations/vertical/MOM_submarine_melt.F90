! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Distributes coupled submarine melt (srunoff) over the upper part of the water column.
module MOM_submarine_melt

use MOM_constants,     only : hlf
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public distribute_submarine_melt, submarine_melt_init

contains

!> Distribute the mass and (optionally) heat content of submarine melt (srunoff) uniformly by
!! thickness over the upper fluxes%srunoff_depth_max of a column. If fluxes%srunoff_depth_max is 0,
!! the distribution depth for a column instead falls back to that column's own top layer thickness,
!! so that srunoff is deposited entirely within the top layer, unless that layer is thinner than
!! minimum_forcing_depth, in which case srunoff is spread over the upper minimum_forcing_depth.
!! This routine is called from both:
!! - applyBoundaryFluxesInOut, to update thickness and temperature and dilute salinity.
!! - applyTracerBoundaryFluxesInOut, to dilute a tracer concentration.
subroutine distribute_submarine_melt(G, GV, dt, fluxes, j, minimum_forcing_depth, h2d, C2d, T2d, C_p, &
                                     g_Hconv2, cTKE, dSV_dT, dSV_dS, dHeat_total, dTempxPmE_total)
  type(ocean_grid_type),            intent(in)    :: G              !< Grid structure
  type(verticalGrid_type),          intent(in)    :: GV             !< Ocean vertical grid structure
  real,                             intent(in)    :: dt             !< Time step [T ~> s]
  type(forcing),                    intent(in)    :: fluxes         !< Surface fluxes, including srunoff
  integer,                          intent(in)    :: j              !< j-index to work on
  real,                             intent(in)    :: minimum_forcing_depth !< The smallest depth over
                                                                    !! which srunoff is distributed when
                                                                    !! fluxes%srunoff_depth_max is 0
                                                                    !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), &
                                    intent(inout) :: h2d            !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), &
                                    intent(inout) :: C2d            !< A concentration diluted by srunoff's
                                                                    !! fresh water, e.g. salinity or a tracer's own
                                                                    !! concentration [CU]
  real, dimension(SZI_(G),SZK_(GV)), &
                          optional, intent(inout) :: T2d            !< Layer temperatures [C ~> degC]
  real,                   optional, intent(in)    :: C_p            !< The heat capacity of seawater
                                                                    !! [Q C-1 ~> J kg-1 degC-1]
  real,                   optional, intent(in)    :: g_Hconv2       !< A conversion factor for use in the
                                                                    !! TKE calculation.
                                                                    !! [Z3 R2 T-2 H-2 ~> kg2 m-5 s-2 or m s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                          optional, intent(inout) :: cTKE           !< Turbulent kinetic energy requirement
                                                                    !! to mix forcing through each layer.
                                                                    !! [R Z3 T-2 ~> J m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                          optional, intent(in)    :: dSV_dT         !< Partial derivative of specific
                                                                    !! volume with potential temperature.
                                                                    !! [R-1 C-1 ~> m3 kg-1 degC-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                          optional, intent(in)    :: dSV_dS         !< Partial derivative of specific
                                                                    !! volume with salinity.
                                                                    !! [R-1 S-1 ~> m3 kg-1 ppt-1]
  real, dimension(SZI_(G)), &
                          optional, intent(out)   :: dHeat_total    !< The total heat added to each
                                                                    !! column by srunoff.
                                                                    !! [Q H ~> J m-2]
  real, dimension(SZI_(G)), &
                          optional, intent(out)   :: dTempxPmE_total !< The total srunoff mass
                                                                    !! added to each column,
                                                                    !! weighted by ambient
                                                                    !! temperature.
                                                                    !! [C H ~> degC m or degC kg m-2]

  ! Local variables
  real :: dM_tot      ! Total mass of srunoff to be distributed this timestep [H ~> m or kg m-2]
  real :: dM_k        ! The portion of dM_tot applied to a layer [H ~> m or kg m-2]
  real :: dHeat       ! The heat added to the column so far, accumulated for diagnostics [Q H ~> J m-2]
  real :: dTempxPmE   ! The ambient-temperature-weighted mass added to the column so far,
                      ! accumulated for diagnostics [C H ~> degC m or degC kg m-2]
  real :: depth       ! The requested depth over which srunoff is distributed [H ~> m or kg m-2]
  real :: total_h     ! The full thickness of the water column [H ~> m or kg m-2]
  real :: eff_depth   ! The depth actually spanned in this column, capped at total_h so that the
                      ! full amount of srunoff is always distributed even in a column shallower
                      ! than depth [H ~> m or kg m-2]
  real :: h_above     ! Cumulative thickness above the current layer, before adding srunoff [H ~> m or kg m-2]
  real :: h_in_range  ! The portion of a layer's thickness within eff_depth [H ~> m or kg m-2]
  real :: T_eff       ! The effective temperature carried by srunoff mass added to a layer [C ~> degC]
  real :: hOld        ! The layer thickness before adding srunoff mass [H ~> m or kg m-2]
  real :: I_eff_depth ! The inverse of eff_depth [H-1 ~> m-1 or m2 kg-1]
  logical :: do_heat       ! If true, T2d and C_p are present and get updated
  logical :: do_energetics ! If true, cTKE and its supporting arguments are also present
  integer :: i, k, is, ie, nz

  is = G%isc ; ie = G%iec ; nz = GV%ke
  do_heat = present(T2d) .and. present(C_p)
  do_energetics = do_heat .and. present(cTKE) .and. present(dSV_dT) .and. present(dSV_dS) &
                  .and. present(g_Hconv2)

  do i=is,ie
    if (present(dHeat_total)) dHeat_total(i) = 0.0
    if (present(dTempxPmE_total)) dTempxPmE_total(i) = 0.0

    if ((G%mask2dT(i,j) > 0.0) .and. (fluxes%srunoff(i,j) > 0.0)) then

      dM_tot = GV%RZ_to_H * dt * fluxes%srunoff(i,j)

      ! When fluxes%srunoff_depth_max is 0 (the default), fall back to this column's own top layer
      ! thickness as the distribution depth, so srunoff is deposited within the top layer. If the
      ! top layer is thinner than minimum_forcing_depth, srunoff is instead spread over the upper
      ! minimum_forcing_depth, as for other surface fluxes, so a near-vanished layer is not
      ! swamped by the melt.
      depth = fluxes%srunoff_depth_max
      if (depth <= 0.0) depth = max(h2d(i,1), minimum_forcing_depth, GV%H_subroundoff)

      ! Cap the depth actually used at the full column thickness, so that a column shallower
      ! than the requested depth still receives all of dM_tot.
      ! Note, with the current implementation, if the column is zero thickness any srunoff
      ! is just dropped.
      total_h = 0.0
      do k=1,nz ; total_h = total_h + h2d(i,k) ; enddo
      eff_depth = min(depth, total_h)
      I_eff_depth = 0.0 ; if (eff_depth > 0.0) I_eff_depth = 1.0 / eff_depth

      h_above = 0.0
      dHeat = 0.0
      dTempxPmE = 0.0

      ! Spread dM_tot over the upper eff_depth of the column, in proportion to how much of each
      ! layer's thickness falls within that range.
      do k=1,nz
        h_in_range = max(0.0, min(h_above + h2d(i,k), eff_depth) - h_above)
        h_above = h_above + h2d(i,k)
        if (h_in_range > 0.0) then
          dM_k = dM_tot * (h_in_range * I_eff_depth)
          hOld = h2d(i,k)

          if (do_heat) then
            T_eff = T2d(i,k)
            if (fluxes%srunoff_latent_heat) &
              T_eff = T_eff - fluxes%latent_heat_fusion/C_p
            dHeat = dHeat + C_p*dM_k*T_eff
            dTempxPmE = dTempxPmE + dM_k*T2d(i,k)

            if (do_energetics) then
              ! Calculate the energy required to mix the srunoff added to this layer.
              cTKE(i,j,k) = cTKE(i,j,k) + 0.5*g_Hconv2*(hOld*dM_k) * &
                 ((T2d(i,k) - T_eff) * dSV_dT(i,j,k) + C2d(i,k) * dSV_dS(i,j,k))
            endif
          endif

          h2d(i,k) = hOld + dM_k
          if (do_heat) T2d(i,k) = (hOld*T2d(i,k) + dM_k*T_eff) / h2d(i,k)
          C2d(i,k) = (hOld*C2d(i,k)) / h2d(i,k)
        endif
      enddo ! k-loop for layer distribution

      if (present(dHeat_total)) dHeat_total(i) = dHeat
      if (present(dTempxPmE_total)) dTempxPmE_total(i) = dTempxPmE

    endif
  enddo

end subroutine distribute_submarine_melt

!> Initialize the parameters that control the distribution of submarine melt, and store them in
!! the forcing type, where they are used by distribute_submarine_melt.
subroutine submarine_melt_init(param_file, GV, US, fluxes)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(forcing),           intent(inout) :: fluxes     !< A structure containing pointers to
                                                       !! thermodynamic forcing fields, in which
                                                       !! the submarine melt parameters are set.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40) :: mdl = "MOM_submarine_melt" ! This module's name.

  ! These are stored in the forcing type so that every routine that applies srunoff, including
  ! the tracer column physics, distributes it in the same way.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SRUNOFF_LATENT_HEAT", fluxes%srunoff_latent_heat, &
                 "If true, submarine melt (srunoff) enters each layer at an effective temperature "//&
                 "equal to the layer temperature minus LATENT_HEAT_FUSION/C_p, so that melting the ice "//&
                 "takes latent heat from the layers the meltwater enters. This is a simplified form "//&
                 "of the Gade (1979) meltwater mixing line, neglecting the heat needed to warm the "//&
                 "ice to the freezing point and the meltwater from the freezing point to the ambient "//&
                 "temperature. If false, srunoff enters the ocean at the ambient temperature.", &
                 default=.false.)
  call get_param(param_file, mdl, "SRUNOFF_DEPTH_MAX", fluxes%srunoff_depth_max, &
                 "The maximum depth over which the mass and heat from coupled submarine melt are spread, "//&
                 "distributed uniformly by thickness. If 0, srunoff is applied within the top "//&
                 "layer, or over MINIMUM_FORCING_DEPTH if the top layer is thinner than that.", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", fluxes%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf, scale=US%J_kg_to_Q, &
                 do_not_log=.not.fluxes%srunoff_latent_heat)

end subroutine submarine_melt_init

end module MOM_submarine_melt
