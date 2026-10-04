!> Contains routines for handling the FMS coupler_bc_type structures that carry the air-sea
!! tracer fluxes and the fields they are calculated from.  These are the fields, e.g., that the
!! generic tracer packages are designed to exchange through the FMS coupler, as distinct from
!! everything else in this driver, which is exchanged through NUOPC.

module MOM_cap_fms_coupler_bcs

use MOM_constants,             only: wtmair, rdgas, vonkarm
use MOM_coupler_types,         only: coupler_1d_bc_type, coupler_2d_bc_type, coupler_type_spawn
use MOM_coupler_types,         only: coupler_type_num_bcs
use MOM_coupler_types,         only: ind_flux, ind_deltap, ind_kw, ind_flux0
use MOM_coupler_types,         only: ind_pcair, ind_u10, ind_psurf
use MOM_coupler_types,         only: ind_alpha, ind_csurf, ind_sc_no
use MOM_coupler_types,         only: ind_runoff, ind_deposition
use MOM_coupler_types,         only: coupler_type_register_restart_fields
use MOM_coupler_types,         only: coupler_type_get_bc, coupler_type_get_field
use MOM_coupler_types,         only: coupler_type_set_field
use MOM_coupler_types,         only: coupler_type_set_diags
use MOM_coupler_types,         only: coupler_type_send_data
use MOM_data_override,         only: data_override
use MOM_domains,               only: domain2d, get_domain_extent
use MOM_error_handler,         only: MOM_error, FATAL
use MOM_file_parser,           only: param_file_type
use MOM_grid,                  only: ocean_grid_type
use MOM_restart,               only: MOM_restart_CS, restart_init, restart_init_end, restart_end
use MOM_restart,               only: restore_state, save_restart, determine_is_new_run
use MOM_time_manager,          only: time_type
use MOM_tracer_flow_control,   only: call_tracer_flux_init

! This FMS module has no MOM6 infrastructure wrapper.  The solo drivers supply a stub of it that
! provides only aof_set_coupler_flux, and adding stubs of the routines used here for the sake of
! the NUOPC cap seems backwards.
use atmos_ocean_fluxes_mod,    only: atmos_ocean_type_fluxes_init, atmos_ocean_fluxes_init

implicit none; private

! Public member functions
public :: coupler_bcs_init, coupler_bcs_end, coupler_bcs_spawn
public :: coupler_bcs_data_override, coupler_bcs_update_fluxes
public :: coupler_bcs_register_restarts, coupler_bcs_save_restart, coupler_bcs_restore
public :: coupler_bcs_get_cmeps_name

character(len=*), parameter      :: mod_name = 'mom_cap_fms_coupler_bcs'
real, parameter                  :: epsln=1.0e-30

!> The name root of the restart files that hold the FMS coupler_bc_type ocean surface fields.
!! Note, the ocean_restart_file entries in the field table are not used by this module.
character(len=*), parameter      :: restart_file_root = 'coupler_bc'

!> The control structure for the air-sea tracer fluxes, holding the FMS coupler_bc_types that
!! describe them and the restarts of the ocean surface fields that the flux calculation consumes.
type, public :: coupler_bcs_CS ; private
  type(coupler_1d_bc_type) :: gas_fields_atm ! tracer fields in atm
      !< Structure containing atmospheric surface variables that are used in the
      !! calculation of the atmosphere-ocean gas fluxes, as well as parameters
      !! regulating these fluxes. The fields in this structure are never actually
      !! set, but the structure is used for initialisation of components and to
      !! spawn other structure whose fields are set.
  type(coupler_1d_bc_type), public :: gas_fields_ocn ! tracer fields atop the ocean
      !< Structure containing ocean surface variables that are used in the
      !! calculation of the atmosphere-ocean gas fluxes, as well as parameters
      !! regulating these fluxes. The fields in this structure are never actually
      !! set, but the structure is used for initialisation of components and to
      !! spawn other structure whose fields are set.
  type(coupler_1d_bc_type) :: gas_fluxes ! tracer fluxes between the atm and ocean
      !< A structure for exchanging gas or tracer fluxes between the atmosphere and
      !! ocean, defined by the field table, as well as a place holder of
      !! intermediate calculations, such as piston velocities, and parameters that
      !! impact the fluxes. The fields in this structure are never actually set,
      !! but the structure is used for initialisation of components and to spawn
      !! other structure whose fields are set.

  type(MOM_restart_CS), pointer :: restart_CSp => NULL()
      !< A pointer to the control structure that the ocean surface fields used in the flux
      !! calculation are registered with, or NULL if there are none to restart.
end type coupler_bcs_CS

contains

!> \brief Define the FMS coupler boundary conditions from the field table.  This must be called
!! before the tracer packages are registered, because a package that registers a coupler flux of
!! its own needs the flux types that this defines to exist already.
!!
!! Copied and adapted slightly from the FMScoupler gas_exchange_init routine at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/7761886/full/flux_exchange.F90#L626.
subroutine coupler_bcs_init(CS)
  type(coupler_bcs_CS), pointer :: CS !< A pointer to the control structure for the FMS coupler
      !! boundary conditions, which is allocated here if any are in use.

  if (associated(CS)) return
  allocate(CS)

  call atmos_ocean_type_fluxes_init()
  ! FMScoupler's gas_exchange_init calls ocean_model_flux_init, but that is just a wrapper on
  ! call_tracer_flux_init and calling the latter directly avoids a circular dependency.
  call call_tracer_flux_init()
  call atmos_ocean_fluxes_init(CS%gas_fluxes, CS%gas_fields_atm, CS%gas_fields_ocn)

  ! No enabled tracer package registered a coupler flux, so there is nothing for this module to
  ! do. Return an unassociated control structure.
  if (coupler_type_num_bcs(CS%gas_fluxes) <= 0) then
    deallocate(CS) ; CS => NULL()
  endif

end subroutine coupler_bcs_init

!> Release the memory that is associated with the FMS coupler boundary conditions.
subroutine coupler_bcs_end(CS)
  type(coupler_bcs_CS), pointer :: CS !< A pointer to the control structure for the FMS coupler
      !! boundary conditions, which is deallocated here

  if (.not.associated(CS)) return

  ! The value arrays for the CS%gas_* coupler types are never allocated so calling the FMS
  ! destructor on them here gives an error.
  if (associated(CS%restart_CSp)) call restart_end(CS%restart_CSp)

  deallocate(CS) ; CS => NULL()

end subroutine coupler_bcs_end

!> Spawn the FMS coupler_bc_types that hold the air-sea tracer fluxes and the atmospheric fields
!! that they are calculated from, and register their diagnostics. Both are left unset when the
!! FMS coupler boundary conditions are not in use.
subroutine coupler_bcs_spawn(CS, fluxes, atm_fields, axes, Time, isc, iec, jsc, jec)
  type(coupler_bcs_CS),     pointer       :: CS !< The control structure for the FMS coupler
                                                !! boundary conditions
  type(coupler_2d_bc_type), intent(inout) :: fluxes !< The structure that is spawned to hold the
                                                !! air-sea tracer fluxes
  type(coupler_2d_bc_type), intent(inout) :: atm_fields !< The structure that is spawned to hold
                                                !! the atmospheric fields used to calculate them
  integer, dimension(2),    intent(in)    :: axes !< The handles of the horizontal axes that the
                                                !! diagnostics of the spawned structures use
  type(time_type),          intent(in)    :: Time !< The model time at which the diagnostics start
  integer,                  intent(in)    :: isc !< The start i-index of the computational domain
  integer,                  intent(in)    :: iec !< The end i-index of the computational domain
  integer,                  intent(in)    :: jsc !< The start j-index of the computational domain
  integer,                  intent(in)    :: jec !< The end j-index of the computational domain

  if (.not.associated(CS)) return

  ! The param arrays are read by the flux calculation, but are not transferred by default when
  ! spawning.
  call coupler_type_spawn(CS%gas_fluxes, fluxes, (/isc,isc,iec,iec/), (/jsc,jsc,jec,jec/), &
                          suffix='_ice_ocn', copy_param=.true.)
  call coupler_type_spawn(CS%gas_fields_atm, atm_fields, (/isc,isc,iec,iec/), &
                          (/jsc,jsc,jec,jec/), suffix='_atm')

  call coupler_type_set_diags(fluxes, "ocean_flux", axes, Time)
  call coupler_type_set_diags(atm_fields, "atmos_sfc", axes, Time)

end subroutine coupler_bcs_spawn

!> Potentially override the FMS coupler_bc_type air-sea tracer fluxes and the atmospheric
!! fields that they are calculated from, using the component name 'OCN'.
subroutine coupler_bcs_data_override(fluxes, atm_fields, Time)
  type(coupler_2d_bc_type), intent(inout) :: fluxes !< The air-sea tracer fluxes
  type(coupler_2d_bc_type), intent(inout) :: atm_fields !< The atmospheric fields that the fluxes
                                                !! are calculated from
  type(time_type),          intent(in)    :: Time !< The model time at which the fields apply

  ! Local variables
  real, dimension(:,:), pointer :: values ! The data of the field being overridden [various]
  character(len=128) :: name  ! The name of the field being overridden
  logical :: overridden       ! True if the field was overridden
  integer :: nfields          ! The number of fields in a boundary condition
  integer :: m, n             ! The indices of a boundary condition and of a field within it

  ! coupler_type_data_override is not used here because it does not set the override flag on the
  ! fields that it overrides, and mom_import and the flux calculation both use that flag to tell
  ! which of these fields the data_table has already provided.
  do n=1,coupler_type_num_bcs(atm_fields)
    call coupler_type_get_bc(atm_fields, n, num_fields=nfields)
    do m=1,nfields
      call coupler_type_get_field(atm_fields, n, m, values=values, name=name)
      call data_override('OCN', trim(name), values, Time, override=overridden)
      call coupler_type_set_field(atm_fields, n, m, override=overridden)
    enddo
  enddo ! n- and m-loops over the atmospheric fields and their components
  do n=1,coupler_type_num_bcs(fluxes)
    call coupler_type_get_bc(fluxes, n, num_fields=nfields)
    do m=1,nfields
      call coupler_type_get_field(fluxes, n, m, values=values, name=name)
      call data_override('OCN', trim(name), values, Time, override=overridden)
      call coupler_type_set_field(fluxes, n, m, override=overridden)
    enddo
  enddo ! n- and m-loops over the tracer fluxes and their components

end subroutine coupler_bcs_data_override

!> Calculate the FMS coupler_bc_type air-sea tracer fluxes for the current coupling step and send
!! the diagnostics of the fluxes and of the atmospheric fields that they are calculated from.
subroutine coupler_bcs_update_fluxes(fluxes, atm_fields, sfc_fields, ice_fraction, domain, Time)
  type(coupler_2d_bc_type), intent(inout) :: fluxes !< The air-sea tracer fluxes that are calculated
  type(coupler_2d_bc_type), intent(inout) :: atm_fields !< The atmospheric fields that the fluxes
                                                !! are calculated from
  type(coupler_2d_bc_type), intent(in)    :: sfc_fields !< The ocean surface fields that the fluxes
                                                !! are calculated from
  real, dimension(:,:),     intent(in)    :: ice_fraction !< The fraction of each cell that is
                                                !! covered by sea ice [nondim]
  type(domain2d),           intent(in)    :: domain !< The domain of the ocean surface fields
  type(time_type),          intent(in)    :: Time !< The model time at which the fluxes apply

  ! Local variables
  integer :: isc, iec, jsc, jec  ! The computational domain index bounds

  call get_domain_extent(domain, isc, iec, jsc, jec)
  call atmos_ocean_fluxes_calc(atm_fields, sfc_fields, fluxes, ice_fraction, isc, iec, jsc, jec)

  call coupler_type_send_data(atm_fields, Time)
  call coupler_type_send_data(fluxes, Time)

end subroutine coupler_bcs_update_fluxes

!> Register for restart the ocean surface fields of the FMS coupler boundary conditions that the
!! air-sea tracer flux calculation consumes.
subroutine coupler_bcs_register_restarts(CS, sfc_fields, param_file)
  type(coupler_bcs_CS),     pointer       :: CS !< The control structure for the FMS coupler
                                                !! boundary conditions
  type(coupler_2d_bc_type), intent(inout) :: sfc_fields !< The ocean surface fields to register
  type(param_file_type),    intent(in)    :: param_file !< A structure to parse for run-time
                                                !! parameters

  if (.not.associated(CS)) return

  ! The fields in a coupler type are on the input grid, so this control structure applies no
  ! index rotation.
  call restart_init(param_file, CS%restart_CSp, restart_root=restart_file_root, turns=0)
  call coupler_type_register_restart_fields(sfc_fields, CS%restart_CSp)
  call restart_init_end(CS%restart_CSp)

end subroutine coupler_bcs_register_restarts

!> Read the ocean surface fields of the FMS coupler boundary conditions from their restart file.
!! Nothing is read on a new run.
subroutine coupler_bcs_restore(CS, G, input_filename, restart_input_dir)
  type(coupler_bcs_CS),  pointer    :: CS !< The control structure for the FMS coupler boundary
                                          !! conditions
  type(ocean_grid_type), intent(in) :: G  !< The ocean's grid structure
  character(len=*),      intent(in) :: input_filename !< The list of ocean restart file names, or
                                          !! a single character indicating how they are named
  character(len=*),      intent(in) :: restart_input_dir !< The directory holding the restart files

  ! Local variables
  type(time_type) :: restart_time  ! The time recorded in the restart file

  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return

  ! Do not restore on a new run, following the test that MOM_initialize_state uses.
  if (determine_is_new_run(input_filename, restart_input_dir, G, CS%restart_CSp)) return

  call restore_state(input_filename, restart_input_dir, restart_time, G, CS%restart_CSp)

end subroutine coupler_bcs_restore

!> Write the ocean surface fields of the FMS coupler boundary conditions to their restart file.
subroutine coupler_bcs_save_restart(CS, G, Time, directory, name_prefix, restartname, num_rest_files)
  type(coupler_bcs_CS),  pointer          :: CS !< The control structure for the FMS coupler
                                             !! boundary conditions
  type(ocean_grid_type), intent(inout)    :: G !< The ocean's grid structure
  type(time_type),       intent(in)       :: Time !< The current model time
  character(len=*),      intent(in)       :: directory !< The directory into which to write the
                                             !! restart files
  character(len=*), optional, intent(in)  :: name_prefix !< If present, a prefix that is prepended
                                             !! to the name of the restart files
  character(len=*), optional, intent(out) :: restartname !< The name root shared by the restart
                                             !! files that are written
  integer,          optional, intent(out) :: num_rest_files !< The number of restart files written

  ! Local variables
  character(len=240) :: filename  ! The name root of the restart files that are written

  if (present(num_rest_files)) num_rest_files = 0
  if (present(restartname)) restartname = ""
  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return

  filename = restart_file_root
  if (present(name_prefix)) filename = trim(name_prefix)//"."//trim(filename)
  if (present(restartname)) restartname = trim(filename)

  call save_restart(directory, Time, G, CS%restart_CSp, filename=filename, &
                    num_rest_files=num_rest_files)

end subroutine coupler_bcs_save_restart

!> \brief Return the CMEPS standard_name of the field that provides the input of an FMS coupler
!! boundary condition, or an empty string if the coupler does not provide it.  Depending on the
!! flux type, that input is an atmospheric concentration, an atmospheric deposition flux or a
!! runoff flux.
!!
!! The CMEPS standard names are a fixed external vocabulary, so a boundary condition can only be
!! driven from the coupler if the mediator has been set up to provide its field.  Any other
!! boundary condition has to have its input field, or its flux, supplied from the data_table, and
!! mom_import issues a fatal error if neither is the case.
function coupler_bcs_get_cmeps_name(name)
  character(len=64)            :: coupler_bcs_get_cmeps_name !< CMEPS standard_name
  character(len=*), intent(in) :: name !< FMS coupler boundary condition name

  ! Add other FMS coupler boundary conditions that the coupler can drive here.
  select case (trim(name))
    case ('co2_flux') ; coupler_bcs_get_cmeps_name = "Sa_co2prog"
    case default      ; coupler_bcs_get_cmeps_name = ""
  end select
end function coupler_bcs_get_cmeps_name

!> \brief Calculate the FMS coupler_bc_type ocean tracer fluxes. Units should be mol/m^2/s.
!! Upward flux is positive.
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
!! and subsequently modified in the following ways:
!! - Operate on 2D inputs, rather than 1D
!! - Add calculation for 'air_sea_deposition' taken from
!!   https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_dep_fluxes_calc.F90
!! - Multiply fluxes by ice_fraction input, rather than masking based on seawater input
!! - Use MOM over FMS modules where easy to do so
!! - Make tsurf input optional, as it is only used by a few implementations
!! - Use ind_runoff rather than ind_deposition in runoff flux calculation (note, their
!!   values are equal)
!! - Rename gas_fields_ice to gas_fields_ocn
subroutine atmos_ocean_fluxes_calc(gas_fields_atm, gas_fields_ocn, gas_fluxes,&
  ice_fraction, isc, iec, jsc, jec, tsurf, ustar, cd_m)
  type(coupler_2d_bc_type), intent(in)     :: gas_fields_atm ! fields in atm
      !< Structure containing atmospheric surface variables that are used in the calculation
      !! of the atmosphere-ocean tracer fluxes.
  type(coupler_2d_bc_type), intent(in)     :: gas_fields_ocn ! fields atop the ocean
      !< Structure containing ocean surface variables that are used in the calculation of the
      !! atmosphere-ocean tracer fluxes.
  type(coupler_2d_bc_type), intent(inout)  :: gas_fluxes ! fluxes between the atm and ocean
      !< Structure containing the gas fluxes between the atmosphere and the ocean and
      !! parameters related to the calculation of these fluxes.
  real, intent(in)                         :: ice_fraction(isc:iec,jsc:jec) !< sea ice fraction
  integer, intent(in)                      :: isc !< The start i-index of cell centers within
                                                  !! the computational domain
  integer, intent(in)                      :: iec !< The end i-index of cell centers within the
                                                  !! computational domain
  integer, intent(in)                      :: jsc !< The start j-index of cell centers within
                                                  !! the computational domain
  integer, intent(in)                      :: jec !< The end j-index of cell centers within the
                                                  !! computational domain
  real, intent(in), optional               :: tsurf(isc:iec,jsc:jec) !< surface temperature
  real, intent(in), optional               :: ustar(isc:iec,jsc:jec) !< friction velocity, not
                                                                      !! used
  real, intent(in), optional               :: cd_m (isc:iec,jsc:jec) !< drag coefficient, not
                                                                      !! used

  ! local variables
  character(len=*), parameter   :: sub_name = 'atmos_ocean_fluxes_calc'
  character(len=*), parameter   :: error_header =&
      & '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
  real, parameter                         :: permeg=1.0e-6

  integer                                 :: n
  integer                                 :: i
  integer                                 :: j
  real, dimension(:,:), allocatable       :: kw
  real, dimension(:,:), allocatable       :: cair
  character(len=128)                      :: error_string

  ! Return if no fluxes to be calculated
  if (gas_fluxes%num_bcs .le. 0) return

  if (.not. associated(gas_fluxes%bc)) then
    if (gas_fluxes%num_bcs .ne. 0) then
      call MOM_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
    else
      return
    endif
  endif

  do n = 1, gas_fluxes%num_bcs
    ! only do calculations if the flux has not been overridden
    if ( .not. gas_fluxes%bc(n)%field(ind_flux)%override) then
      if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux_generic') then
        if (.not. allocated(kw)) then
          allocate( kw(isc:iec,jsc:jec) )
          allocate ( cair(isc:iec,jsc:jec) )
        elseif ((size(kw(:,:), dim=1) .ne. iec-isc+1) .or. (size(kw(:,:), dim=2) .ne. jec-jsc+1)) then
          call MOM_error(FATAL, trim(error_header) // ' Sizes of flux fields do not match')
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_kw)%values(i,j) =&
                  & (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) * &
                  & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)**2
              cair(i,j) = &
                  gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * &
                  gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) * &
                  gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & sqrt(660. / (gas_fields_ocn%bc(n)%field(ind_sc_no)%values(i,j) + epsln)) *&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
              gas_fluxes%bc(n)%field(ind_flux0)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & sqrt(660. / (gas_fields_ocn%bc(n)%field(ind_sc_no)%values(i,j) + epsln)) *&
                  & gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j)
              gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j)) / &
                (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'duce') then
          if (.not. present(tsurf)) then
            call MOM_error(FATAL, trim(error_header) // ' Implementation ' //&
                trim(gas_fluxes%bc(n)%implementation) // ' for ' // trim(gas_fluxes%bc(n)%name) //&
                ' requires input tsurf')
          endif
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_kw)%values(i,j) = &
                  & (1 - ice_fraction(i,j)) * gas_fields_atm%bc(n)%field(ind_u10)%values(i,j) /&
                  & (770.+45.*gas_fluxes%bc(n)%param(1)**(1./3.)) *&
                  & 101325./(rdgas*wtmair*1e-3*tsurf(i,j) *&
                  & max(gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j),epsln))
              !alpha: mol/m3/atm
              cair(i,j) = &
                  gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * &
                  gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) * &
                  gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * 9.86923e-6
              cair(i,j) = max(cair(i,j),0.)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j))
              gas_fluxes%bc(n)%field(ind_flux0)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.)
              gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                  & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j)) /&
                  & (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'johnson') then
          if (.not. present(tsurf)) then
            call MOM_error(FATAL, trim(error_header) // ' Implementation ' //&
                trim(gas_fluxes%bc(n)%implementation) // ' for ' // trim(gas_fluxes%bc(n)%name) //&
                ' requires input tsurf')
          endif
          !f1p: not sure how to pass salinity. For now, just force at 35.
          do j = jsc,jec
            do i = isc,iec
              !calc_kw(tk,p,u10,h,vb,mw,sc_w,ustar,cd_m)
              gas_fluxes%bc(n)%field(ind_kw)%values(i,j) =&
                  & (1 - ice_fraction(i,j)) * calc_kw(tsurf(i,j),&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j),&
                  & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j),&
                  & 101325./(rdgas*wtmair*1e-3*tsurf(i,j)*&
                  & max(gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j),epsln)),&
                  & gas_fluxes%bc(n)%param(2),&
                  & gas_fluxes%bc(n)%param(1),&
                  & gas_fields_ocn%bc(n)%field(ind_sc_no)%values(i,j))
              cair(i,j) =&
                  & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * 9.86923e-6
              cair(i,j) = max(cair(i,j),0.)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j))
              gas_fluxes%bc(n)%field(ind_flux0)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.)
              gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                  & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j)) /&
                  & (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
            enddo
          enddo
        else
          call MOM_error(FATAL, ' Unknown implementation (' //&
              & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux') then
        if (.not. allocated(kw)) then
          allocate( kw(isc:iec,jsc:jec) )
          allocate ( cair(isc:iec,jsc:jec) )
        elseif ((size(kw(:,:), dim=1) .ne. iec-isc+1) .or. (size(kw(:,:), dim=2) .ne. jec-jsc+1)) then
          call MOM_error(FATAL, trim(error_header) // ' Sizes of flux fields do not match')
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2_data') then
          do j = jsc,jec
            do i = isc,iec
              kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                  & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)
              cair(i,j) =&
                  & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
          do j = jsc,jec
            do i = isc,iec
              kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                  & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)**2
              cair(i,j) =&
                  & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'linear') then
          do j = jsc,jec
            do i = isc,iec
              kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                  & max(0.0, gas_fields_atm%bc(n)%field(ind_u10)%values(i,j) - gas_fluxes%bc(n)%param(2))
              cair(i,j) =&
                  & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(3)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
            enddo
          enddo
        else
          call MOM_error(FATAL, ' Unknown implementation (' //&
              & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then
        if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
          write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
          call MOM_error(FATAL, 'Bad parameter (' // trim(error_string) //&
              & ') for air_sea_deposition for ' // trim(gas_fluxes%bc(n)%name))
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'dry') then
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                  gas_fields_atm%bc(n)%field(ind_deposition)%values(i,j) / gas_fluxes%bc(n)%param(1)
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'wet') then
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                  gas_fields_atm%bc(n)%field(ind_deposition)%values(i,j) / gas_fluxes%bc(n)%param(1)
            enddo
          enddo
        else
          call MOM_error(FATAL, 'Unknown implementation (' //&
              & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      elseif (gas_fluxes%bc(n)%flux_type .eq. 'land_sea_runoff') then
        if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
          write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
          call MOM_error(FATAL, ' Bad parameter (' // trim(error_string) //&
              & ') for land_sea_runoff for ' // trim(gas_fluxes%bc(n)%name))
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'river') then
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                  & gas_fields_atm%bc(n)%field(ind_runoff)%values(i,j) /&
                  & gas_fluxes%bc(n)%param(1)
            enddo
          enddo
        else
          call MOM_error(FATAL, ' Unknown implementation (' //&
              & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      else
        call MOM_error(FATAL, ' Unknown flux_type (' // trim(gas_fluxes%bc(n)%flux_type) //&
            & ') for ' // trim(gas_fluxes%bc(n)%name))
      endif
    endif
  enddo

  if (allocated(kw)) then
    deallocate(kw)
    deallocate(cair)
  endif
end subroutine  atmos_ocean_fluxes_calc

!> Calculate \f$k_w\f$
!!
!! Taken from Johnson, Ocean Science, 2010. (http://doi.org/10.5194/os-6-913-2010)
!!
!! Uses equations defined in Liss[1974],
!! \f[
!!  F = K_g(c_g - H C_l) = K_l(c_g/H - C_l)
!! \f]
!! where \f$c_g\f$ and \f$C_l\f$ are the bulk gas and liquid concentrations, \f$H\f$
!! is the Henry's law constant (\f$H = c_{sg}/C_{sl}\f$, where \f$c_{sg}\f$ is the
!! equilibrium concentration in gas phase (\f$g/cm^3\f$ of air) and \f$C_{sl}\f$ is the
!! equilibrium concentration of unionised dissolved gas in liquid phase (\f$g/cm^3\f$
!! of water)),
!! \f[
!!    1/K_g = 1/k_g + H/k_l
!! \f]
!! and
!! \f[
!!    1/K_l = 1/k_l + 1/{Hk_g}
!! \f]
!! where \f$k_g\f$ and \f$k_l\f$ are the exchange constants for the gas and liquid
!! phases, respectively.
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
real function calc_kw(tk, p, u10, h, vb, mw, sc_w, ustar, cd_m)
  real, intent(in) :: tk !< temperature at surface in kelvin
  real, intent(in) :: p !< pressure at surface in pa
  real, intent(in) :: u10 !< wind speed at 10m above the surface in m/s
  real, intent(in) :: h !< Henry's law constant (\f$H=c_sg/C_sl\f$) (unitless)
  real, intent(in) :: vb !< Molar volume
  real, intent(in) :: mw !< molecular weight (g/mol)
  real, intent(in) :: sc_w
  real, intent(in), optional :: ustar !< Friction velocity (m/s).  If not provided,
                                      !! ustar = \f$u_{10} \sqrt{C_D}\f$.
  real, intent(in), optional :: cd_m !< Drag coefficient (\f$C_D\f$).  Used only if
                                      !! ustar is provided.
                                      !! If ustar is not provided,
                                      !! cd_m = \f$6.1 \times 10^{-4} + 0.63 \times 10^{-4} *u_10\f$

  real :: ra,rl,tc

  tc = tk-273.15
  ra = 1./max(h*calc_ka(tc,p,mw,vb,u10,ustar,cd_m),epsln)
  rl = 1./max(calc_kl(tc,u10,sc_w),epsln)
  calc_kw = 1./max(ra+rl,epsln)
end function calc_kw

!> Calculate \f$k_a\f$
!!
!! See calc_kw
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
real function calc_ka(t, p, mw, vb, u10, ustar, cd_m)
  real, intent(in) :: t !< temperature at surface in C
  real, intent(in) :: p !< pressure at surface in pa
  real, intent(in) :: mw !< molecular weight (g/mol)
  real, intent(in) :: vb !< molar volume
  real, intent(in) :: u10 !< wind speed at 10m above the surface in m/s
  real, intent(in), optional :: ustar !< Friction velocity (m/s).  If not provided,
                                      !! ustar = \f$u_{10} \sqrt{C_D}\f$.
  real, intent(in), optional :: cd_m !< Drag coefficient (\f$C_D\f$).  Used only if
                                      !! ustar is provided.
                                      !! If ustar is not provided,
                                      !! cd_m = \f$6.1 \times 10^{-4} + 0.63 \times 10^{-4} *u_10\f$

  real             :: sc
  real             :: ustar_t, cd_m_t

  if (.not. present(ustar)) then
    !drag coefficient
    cd_m_t = 6.1e-4 +0.63e-4*u10
    !friction velocity
    ustar_t = u10*sqrt(cd_m_t)
  else
    cd_m_t = cd_m
    ustar_t = ustar
  end if
  sc = schmidt_g(t,p,mw,vb)
  calc_ka = 1e-3+ustar_t/(13.3*sqrt(sc)+1/sqrt(cd_m_t)-5.+log(sc)/(2.*vonkarm))
end function calc_ka

!> Calculate \f$k_l\f$
!!
!! See calc_kw, and Nightingale, Global Biogeochemical Cycles, 2000
!! (https://doi.org/10.1029/1999GB900091)
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
real function calc_kl(t, v, sc)
  real, intent(in) :: t !< temperature at surface in C
  real, intent(in) :: v !< wind speed at surface in m/s
  real, intent(in) :: sc

  calc_kl = (((0.222*v**2)+0.333*v)*(max(sc,epsln)/600.)**(-0.5))/(100.*3600.)
end function calc_kl

!> Schmidt number of the gas in air
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
real function schmidt_g(t, p, mw, vb)
  real, intent(in) :: t !< temperature at surface in C
  real, intent(in) :: p !< pressure at surface in pa
  real, intent(in) :: mw !< molecular weight (g/mol)
  real, intent(in) :: vb !< molar volume

  real :: d,v

  d = d_air(t,p,mw,vb)
  v = v_air(t)
  schmidt_g = v / d
end function schmidt_g

!> From Fuller, Industrial & Engineering Chemistry (https://doi.org/10.1021/ie50677a007)
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
real function d_air(t, p, mw, vb)
  real, intent(in) :: t  !< temperature in c
  real, intent(in) :: p  !< pressure in pa
  real, intent(in) :: mw !< molecular weight (g/mol)
  real, intent(in) :: vb !< diffusion coefficient (\f$cm3/mol\f$)

  real, parameter :: ma = 28.97d0 !< molecular weight air in g/mol
  real, parameter :: va = 20.1d0  !< diffusion volume for air (\f$cm^3/mol\f$)

  real            :: pa

  ! convert p to atm
  pa = 9.8692d-6*p
  d_air = 1d-3 *&
      & (t+273.15d0)**(1.75d0)*sqrt(1d0/ma + 1d0/mw)/(pa*(va**(1d0/3d0)+vb**(1d0/3d0))**2d0)
  ! d_air is in cm2/s convert to m2/s
  d_air = d_air * 1d-4
end function d_air

!> kinematic viscosity in air
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
real function p_air(t)
  real, intent(in) :: t

  real, parameter :: sd_0 = 1.293393662d0,&
      & sd_1 = -5.538444326d-3,&
      & sd_2 = 3.860201577d-5,&
      & sd_3 = -5.2536065d-7
  p_air = sd_0+(sd_1*t)+(sd_2*t**2)+(sd_3*t**3)
end function p_air

!> Kinematic viscosity in air (\f$m^2/s\f$
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
real function v_air(t)
  real, intent(in) :: t !< temperature in C
  v_air = n_air(t)/p_air(t)
end function v_air

!> dynamic viscosity in air
!!
!! This routine was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d38/full/atmos_ocean_fluxes_calc.F90
real function n_air(t)
  real, intent(in) :: t !< temperature in C

  real, parameter :: sv_0 = 1.715747771d-5,&
      & sv_1 = 4.722402075d-8,&
      & sv_2 = -3.663027156d-10,&
      & sv_3 = 1.873236686d-12,&
      & sv_4 = -8.050218737d-14
  ! in n.s/m^2 (pa.s)
  n_air = sv_0+(sv_1*t)+(sv_2*t**2)+(sv_3*t**3)+(sv_4*t**4)
end function n_air

end module MOM_cap_fms_coupler_bcs