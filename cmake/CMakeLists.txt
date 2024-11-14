cmake_minimum_required(VERSION 3.18)

function(add_fortran_library LIB MOD_DIR)
  add_library(${LIB} ${ARGN})

  get_target_property(LIB_DIR ${LIB} BINARY_DIR)
  set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/${MOD_DIR})

  target_include_directories(${LIB} INTERFACE "$<BUILD_INTERFACE:${LIB_DIR}/${MOD_DIR}>")
endfunction(add_fortran_library)

project(MOM6 DESCRIPTION "MOM6" LANGUAGES C Fortran)

function(add_fortran_library LIB MOD_DIR)
    add_library(${LIB} ${ARGN})

    get_target_property(LIB_DIR ${LIB} BINARY_DIR)
    set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/${MOD_DIR})

    target_include_directories(${LIB} INTERFACE "$<BUILD_INTERFACE:${LIB_DIR}/${MOD_DIR}>")
  endfunction(add_fortran_library)
  
# Common compiler flags and definitions
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Og -ffpe-trap=zero,overflow -fcheck=bounds")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
  set(CMAKE_C_FLAGS_RELEASE "-O")
  set(CMAKE_C_FLAGS_DEBUG "-g -Wall -Og -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -traceback -qno-opt-dynamic-align -fp-model precise -std=gnu99")
  set(CMAKE_C_FLAGS_RELEASE "-O2 -debug minimal")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
else()
  message(WARNING "C compiler with ID ${CMAKE_C_COMPILER_ID} will be used with CMake default options")
endif()

add_compile_definitions(
  _USE_GENERIC_TRACER
  _USE_MOM6_DIAG
)

### Targets

set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/../src")
set(CONFIG_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../config_src")

## MOM6 library
add_fortran_library(OM3_mom6 mod STATIC)
add_library(AccessOM3::mom6 ALIAS OM3_mom6)
target_include_directories(OM3_mom6 PRIVATE $<BUILD_INTERFACE:${SRC}/framework>)
if(OM3_MOM_SYMMETRIC)
  target_include_directories(OM3_mom6 PRIVATE $<BUILD_INTERFACE:${CONFIG_SRC}/memory/dynamic_symmetric>)
else()
  target_include_directories(OM3_mom6 PRIVATE $<BUILD_INTERFACE:${CONFIG_SRC}/memory/dynamic_nonsymmetric>)
endif()
target_compile_options(OM3_mom6 PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:${fortran_compile_flags}>")
target_link_libraries(OM3_mom6
  PUBLIC esmf
  PRIVATE AccessOM3::nuopc_cap_share AccessOM3::share FMS::fms_r8
)


target_sources(OM3_mom6 PRIVATE
  ${SRC}/ALE/coord_adapt.F90
  ${SRC}/ALE/coord_hycom.F90
  ${SRC}/ALE/coord_rho.F90
  ${SRC}/ALE/coord_sigma.F90
  ${SRC}/ALE/coord_zlike.F90
  ${SRC}/ALE/MOM_ALE.F90
  ${SRC}/ALE/MOM_hybgen_regrid.F90
  ${SRC}/ALE/MOM_hybgen_remap.F90
  ${SRC}/ALE/MOM_hybgen_unmix.F90
  ${SRC}/ALE/MOM_regridding.F90
  ${SRC}/ALE/MOM_remapping.F90
  ${SRC}/ALE/P1M_functions.F90
  ${SRC}/ALE/P3M_functions.F90
  ${SRC}/ALE/PCM_functions.F90
  ${SRC}/ALE/PLM_functions.F90
  ${SRC}/ALE/polynomial_functions.F90
  ${SRC}/ALE/PPM_functions.F90
  ${SRC}/ALE/PQM_functions.F90
  ${SRC}/ALE/regrid_consts.F90
  ${SRC}/ALE/regrid_edge_values.F90
  ${SRC}/ALE/regrid_interp.F90
  ${SRC}/ALE/regrid_solvers.F90
  ${SRC}/ALE/remapping_attic.F90

  ${SRC}/core/MOM_barotropic.F90
  ${SRC}/core/MOM_boundary_update.F90
  ${SRC}/core/MOM_check_scaling.F90
  ${SRC}/core/MOM_checksum_packages.F90
  ${SRC}/core/MOM_continuity.F90
  ${SRC}/core/MOM_continuity_PPM.F90
  ${SRC}/core/MOM_CoriolisAdv.F90
  ${SRC}/core/MOM_density_integrals.F90
  ${SRC}/core/MOM_dynamics_split_RK2.F90
  ${SRC}/core/MOM_dynamics_unsplit.F90
  ${SRC}/core/MOM_dynamics_unsplit_RK2.F90
  ${SRC}/core/MOM.F90
  ${SRC}/core/MOM_grid.F90
  ${SRC}/core/MOM_interface_heights.F90
  ${SRC}/core/MOM_isopycnal_slopes.F90
  ${SRC}/core/MOM_open_boundary.F90
  ${SRC}/core/MOM_porous_barriers.F90
  ${SRC}/core/MOM_PressureForce.F90
  ${SRC}/core/MOM_PressureForce_FV.F90
  ${SRC}/core/MOM_PressureForce_Montgomery.F90
  ${SRC}/core/MOM_stoch_eos.F90
  ${SRC}/core/MOM_transcribe_grid.F90
  ${SRC}/core/MOM_unit_tests.F90
  ${SRC}/core/MOM_variables.F90
  ${SRC}/core/MOM_verticalGrid.F90

  ${SRC}/diagnostics/MOM_debugging.F90
  ${SRC}/diagnostics/MOM_diagnostics.F90
  ${SRC}/diagnostics/MOM_obsolete_diagnostics.F90
  ${SRC}/diagnostics/MOM_obsolete_params.F90
  ${SRC}/diagnostics/MOM_PointAccel.F90
  ${SRC}/diagnostics/MOM_spatial_means.F90
  ${SRC}/diagnostics/MOM_sum_output.F90
  ${SRC}/diagnostics/MOM_wave_speed.F90

  ${SRC}/equation_of_state/MOM_EOS.F90
  ${SRC}/equation_of_state/MOM_EOS_Jackett06.F90
  ${SRC}/equation_of_state/MOM_EOS_linear.F90
  ${SRC}/equation_of_state/MOM_EOS_Roquet_rho.F90
  ${SRC}/equation_of_state/MOM_EOS_Roquet_SpV.F90
  ${SRC}/equation_of_state/MOM_EOS_TEOS10.F90
  ${SRC}/equation_of_state/MOM_EOS_UNESCO.F90
  ${SRC}/equation_of_state/MOM_EOS_Wright.F90
  ${SRC}/equation_of_state/MOM_EOS_Wright_full.F90
  ${SRC}/equation_of_state/MOM_EOS_Wright_red.F90
  ${SRC}/equation_of_state/MOM_temperature_convert.F90
  ${SRC}/equation_of_state/MOM_TFreeze.F90

  ${SRC}/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90
  ${SRC}/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90
  ${SRC}/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90
  ${SRC}/equation_of_state/TEOS10/gsw_ct_from_pt.f90
  ${SRC}/equation_of_state/TEOS10/gsw_ct_from_t.f90
  ${SRC}/equation_of_state/TEOS10/gsw_entropy_part.f90
  ${SRC}/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90
  ${SRC}/equation_of_state/TEOS10/gsw_gibbs.f90
  ${SRC}/equation_of_state/TEOS10/gsw_gibbs_ice.f90
  ${SRC}/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90
  ${SRC}/equation_of_state/TEOS10/gsw_mod_error_functions.f90
  ${SRC}/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90
  ${SRC}/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90
  ${SRC}/equation_of_state/TEOS10/gsw_mod_kinds.f90
  ${SRC}/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90
  ${SRC}/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90
  ${SRC}/equation_of_state/TEOS10/gsw_mod_toolbox.f90
  ${SRC}/equation_of_state/TEOS10/gsw_pt0_from_t.f90
  ${SRC}/equation_of_state/TEOS10/gsw_pt_from_ct.f90
  ${SRC}/equation_of_state/TEOS10/gsw_pt_from_t.f90
  ${SRC}/equation_of_state/TEOS10/gsw_rho.f90
  ${SRC}/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90
  ${SRC}/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90
  ${SRC}/equation_of_state/TEOS10/gsw_sp_from_sr.f90
  ${SRC}/equation_of_state/TEOS10/gsw_sr_from_sp.f90
  ${SRC}/equation_of_state/TEOS10/gsw_specvol.f90
  ${SRC}/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90
  ${SRC}/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90
  ${SRC}/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90
  ${SRC}/equation_of_state/TEOS10/gsw_t_freezing_exact.f90
  ${SRC}/equation_of_state/TEOS10/gsw_t_freezing_poly.f90
  ${SRC}/equation_of_state/TEOS10/gsw_t_from_ct.f90

  ${SRC}/framework/MOM_array_transform.F90
  ${SRC}/framework/MOM_checksums.F90
  ${SRC}/framework/MOM_coms.F90
  ${SRC}/framework/MOM_cpu_clock.F90
  ${SRC}/framework/MOM_data_override.F90
  ${SRC}/framework/MOM_diag_mediator.F90
  ${SRC}/framework/MOM_diag_remap.F90
  ${SRC}/framework/MOM_document.F90
  ${SRC}/framework/MOM_domains.F90
  ${SRC}/framework/MOM_dyn_horgrid.F90
  ${SRC}/framework/MOM_ensemble_manager.F90
  ${SRC}/framework/MOM_error_handler.F90
  ${SRC}/framework/MOM_file_parser.F90
  ${SRC}/framework/MOM_get_input.F90
  ${SRC}/framework/MOM_hor_index.F90
  ${SRC}/framework/MOM_horizontal_regridding.F90
  ${SRC}/framework/MOM_interpolate.F90
  ${SRC}/framework/MOM_intrinsic_functions.F90
  ${SRC}/framework/MOM_io.F90
  ${SRC}/framework/MOM_io_file.F90
  ${SRC}/framework/MOM_memory_macros.h
  ${SRC}/framework/MOM_netcdf.F90
  ${SRC}/framework/MOM_random.F90
  ${SRC}/framework/MOM_restart.F90
  ${SRC}/framework/MOM_safe_alloc.F90
  ${SRC}/framework/MOM_string_functions.F90
  ${SRC}/framework/MOM_unique_scales.F90
  ${SRC}/framework/MOM_unit_scaling.F90
  ${SRC}/framework/MOM_write_cputime.F90
  ${SRC}/framework/posix.F90

  ${SRC}/ice_shelf/MOM_ice_shelf_diag_mediator.F90
  ${SRC}/ice_shelf/MOM_ice_shelf_dynamics.F90
  ${SRC}/ice_shelf/MOM_ice_shelf.F90
  ${SRC}/ice_shelf/MOM_ice_shelf_initialize.F90
  ${SRC}/ice_shelf/MOM_ice_shelf_state.F90
  ${SRC}/ice_shelf/MOM_marine_ice.F90
  ${SRC}/ice_shelf/user_shelf_init.F90

  ${SRC}/initialization/MOM_coord_initialization.F90
  ${SRC}/initialization/MOM_fixed_initialization.F90
  ${SRC}/initialization/MOM_grid_initialize.F90
  ${SRC}/initialization/MOM_shared_initialization.F90
  ${SRC}/initialization/MOM_state_initialization.F90
  ${SRC}/initialization/MOM_tracer_initialization_from_Z.F90

  ${SRC}/ocean_data_assim/MOM_oda_driver.F90
  ${SRC}/ocean_data_assim/MOM_oda_incupd.F90

  ${SRC}/parameterizations/CVmix/cvmix_background.F90
  ${SRC}/parameterizations/CVmix/cvmix_convection.F90
  ${SRC}/parameterizations/CVmix/cvmix_ddiff.F90
  ${SRC}/parameterizations/CVmix/cvmix_kinds_and_types.F90
  ${SRC}/parameterizations/CVmix/cvmix_kpp.F90
  ${SRC}/parameterizations/CVmix/cvmix_math.F90
  ${SRC}/parameterizations/CVmix/cvmix_put_get.F90
  ${SRC}/parameterizations/CVmix/cvmix_shear.F90
  ${SRC}/parameterizations/CVmix/cvmix_tidal.F90
  ${SRC}/parameterizations/CVmix/cvmix_utils.F90

  ${SRC}/parameterizations/lateral/MOM_hor_visc.F90
  ${SRC}/parameterizations/lateral/MOM_interface_filter.F90
  ${SRC}/parameterizations/lateral/MOM_internal_tides.F90
  ${SRC}/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90
  ${SRC}/parameterizations/lateral/MOM_load_love_numbers.F90
  ${SRC}/parameterizations/lateral/MOM_MEKE.F90
  ${SRC}/parameterizations/lateral/MOM_MEKE_types.F90
  ${SRC}/parameterizations/lateral/MOM_mixed_layer_restrat.F90
  ${SRC}/parameterizations/lateral/MOM_spherical_harmonics.F90
  ${SRC}/parameterizations/lateral/MOM_thickness_diffuse.F90
  ${SRC}/parameterizations/lateral/MOM_tidal_forcing.F90
  ${SRC}/parameterizations/lateral/MOM_Zanna_Bolton.F90

  ${SRC}/parameterizations/stochastic/MOM_stochastics.F90

  ${SRC}/parameterizations/vertical/MOM_ALE_sponge.F90
  ${SRC}/parameterizations/vertical/MOM_bkgnd_mixing.F90
  ${SRC}/parameterizations/vertical/MOM_bulk_mixed_layer.F90
  ${SRC}/parameterizations/vertical/MOM_CVMix_conv.F90
  ${SRC}/parameterizations/vertical/MOM_CVMix_ddiff.F90
  ${SRC}/parameterizations/vertical/MOM_CVMix_KPP.F90
  ${SRC}/parameterizations/vertical/MOM_CVMix_shear.F90
  ${SRC}/parameterizations/vertical/MOM_diabatic_aux.F90
  ${SRC}/parameterizations/vertical/MOM_diabatic_driver.F90
  ${SRC}/parameterizations/vertical/MOM_diapyc_energy_req.F90
  ${SRC}/parameterizations/vertical/MOM_energetic_PBL.F90
  ${SRC}/parameterizations/vertical/MOM_entrain_diffusive.F90
  ${SRC}/parameterizations/vertical/MOM_full_convection.F90
  ${SRC}/parameterizations/vertical/MOM_geothermal.F90
  ${SRC}/parameterizations/vertical/MOM_internal_tide_input.F90
  ${SRC}/parameterizations/vertical/MOM_kappa_shear.F90
  ${SRC}/parameterizations/vertical/MOM_opacity.F90
  ${SRC}/parameterizations/vertical/MOM_regularize_layers.F90
  ${SRC}/parameterizations/vertical/MOM_set_diffusivity.F90
  ${SRC}/parameterizations/vertical/MOM_set_viscosity.F90
  ${SRC}/parameterizations/vertical/MOM_sponge.F90
  ${SRC}/parameterizations/vertical/MOM_tidal_mixing.F90
  ${SRC}/parameterizations/vertical/MOM_vert_friction.F90

  ${SRC}/tracer/advection_test_tracer.F90
  ${SRC}/tracer/boundary_impulse_tracer.F90
  ${SRC}/tracer/DOME_tracer.F90
  ${SRC}/tracer/dyed_obc_tracer.F90
  ${SRC}/tracer/dye_example.F90
  ${SRC}/tracer/ideal_age_example.F90
  ${SRC}/tracer/ISOMIP_tracer.F90
  ${SRC}/tracer/MOM_CFC_cap.F90
  ${SRC}/tracer/MOM_hor_bnd_diffusion.F90
  ${SRC}/tracer/MOM_neutral_diffusion.F90
  ${SRC}/tracer/MOM_OCMIP2_CFC.F90
  ${SRC}/tracer/MOM_offline_aux.F90
  ${SRC}/tracer/MOM_offline_main.F90
  ${SRC}/tracer/MOM_tracer_advect.F90
  ${SRC}/tracer/MOM_tracer_diabatic.F90
  ${SRC}/tracer/MOM_tracer_flow_control.F90
  ${SRC}/tracer/MOM_tracer_hor_diff.F90
  ${SRC}/tracer/MOM_tracer_registry.F90
  ${SRC}/tracer/MOM_tracer_types.F90
  ${SRC}/tracer/MOM_tracer_Z_init.F90
  ${SRC}/tracer/nw2_tracers.F90
  ${SRC}/tracer/oil_tracer.F90
  ${SRC}/tracer/pseudo_salt_tracer.F90
  ${SRC}/tracer/RGC_tracer.F90
  ${SRC}/tracer/tracer_example.F90

  ${SRC}/user/adjustment_initialization.F90
  ${SRC}/user/baroclinic_zone_initialization.F90
  ${SRC}/user/basin_builder.F90
  ${SRC}/user/benchmark_initialization.F90
  ${SRC}/user/BFB_initialization.F90
  ${SRC}/user/BFB_surface_forcing.F90
  ${SRC}/user/circle_obcs_initialization.F90
  ${SRC}/user/dense_water_initialization.F90
  ${SRC}/user/DOME2d_initialization.F90
  ${SRC}/user/DOME_initialization.F90
  ${SRC}/user/dumbbell_initialization.F90
  ${SRC}/user/dumbbell_surface_forcing.F90
  ${SRC}/user/dyed_channel_initialization.F90
  ${SRC}/user/dyed_obcs_initialization.F90
  ${SRC}/user/external_gwave_initialization.F90
  ${SRC}/user/Idealized_Hurricane.F90
  ${SRC}/user/ISOMIP_initialization.F90
  ${SRC}/user/Kelvin_initialization.F90
  ${SRC}/user/lock_exchange_initialization.F90
  ${SRC}/user/MOM_controlled_forcing.F90
  ${SRC}/user/MOM_wave_interface.F90
  ${SRC}/user/Neverworld_initialization.F90
  ${SRC}/user/Phillips_initialization.F90
  ${SRC}/user/RGC_initialization.F90
  ${SRC}/user/Rossby_front_2d_initialization.F90
  ${SRC}/user/SCM_CVMix_tests.F90
  ${SRC}/user/seamount_initialization.F90
  ${SRC}/user/shelfwave_initialization.F90
  ${SRC}/user/sloshing_initialization.F90
  ${SRC}/user/soliton_initialization.F90
  ${SRC}/user/supercritical_initialization.F90
  ${SRC}/user/tidal_bay_initialization.F90
  ${SRC}/user/user_change_diffusivity.F90
  ${SRC}/user/user_initialization.F90
  ${SRC}/user/user_revise_forcing.F90

  ${CONFIG_SRC}/external/database_comms/MOM_database_comms.F90
  ${CONFIG_SRC}/external/database_comms/database_client_interface.F90

  ${CONFIG_SRC}/external/drifters/MOM_particles.F90
  ${CONFIG_SRC}/external/drifters/MOM_particles_types.F90

  ${CONFIG_SRC}/external/ODA_hooks/kdtree.f90
  ${CONFIG_SRC}/external/ODA_hooks/ocean_da_core.F90
  ${CONFIG_SRC}/external/ODA_hooks/ocean_da_types.F90
  ${CONFIG_SRC}/external/ODA_hooks/write_ocean_obs.F90

  ${CONFIG_SRC}/external/stochastic_physics/get_stochy_pattern.F90
  ${CONFIG_SRC}/external/stochastic_physics/stochastic_physics.F90

  ${CONFIG_SRC}/infra/FMS2/MOM_coms_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_constants.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_cpu_clock_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_data_override_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_diag_manager_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_domain_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_ensemble_manager_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_error_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_interp_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_io_infra.F90
  ${CONFIG_SRC}/infra/FMS2/MOM_time_manager.F90

  ${CONFIG_SRC}/drivers/nuopc_cap/mom_cap_time.F90
  ${CONFIG_SRC}/drivers/nuopc_cap/mom_surface_forcing_nuopc.F90
  ${CONFIG_SRC}/drivers/nuopc_cap/ocn_comp_NUOPC.F90
  ${CONFIG_SRC}/drivers/nuopc_cap/time_utils.F90

  #extra_sources/mom_cap_gtracer_flux.F90
)

### Install and Export

## Library
if(OM3_LIB_INSTALL)
  set_target_properties(OM3_mom6 PROPERTIES
    OUTPUT_NAME access-mom6
    EXPORT_NAME mom6
  )
  install(TARGETS OM3_mom6
    EXPORT AccessOM3mom6_Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT AccessOM3_RunTime NAMELINK_COMPONENT AccessOM3_Development
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT AccessOM3_Development
  )
  # Fortran module files are a special case, as currently there is no standard
  # way of handling them in CMake
  target_include_directories(OM3_mom6 PUBLIC "$<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}/access-mom6>")
  get_target_property(mom_moddir OM3_mom6 Fortran_MODULE_DIRECTORY)
  install(FILES ${mom_moddir}/ocn_comp_nuopc.mod ${mom_moddir}/mom_cap_mod.mod
    DESTINATION ${CMAKE_INSTALL_MODULEDIR}/access-mom6
    COMPONENT AccessOM3_Development
  )
  install(EXPORT AccessOM3mom6_Targets
    FILE AccessOM3mom6Targets.cmake
    NAMESPACE AccessOM3::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AccessOM3
  )
endif()
