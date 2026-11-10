# Copyright ACCESS-NRI and contributors. See the top-level LICENSE file for details.

cmake_minimum_required(VERSION 3.18)

#[==============================================================================[
#                         Basic project definition                              #
#]==============================================================================]
project(MOM6
        DESCRIPTION "Modular Ocean Model 6"
        HOMEPAGE_URL https://github.com/ACCESS-NRI/MOM6
        LANGUAGES C Fortran)

#[==============================================================================[
#                           Project configuration                               #
#]==============================================================================]
message(STATUS "Build options")

option(MOM6_ACCESS3    "Build MOM6 library for use with ACCESS3"  OFF)
option(MOM6_SOLO       "Build MOM6 solo executable"               ON)
option(MOM6_OPENMP     "Enable OpenMP threading"                  OFF)
option(MOM6_ASYMMETRIC "Use MOM asymmetric memory"                OFF)
option(MOM6_FAST_RINT  "Use the fast round-to-nearest-integer in the reproducible exp and log" ON)
option(MOM6_OPENMP_OFFLOAD "Enable NVHPC OpenMP/OpenACC/do-concurrent GPU offload flags (requires NVHPC)" OFF)
set(MOM6_OFFLOAD_MEM_MODEL "separate" CACHE STRING "GPU memory model for MOM6_OPENMP_OFFLOAD: unified, managed, or separate (maps to -gpu=mem:<model>)")
set_property(CACHE MOM6_OFFLOAD_MEM_MODEL PROPERTY STRINGS unified managed separate)
set(MOM6_OFFLOAD_CC "" CACHE STRING "GPU compute capabilities for MOM6_OPENMP_OFFLOAD, as a semicolon-separated list (e.g. 70;80;90); each entry NN maps to ccNN in -gpu=. Empty (default) omits -gpu=cc... entirely, so NVHPC targets the native/local GPU.")
message(STATUS "  - MOM6_ACCESS3     ${MOM6_ACCESS3}")
message(STATUS "  - MOM6_SOLO        ${MOM6_SOLO}")
message(STATUS "  - MOM6_OPENMP      ${MOM6_OPENMP}")
message(STATUS "  - MOM6_ASYMMETRIC  ${MOM6_ASYMMETRIC}")
message(STATUS "  - MOM6_FAST_RINT   ${MOM6_FAST_RINT}")
message(STATUS "  - MOM6_OPENMP_OFFLOAD ${MOM6_OPENMP_OFFLOAD}")
if(MOM6_OPENMP_OFFLOAD)
  message(STATUS "  - MOM6_OFFLOAD_MEM_MODEL ${MOM6_OFFLOAD_MEM_MODEL}")
  if(MOM6_OFFLOAD_CC)
    message(STATUS "  - MOM6_OFFLOAD_CC  ${MOM6_OFFLOAD_CC}")
  else()
    message(STATUS "  - MOM6_OFFLOAD_CC  (unset - compiler defaults to native GPU)")
  endif()
endif()

if(NOT MOM6_ACCESS3 AND NOT MOM6_SOLO)
  message(FATAL_ERROR "One of MOM6_ACCESS3 or MOM6_SOLO must be enabled")
endif()

if(MOM6_OPENMP_OFFLOAD AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  message(FATAL_ERROR "MOM6_OPENMP_OFFLOAD requires the NVHPC Fortran compiler (found ${CMAKE_Fortran_COMPILER_ID})")
endif()

if(MOM6_OPENMP_OFFLOAD AND NOT MOM6_OFFLOAD_MEM_MODEL MATCHES "^(unified|managed|separate)$")
  message(FATAL_ERROR "MOM6_OFFLOAD_MEM_MODEL must be one of unified, managed, separate (got '${MOM6_OFFLOAD_MEM_MODEL}')")
endif()

#[==============================================================================[
#                           Project configuration                               #
#]==============================================================================]

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(GNUInstallDirs)
include(FortranLib)
include(CMakePackageConfigHelpers)

# Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -fdefault-real-8 -fdefault-double-8 -g -grecord-gcc-switches")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O")
  set(CMAKE_Fortran_FLAGS_DEBUG "-Wall -Og -ffpe-trap=zero,overflow -fcheck=bounds")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model precise -r8 -g -grecord-gcc-switches")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check noarg_temp_created -check all") #  -fpe0 see https://github.com/ACCESS-NRI/MOM6/issues/39
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  # -r8 promotes default REAL to 8 bytes, matching the other compiler branches
  # above. Without it, this library's plain `real` declarations stay 4 bytes
  # while FMS (R8 component) and the GSW-Fortran submodule (explicit real(r8),
  # selected_real_kind(14,30)) are double precision, causing kind-mismatch
  # compile errors (e.g. gsw_ct_freezing_exact in MOM_TFreeze.F90).
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
  if(MOM6_OPENMP_OFFLOAD)
    # NVHPC accepts multiple -gpu= flags and accumulates their suboptions, so
    # the memory model and compute-capability list are kept as separate flags
    # rather than combined into one comma-joined string. The compute-capability
    # list is optional: if MOM6_OFFLOAD_CC is unset, -gpu=cc... is omitted
    # entirely and NVHPC defaults to targeting the native/local GPU.
    set(_mom6_offload_gpu "-gpu=mem:${MOM6_OFFLOAD_MEM_MODEL}")
    if(MOM6_OFFLOAD_CC)
      list(TRANSFORM MOM6_OFFLOAD_CC PREPEND "cc" OUTPUT_VARIABLE _mom6_offload_cc)
      list(JOIN _mom6_offload_cc "," _mom6_offload_cc)
      string(APPEND _mom6_offload_gpu " -gpu=${_mom6_offload_cc}")
      unset(_mom6_offload_cc)
    endif()
    # -Mnofma/-mno-fma/-Mnovect preserve bit-reproducibility with the
    # reference CPU build; see SPACK-NVHPC-FMA-ANSWER-CHANGE.md. Requires an
    # optimized build type (e.g. -DCMAKE_BUILD_TYPE=Release): compiling with
    # no -O level at all has been observed to segfault inside NVHPC's runtime
    # array-gather while writing the final restart file (MOM_restart.F90).
    #
    # These must land AFTER the build type's -O.../-g flags on the actual
    # compiler command line, not in CMAKE_Fortran_FLAGS: CMake always emits
    # CMAKE_Fortran_FLAGS before CMAKE_Fortran_FLAGS_<CONFIG>, regardless of
    # the order these set() calls run in. NVHPC applies -M/-m options
    # left-to-right with last-wins semantics, and -O<n> re-asserts its own
    # bundled vectorization -- so with -Mnovect placed before -O2, -O2 was
    # silently re-enabling vectorization, defeating the bit-reproducibility
    # goal for every optimized build type. Confirmed empirically on sm_70:
    # with -O2 before these flags, a representative do-concurrent kernel
    # compiles to scalar host code; with -O2 after them (the old order),
    # the same kernel vectorizes (vmulpd/vaddpd), changing the summation
    # order and hence the answers. GPU device code was unaffected either
    # way in that test, but keep these ordered after -O for safety on the
    # host side regardless.
    set(_mom6_offload_flags "-mp=gpu -stdpar=gpu -acc=gpu ${_mom6_offload_gpu} -Mnofma -Mnovect -mno-fma -Minline=name:flux_elem,name:flux_elem_OBC,name:ratio_max")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${_mom6_offload_flags}")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_mom6_offload_flags}")
    # -gpu=... must also be passed at link (not just compile), since it selects
    # which GPU runtime/memory-model/compute-capability support gets linked in.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mp=gpu ${_mom6_offload_gpu}")
    # Keep the offload flags for the OpenMP feature probes below.
    set(MOM6_OFFLOAD_FC_FLAGS "${_mom6_offload_flags}")
    unset(_mom6_offload_gpu)
    unset(_mom6_offload_cc)
    unset(_mom6_offload_flags)
  endif()
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

# Check whether the Fortran compiler supports do-concurrent locality
# specifiers (local/local_init/reduce, Fortran 2018). Mirrors the autoconf
# MOM6_FC_DO_CONCURRENT_LOCALITY check (ac/m4/mom6_fc_do_concurrent_locality.m4),
# and uses the same test program.
# Without this define, DO_LOCALITY() in do_concurrent_compat.h strips all
# locality clauses from `do concurrent` loops, silently turning per-iteration
# scratch variables and reductions (e.g. MOM_coms.F90, MOM_sum_output.F90)
# into shared/racy state when the loops are actually parallelized (-stdpar=gpu).
include(CheckFortranSourceCompiles)
check_fortran_source_compiles(
"program test
  implicit none
  integer :: i, a, b, total, bins(2), largest
  logical :: found
  total = 0
  bins(:) = 0
  largest = 0
  found = .false.
  b = 2
  do concurrent(i=1:2) local(a) local_init(b) reduce(+:total,bins)
    a = i
    b = b * i
    total = total + a
    bins(:) = bins(:) + b
  end do
  do concurrent(i=1:2) reduce(max:largest) reduce(.or.:found)
    largest = max(largest, i)
    found = found .or. (i == 2)
  end do
end program test"
  HAVE_FC_DO_CONCURRENT_LOCALITY
  SRC_EXT F90
)
if(HAVE_FC_DO_CONCURRENT_LOCALITY)
  add_compile_definitions(HAVE_FC_DO_CONCURRENT_LOCALITY=1)
endif()

# Find out whether the Fortran compiler honors parentheses in floating-point
# expressions, which the fast round-to-nearest-integer in the reproducible exp
# and log (src/framework/intrinsics/MOM_exp.h) relies on: (x + bias) - bias must
# not be reassociated into x.  Mirrors the autoconf MOM6_FC_FAST_RINT check
# (ac/m4/mom6_fc_fast_rint.m4), including its --disable-fast-rint (here
# MOM6_FAST_RINT=OFF), the candidates it tries and the test program.  The first
# candidate whose program runs and gives the right answer is used: ENABLE_FAST_RINT
# is defined and its flag, if it needs one, is added to the Fortran flags.  If none
# works, the ieee_rint fallback of MOM_exp.h is used.
if(MOM6_FAST_RINT)
  include(CheckFortranSourceRuns)
  set(_mom6_fast_rint_found FALSE)
  if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    message(WARNING "Cross-compiling without an emulator: cannot test whether the Fortran "
                    "compiler honors parentheses; using the ieee_rint fallback")
  else()
    set(_mom6_fast_rint_n 0)
    foreach(_mom6_fast_rint_flag IN ITEMS "none" "-fprotect-parens" "-assume protect_parens")
      # The program is run with the flags of the Release build, since that is where an
      # optimizer could reassociate the expression.
      set(CMAKE_REQUIRED_FLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
      if(NOT _mom6_fast_rint_flag STREQUAL "none")
        string(APPEND CMAKE_REQUIRED_FLAGS " ${_mom6_fast_rint_flag}")
      endif()
      check_fortran_source_runs(
"program test
  implicit none
  real, parameter :: bias = 1.5 * 2.**(digits(1.) - 1)
  if (fast_rint(1.23) == 1.) stop 0
  stop 1
contains
  function fast_rint(x) result(y)
    real, intent(in) :: x
    real :: y
    y = (x + bias) - bias
  end function
end program test"
        MOM6_FC_FAST_RINT_OK_${_mom6_fast_rint_n}
        SRC_EXT F90
      )
      if(MOM6_FC_FAST_RINT_OK_${_mom6_fast_rint_n})
        set(_mom6_fast_rint_found TRUE)
        set(_mom6_fast_rint_chosen "${_mom6_fast_rint_flag}")
        break()
      endif()
      math(EXPR _mom6_fast_rint_n "${_mom6_fast_rint_n} + 1")
    endforeach()
    unset(CMAKE_REQUIRED_FLAGS)
  endif()
  if(_mom6_fast_rint_found)
    if(_mom6_fast_rint_chosen STREQUAL "none")
      message(STATUS "Parentheses are honored by the Fortran compiler: using the fast rint")
    else()
      message(STATUS "Using the fast rint with the Fortran flag ${_mom6_fast_rint_chosen}")
      # Like the offload flags, added to the build type flags so that nothing overrides it.
      string(APPEND CMAKE_Fortran_FLAGS_RELEASE " ${_mom6_fast_rint_chosen}")
      string(APPEND CMAKE_Fortran_FLAGS_DEBUG " ${_mom6_fast_rint_chosen}")
    endif()
    add_compile_definitions(ENABLE_FAST_RINT=1)
  else()
    message(WARNING "No known flag found to protect parentheses; using the ieee_rint fallback")
  endif()
  unset(_mom6_fast_rint_n)
  unset(_mom6_fast_rint_flag)
  unset(_mom6_fast_rint_chosen)
  unset(_mom6_fast_rint_found)
endif()

# C compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -g -grecord-gcc-switches")
  set(CMAKE_C_FLAGS_RELEASE "-O")
  set(CMAKE_C_FLAGS_DEBUG "-Wall -Og -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -traceback -qno-opt-dynamic-align -fp-model precise -std=gnu99 -g -grecord-gcc-switches")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
else()
  message(WARNING "C compiler with ID ${CMAKE_C_COMPILER_ID} will be used with CMake default options")
endif()

# See https://github.com/ACCESS-NRI/GFDL-generic-tracers/issues/29
# add_compile_definitions(_USE_MOM6_DIAG)

## Fortran modules path; currently this is simply set to be the include dir
set(CMAKE_INSTALL_MODULEDIR ${CMAKE_INSTALL_INCLUDEDIR}
  CACHE STRING
  "Fortran module installation path (Not a cmake native variable)"
)

#[==============================================================================[
#                              External packages                                #
#]==============================================================================]

find_package(fms COMPONENTS R8 REQUIRED)
find_package(GFDLGTracers REQUIRED)

if(NOT TARGET NetCDF::NetCDF_Fortran OR NOT TARGET NetCDF::NetCDF_C) #fms provides NETCDF inteface
  message(FATAL_ERROR "NetCDF interface missing from FMS package")
endif()

if(MOM6_ACCESS3)
  find_package(Access3Share REQUIRED cdeps timing share nuopc_cap_share)
  if(NOT TARGET ESMF::ESMF)
    message(FATAL_ERROR "ESMF interface missing from Access3Share package")
  endif()
endif()

if(MOM6_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# Check whether the OpenMP loop construct supports the bind clause, which the
# GPU-offloaded loops request as `!$omp target teams loop LOOP_BIND_TEAMS_PARALLEL`
# (see src/framework/omp_loop_bind_compat.h).  Mirrors the autoconf
# MOM6_FC_OMP_LOOP_BIND check (ac/m4/mom6_fc_omp_loop_bind.m4) and uses the same
# test program.  Without this define the macro expands to nothing, and NVHPC then
# launches the loop with one thread per team.  The directives are comments unless
# OpenMP is enabled, so the check is only done then, and it has to be run with the
# OpenMP flags that the sources are compiled with.
if(MOM6_OPENMP OR MOM6_OPENMP_OFFLOAD)
  set(CMAKE_REQUIRED_FLAGS "")
  if(MOM6_OPENMP)
    string(APPEND CMAKE_REQUIRED_FLAGS " ${OpenMP_Fortran_FLAGS}")
  endif()
  if(MOM6_OPENMP_OFFLOAD)
    string(APPEND CMAKE_REQUIRED_FLAGS " ${MOM6_OFFLOAD_FC_FLAGS}")
  endif()
  # Like the autoconf check, only compile the test program.
  if(DEFINED CMAKE_TRY_COMPILE_TARGET_TYPE)
    set(_mom6_try_compile_target_type ${CMAKE_TRY_COMPILE_TARGET_TYPE})
  endif()
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
  include(CheckFortranSourceCompiles)
  check_fortran_source_compiles(
"program test
  implicit none
  integer :: i
  !$omp target teams
  !$omp loop bind(teams,parallel)
  do i=1,2
  end do
  !$omp end target teams
end program test"
    HAVE_FC_OMP_LOOP_BIND
    SRC_EXT F90
  )
  if(DEFINED _mom6_try_compile_target_type)
    set(CMAKE_TRY_COMPILE_TARGET_TYPE ${_mom6_try_compile_target_type})
    unset(_mom6_try_compile_target_type)
  else()
    unset(CMAKE_TRY_COMPILE_TARGET_TYPE)
  endif()
  unset(CMAKE_REQUIRED_FLAGS)
  if(HAVE_FC_OMP_LOOP_BIND)
    add_compile_definitions(HAVE_FC_OMP_LOOP_BIND=1)
  endif()
endif()

#[==============================================================================[
#                             Main definitions                                  #
#]==============================================================================]

## Build MOM6 shared library
add_subdirectory(src)

## Build MOM6 ACCESS3 library and/or MOM6 solo executable
add_subdirectory(config_src)
