cmake_minimum_required(VERSION 3.18)
project(osqp_matlab LANGUAGES C CXX)

include(FetchContent)

# -----------------------------------------------
# Options
# -----------------------------------------------
set(OSQP_ALGEBRA_BACKEND "builtin" CACHE STRING "Algebra backend (builtin/mkl/cuda)")
option(OSQP_USE_FLOAT    "Use single-precision floats" OFF)
option(OSQP_USE_LONG     "Use 64-bit integers"         ON)

# Pass through to OSQP sub-project
set(OSQP_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
set(OSQP_BUILD_STATIC_LIB ON  CACHE BOOL "" FORCE)
set(OSQP_BUILD_DEMO_EXE   OFF CACHE BOOL "" FORCE)
set(OSQP_BUILD_UNITTESTS  OFF CACHE BOOL "" FORCE)
set(OSQP_ENABLE_PRINTING    ON  CACHE BOOL "" FORCE)
set(OSQP_ENABLE_PROFILING  ON  CACHE BOOL "" FORCE)
set(OSQP_ENABLE_INTERRUPT  ON  CACHE BOOL "" FORCE)
set(OSQP_ENABLE_DERIVATIVES ON CACHE BOOL "" FORCE)
set(OSQP_CODEGEN           ON  CACHE BOOL "" FORCE)

# -----------------------------------------------
# Fetch OSQP v1.0.0 (brings in QDLDL v0.1.8 automatically)
# -----------------------------------------------
FetchContent_Declare(
    osqp
    GIT_REPOSITORY https://github.com/osqp/osqp.git
    GIT_TAG        v1.0.0
)
FetchContent_MakeAvailable(osqp)

# -----------------------------------------------
# Export paths for make_osqp.m to discover
# -----------------------------------------------
# After building, make_osqp.m reads these paths from the CMake cache.
set(OSQP_INC_DIR "${osqp_SOURCE_DIR}/include/public;${osqp_BINARY_DIR}/include/public"
    CACHE STRING "OSQP include directories" FORCE)
set(OSQP_LIB_DIR "${osqp_BINARY_DIR}/out"
    CACHE STRING "OSQP library output directory" FORCE)
