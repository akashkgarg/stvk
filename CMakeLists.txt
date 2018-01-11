cmake_minimum_required (VERSION 3.9.4)
project (stvk)

# Default values for CMAKE_BUILD_TYPE
#-------------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: Debug, Release, RelWithDebInfo" FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if ((NOT CMAKE_BUILD_TYPE MATCHES Debug) AND
    (NOT CMAKE_BUILD_TYPE MATCHES Release) AND
    (NOT CMAKE_BUILD_TYPE MATCHES RelWithDebInfo))
  message (FATAL_ERROR "Invalid CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif ((NOT CMAKE_BUILD_TYPE MATCHES Debug) AND
    (NOT CMAKE_BUILD_TYPE MATCHES Release) AND
    (NOT CMAKE_BUILD_TYPE MATCHES RelWithDebInfo))

# Compiler Flags
#-------------------------------------------------------------------------------

option(USE_CLANG "build application with clang" OFF) # OFF is the default

if (USE_CLANG)
  SET (CMAKE_C_COMPILER "clang")
  SET (CMAKE_CXX_COMPILER "clang++")
endif (USE_CLANG)

IF (NOT WIN32 AND NOT ANDROID)
  # Add common flags.
  set (CMAKE_CXX_FLAGS
      "-Wall -Wextra -fno-common -mfix-and-continue -stdlib=libc++ -std=c++11"
      CACHE STRING "" FORCE)

  # set debug flags
  set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG" CACHE STRING "" FORCE)

  # set release flags
  set (CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG" CACHE STRING "" FORCE)
ENDIF ()

IF (ANDROID)
  # common flags
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fno-rtti -fexceptions")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fexceptions")

  # set debug flags
  set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG" CACHE STRING "" FORCE)

  # set release flags
  set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG" CACHE STRING "" FORCE)
ENDIF ()

# Extra modules for cmake
#-------------------------------------------------------------------------------

# Find Necessary Packages
#-------------------------------------------------------------------------------

# Find Eigen3
set (EIGEN3_INCLUDE_DIR "/Users/akashgarg/projects/loomai/third-party/eigen/3.3.4/include")
include_directories (SYSTEM ${EIGEN3_INCLUDE_DIR})
message (STATUS "EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR}")

set (SRC_DIRECTORIES
  ${PROJECT_SOURCE_DIR}
)

include_directories (${PROJECT_SOURCE_DIR}/)

add_executable (stvk main.cc tri_mesh.cc obj_file.cc StVk.h tri_mesh.h types.h)
