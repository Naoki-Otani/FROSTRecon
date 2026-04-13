################ Generic settings ################

# Be verbose
# set(CMAKE_VERBOSE_MAKEFILE ON)

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "DO NOT BUILD in-tree.")
endif ()

# CMAKE 3.10.2 is the default version in the Ubuntu 18.04 repositories
cmake_minimum_required(VERSION 3.10.2)

# Name of the project
project(NinjaReconstruction
        DESCRIPTION "Track reconstruction and track matching software for the NINJA Physics Run A"
        LANGUAGES CXX)

# Summary of Build types :
# RelWithDebInfo [-O2 –g] (default)
# Debug [-g -O0 -ggdb]
# Release [-O3]

# Build type
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
            "MinSizeRel" "RelWithDebInfo")
endif ()

# The version number
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 1)

set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
add_compile_definitions(VERSION=${PROJECT_VERSION})

# Only UNIX OSes
if (NOT UNIX)
    message(FATAL_ERROR "Non UNIX OSes are not supported")
endif ()

# Additional include directories for MacOS
if (APPLE)
    include_directories(/usr/local/include)
    include_directories(/opt/local/include)
    link_directories(/usr/local/lib)
    link_directories(/opt/local/lib)
endif ()

# Set default install path to the same path of the source code
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH "default install path" FORCE)
endif ()

# Set install paths for binaries, libraries and headers
set(CMAKE_INSTALL_BINDIR ${CMAKE_INSTALL_PREFIX}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib/ninja/recon)
set(CMAKE_INSTALL_INCDIR ${CMAKE_INSTALL_PREFIX}/include/ninja/recon)

################ Find packages ################

# where to look first for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

if (${CANVAS})
add_definitions(-DCANVAS=ON)
endif()

# ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT 6.18 CONFIG REQUIRED)
include(${ROOT_USE_FILE})
message(STATUS "Found ROOT: ${ROOT_INCLUDE_DIRS}")

# Geant4
# for a nice description of how to include Geant4 in a CMake project
# refer to this web page: http://www.sixiangguo.net/code/geant4/AppDevelop/apas02.html
find_package(Geant4 REQUIRED)
find_package(Geant4 REQUIRED COMPONENTS ui_tcsh)
include(${Geant4_USE_FILE})
message(STATUS "Found Geant4: ${Geant4_INCLUDE_DIRS}")
if (${Geant4_ui_tcsh_FOUND})
    message(STATUS "Found Geant4 UI tcsh: ${Geant4_ui_tcsh_FOUND}")
endif ()

# boost program_options, system, filesystem and log are required
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost 1.73 COMPONENTS system filesystem program_options unit_test_framework log_setup log REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost library -- found")
    message(STATUS "Boost library include dir : ${Boost_INCLUDE_DIRS}")
    include_directories("${Boost_INCLUDE_DIRS}")
endif()

# B2MC library
find_package(B2MC 0.1.10 REQUIRED)
if (B2MC_FOUND)
    message(STATUS "B2MC library -- found")
    message(STATUS "B2MC library include dir : ${B2MC_INCLUDE_DIR}")
    include_directories("${B2MC_INCLUDE_DIR}")
endif ()

# Event display library
#find_package(EventDisplay 0.0.5 REQUIRED)
#if (EVENT_DISPLAY_FOUND)
#    message(STATUS "Event display library -- found")
#    message(STATUS "Event display library include dir : ${EVENT_DISPLAY_INCLUDE_DIR}")
#    include_directories("${EVENT_DISPLAY_INCLUDE_DIR}")
#endif ()

# Wagasci reconstruction library
#find_package(WagasciRecon 0.0.1 REQUIRED)
#if(WAGASCI_RECON_FOUND)
#    message(STATUS "Wagasci reconstruction library -- found")
#    message(STATUS "Wagasci reconstruction library include dir : ${WAGASCI_RECON_INCLUDE_DIR}")
#    include_directories("${WAGASCI_EVENT_INCLUDE_DIR}")
#endif()

################ Libraries ################

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(NINJA_RECON_LIBRARIES "")
if (NOT APPLE)
    list(APPEND NINJA_RECON_LIBRARIES
            ${ROOT_LIBRARIES}
            ${Boost_LIBRARIES}
            ${B2MC_LIBRARY}
            ${EVENT_DISPLAY_LIBRARY}
#	    ${WAGASCI_RECON_LIBRARY}
            )
else ()
    list(APPEND NINJA_RECON_LIBRARIES
            ${ROOT_LIBRARIES}
            Boost::program_options
            Boost::system
            Boost::filesystem
            Boost::log_setup
            Boost::log
            ${B2MC_LIBRARY}
            ${EVENT_DISPLAY_LIBRARY}
#	    ${WAGASCI_RECON_LIBRARY}
            )
endif ()

################ Compiler flags ################

# Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)
add_definitions(-DBOOST_TEST_DYN_LINK)
add_definitions(-DBOOST_LOG_DYN_LINK)

message(STATUS "Compiler ID : ${CMAKE_CXX_COMPILER_ID}")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wno-overloaded-virtual -Wno-shadow")
elseif (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    # Add -O0 to remove optimizations when using gcc and debugging
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -ggdb -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0")
    # show all warnings -Wall
    # be pedantic -pedantic
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wno-overloaded-virtual -Wno-shadow \
     -Wno-vla -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-sign-compare -Wno-deprecated-copy \
     -Wno-unused-variable -Wno-return-type -Wno-comment")
endif ()

################ Linker flags ################

set(CMAKE_INSTALL_PPATH "${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(NINJA_RECON_LIBRARIES ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${B2MC_LIBRARY} ${EVENT_DISPLAY_LIBRARY})

################ Sub Directories ################

add_subdirectory(lib)
add_subdirectory(src/HitConverter)
add_subdirectory(src/TrackMatch)
