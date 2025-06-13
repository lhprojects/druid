##########################################################
# Date: 2023-05-28
# Yuzhi CHE
#
# The Original GNUMakefile does not work in newer ROOT and
# LCIO version, since 2020.
# This CMakeLists.txt builds Druid under ROOT-6.24 and 
# LCIO-2.16 in CentOS7 and MacOS 12.6.6
# Note that one should go into ROOT & LCIO surrounding (
# by `source $LCIO/setup.sh` and `source $ROOTSYS/bin/thisroot.sh`
# ) or LCG releases (version > 100) before compiling. 
##########################################################

cmake_minimum_required(VERSION 2.0 FATAL_ERROR)
project(Druid)

# Set up compiler
set(CMAKE_CXX_COMPILER "g++")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LIBRARY_PATH_VAR "LD_LIBRARY_PATH")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(LIBRARY_PATH_VAR "DYLD_LIBRARY_PATH")
endif()

# find ROOT and LCIO
find_package(ROOT COMPONENTS Gui Core Imt RIO Net Hist Graf Graf3d Gpad ROOTVecOps Tree TreePlayer Rint Postscript Matrix Physics MathCore Thread MultiProc ROOTDataFrame Geom Ged RGL EG Eve REQUIRED)
find_package(LCIO REQUIRED)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})

# Execute rootcint to generate Class Dictionary.
execute_process(COMMAND ${ROOT_rootcint_CMD} -f ${CMAKE_CURRENT_SOURCE_DIR}/src/EventNavigatorDict.cc -c EventNavigator.hh
								WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/EventNavigatorDict_rdict.pcm
								${CMAKE_INSTALL_PREFIX}/bin/EventNavigatorDict_rdict.pcm
								COPYONLY)

# Generate setup.sh
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in
								${CMAKE_INSTALL_PREFIX}/bin/setup.sh
								@ONLY)

# add include PATH
include_directories(${ROOT_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS} include)

# add source code
file(GLOB SOURCES "src/*.cc")

# add dependency
set(LIBRARIES
  ${ROOT_LIBRARIES}
  ${LCIO_LIBRARIES}
)

# Generate executable file
add_executable(Druid ${SOURCES})

# linke library
target_link_libraries(Druid ${LIBRARIES})

# install
install(TARGETS Druid DESTINATION bin)
