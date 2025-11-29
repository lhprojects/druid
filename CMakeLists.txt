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
set(CMAKE_CXX_STANDARD_REQUIRED ON)



# find ROOT and LCIO
find_package(ROOT COMPONENTS Gui Core Imt RIO Net Hist Graf Graf3d Gpad ROOTVecOps Tree TreePlayer Rint Postscript Matrix Physics MathCore Thread MultiProc ROOTDataFrame Geom Ged RGL EG Eve REQUIRED)

# Find GSL
find_package(GSL REQUIRED)

# Try to find LCIO
find_path(LCIO_INCLUDE_DIRS lcio.h)
find_library(LCIO_LIBRARIES NAMES lcio)

if(LCIO_INCLUDE_DIRS AND LCIO_LIBRARIES)
    message(STATUS "Found system LCIO")
    add_library(lcio INTERFACE)
    target_link_libraries(lcio INTERFACE ${LCIO_LIBRARIES})
else()
    message(STATUS "LCIO not found, downloading and building...")

    include(ExternalProject)
    ExternalProject_Add(
        LCIO_ext
        PREFIX ${CMAKE_BINARY_DIR}/_deps/lcio
        GIT_REPOSITORY https://github.com/iLCSoft/LCIO.git
        GIT_TAG master
        UPDATE_COMMAND ""
        INSTALL_DIR ${CMAKE_BINARY_DIR}/lcio_install
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_CXX_STANDARD=17
	  -DCMAKE_C_FLAGS=-Wno-deprecated-non-prototype
)

    set(LCIO_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/lcio_install/include)

    if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        set(LCIO_LIBRARY_PATH ${CMAKE_BINARY_DIR}/lcio_install/lib/liblcio.dylib)
    elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        set(LCIO_LIBRARY_PATH ${CMAKE_BINARY_DIR}/lcio_install/lib/liblcio.so)
    else()
        message(FATAL_ERROR "Unsupported platform: ${CMAKE_SYSTEM_NAME}")
    endif()

    add_library(lcio SHARED IMPORTED GLOBAL)
    set_target_properties(lcio PROPERTIES
    IMPORTED_LOCATION ${LCIO_LIBRARY_PATH}
)

    add_dependencies(lcio LCIO_ext)


endif()


set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})

# Execute rootcint to generate Class Dictionary.
execute_process(COMMAND ${ROOT_rootcint_CMD} -f ${CMAKE_CURRENT_SOURCE_DIR}/src/EventNavigatorDict.cc -c EventNavigator.hh
								WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/EventNavigatorDict_rdict.pcm
								${CMAKE_INSTALL_PREFIX}/bin/EventNavigatorDict_rdict.pcm
								COPYONLY)


# add include PATH
include_directories(${ROOT_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} include thirdparty/MLPFA/core/include thirdparty/MLPFA/thirdparty/PFAA/include)
include_directories(${ROOT_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS} include thirdparty/MLPFA/thirdparty/PFATime/include )


# add source code
file(GLOB SOURCES "src/*.cc")

# add dependency
set(LIBRARIES
  ${ROOT_LIBRARIES}
  ${GSL_LIBRARIES}
)


# Use RPATH in the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# Add paths where shared libraries are installed
list(APPEND CMAKE_INSTALL_RPATH
    "${CMAKE_BINARY_DIR}/lcio_install/lib"
)

# On Linux, add origin-relative paths
if(UNIX AND NOT APPLE)
    list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

# On macOS, add @loader_path-relative fallback
if(APPLE)
    list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../lib")
endif()



# Generate executable file
add_executable(Druid ${SOURCES})

# linke library
target_link_libraries(Druid PRIVATE ${LIBRARIES} lcio)


# install
install(TARGETS Druid DESTINATION bin)
install(DIRECTORY icons DESTINATION bin)
