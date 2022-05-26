cmake_minimum_required(VERSION 2.8)
project(MEX)

if(CMAKE_VERSION VERSION_LESS "3.18")
  message(WARNING "CMake version less than 3.18 detected"
    " -- using local FindMatlab module extracted from v3.18")
  # Where we keep the FindMatlab module
  set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-compat ${CMAKE_MODULE_PATH})
endif()

if(APPLE)
    set(CC "/usr/local/bin/gcc-9")
    set(CXX "/usr/local/bin/g++-9")

    set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")

    set(Matlab_ROOT_DIR "/Applications/MATLAB_R2021a.app")
    set(CMAKE_MACOSX_RPATH 1)

endif(APPLE)

find_package(Matlab REQUIRED MAIN_PROGRAM MX_LIBRARY)
if(NOT Matlab_FOUND)
  message(FATAL_ERROR "MATLAB not found.")
endif()

# NOTE: Ewald code for K0 requires GSL

#find_package(GSL REQUIRED)
#if(NOT GSL_FOUND)
#  message(FATAL_ERROR "GSL not found.")
#else()
#  message("GSL found.")
#endif()


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")

# Assuming gcc for now
# with parallelization
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O0 -fopenmp -msse3 -falign-loops=16 -DMEX_DOUBLE_HANDLE")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O0 -fopenmp -falign-loops=16 -msse4.1 -DMEX_DOUBLE_HANDLE")

# without paralellization
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O0 -msse3 -falign-loops=16 -DMEX_DOUBLE_HANDLE")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O0 -falign-loops=16 -msse4.1 -DMEX_DOUBLE_HANDLE")

include_directories(
  ${Matlab_INCLUDE_DIRS}
)

# Add modules with MEX to be built
add_subdirectory("${PROJECT_SOURCE_DIR}/mex/StokesSLP")
add_subdirectory("${PROJECT_SOURCE_DIR}/mex/StokesDLP")
