cmake_minimum_required(VERSION 2.6)
include(FindPkgConfig)

project(sasir)

#
# Find libcfitsio
#
# Locate cfitsio using pkg-config
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)
find_package(CFITSIO REQUIRED)

set(BIN_DIR  bin)
if (DEFINED ENV{ISAP})  
set(BIN_DIR  $ENV{ISAP}/bin)
endif  (DEFINED ENV{ISAP}) 

#
# Compilation flags
#
if(APPLE)
  set(CMAKE_CXX_FLAGS "-DMACOS -DNO_DISP_IO -g0 -O2 -ffast-math -fomit-frame-pointer -fopenmp -Wno-write-strings -DNDEBUG")
#  set(CMAKE_CXX_FLAGS "-DMACOS -g -Wno-write-strings -DNDEBUG")
else(APPLE)
  set(CMAKE_CXX_FLAGS "-DNO_DISP_IO  -g0 -O2 -ffast-math -fomit-frame-pointer -fopenmp -Wno-write-strings -DNDEBUG")
endif(APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -lcfitsio")

#
# Build tools library
#
FILE(GLOB src_lib1  "${PROJECT_SOURCE_DIR}/src/libtools/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/libtools")
add_library(tools STATIC ${src_lib1})
target_link_libraries(tools ${CFITSIO_LIBRARIES})

#
# Build LOFAR library
#
FILE(GLOB src_lib2  "${PROJECT_SOURCE_DIR}/src/liblofar/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/liblofar")
add_library(lofar STATIC ${src_lib2})
target_link_libraries(lofar ${CFITSIO_LIBRARIES})


########################### SASIR INFO  ########################### 
# SASIR PACKAGE
FILE(GLOB sasir_src_lib  "${PROJECT_SOURCE_DIR}/src/libsasir/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/libsasir")
add_library(MRsasir STATIC ${sasir_src_lib})


set(SASIR_BIN  sasir)
foreach(program ${SASIR_BIN})
  add_executable(${program} ${PROJECT_SOURCE_DIR}/src/${program}.cc)
  target_link_libraries(${program} MRsasir tools)
endforeach(program)

set(LOFAR_BIN   lofar_test) 
foreach(program ${LOFAR_BIN})
  add_executable(${program} ${PROJECT_SOURCE_DIR}/src/${program}.cc)
  target_link_libraries(${program} lofar  MRsasir tools)
endforeach(program)

#
# Install (by default in the project directory)
#

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

#install sparse1d headers
FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libtools/*.h")
INSTALL(FILES ${inc_lib} DESTINATION include)
FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libsasir/*.h")
INSTALL(FILES ${inc_lib} DESTINATION include)
FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/liblofar/*.h")
INSTALL(FILES ${inc_lib} DESTINATION include)

# Install sparse1d library
INSTALL(TARGETS MRsasir tools DESTINATION lib)

# install sparse1d executables
INSTALL(TARGETS ${SASIR_BIN} ${LOFAR_BIN} DESTINATION ${BIN_DIR})

INCLUDE_DIRECTORIES(${CFITSIO_INCLUDE_DIR})

set(OTHER_LIBRARIES ${CFITSIO_LIBRARIES})

