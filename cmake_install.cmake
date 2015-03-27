# Install script for directory: /home/garsden/SASIR

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/home/garsden/SASIR")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/garsden/SASIR/src/libtools/IM_IO.h"
    "/home/garsden/SASIR/src/libtools/FFTN_1D.h"
    "/home/garsden/SASIR/src/libtools/IM_Color.h"
    "/home/garsden/SASIR/src/libtools/VMS.h"
    "/home/garsden/SASIR/src/libtools/IM_Sigma.h"
    "/home/garsden/SASIR/src/libtools/DefFunc.h"
    "/home/garsden/SASIR/src/libtools/NR_util.h"
    "/home/garsden/SASIR/src/libtools/IM_Obj.h"
    "/home/garsden/SASIR/src/libtools/OptMedian.h"
    "/home/garsden/SASIR/src/libtools/IM_Lut.h"
    "/home/garsden/SASIR/src/libtools/SoftInfo.h"
    "/home/garsden/SASIR/src/libtools/GenMedian.h"
    "/home/garsden/SASIR/src/libtools/IM_IOCol.h"
    "/home/garsden/SASIR/src/libtools/IM_Prob.h"
    "/home/garsden/SASIR/src/libtools/IM_Noise.h"
    "/home/garsden/SASIR/src/libtools/IM_Rot.h"
    "/home/garsden/SASIR/src/libtools/IM1D_IO.h"
    "/home/garsden/SASIR/src/libtools/GlobalInc.h"
    "/home/garsden/SASIR/src/libtools/Usage.h"
    "/home/garsden/SASIR/src/libtools/IM3D_IO.h"
    "/home/garsden/SASIR/src/libtools/TempMemory.h"
    "/home/garsden/SASIR/src/libtools/CErf.h"
    "/home/garsden/SASIR/src/libtools/FFTN_2D.h"
    "/home/garsden/SASIR/src/libtools/TempArray.h"
    "/home/garsden/SASIR/src/libtools/macro.h"
    "/home/garsden/SASIR/src/libtools/NR.h"
    "/home/garsden/SASIR/src/libtools/IM_IOTools.h"
    "/home/garsden/SASIR/src/libtools/Border.h"
    "/home/garsden/SASIR/src/libtools/IM_Deconv.h"
    "/home/garsden/SASIR/src/libtools/GetLongOptions.h"
    "/home/garsden/SASIR/src/libtools/DefMath.h"
    "/home/garsden/SASIR/src/libtools/writefits3d.h"
    "/home/garsden/SASIR/src/libtools/Memory.h"
    "/home/garsden/SASIR/src/libtools/IM_Math.h"
    "/home/garsden/SASIR/src/libtools/FFTN.h"
    "/home/garsden/SASIR/src/libtools/MatrixOper.h"
    "/home/garsden/SASIR/src/libtools/Array.h"
    "/home/garsden/SASIR/src/libtools/Licence.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/garsden/SASIR/src/libsasir/MCA.h"
    "/home/garsden/SASIR/src/libsasir/SB_Filter.h"
    "/home/garsden/SASIR/src/libsasir/Fista.h"
    "/home/garsden/SASIR/src/libsasir/LineCol.h"
    "/home/garsden/SASIR/src/libsasir/Atrou1D.h"
    "/home/garsden/SASIR/src/libsasir/Nesterov.h"
    "/home/garsden/SASIR/src/libsasir/FloatTrans.h"
    "/home/garsden/SASIR/src/libsasir/Pyr1D.h"
    "/home/garsden/SASIR/src/libsasir/SB_Filter1D.h"
    "/home/garsden/SASIR/src/libsasir/MeyerWT.h"
    "/home/garsden/SASIR/src/libsasir/Filter.h"
    "/home/garsden/SASIR/src/libsasir/FCur.h"
    "/home/garsden/SASIR/src/libsasir/SB_Filter_float.h"
    "/home/garsden/SASIR/src/libsasir/FCur_float.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/garsden/SASIR/src/liblofar/CEA_comp_sens.h"
    "/home/garsden/SASIR/src/liblofar/MCA_LOFAR.h"
    "/home/garsden/SASIR/src/liblofar/Fista_LOFAR.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/garsden/SASIR/libMRsasir.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/garsden/SASIR/libtools.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sasir")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sasir"
         RPATH "")
  ENDIF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sasir")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/garsden/SASIR/sasir")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sasir")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sasir")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sasir")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/lofar_test")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/lofar_test"
         RPATH "")
  ENDIF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/lofar_test")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/garsden/SASIR/lofar_test")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/lofar_test")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/lofar_test")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/lofar_test")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/garsden/SASIR/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/garsden/SASIR/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
