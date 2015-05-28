# Install script for directory: D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/seqan")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
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

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/eigth_app/cmake_install.cmake")
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/first_app/cmake_install.cmake")
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/five_app/cmake_install.cmake")
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/fourth_app/cmake_install.cmake")
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/nine.app/cmake_install.cmake")
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/second_app/cmake_install.cmake")
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/seven_app/cmake_install.cmake")
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/six_app/cmake_install.cmake")
  INCLUDE("D:/Informatik/Development/SeqAn/sandbox/mysandbox/apps/third_app/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

