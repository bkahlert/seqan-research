# Install script for directory: /home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
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

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/sequence_1/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/Indices2/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/qgram1/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/PM1/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/sequence_2/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/align1/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/graphs1/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/the_example/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/sequence_workshop_1/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/align2/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/IndexIt2/cmake_install.cmake")
  INCLUDE("/home/dkersting/Development/seqan-trunk/sandbox/dkersting/apps/fabi/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

