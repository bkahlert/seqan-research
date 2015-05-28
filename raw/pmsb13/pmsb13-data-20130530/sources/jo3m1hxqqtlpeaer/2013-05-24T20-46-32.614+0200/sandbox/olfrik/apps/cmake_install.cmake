# Install script for directory: /home/olfrik/Projects/seqan/src/sandbox/olfrik/apps

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
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t2SeqA1/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t3Ala2/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/first_app2/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t5Graphs/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t1SeqA2/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t3Io/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t3BuildFai/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/first_app/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/app1/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/mini_bowtie/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/workshop/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t3AlGloA1/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t3Ala1/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t3Align/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t4Patterns/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/t4Index/cmake_install.cmake")
  INCLUDE("/home/olfrik/Projects/seqan/src/sandbox/olfrik/apps/the_example/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

