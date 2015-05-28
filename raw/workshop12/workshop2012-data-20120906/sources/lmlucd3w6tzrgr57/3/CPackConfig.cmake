# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


SET(CPACK_BINARY_BUNDLE "")
SET(CPACK_BINARY_CYGWIN "")
SET(CPACK_BINARY_DEB "")
SET(CPACK_BINARY_DRAGNDROP "")
SET(CPACK_BINARY_NSIS "")
SET(CPACK_BINARY_OSXX11 "")
SET(CPACK_BINARY_PACKAGEMAKER "")
SET(CPACK_BINARY_RPM "")
SET(CPACK_BINARY_STGZ "")
SET(CPACK_BINARY_TBZ2 "")
SET(CPACK_BINARY_TGZ "")
SET(CPACK_BINARY_TZ "")
SET(CPACK_BINARY_ZIP "")
SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
SET(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
SET(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
SET(CPACK_GENERATOR "ZIP;TBZ2;DEB")
SET(CPACK_INSTALL_CMAKE_PROJECTS "/Users/grunert/Development/seqan-trunk;seqan;ALL;/")
SET(CPACK_INSTALL_PREFIX "/usr/local")
SET(CPACK_MODULE_PATH "/Users/grunert/Development/seqan-trunk/util/cmake/FindTBB;/Users/grunert/Development/seqan-trunk/util/cmake/cuda;/Users/grunert/Development/seqan-trunk/util/cmake")
SET(CPACK_NSIS_DISPLAY_NAME "CMake .")
SET(CPACK_NSIS_INSTALLER_ICON_CODE "")
SET(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
SET(CPACK_NSIS_PACKAGE_NAME "CMake .")
SET(CPACK_OUTPUT_CONFIG_FILE "/Users/grunert/Development/seqan-trunk/CPackConfig.cmake")
SET(CPACK_PACKAGE_DEFAULT_LOCATION "/")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "/Users/grunert/Development/seqan-trunk/README")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeqAn - The C++ library for sequence analysis.")
SET(CPACK_PACKAGE_FILE_NAME "seqan-pre1.4.0-Darwin")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake .")
SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "CMake .")
SET(CPACK_PACKAGE_NAME "seqan")
SET(CPACK_PACKAGE_RELOCATABLE "true")
SET(CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
SET(CPACK_PACKAGE_VERSION "pre1.4.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "4")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_RESOURCE_FILE_LICENSE "/Users/grunert/Development/seqan-trunk/LICENSE")
SET(CPACK_RESOURCE_FILE_README "/usr/usr/share/cmake-2.8/Templates/CPack.GenericDescription.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "/usr/usr/share/cmake-2.8/Templates/CPack.GenericWelcome.txt")
SET(CPACK_SET_DESTDIR "OFF")
SET(CPACK_SOURCE_CYGWIN "")
SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2;TZ")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/Users/grunert/Development/seqan-trunk/CPackSourceConfig.cmake")
SET(CPACK_SOURCE_TBZ2 "ON")
SET(CPACK_SOURCE_TGZ "ON")
SET(CPACK_SOURCE_TZ "ON")
SET(CPACK_SOURCE_ZIP "OFF")
SET(CPACK_SYSTEM_NAME "Darwin")
SET(CPACK_TOPLEVEL_TAG "Darwin")