cmake_minimum_required(VERSION 2.6)
PROJECT(FLEXBAR)

#file(MAKE_DIRECTORY build)
set(EXECUTABLE_OUTPUT_PATH ${FLEXBAR_BINARY_DIR})

ADD_SUBDIRECTORY(src)

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

