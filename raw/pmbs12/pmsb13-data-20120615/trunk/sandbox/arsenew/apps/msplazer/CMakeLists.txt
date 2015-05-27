cmake_minimum_required (VERSION 2.6)
project (sandbox_arsenew_apps_msplazer)

add_definitions(-pg -g)
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# Update the list of file names below if you add source files to your application.
seqan_add_executable(msplazer msplazer_rna.cpp msplazer.h)
# app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_msplazer COMMAND ${PYTHON_EXECUTABLE}
	${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)
