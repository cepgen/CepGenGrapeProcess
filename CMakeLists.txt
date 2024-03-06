cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(CepGenGrapeProcess)
enable_language(Fortran)

set(CEPGEN_DIRS $ENV{CEPGEN_PATH} /usr /usr/local)
find_library(CEPGEN CepGen HINTS ${CEPGEN_DIRS} PATH_SUFFIXES build)
find_path(CEPGEN_INCLUDE CepGen/Generator.h HINTS ${CEPGEN_DIRS})

if(NOT CEPGEN OR NOT CEPGEN_INCLUDE)
  message(FATAL_ERROR "CepGen path was not found. Did you set the $CEPGEN_DIR environment variable?")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CEPGEN_INCLUDE}/cmake)
include(BuildUtils)

file(GLOB SOURCES *.cpp)
file(GLOB_RECURSE F77_SOURCES src/*.f)

add_library(CepGenGrape SHARED ${SOURCES} ${F77_SOURCES})
target_link_libraries(CepGenGrape PRIVATE ${CEPGEN})
target_include_directories(CepGenGrape PRIVATE include ${CEPGEN_INCLUDE})
target_compile_options(CepGenGrape PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-132 -Wno-line-truncation -Wno-conversion -Wno-unused-label -Wno-unused-dummy-argument -Wno-maybe-uninitialized>)
