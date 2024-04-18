cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(CepGenGrapeProcess)
enable_language(Fortran)

set(CEPGEN_PATH ${CEPGEN_PATH} $ENV{CEPGEN_PATH})
set(CMAKE_MODULE_PATH ${CEPGEN_PATH}/cmake)
find_package(CepGen REQUIRED)

file(GLOB SOURCES *.cpp)
file(GLOB_RECURSE F77_SOURCES src/*.f)

add_library(CepGenGrape SHARED ${SOURCES} ${F77_SOURCES})
target_link_libraries(CepGenGrape PRIVATE ${CepGen_LIBRARIES})
target_include_directories(CepGenGrape PRIVATE include ${CepGen_INCLUDE_DIRS})
target_compile_options(CepGenGrape
                       PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-132
                               -Wno-line-truncation -Wno-conversion -Wno-unused-label -Wno-unused-dummy-argument -Wno-maybe-uninitialized>)

# copy the input cards
file(GLOB_RECURSE input_cards RELATIVE ${PROJECT_SOURCE_DIR} cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()

