cmake_minimum_required(VERSION 3.19)
project(packing_solution)

set(platform x86)
set(CMAKE_CXX_STANDARD 17)

include_directories(./)
include_directories(./include)
include_directories(./CarvingMachine)
#include_directories(./heuristic)

link_directories(./lib)

aux_source_directory(./ SRC_FILES)
#aux_source_directory(./heuristic HEURISTIC_FILES)
aux_source_directory(./CarvingMachine CARVINGMACHINE_FILES)
#aux_source_directory(./marching MARCHING_FILES)

INCLUDE(FindPkgConfig)
pkg_check_modules(JSONCPP REQUIRED jsoncpp)
include_directories(${JSONCPP_INCLUDE_DIRS})
link_directories(${JSONCPP_LIBRARY_DIRS})

add_executable(test_carvingMarchine ${SRC_FILES} ${CARVINGMACHINE_FILES} ${HEURISTIC_FILES} ${MARCHING_FILES} amopt_pack_base.cpp amopt_pack_base.h)
target_link_libraries(test_carvingMarchine jsoncpp)