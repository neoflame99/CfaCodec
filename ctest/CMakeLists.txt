cmake_minimum_required(VERSION 3.5)
set ( PROJECT_NAME "BayerJpegCodec" )
PROJECT( ${PROJECT_NAME} )

set( CMAKE_EXPORT_COMPILE_COMMANDS 1)
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O0 -g3 ")

#message(STATUS "Systemc_home: $ENV{SYSTEMC_HOME}")
include_directories(
    ${PROJECT_SOURCE_DIR}/inc
    )
    
link_directories(
#    $ENV{SYSTEMC_HOME}/lib
    )
add_executable(
    ${PROJECT_NAME}  # Executable file name
    ${PROJECT_SOURCE_DIR}/main.cpp
    ${PROJECT_SOURCE_DIR}/src/dwt.cpp
    ${PROJECT_SOURCE_DIR}/src/imgloader.cpp
    ${PROJECT_SOURCE_DIR}/src/bayer_dwt.cpp
    ${PROJECT_SOURCE_DIR}/src/loadparam.cpp
    ${PROJECT_SOURCE_DIR}/src/monitor.cpp
    ${PROJECT_SOURCE_DIR}/src/quant.cpp
    ${PROJECT_SOURCE_DIR}/src/entropy.cpp
    # other source file list
    )
#target_link_libraries(
#    ${PROJECT_NAME}
#    systemc
#    uvm-systemc
#    )
    
