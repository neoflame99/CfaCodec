cmake_minimum_required(VERSION 3.16)

set ( PROJECT_NAME "CfaCodec" )
PROJECT( ${PROJECT_NAME} )

set( CMAKE_EXPORT_COMPILE_COMMANDS 1)
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O0 -g3 ")

include_directories(
    ${PROJECT_SOURCE_DIR}/../inc
    ${PROJECT_SOURCE_DIR}
    )
    
link_directories(
    #$ENV{SYSTEMC_HOME}/lib
    )
add_executable(
    ${PROJECT_NAME}  # Executable file name
    ${PROJECT_SOURCE_DIR}/main.cpp
    ${PROJECT_SOURCE_DIR}/../src/loadparam.cpp
    ${PROJECT_SOURCE_DIR}/../src/imgloader.cpp
    ${PROJECT_SOURCE_DIR}/../src/msst.cpp
    ${PROJECT_SOURCE_DIR}/../src/bayer_dwt.cpp
    ${PROJECT_SOURCE_DIR}/../src/monitor.cpp
    ${PROJECT_SOURCE_DIR}/../src/dwt.cpp
    ${PROJECT_SOURCE_DIR}/../src/quantize.cpp
    ${PROJECT_SOURCE_DIR}/../src/entropy.cpp
    # other source file list
    )

target_include_directories(${PROJECT_NAME} PRIVATE ../src)

#target_link_libraries(
    #${PROJECT_NAME}
    #systemc
#    )
    
