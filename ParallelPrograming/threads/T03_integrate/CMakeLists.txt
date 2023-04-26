cmake_minimum_required( VERSION 3.11 )

project( Integrate VERSION 0.1 LANGUAGES CXX )

set( INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include" )
set( SOURCE_DIR "${CMAKE_SOURCE_DIR}/src" )

set( EXEC_NAME integrate )
add_executable( ${EXEC_NAME} )

target_include_directories( ${EXEC_NAME} PRIVATE ${INCLUDE_DIR} )
target_sources( ${EXEC_NAME} PRIVATE "${SOURCE_DIR}/main.cpp" )

target_compile_features( ${EXEC_NAME} PRIVATE cxx_std_20 )

set( COMMON_FLAGS -Wall -Wextra -Wpedantic -Werror )

target_compile_options( ${EXEC_NAME} PRIVATE ${COMMON_FLAGS} )
