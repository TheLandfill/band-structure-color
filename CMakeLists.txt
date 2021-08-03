cmake_minimum_required(VERSION 3.0.0)
project(energy-gap-color)

include_directories(${CMAKE_SOURCE_DIR}/includes)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

find_package(OpenMP REQUIRED)

file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/common/*.c*")
add_library(color-from-absorption-spectrum SHARED ${SOURCES})
target_link_libraries(color-from-absorption-spectrum ${PNG_LIBRARY} OpenMP::OpenMP_CXX)

add_executable(energy-gap-color "${CMAKE_SOURCE_DIR}/src/band-gap-color.cpp")
target_link_libraries(energy-gap-color color-from-absorption-spectrum)

add_executable(acceptor-color-2D "${CMAKE_SOURCE_DIR}/src/2D-acceptor-color.cpp")
target_link_libraries(acceptor-color-2D color-from-absorption-spectrum)

add_executable(visible-spectrum "${CMAKE_SOURCE_DIR}/src/visible-spectrum.cpp")
target_link_libraries(visible-spectrum color-from-absorption-spectrum)
