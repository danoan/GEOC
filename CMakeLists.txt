cmake_minimum_required(VERSION 3.2)
project(geoc)

#CMAKE_INSTALL_PREFIX: install location

message("Project Root Dir: ${CMAKE_SOURCE_DIR}")
message("Build files at: ${CMAKE_BINARY_DIR}")
message("Installation at: ${CMAKE_INSTALL_PREFIX}")

find_package(DGtal REQUIRED)

set(BOOST_INCLUDE_DIRS /usr/local/include/boost CACHE PATH "Boost include location path." )

set(DIPACUS_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/ext/include CACHE PATH "DIPaCUS include location path.")
set(DIPACUS_LIBS_DIR /usr/local/lib CACHE PATH "DIPaCUS libraries location path.")

set(GCURVE_INCLUDE_DIRS /usr/local/include/gcurve CACHE PATH "gcurve include location path.")
set(GCURVE_LIBS_DIR /usr/local/lib CACHE PATH "gcurve librarie location path.")

set(BUILD_TESTS OFF CACHE BOOL "Build tests.")

link_directories(${CMAKE_BINARY_DIR}/modules/estimator
        ${CMAKE_BINARY_DIR}/modules/adapter
        ${CMAKE_BINARY_DIR}/modules/energy
        ${DIPACUS_LIBS_DIR}
        ${GCURVE_LIBS_DIR})

add_subdirectory(modules/estimator)
add_subdirectory(modules/adapter)

if (${BUILD_TESTS})
    add_subdirectory(lab/test)
endif()