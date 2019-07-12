cmake_minimum_required(VERSION 3.2)

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(BOOST_INCLUDE_DIRS /usr/local/include/boost CACHE PATH "Boost include location path." )
set(BOOST_LIB_DIR /usr/local/lib CACHE PATH "Boost libraries location path." )

set(DIPACUS_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/ext-projects/build/include CACHE PATH "DIPaCUS include location path.")
set(DIPACUS_LIBS_DIR ${CMAKE_SOURCE_DIR}/ext-projects/build/lib CACHE PATH "DIPaCUS libraries location path.")

set(GCURVE_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/ext-projects/build/include CACHE PATH "DIPaCUS include location path.")
set(GCURVE_LIBS_DIR ${CMAKE_SOURCE_DIR}/ext-projects/build/lib CACHE PATH "DIPaCUS libraries location path.")

set(BUILD_TESTS OFF CACHE BOOL "Build tests.")
set(USE_REMOTE_REPOSITORIES ON CACHE BOOL "Pull external repositories from GitHub.")


project(geoc)

message("Project Root Dir: ${CMAKE_SOURCE_DIR}")
message("Build files at: ${CMAKE_BINARY_DIR}")
message("Installation at: ${CMAKE_INSTALL_PREFIX}")

if(${USE_REMOTE_REPOSITORIES})
    include(ext-projects/ExtProjects.cmake)
endif()

find_package(DGtal REQUIRED)


link_directories(${CMAKE_BINARY_DIR}/modules/estimator
        ${CMAKE_BINARY_DIR}/modules/adapter
        ${CMAKE_BINARY_DIR}/modules/energy
        ${CMAKE_SOURCE_DIR}/ext/lib
        ${DIPACUS_LIBS_DIR}
        ${GCURVE_LIBS_DIR}
        ${BOOST_LIBS_DIR})

add_subdirectory(modules/estimator)
add_subdirectory(modules/adapter)
add_subdirectory(modules/api)

if (${BUILD_TESTS})
    add_subdirectory(lab/test)
endif()
