cmake_minimum_required(VERSION 3.2)
project(lib-geoc-estimators)

find_package(DGtal REQUIRED)

add_library(geoc-estimator INTERFACE)

if(${USE_REMOTE_REPOSITORIES})
    add_dependencies(geoc-estimator dipacus gcurve)
endif()

target_include_directories(geoc-estimator INTERFACE
        include
        ${DIPACUS_INCLUDE_DIRS}
        ${GCURVE_INCLUDE_DIRS}
        ${DGTAL_INCLUDE_DIRS})


install(DIRECTORY include
        DESTINATION ${CMAKE_INSTALL_PREFIX})