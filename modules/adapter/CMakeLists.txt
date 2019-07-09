cmake_minimum_required(VERSION 3.2)
project(lib-geoc-adapter)

find_package(DGtal REQUIRED)

add_library(geoc-adapter INTERFACE)

if(${USE_LOCAL_REPOSITORIES})
    add_dependencies(geoc-adapter dipacus gcurve)
endif()

target_include_directories(geoc-adapter INTERFACE
        include
        ${CMAKE_SOURCE_DIR}/estimator/include
        ${DIPACUS_INCLUDE_DIRS}
        ${GCURVE_INCLUDE_DIRS}
        ${DGTAL_INCLUDE_DIRS})


install(DIRECTORY include
        DESTINATION ${CMAKE_INSTALL_PREFIX})

