cmake_minimum_required(VERSION 3.2)
project(geoc_api)

find_package(DGtal REQUIRED)

add_library(geoc-api INTERFACE)

add_dependencies(geoc-api geoc-adapter geoc-estimator)

target_include_directories(geoc-api INTERFACE
        include
        ${CMAKE_SOURCE_DIR}/modules/adapter/include
        ${CMAKE_SOURCE_DIR}/modules/estimator/include
        ${DGTAL_INCLUDE_DIRS}
        )

install(DIRECTORY include
        DESTINATION ${CMAKE_INSTALL_PREFIX})