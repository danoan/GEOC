cmake_minimum_required(VERSION 3.2)
project(lib-geoc-adapter)

find_package(DGtal REQUIRED)

add_library(geoc-adapter INTERFACE)

target_include_directories(geoc-adapter INTERFACE
        include
        ../estimator/include
        ../../ext/include
        ${DGTAL_INCLUDE_DIRS})


install(DIRECTORY include
        DESTINATION ${CMAKE_INSTALL_PREFIX})

