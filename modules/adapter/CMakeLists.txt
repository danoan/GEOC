cmake_minimum_required(VERSION 3.2)
project(lib-geoc-adapter)

add_library(geoc-adapter SHARED src/temp.cpp)

target_include_directories(geoc-adapter PRIVATE
        include
        ../estimator/include
        ../../ext/include
        ${DGTAL_INCLUDE_DIRS})

target_link_libraries(geoc-adapter
        ${DGTAL_LIBRARIES}
        -lgeoc-estimator
        -ldipacus-derivates
        -ldipacus-components
        -ldipacus-base)


install(DIRECTORY include
        DESTINATION ${CMAKE_INSTALL_PREFIX}})

