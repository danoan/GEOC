cmake_minimum_required(VERSION 3.2)
project(lib-geoc-estimators)

add_library(geoc-estimator SHARED src/temp.cpp)

target_include_directories(geoc-estimator PRIVATE
        include
        ../../ext/include
        ${DGTAL_INCLUDE_DIRS})

target_link_libraries(geoc-estimator
        ${DGTAL_LIBRARIES}
        -lgcurve)

install(DIRECTORY include
        DESTINATION ${CMAKE_INSTALL_PREFIX})