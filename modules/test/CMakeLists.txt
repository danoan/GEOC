cmake_minimum_required(VERSION 3.2)
project(test-geoc-estimator)

find_package(DGtal REQUIRED)

add_executable(testEstimator testEstimator.cpp)

target_include_directories(testEstimator PRIVATE
        ../adapter/include
        ../estimator/include
        ../../ext/include
        ${DGTAL_INCLUDE_DIRS})

target_link_libraries(testEstimator
        -ldipacus-derivates
        -ldipacus-components
        -ldipacus-base
        ${DGTAL_LIBRARIES})