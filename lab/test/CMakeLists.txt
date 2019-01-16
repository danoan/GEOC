cmake_minimum_required(VERSION 3.2)
project(test-geoc)

find_package(DGtal REQUIRED)

add_executable(testGeoc testGeoc.cpp src/testConvergence.cpp)

target_include_directories(testGeoc PRIVATE
        include
        ${CMAKE_SOURCE_DIR}/modules/adapter/include
        ${CMAKE_SOURCE_DIR}/modules/estimator/include
        ${CMAKE_SOURCE_DIR}/modules/api/include
        ${DIPACUS_INCLUDE_DIRS}
        ${GCURVE_INCLUDE_DIRS}
        ${DGTAL_INCLUDE_DIRS})

target_link_libraries(testGeoc
        -ldipacus-derivates
        -ldipacus-components
        -ldipacus-base
        ${DGTAL_LIBRARIES})