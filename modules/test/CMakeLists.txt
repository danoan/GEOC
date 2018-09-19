cmake_minimum_required(VERSION 3.2)
project(test-geoc)

find_package(DGtal REQUIRED)

add_executable(testGeoc testGeoc.cpp src/testConvergence.cpp)

target_include_directories(testGeoc PRIVATE
        include
        ../adapter/include
        ../estimator/include
        ../../ext/include
        ${DGTAL_INCLUDE_DIRS})

target_link_libraries(testGeoc
        -ldipacus-derivates
        -ldipacus-components
        -ldipacus-base
        ${DGTAL_LIBRARIES})