cmake_minimum_required(VERSION 3.2)
project(test-geoc)

find_package(DGtal REQUIRED)
find_package(OpenCV REQUIRED)

add_executable(testGeoc testGeoc.cpp src/testConvergence.cpp)

if(${USE_LOCAL_REPOSITORIES})
    add_dependencies(testGeoc dipacus gcurve geoc-api)
endif()

target_include_directories(testGeoc PUBLIC
        include
        ${CMAKE_SOURCE_DIR}/modules/adapter/include
        ${CMAKE_SOURCE_DIR}/modules/estimator/include
        ${CMAKE_SOURCE_DIR}/modules/api/include
        ${DIPACUS_INCLUDE_DIRS}
        ${GCURVE_INCLUDE_DIRS}
        ${OpenCV_INCLUDE_DIRS}
        ${DGTAL_INCLUDE_DIRS})

target_link_libraries(testGeoc
        -ldipacus-derivates
        ${DGTAL_LIBRARIES})
