cmake_minimum_required(VERSION 3.2)
project(app-image-estimator)

add_executable(app-image-estimator image-estimator.cpp)

if(${USE_REMOTE_REPOSITORIES})
    add_dependencies(app-image-estimator dipacus gcurve)
endif()

find_package(OpenCV REQUIRED)

target_include_directories(app-image-estimator PRIVATE
        include
        ${CMAKE_SOURCE_DIR}/modules/api/include
        ${CMAKE_SOURCE_DIR}/modules/adapter/include
        ${CMAKE_SOURCE_DIR}/modules/estimator/include
        ${DIPACUS_INCLUDE_DIRS}
        ${BOOST_INCLUDE_DIRS}
        ${DGTAL_INCLUDE_DIRS})

target_link_libraries(app-image-estimator
        -ldipacus-derivates
        -ldipacus-components
        -ldipacus-base
        -lboost_filesystem
        -lboost_system
        ${DGTAL_LIBRARIES}
        ${OpenCV_LIBS})

install(TARGETS app-image-estimator
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )