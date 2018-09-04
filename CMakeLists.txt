cmake_minimum_required(VERSION 3.2)
project(geoc)

#It expects variables
#CMAKE_PREFIX_PATH: find_[package] macro directories
#CMAKE_INSTALL_PREFIX: install location
#EXTRA_LINK_DIR: extra link directories
#BUILD_PATH: location of build files

link_directories(${EXTRA_LINK_DIR}
        ${BUILD_PATH}/modules/estimator
        ${BUILD_PATH}/modules/adapter
        ${BUILD_PATH}/modules/energy)

add_subdirectory(modules/estimator)
add_subdirectory(modules/adapter)