cmake_minimum_required(VERSION 3.0)
project(geoc)


add_library(geoc INTERFACE)

target_include_directories(geoc INTERFACE
        include
        ext/include )

target_link_libraries(geoc INTERFACE
        ext/lib)