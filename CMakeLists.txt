cmake_minimum_required(VERSION 3.5)  

project(Ueb2 CXX)  

set(SOURCES
    main.cpp
)
file(GLOB_RECURSE SOURCES 
    "*.cpp"
    "*.h"
)


add_executable(${PROJECT_NAME} ${SOURCES})

