cmake_minimum_required(VERSION 3.5)  

project(Ueb2 CXX)  

message(matlib.h and matlib.cpp are taken from work for my bachelors thesis: https://github.com/kenscl/Bachelorarbeit-Code)

set(SOURCES
    main.cpp
    helper/mathhelper.cpp
    helper/mathhelper.h
    tle/tle.cpp
    tle/tle.h
    tle/tlereader.cpp
    tle/tlereader.h
    sgp4/coordinates.h
    sgp4/SGP4Propagator.cpp
    sgp4/SGP4Propagator.h
    sgp4/timeDate.h
    math/matlib.h
    math/matlib.cpp
)


add_executable(${PROJECT_NAME} ${SOURCES})
