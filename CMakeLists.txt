cmake_minimum_required(VERSION 3.5)  

project(Ueb2 CXX)  

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

)


add_executable(${PROJECT_NAME} ${SOURCES})
