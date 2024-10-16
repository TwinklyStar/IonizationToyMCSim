cmake_minimum_required(VERSION 3.26)
project(IonizationToyMCSim)

set(CMAKE_CXX_STANDARD 17)

# boost library
set(BOOST_ROOT /opt/homebrew/Cellar/boost/1.86.0)
set(BOOST_INCLUDEDIR /opt/homebrew/Cellar/boost/1.86.0/include)
set(BOOST_LIBRARY_DIR /opt/homebrew/Cellar/boost/1.86.0/lib)

find_package(Boost 1.86.0 COMPONENTS system filesystem REQUIRED)

include_directories(${Boost_INCLUDE_DIR})

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()

# ROOT library
FIND_PACKAGE(ROOT COMPONENTS Geom TMVA Eve Gui RGL EG)
IF (ROOT_FOUND)
    INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
    LINK_LIBRARIES(${ROOT_LIBRARIES})
    LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})
ENDIF ()

include_directories(include)

file(GLOB SOURCES src/*.cpp main.cpp)
file(GLOB HEADERS include/*.h)

add_executable(IonizationToyMCSim ${SOURCES} ${HEADERS})
add_executable(TestHarmonic test/harmonic_oscillator.cpp)

