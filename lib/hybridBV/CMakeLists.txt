cmake_minimum_required(VERSION 3.5)
set(CMAKE_OSX_ARCHITECTURES x86_64)
project(HybridBV)
set(CMAKE_PROJECT_NAME "HybridBV")
#set(CMAKE_COMPILER_IS_GNUCXX)

include_directories(include)

set(SOURCES
        src/hybridBV/basics.c
        include/hybridBV/basics.h
        src/hybridBV/hybridBV.c
        include/hybridBV/hybridBV.h
        src/hybridBV/hybridId.c
        include/hybridBV/hybridId.h
        src/hybridBV/leafBV.c
        include/hybridBV/leafBV.h
        src/hybridBV/leafId.c
        include/hybridBV/leafId.h
        src/hybridBV/staticBV.c
        include/hybridBV/staticBV.h
        include/hybridBV/hybridBVId.h
        include/hybridBV/leafBVId.h
        include/hybridBV/staticBVId.h
        src/hybridBV/hybridBVId.c
        src/hybridBV/leafBVId.c
        src/hybridBV/staticBVId.c

)

# Flags
set(CMAKE_CXX_FLAGS " -Wall")
set(CMAKE_CXX_FLAGS_RELEASE " -std=c++0x -Wall -O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -O9 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -g3  -DDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -O3 -pg -DDEBUG")

# set default compilation and show details for debug build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif(NOT CMAKE_BUILD_TYPE)

# Create the library
add_library(hybridbv_gn STATIC ${SOURCES})
