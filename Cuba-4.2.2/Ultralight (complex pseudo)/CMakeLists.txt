CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

PROJECT(ZprimeFIFO)

include_directories(/home/parracho/Desktop/Project/Cuba-4.2.2 include)
link_directories(/home/parracho/Desktop/Project/Cuba-4.2.2)

find_package(GSL 2.1 REQUIRED)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/include/boost")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/lib/x86_64-linux-gnu")

find_package(Boost 1.40 REQUIRED)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    add_definitions( "-DHAS_BOOST" )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=gnu++11")

add_subdirectory(src)

