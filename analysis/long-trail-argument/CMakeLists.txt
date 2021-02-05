# Time-stamp: <2019-02-04 16:30:13 lperrin>

# LICENSE

cmake_minimum_required(VERSION 0.4)
PROJECT(Trail-Search)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -std=c++11 -g")

ADD_EXECUTABLE(Trail-Search 
trail-search.cpp
general_utils.cpp
linear_layer.cpp
truncated-trail.cpp)

