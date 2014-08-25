CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -std=c++11 -Wall")

FIND_PACKAGE(OpenCV COMPONENTS core highgui imgproc REQUIRED)
#FIND_PACKAGE(Boost 1.55.0 REQUIRED)

set(FlameSRC main.cpp)
add_executable(fflame ${FlameSRC})
target_link_libraries(fflame ${OpenCV_LIBS})
