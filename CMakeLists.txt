CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -std=c++1y -Wall")

FIND_PACKAGE(OpenCV COMPONENTS core highgui imgproc REQUIRED)
#FIND_PACKAGE(Boost 1.55.0 REQUIRED)

#set(FFSRC fractal_flame.cpp)
#add_library(fflame SHARED ${FFSRC})

set(FlameSRC main.cpp)
add_executable(fflame_gen ${FlameSRC})
target_link_libraries(fflame_gen ${OpenCV_LIBS})
