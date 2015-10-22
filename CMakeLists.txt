CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -std=c++11 -Wall")

FIND_PACKAGE(OpenCV COMPONENTS core highgui imgproc REQUIRED)

set(boost_components system filesystem program_options)
add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_FILESYSTEM_VERSION=3 -DBOOST_FILESYSTEM_NO_DEPRECATED)
find_package(Boost REQUIRED COMPONENTS ${boost_components})
include_directories(${Boost_INCLUDE_DIRS})

set(FlameSRC main.cpp)
add_executable(fflame_gen ${FlameSRC})
target_link_libraries(fflame_gen ${OpenCV_LIBS} ${Boost_LIBRARIES})
