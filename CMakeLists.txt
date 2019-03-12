CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2)

include (CTest)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -std=c++17 -Wall -pthread -Wextra")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
include_directories(${CMAKE_SOURCE_DIR}/.)
include_directories(${CMAKE_SOURCE_DIR}/fflames)

FIND_PACKAGE(OpenCV COMPONENTS core highgui imgproc imgcodecs REQUIRED)

set(GTEST_DIR "${CMAKE_SOURCE_DIR}/gtest")
add_subdirectory(${GTEST_DIR}) #${googletest_BINARY_DIR})

set(boost_components system filesystem program_options)
add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_FILESYSTEM_VERSION=3 -DBOOST_FILESYSTEM_NO_DEPRECATED)
find_package(Boost REQUIRED COMPONENTS ${boost_components})
include_directories(${Boost_INCLUDE_DIRS})

set(FlameSRC main.cpp)
add_executable(fflame_gen ${FlameSRC})
target_link_libraries(fflame_gen ${OpenCV_LIBS} ${Boost_LIBRARIES} gtest_main)

set(TestFlameSRC tests/variant_tests.cpp)
add_executable(TestFlames ${TestFlameSRC})
target_link_libraries(TestFlames gtest gtest_main)
add_test(NAME test_fflames COMMAND TestFlames)
