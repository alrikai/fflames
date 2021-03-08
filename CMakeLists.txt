CMAKE_MINIMUM_REQUIRED(VERSION 2.13)
project(FractalFlames)

include (CTest)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -std=c++17 -Wall -pthread -Wextra")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
include_directories(${CMAKE_SOURCE_DIR}/.)
include_directories(${CMAKE_SOURCE_DIR}/fflames)

FIND_PACKAGE(OpenCV COMPONENTS core imgproc imgcodecs REQUIRED)
MESSAGE("opencv version: ${OpenCV_VERSION}")


# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

set(boost_components system filesystem program_options)
add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_FILESYSTEM_VERSION=3 -DBOOST_FILESYSTEM_NO_DEPRECATED)
find_package(Boost REQUIRED COMPONENTS ${boost_components})
include_directories(${Boost_INCLUDE_DIRS})

set(FlameSRC main.cpp)
add_executable(fflame_gen ${FlameSRC})
target_link_libraries(fflame_gen ${OpenCV_LIBS} ${Boost_LIBRARIES} gtest_main)

add_subdirectory(tests)
