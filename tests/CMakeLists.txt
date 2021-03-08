set(TestFlameSRC variant_tests.cpp)
add_executable(TestFlames ${TestFlameSRC})
target_link_libraries(TestFlames gtest gtest_main)
add_test(NAME test_fflames COMMAND TestFlames)
