#include "gtest/gtest.h"

#include "../fflames/ifs_types.hpp"
#include "./fflames/ifs.hpp"

#include <memory>
#include <vector>
#include <string>

namespace {

//breakpoint on failure:
//gdb --args ./TestFlames --gtest_break_on_failure	
TEST (VariantTest, LinearDefault) {
    using data_t = double;
	std::vector<std::string> test_variant {"linear"};
    auto flamer = std::make_unique<FFlames::affine_fcns::invoker<data_t>>(std::move(test_variant));
	flamer->seed_rng(0, 0);
    flamer->randomize_parameters(0, 0);
    FFlames::flame_point<data_t> flame_pt(0.0, 0.0);
    flamer->invoke(0, flame_pt);

    EXPECT_EQ(0, flame_pt.x);
    EXPECT_EQ(0, flame_pt.y);
}

TEST (VariantTest, LinearParam) {
    using data_t = double;
	std::vector<std::string> test_variant {"linear"};
    auto flamer = std::make_unique<FFlames::affine_fcns::invoker<data_t>>(std::move(test_variant));
	flamer->seed_rng(0, 0);
    flamer->randomize_parameters(1, 1);
    FFlames::flame_point<data_t> flame_pt(1.0, 1.0);
    flamer->invoke(0, flame_pt);

    EXPECT_EQ(9, flame_pt.x);
    EXPECT_EQ(15, flame_pt.y);
}

TEST (VariantTest, SinusoidDefault) {
    using data_t = double;
	std::vector<std::string> test_variant {"sinusoidal"};
    auto flamer = std::make_unique<FFlames::affine_fcns::invoker<data_t>>(std::move(test_variant));
	flamer->seed_rng(0, 0);
    flamer->randomize_parameters(0, 0);
    FFlames::flame_point<data_t> flame_pt(0.0, 0.0);
    flamer->invoke(0, flame_pt);

    EXPECT_EQ(0, flame_pt.x);
    EXPECT_EQ(0, flame_pt.y);
}

TEST (VariantTest, SinusoidParam) {
    using data_t = double;
	std::vector<std::string> test_variant {"sinusoidal"};
    auto flamer = std::make_unique<FFlames::affine_fcns::invoker<data_t>>(std::move(test_variant));
	flamer->seed_rng(0, 0);
    flamer->randomize_parameters(1, 1);
    FFlames::flame_point<data_t> flame_pt(1.0, 1.0);
    flamer->invoke(0, flame_pt);

    ASSERT_DOUBLE_EQ(0.18219573339672879, flame_pt.x);
    ASSERT_DOUBLE_EQ(0.22327145873359033, flame_pt.y);
}

TEST (VariantTest, SphericalDefault) {
    using data_t = double;
	std::vector<std::string> test_variant {"spherical"};
    auto flamer = std::make_unique<FFlames::affine_fcns::invoker<data_t>>(std::move(test_variant));
	flamer->seed_rng(0, 0);
    flamer->randomize_parameters(0, 0);
    FFlames::flame_point<data_t> flame_pt(0.0, 0.0);
    flamer->invoke(0, flame_pt);

    EXPECT_EQ(0, flame_pt.x);
    EXPECT_EQ(0, flame_pt.y);
}

TEST (VariantTest, SphericalParam) {
    using data_t = double;
	std::vector<std::string> test_variant {"spherical"};
    auto flamer = std::make_unique<FFlames::affine_fcns::invoker<data_t>>(std::move(test_variant));
	flamer->seed_rng(0, 0);
    flamer->randomize_parameters(1, 1);
    FFlames::flame_point<data_t> flame_pt(1.0, 1.0);
    flamer->invoke(0, flame_pt);

    EXPECT_EQ(0, flame_pt.x);
    EXPECT_EQ(0, flame_pt.y);
}



}

