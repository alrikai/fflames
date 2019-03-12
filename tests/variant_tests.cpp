#include "gtest/gtest.h"

#include "../fflames/ifs_types.hpp"
#include "./fflames/ifs.hpp"

#include <memory>
#include <vector>
#include <string>

namespace {

TEST (VariantTest, Linear) {
    using data_t = double;
	std::vector<std::string> test_variant {"linear"};
    auto flamer = std::make_unique<affine_fcns::invoker<data_t>>(std::move(test_variant));
	flamer->seed_rng(0, 0);
    flamer->randomize_parameters(1, 1);
    flame_point<data_t> flame_pt(0.0, 0.0);
    flamer->invoke(0, flame_pt);

    EXPECT_EQ(0, flame_pt.x);
    EXPECT_EQ(0, flame_pt.y);
}

}

