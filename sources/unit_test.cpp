#define BOOST_TEST_MODULE test_module_name

#include <boost/test/included/unit_test.hpp>

#include "include/Matrix4D.hpp"
#include "include/graph_weight.hpp"

static constexpr int NI_b = 18;
static constexpr int NJ_b = 18;
static constexpr int NK_b = 18;

static constexpr std::array<size_t, 3> shape3D_b = { {NI_b, NJ_b, NK_b} };

BOOST_AUTO_TEST_CASE(test_no1) {
    _spatial::Matrix3D<std::complex<double>> AC3(shape3D_b);
    _spatial::Matrix3D<std::complex<double>> AC3_(shape3D_b);
    AC3_.Random(0.2, 0.8);
    BOOST_CHECK_NE(AC3, AC3_);
}
