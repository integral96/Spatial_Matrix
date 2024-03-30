#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <chrono>
#include <thread>

#include <boost/timer/timer.hpp>
#include <boost/system/system_error.hpp>

#include "include/Matrix4D.hpp"

static constexpr int NI = 9;
static constexpr int NJ = 9;
static constexpr int NK = 9;
static constexpr int NL = 9;

static constexpr int NI_b = 18;
static constexpr int NJ_b = 18;
static constexpr int NK_b = 18;


static constexpr std::array<size_t, 2> shape2D = { {NI, NJ} };
static constexpr std::array<size_t, 3> shape3D = { {NI, NJ, NK} };
static constexpr std::array<size_t, 3> shape3D_b = { {NI_b, NJ_b, NK_b} };
static constexpr std::array<size_t, 4> shape4D = { {NI, NJ, NK, NL} };


using C_T = int;

int main()
{
    using namespace _spatial;
    using namespace std::complex_literals;
    std::cout << std::fixed << std::setprecision(1);

    const auto int_perm = _my::permutation_inv<0>(1, 2, 3, 4);
    std::cout << int_perm << std::endl;
    const auto chr_perm = _my::permutation_inv<0>('a', 'b', 'c', 'd', 'e');
    std::cout << chr_perm << std::endl;
    try {



        Matrix2D<double> A2(shape2D);
        Matrix2D<std::complex<double>> AC2(shape2D);
        A2.Random(1.1, 2.2);
        AC2.Random(0.01, 0.72);
        Matrix3D<std::complex<double>> AC3(shape3D);
        Matrix3D<std::complex<double>> AC3_(shape3D);
        Matrix3D<std::complex<double>> AC3__(shape3D);
        AC3.Random(.1, .3);
        AC3_.Random(0.2, 0.8);
        AC3__.Random(0.1, 0.6);
        Matrix4D<std::complex<double>> AC4(shape4D);
        AC4.Random(1.1, 2.2);

        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_2D ..." << std::endl;

            std::cout << "DETERMINAT 2D = " << A2.DET() << std::endl;
            std::cout << "DETERMINAT 2D complex = " << AC2.DET() << std::endl;
            std::cout << "END TEST MATRIX_2D " << tmr.format() << std::endl;
        }
        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_3D. Dim:" << NI << "x" << NJ << "x" << NK << " ..." << std::endl;

            std::cout << "RULE DISTRIBUTED 3D + 2D (complex, index j) = " << std::endl;
            auto tmp2D = AC3.template product<'j'>(AC2) + AC3_.template product<'j'>(AC2);
            auto tmp2D_ = (AC3 + AC3_).template product<'j'>(AC2);
            std::cout << "RETURNT 3D (complex) = " << std::endl;
            std::cout << std::boolalpha << (tmp2D == tmp2D_)  << std::endl;
            std::cout << "RULE DISTRIBUTED 3D + 3D (complex, index j) = " << std::endl;
            auto tmp3D = AC3.template product<'j'>(AC3__) + AC3_.template product<'j'>(AC3__);
            auto tmp3D_ = (AC3 + AC3_).template product<'j'>(AC3__);
            std::cout << "RETURNT 4D (complex) = " << std::endl;
            std::cout << std::boolalpha << (tmp3D  == tmp3D_) << std::endl;
            for(size_t i = 0; i < tmp2D_.size(0); ++i) {
                std::cout << "MATRIX_3D cross-sections of orientation i = \n" << tmp2D_.transversal_matrix<'i'>(i) << std::endl;
            }
            for(size_t i = 0; i < tmp2D_.size(0) ; ++i) {
                std::cout << "Determinant of matrix 3D orientation i = " << i << "; " << tmp2D_.DET_orient<'i'>(i) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t j = 0; j < tmp2D_.size(1); ++j) {
                std::cout << "Determinant of matrix 3D orientation j = " << j << "; " << tmp2D_.DET_orient<'j'>(j) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t k = 0; k < tmp2D_.size(2); ++k) {
                std::cout << "Determinant of matrix 3D orientation k = " << k << "; " << tmp2D_.DET_orient<'k'>(k) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            std::cout << "The complete determinant of the matrix 3D = " << tmp2D_.DET_FULL() << std::endl;
            std::cout << "END TEST MATRIX_3D, ITERATION = " << ", " << tmr.format() << std::endl;

        }
        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_4D. Dim:" << NI << "x" << NJ << "x" << NK << "x" << NL << " ..." << std::endl;
            std::cout << "PRODUCT 4D * 3D (complex, index l) = " << std::endl;
            auto tmp3D_ = AC4.template product<'l'>(AC3);
            std::cout << tmp3D_ << std::endl;
            for(size_t i = 0; i < tmp3D_.size(0); ++i) {
                std::cout << "MATRIX_4D cross-sections of orientation i = \n" << tmp3D_.transversal_matrix<'i'>(i);
            }
            for(size_t i = 0; i < tmp3D_.size(0); ++i) {
                std::cout << "Determinant of matrix 4D orientation i = " << i << "; " << tmp3D_.DET_orient<'i'>(i) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t j = 0; j < tmp3D_.size(1); ++j) {
                std::cout << "Determinant of matrix 4D orientation j = " << j << "; " << tmp3D_.DET_orient<'j'>(j) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t k = 0; k < tmp3D_.size(2); ++k) {
                std::cout << "Determinant of matrix 4D orientation k = " << k << "; " << tmp3D_.DET_orient<'k'>(k) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t l = 0; l < tmp3D_.size(3); ++l) {
                std::cout << "Determinant of matrix 4D orientation l = " << l << "; " << tmp3D_.DET_orient<'l'>(l) << std::endl;
            }
            std::cout << "======================================" << std::endl;
//            std::cout << "The complete determinant of the matrix 4D = " << A4.DET_FULL() << std::endl;
            std::cout << "END TEST MATRIX_4D, ITERATION = " << ", " << tmr.format() << std::endl;
        }

    }
    catch (boost::system::system_error& ec) {
        std::cerr << "Error ocured! # " << ec.code() << " Messaga "
            << ec.what() << std::endl;
    }


    return 0;
}
