#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <thread>

#include <boost/timer/timer.hpp>
#include <boost/system/system_error.hpp>


#include "../include/generate.hpp"


static constexpr int NI = 9;
static constexpr int NJ = 9;
static constexpr int NK = 9;
static constexpr int NL = 9;

static constexpr int NI1 = 64;
static constexpr int NJ1 = 64;
static constexpr int NK1 = 64;

static constexpr std::array<size_t, 2> shape2D = { {NI, NJ} };
static constexpr std::array<size_t, 3> shape3D = { {NI, NJ, NK} };
static constexpr std::array<size_t, 4> shape4D = { {NI, NJ, NK, NL} };

static constexpr std::array<size_t, 2> shape2D_big = { {NI1, NJ1} };
static constexpr std::array<size_t, 3> shape3D_big = { {NI1, NJ1, NK1} };


int main()
{
    using namespace std::complex_literals;
    std::cout << std::fixed << std::setprecision(1);
    MATRIX<2, double> A2(shape2D);
    MATRIX<2, double> A2_big(shape2D_big);
    MATRIX<2, std::complex<double>> AC2(shape2D);
    A2.Random(1.1, 2.2);
    AC2.Random(0.01, 0.72);
    MATRIX<3, int> A_maze(shape3D);
    MATRIX<3, double> A3(shape3D);
    MATRIX<3, double> A3_(shape3D);
    MATRIX<3, int> A3_big(shape3D_big);
    MATRIX<3, std::complex<double>> AC3(shape3D);
    MATRIX<3, std::complex<double>> AC3_(shape3D);
    MATRIX<3, std::complex<double>> AC3__(shape3D);
    std::cout << "CTOR MATRIX_3D ..." << std::endl;
//    A3.Random(1.1, 2.2);
//    A3_.Random(1.1, 2.2);
    AC3.Random(0.1, 0.5);
    AC3_.Random(0.2, 0.8);
    AC3__.Random(0.1, 0.6);
    A2_big.Random(1.2, 2.2);
    MATRIX<4, double> A4(shape4D);
    MATRIX<4, std::complex<double>> AC4(shape4D);
    std::cout << "CTOR MATRIX_3D BIG ..." << std::endl;
    A4.Random(0.1, 1.2);
    AC4.Random(1.1, 2.2);
    try {
        maze_weight<std::complex<int>, _spatial::Matrix3D> maze_3D(shape3D);
        maze_weight<std::complex<int>, _spatial::Matrix4D> maze_4D(shape4D);
        maze_4D.calc_maze(.5, .5);
        std::cout << maze_4D.get_maze() << std::endl;


//        {
//            boost::timer::cpu_timer tmr;
//            std::cout << "START TEST MATRIX_2D ..." << std::endl;
//
//            std::cout << "DETERMINAT 2D = " << A2.DET() << std::endl;
//            std::cout << "DETERMINAT 2D complex = " << AC2.DET() << std::endl;
//            std::cout << "END TEST MATRIX_2D " << tmr.format() << std::endl;
//        }
//        {
//            boost::timer::cpu_timer tmr;
//            std::cout << "START TEST MATRIX_3D. Dim:" << NI << "x" << NJ << "x" << NK << " ..." << std::endl;
//
//            std::cout << "RULE DISTRIBUTED 3D + 2D (complex, index j) = " << std::endl;
//            auto tmp2D = AC3.template product<'j'>(AC2) + AC3_.template product<'j'>(AC2);
//            auto tmp2D_ = (AC3 + AC3_).template product<'j'>(AC2);
//            std::cout << "RETURNT 3D (complex) = " << std::endl;
//            std::cout << std::boolalpha << (tmp2D == tmp2D_)  << std::endl;
//            std::cout << "RULE DISTRIBUTED 3D + 3D (complex, index j) = " << std::endl;
//            auto tmp3D = AC3.template product<'j'>(AC3__) + AC3_.template product<'j'>(AC3__);
//            auto tmp3D_ = (AC3 + AC3_).template product<'j'>(AC3__);
//            std::cout << "RETURNT 4D (complex) = " << std::endl;
//            std::cout << std::boolalpha << (tmp3D  == tmp3D_) << std::endl;
//            for(size_t i = 0; i < tmp2D_.size(0); ++i) {
//                std::cout << "MATRIX_3D cross-sections of orientation i = \n" << tmp2D_.transversal_matrix<'i'>(i) << std::endl;
//            }
//            for(size_t i = 0; i < tmp2D_.size(0) ; ++i) {
//                std::cout << "Determinant of matrix 3D orientation i = " << i << "; " << tmp2D_.DET_orient<'i'>(i) << std::endl;
//            }
//            std::cout << "======================================" << std::endl;
//            for(size_t j = 0; j < tmp2D_.size(1); ++j) {
//                std::cout << "Determinant of matrix 3D orientation j = " << j << "; " << tmp2D_.DET_orient<'j'>(j) << std::endl;
//            }
//            std::cout << "======================================" << std::endl;
//            for(size_t k = 0; k < tmp2D_.size(2); ++k) {
//                std::cout << "Determinant of matrix 3D orientation k = " << k << "; " << tmp2D_.DET_orient<'k'>(k) << std::endl;
//            }
//            std::cout << "======================================" << std::endl;
//            std::cout << "The complete determinant of the matrix 3D = " << tmp2D_.DET_FULL() << std::endl;
//            std::cout << "END TEST MATRIX_3D, ITERATION = " << ", " << tmr.format() << std::endl;
//
//        }
//        {
//            boost::timer::cpu_timer tmr;
//            std::cout << "START TEST MATRIX_4D. Dim:" << NI << "x" << NJ << "x" << NK << "x" << NL << " ..." << std::endl;
//            std::cout << "PRODUCT 4D * 3D (complex, index l) = " << std::endl;
//            auto tmp3D_ = AC4.template product<'l'>(AC3);
//            std::cout << tmp3D_ << std::endl;
//            for(size_t i = 0; i < tmp3D_.size(0); ++i) {
//                std::cout << "MATRIX_4D cross-sections of orientation i = \n" << tmp3D_.transversal_matrix<'i'>(i);
//            }
//            for(size_t i = 0; i < tmp3D_.size(0); ++i) {
//                std::cout << "Determinant of matrix 4D orientation i = " << i << "; " << tmp3D_.DET_orient<'i'>(i) << std::endl;
//            }
//            std::cout << "======================================" << std::endl;
//            for(size_t j = 0; j < tmp3D_.size(1); ++j) {
//                std::cout << "Determinant of matrix 4D orientation j = " << j << "; " << tmp3D_.DET_orient<'j'>(j) << std::endl;
//            }
//            std::cout << "======================================" << std::endl;
//            for(size_t k = 0; k < tmp3D_.size(2); ++k) {
//                std::cout << "Determinant of matrix 4D orientation k = " << k << "; " << tmp3D_.DET_orient<'k'>(k) << std::endl;
//            }
//            std::cout << "======================================" << std::endl;
//            for(size_t l = 0; l < tmp3D_.size(3); ++l) {
//                std::cout << "Determinant of matrix 4D orientation l = " << l << "; " << tmp3D_.DET_orient<'l'>(l) << std::endl;
//            }
//            std::cout << "======================================" << std::endl;
////            std::cout << "The complete determinant of the matrix 4D = " << A4.DET_FULL() << std::endl;
//            std::cout << "END TEST MATRIX_4D, ITERATION = " << ", " << tmr.format() << std::endl;
//        }

    }
    catch (boost::system::system_error& ec) {
        std::cerr << "Error ocured! # " << ec.code() << " Messaga "
            << ec.what() << std::endl;
    }


    return 0;
}
