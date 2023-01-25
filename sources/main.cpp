#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <thread>

#include <boost/timer/timer.hpp>
#include <boost/system/system_error.hpp>


#include "../include/Matrix_type.hpp"

static constexpr int NI = 9;
static constexpr int NJ = 9;
static constexpr int NK = 9;
static constexpr int NL = 9;

static constexpr std::array<size_t, 2> shape2D = { {NI, NJ} };
static constexpr std::array<size_t, 3> shape3D = { {NI, NJ, NK} };
static constexpr std::array<size_t, 4> shape4D = { {NI, NJ, NK, NL} };


int main()
{
    using namespace std::complex_literals;
    std::cout << std::fixed << std::setprecision(1);
    try {

        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_2D ..." << std::endl;
            MATRIX<2, double> A2(shape2D);
            MATRIX<2, std::complex<double>> AC(shape2D);
            std::cout << "CTOR MATRIX_2D ..." << std::endl;
            A2.Random(1.1, 2.2);
            AC.Random(1.1, 2.2);
            std::cout << "RANDOM MATRIX_2D ..." << std::endl;
//            std::cout << A2 << std::endl;
//            std::cout << A2*A2 << std::endl;
//            std::cout << AC*AC << std::endl;
            std::cout << "DETERMINAT 2D = " << A2.DET() << std::endl;
            std::cout << "DETERMINAT 2D complex = " << AC.DET() << std::endl;
            std::cout << "END TEST MATRIX_2D " << tmr.format() << std::endl;
        }
        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_3D. Dim:" << NI << "x" << NJ << "x" << NK << " ..." << std::endl;
            MATRIX<3, double> A3(shape3D);
            MATRIX<3, std::complex<double>> AC(shape3D);
            std::cout << "CTOR MATRIX_3D ..." << std::endl;
            A3.Random(1.1, 2.2);
            AC.Random(1.1, 2.2);
//            std::cout << A3*A3 << std::endl;
//            std::cout << AC*AC << std::endl;
            for(size_t i = 0; i < AC.size(0); ++i) {
                std::cout << "MATRIX_3D cross-sections of orientation i = \n" << AC.transversal_matrix('i', i) << std::endl;
            }
            for(size_t i = 0; i < AC.size(0) ; ++i) {
                std::cout << "Determinant of matrix 3D orientation i = " << i << "; " << AC.DET_orient('i', i) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t j = 0; j < AC.size(1); ++j) {
                std::cout << "Determinant of matrix 3D orientation j = " << j << "; " << AC.DET_orient('j', j) << std::endl;
            }
            
            std::cout << "======================================" << std::endl;
            for(size_t k = 0; k < AC.size(2); ++k) {
                std::cout << "Determinant of matrix 3D orientation k = " << k << "; " << AC.DET_orient('k', k) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            std::cout << "The complete determinant of the matrix 4D = " << AC.DET_FULL() << std::endl;
            std::cout << "END TEST MATRIX_3D, ITERATION = " << ", " << tmr.format() << std::endl;

        }
        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_4D. Dim:" << NI << "x" << NJ << "x" << NK << "x" << NL << " ..." << std::endl;
            MATRIX<4, double> A4(shape4D);
            MATRIX<4, std::complex<double>> AC(shape4D);
            std::cout << "CTOR MATRIX_4D ..." << std::endl;
            A4.Random(1.1, 2.2);
            AC.Random(1.1, 2.2);
            std::cout << AC*AC << std::endl;
            for(size_t i = 0; i < A4.size(0); ++i) {
                std::cout << "MATRIX_4D cross-sections of orientation i = \n" << AC.transversal_matrix('i', i);
            }
            for(size_t i = 0; i < A4.size(0); ++i) {
                std::cout << "Determinant of matrix 4D orientation i = " << i << "; " << AC.DET_orient('i', i) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t j = 0; j < A4.size(1); ++j) {
                std::cout << "Determinant of matrix 4D orientation j = " << j << "; " << AC.DET_orient('j', j) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t k = 0; k < A4.size(2); ++k) {
                std::cout << "Determinant of matrix 4D orientation k = " << k << "; " << AC.DET_orient('k', k) << std::endl;
            }
            std::cout << "======================================" << std::endl;
            for(size_t l = 0; l < A4.size(3); ++l) {
                std::cout << "Determinant of matrix 4D orientation l = " << l << "; " << AC.DET_orient('l', l) << std::endl;
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
