#pragma once
#define TBB_PREVIEW_BLOCKED_RANGE_ND 1

#include <type_traits>

#include <boost/multi_array.hpp>
#include <boost/mpl/int.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/random.hpp>
#include <boost/noncopyable.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/core/demangle.hpp>
#include <boost/core/typeinfo.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_rangeNd.h>

namespace mpl = boost::mpl;
namespace proto = boost::proto;
namespace big = boost::multiprecision;

namespace _spatial {

    /////////MATRIX1D
    template<typename T>
    struct matrix1_ : boost::multi_array<T, 1> {
        typedef boost::multi_array<T, 1> array_type;
        typedef boost::multi_array_types::index_range range;
        typedef typename array_type::template array_view<1>::type sub_matrix_view1D;
        typedef T value_type;
    public:
        constexpr matrix1_(const std::array<size_t, 1>& shape) : array_type(shape) {}
        constexpr matrix1_(size_t N) : array_type({ boost::extents[N] }) {}

        constexpr void is_Matrix() {}

        inline T& operator () (size_t i) {
            BOOST_ASSERT_MSG((i < this->shape()[0]), "Error i >= size");
            return (*this)[i];
        }
        inline T const& operator () (size_t i) const {
            BOOST_ASSERT_MSG((i < this->shape()[0]), "Error i >= size");
            return (*this)[i];
        }
    };

    /////////MATRIX2D
    template<typename T>
    struct matrix2_ : boost::multi_array<T, 2> {
        typedef boost::multi_array<T, 2> array_type;
        typedef boost::multi_array_types::index_range range;
        typedef typename array_type::template array_view<1>::type sub_matrix_view1D;
        typedef typename array_type::template array_view<2>::type sub_matrix_view2D;
        typedef T value_type;
    public:
        constexpr matrix2_(const std::array<size_t, 2>& shape) : array_type(shape) {}
        constexpr matrix2_(size_t N, size_t M) : array_type({ boost::extents[N][M] }) {}

        constexpr void is_Matrix() {}

        inline T& operator () (size_t i, size_t j) {
            BOOST_ASSERT_MSG((i < this->shape()[0]), "Error i >= size N");
            BOOST_ASSERT_MSG((j < this->shape()[1]), "Error j >= size M");
            return (*this)[i][j];
        }
        inline T const& operator () (size_t i, size_t j) const {
            BOOST_ASSERT_MSG((i < this->shape()[0]), "Error i >= size N");
            BOOST_ASSERT_MSG((j < this->shape()[1]), "Error j >= size M");
            return (*this)[i][j];
        }
    };

    /////////MATRIX3D
    template<typename T>
    struct matrix3_ : boost::multi_array<T, 3> {
        using range_tbb = tbb::blocked_rangeNd<size_t, 3>;
        typedef boost::multi_array<T, 3> array_type;
        typedef boost::multi_array_types::index_range range;
        typedef typename array_type::template array_view<1>::type sub_matrix_view1D;
        typedef typename array_type::template array_view<2>::type sub_matrix_view2D;
        typedef typename array_type::template array_view<3>::type sub_matrix_view3D;
        typedef T value_type;
        array_type MTRX;
    public:
        constexpr matrix3_(const std::array<size_t, 3>& shape) : array_type(shape) {}
        constexpr matrix3_(size_t N, size_t M, size_t K) : array_type({ boost::extents[N][M][K] }) {}

        constexpr void is_Matrix() {}

        inline T& operator () (size_t i, size_t j, size_t k) {
            BOOST_ASSERT_MSG((i < this->shape()[0]), "Error i >= size N");
            BOOST_ASSERT_MSG((j < this->shape()[1]), "Error j >= size M");
            BOOST_ASSERT_MSG((k < this->shape()[2]), "Error k >= size K");
            return (*this)[i][j][k];
        }
        inline T const& operator () (size_t i, size_t j, size_t k) const {
            BOOST_ASSERT_MSG((i < this->shape()[0]), "Error i >= size N");
            BOOST_ASSERT_MSG((j < this->shape()[1]), "Error j >= size M");
            BOOST_ASSERT_MSG((k < this->shape()[2]), "Error k >= size K");
            return (*this)[i][j][k];
        }
    };

    /////////MATRIX4D
    template<typename T>
    struct matrix4_ : boost::multi_array<T, 4> {
        using range_tbb = tbb::blocked_rangeNd<size_t, 4>;
        typedef boost::multi_array<T, 4> array_type;
        typedef boost::multi_array_types::index_range range;
        typedef typename array_type::template array_view<1>::type sub_matrix_view1D;
        typedef typename array_type::template array_view<2>::type sub_matrix_view2D;
        typedef typename array_type::template array_view<3>::type sub_matrix_view3D;
        typedef typename array_type::template array_view<4>::type sub_matrix_view4D;
        typedef T value_type;
        array_type MTRX;
    public:
        constexpr matrix4_(const std::array<size_t, 4>& shape) : array_type(shape), MTRX(shape) {}
        constexpr matrix4_(size_t N, size_t M, size_t K, size_t L) : array_type({ boost::extents[N][M][K][L] }) {}

        constexpr void is_Matrix() {}

        inline T& operator () (size_t i, size_t j, size_t k, size_t l) {
            BOOST_ASSERT_MSG((i < this->shape()[0]), "Error i >= size N");
            BOOST_ASSERT_MSG((j < this->shape()[1]), "Error j >= size M");
            BOOST_ASSERT_MSG((k < this->shape()[2]), "Error k >= size K");
            BOOST_ASSERT_MSG((k < this->shape()[3]), "Error k >= size L");
            return (*this)[i][j][k][l];
        }
        inline T const& operator () (size_t i, size_t j, size_t k, size_t l) const {
            BOOST_ASSERT_MSG((i < this->shape()[0]), "Error i >= size N");
            BOOST_ASSERT_MSG((j < this->shape()[1]), "Error j >= size M");
            BOOST_ASSERT_MSG((k < this->shape()[2]), "Error k >= size K");
            BOOST_ASSERT_MSG((k < this->shape()[3]), "Error k >= size L");
            return (*this)[i][j][k][l];
        }
    };

    template<typename T>
    struct Matrix1D;
    template<typename T>
    struct Matrix2D;
    template<typename T>
    struct Matrix3D;
    template<typename T>
    struct Matrix4D;

}

template<typename Matrix, class = void>
struct IsMatrix : mpl::false_ {};
template<typename Matrix>
struct IsMatrix<Matrix, std::void_t<decltype (std::declval<Matrix&>().is_Matrix())> > : mpl::true_ {};

template<typename T, class Matrix>
struct IsMatrix1D : mpl::false_ {};
template<typename T>
struct IsMatrix1D<T, _spatial::Matrix1D<T>> : mpl::true_ {};

template<typename T, class Matrix>
struct IsMatrix2D : mpl::false_ {};
template<typename T>
struct IsMatrix2D<T, _spatial::Matrix2D<T>> : mpl::true_ {};

template<typename T, class Matrix>
struct IsMatrix3D : mpl::false_ {};
template<typename T>
struct IsMatrix3D<T, _spatial::Matrix3D<T>> : mpl::true_ {};

template<typename T, class Matrix>
struct IsMatrix4D : mpl::false_ {};
template<typename T>
struct IsMatrix4D<T, _spatial::Matrix4D<T>> : mpl::true_ {};

template<typename T>
struct IsBigInt : mpl::false_ {};
template<>
struct IsBigInt<big::int128_t> : mpl::true_ {};

template<typename T>
using type_end = typename mpl::if_c<std::is_same_v<T, size_t>, size_t, big::int128_t>::type;
template<typename T>
using type_two = typename mpl::if_c<std::is_same_v<T, int>, int, type_end<T> >::type;

///ABS
namespace _my {
template<typename T>
struct is_type : mpl::if_c<std::is_same_v<T, double>, double,
                 typename mpl::if_c<std::is_same_v<T, float>, float, type_two<T>>::type
                > {};
template<typename T, typename = std::enable_if_t<std::is_same<T, std::complex<double>>::value>>
inline double abs(const T& x) {
    return std::norm(x);
}
template<typename T, typename = std::enable_if_t<!std::is_same<T, std::complex<double>>::value>>
inline T abs(const T& x) {
    return (x < 0) ? -x : x;
}
template<typename T, template<typename Elem, typename = std::allocator<Elem>> class Array = std::vector>
void swap(Array<T>& a, int i, int j) {
    T s = a[i];
    a[i] = a[j];
    a[j] = s;
}

struct print_type
{
    template <class T>
    void operator() (T) const
    {
        auto const& ti = BOOST_CORE_TYPEID(T);
        std::cout << boost::core::demangled_name(ti) << std::endl;
    }
};
template<size_t DimN>
size_t invers(int I, int J, int K) {
    std::array<int, DimN> a{I, J, K};
    size_t count{};
    for(size_t i = 1; i < DimN - 1; ++i) {
        if(((a[i] > a[i+1]) && (a[i]>a[i-1])) || ((a[i] < a[i+1]) && (a[i] < a[i-1]))) {
            count++;
        }
    }
    if ((a[0] > a[1]) || (a[0] < a[1])) {
        count++;
    }
    return count;
}
template<size_t DimN>
size_t invers(int I, int J, int K, int L) {
    std::array<int, DimN> a{I, J, K, L};
    size_t count{};
    for(size_t i = 1; i < DimN - 1; ++i) {
        if(((a[i] > a[i+1]) && (a[i]>a[i-1])) || ((a[i] < a[i+1]) && (a[i] < a[i-1]))) {
            count++;
        }
    }
    if ((a[0] > a[1]) || (a[0] < a[1])) {
        count++;
    }
    return count;
}
}
