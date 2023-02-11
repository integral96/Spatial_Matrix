#pragma once

#include "../include/Matrix3D.hpp"


static constexpr uint8_t FORW  = 64;  //01 00 00 00;
static constexpr uint8_t BACKW = 128; //10 00 00 00;

namespace _spatial {

template<typename Expr>
struct Matrix4D_expr;

struct Matrix4DGrammar : proto::or_<
                              proto::terminal<matrix4_<proto::_>>,
//                              proto::plus<Matrix4DGrammar, Matrix4DGrammar>,
                              proto::minus<Matrix4DGrammar, Matrix4DGrammar>,
                              proto::negate< Matrix4DGrammar>,
                              proto::less_equal< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::greater_equal< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::less< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::greater< Matrix4DGrammar, Matrix4DGrammar>,
                              proto::not_equal_to< Matrix4DGrammar, Matrix4DGrammar>
//                              proto::equal_to< Matrix4DGrammar, Matrix4DGrammar>
                                    >{};
struct Matrix4D_domain : proto::domain<proto::generator<Matrix4D_expr>, Matrix4DGrammar> {};

struct Matrix4D_context : proto::callable_context< Matrix4D_context const > {
    Matrix4D_context(size_t i, size_t j, size_t k, size_t l) : i(i), j(j), k(k), l(l) {}

    template<typename Expr, typename Tag = typename Expr::proto_tag>
    struct eval : proto::default_eval<Expr, Matrix4D_context> {};

    template<typename Expr>
    struct eval<Expr, proto::tag::terminal> {
        using result_type = typename proto::result_of::value<Expr>::type::value_type;

        result_type operator()(const Expr& expr, Matrix4D_context& ctx) const {
            return proto::value(expr)(ctx.i, ctx.j, ctx.k, ctx.l);
        }
    };

public:
    size_t i;
    size_t j;
    size_t k;
    size_t l;
};

struct SizeMatrix4D_context {
    SizeMatrix4D_context(size_t Ni, size_t Nj, size_t Nk, size_t Nl)
      : NI(Ni), NJ(Nj), NK(Nk), NL(Nl) {}
    template<typename Expr, typename EnableIf = void>
    struct eval : proto::null_eval<Expr, SizeMatrix4D_context const> {};

    template<typename Expr>
    struct eval<Expr, typename boost::enable_if<
            proto::matches<Expr, proto::terminal<matrix4_<proto::_> > >
        >::type
    >
    {
        typedef void result_type;

        result_type operator ()(Expr &expr, SizeMatrix4D_context const &ctx) const
        {
            if(ctx.NI != proto::value(expr).shape()[0]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу i");
            } else if (ctx.NJ != proto::value(expr).shape()[1]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу j");
            } else if (ctx.NK != proto::value(expr).shape()[2]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу k");
            } else if (ctx.NL != proto::value(expr).shape()[3]) {
                throw std::runtime_error("Матрицы не совпадают в размерности по индексу l");
            }
        }
    };

    size_t NI;
    size_t NJ;
    size_t NK;
    size_t NL;
};

template<typename Expr>
struct Matrix4D_expr : proto::extends<Expr, Matrix4D_expr<Expr>, Matrix4D_domain> {
    Matrix4D_expr(const Expr& expr = Expr()) : Matrix4D_expr::proto_extends(expr) {}

    typename proto::result_of::eval<Expr, Matrix4D_context>::type
    operator () (size_t i, size_t j, size_t k, size_t l) const {
        Matrix4D_context ctx(i, j, k, l);
        return proto::eval(*this, ctx);
    }
};
template<typename T>
struct Matrix4D : Matrix4D_expr<typename proto::terminal< matrix4_<T>>::type> {
    using expr_type = typename proto::terminal< matrix4_<T>>::type;
    using range_tbb = tbb::blocked_rangeNd<size_t, 4>;

    using array_type = typename matrix4_<T>::array_type;
    using range = boost::multi_array_types::index_range;
    using sub_matrix_view1D = typename matrix4_<T>::sub_matrix_view1D;

    const std::array<size_t, 4>& shape_;
    static constexpr size_t dim = 4;
    constexpr Matrix4D(const std::array<size_t, 4>& shape) :
        Matrix4D_expr<expr_type>(expr_type::make(matrix4_<T>(shape))), shape_(shape) {

    }
    size_t size(size_t i) const {
        BOOST_ASSERT_MSG((i < 4), "Error i >= 4");
        return proto::value(*this).shape()[i];
    }
    void Random(T min, T max, long shift = 0) {
        std::time_t now = std::time(0);
        boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
            if constexpr(std::is_integral_v<T> || IsBigInt<T>::value) {
                boost::random::uniform_int_distribution<> dist{int(min), int(max)};
                for(size_t i = 0; i < size(0); ++i)
                    for(size_t j = 0; j < size(1); ++j)
                        for(size_t k = 0; k < size(2); ++k)
                            for(size_t l = 0; l < size(3); ++l)
                                proto::value(*this)(i, j, k, l) = dist(gen);
            }
            if constexpr(!std::is_integral_v<T> && std::is_same_v<T, double>) {
                boost::random::uniform_real_distribution<> dist{double(min), double(max)};
                for(size_t i = 0; i < size(0); ++i)
                    for(size_t j = 0; j < size(1); ++j)
                        for(size_t k = 0; k < size(2); ++k)
                            for(size_t l = 0; l < size(3); ++l)
                                proto::value(*this)(i, j, k, l) = dist(gen);
            }
    }
    template<typename F,  typename = std::enable_if_t<std::is_same_v<T, std::complex<typename _my::is_type<F>::type>> &&
                                                      std::is_same_v<F, typename _my::is_type<F>::type>>>
    void Random(F min, F max) {
        std::time_t now = std::time(0);
        boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
        boost::random::uniform_real_distribution<> dist{double(min), double(max)};
        for(size_t i = 0; i < size(0); ++i)
            for(size_t j = 0; j < size(1); ++j)
                for(size_t k = 0; k < size(2); ++k)
                    for(size_t l = 0; l < size(3); ++l)
                    proto::value(*this)(i, j, k, l) = std::complex(dist(gen), dist(gen));
    }
    template< typename Expr >
    Matrix4D<T> &operator = (Expr const & expr) {
        SizeMatrix4D_context const sizes(size(0), size(1), size(2), size(3));
        proto::eval(proto::as_expr<Matrix4D_domain>(expr), sizes);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
            const auto& out_i = out.dim(0);
            const auto& out_j = out.dim(1);
            const auto& out_k = out.dim(2);
            const auto& out_l = out.dim(3);
            for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                    for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                        for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                            proto::value(*this)(i, j, k, l) = expr(i, j, k, l);
        });
        return *this;
    }
    template< typename Expr >
    Matrix4D<T> &operator += (Expr const & expr) {
        SizeMatrix4D_context const sizes(size(0), size(1), size(2), size(3));
        proto::eval(proto::as_expr<Matrix4D_domain>(expr), sizes);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
            const auto& out_i = out.dim(0);
            const auto& out_j = out.dim(1);
            const auto& out_k = out.dim(2);
            const auto& out_l = out.dim(3);
            for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                    for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                        for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                            proto::value(*this)(i, j, k, l) += expr(i, j, k, l);
        });
        return *this;
    }
    template< typename Expr >
    friend Matrix4D<T> operator + (Matrix4D<T> lhs, Expr const& expr) {
        SizeMatrix4D_context const sizes(lhs.size(0), lhs.size(1), lhs.size(2), lhs.size(2));
        proto::eval(proto::as_expr<Matrix4D_domain>(expr), sizes);
        lhs += expr;
        return lhs;
    }
    template< typename Expr >
    bool operator == (Expr const& expr) {
        SizeMatrix4D_context const sizes(size(0), size(1), size(2), size(3));
        proto::eval(proto::as_expr<Matrix4D_domain>(expr), sizes);
        size_t tmp_ = tbb::parallel_reduce(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(2) }), size_t(0),
                [=](const range_tbb& out, size_t tmp) {
                const auto& out_i = out.dim(0);
                const auto& out_j = out.dim(1);
                const auto& out_k = out.dim(2);
                const auto& out_l = out.dim(3);
                for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                            for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                            tmp += proto::value(*this)(i, j, k, l) != expr(i, j, k, l) ? 1 : 0;
                return tmp; }, std::plus<size_t>() );
        size_t tmp2_ = tbb::parallel_reduce(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(2) }), size_t(0),
                [=](const range_tbb& out, size_t tmp) {
                const auto& out_i = out.dim(0);
                const auto& out_j = out.dim(1);
                const auto& out_k = out.dim(2);
                const auto& out_l = out.dim(3);
                for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                            for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                            tmp += proto::value(*this)(i, j, k, l) == expr(i, j, k, l) ? 1 : 0;
                return tmp; }, std::plus<size_t>() );
        std::cout << "EQUAL = " << tmp_ << "; " << tmp2_ << std::endl;
        if(tmp_ == 0) return true;
        else return false;
    }
    template<char Index, typename F, typename =
             std::enable_if_t<std::is_same_v<T, F> && (Index == 'i' || Index == 'j' || Index == 'k' || Index == 'l')> >
    Matrix4D<T> product(Matrix3D<F> const& matr1) {
        if constexpr(Index == 'i') {
            BOOST_ASSERT_MSG(size(0) <= matr1.size(1),
                             "MATRICES 3D * 2D SIZE IS NOT EQUAL, Index == i");
            SizeMatrix3D_context const sizes(size(0), size(0), size(0));
            proto::eval(proto::as_expr<Matrix3D_domain>(matr1), sizes);
            Matrix4D<T> matrix(shape_);
            auto value_element([&](size_t i, size_t j, size_t k, size_t l) {
                for (size_t z = 0; z < size(0); ++z)
                    proto::value(matrix)(i, j, k, l) += proto::value(*this)(z, j, k, l) * proto::value(matr1)(z, i, l);
            });
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(3) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    const auto& out_l = out.dim(3);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                                    value_element(i, j, k, l);
                });

            return matrix;
        }
        if constexpr(Index == 'j') {
            BOOST_ASSERT_MSG(size(0) <= matr1.size(1),
                             "MATRICES 3D * 2D SIZE IS NOT EQUAL, Index == j");
            SizeMatrix3D_context const sizes(size(1), size(1), size(1));
            proto::eval(proto::as_expr<Matrix3D_domain>(matr1), sizes);
            Matrix4D<T> matrix(shape_);
            auto value_element([&](size_t i, size_t j, size_t k, size_t l) {
                for (size_t z = 0; z < size(1); ++z)
                    proto::value(matrix)(i, j, k, l) += proto::value(*this)(i, z, k, l) * proto::value(matr1)(z, j, l);
            });
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(3) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    const auto& out_l = out.dim(3);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                                    value_element(i, j, k, l);
                });

            return matrix;
        }
        if constexpr(Index == 'k') {
            BOOST_ASSERT_MSG(size(2) <= matr1.size(1),
                             "MATRICES 3D * 2D SIZE IS NOT EQUAL, Index == i");
            SizeMatrix3D_context const sizes(size(2), size(2), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(matr1), sizes);
            Matrix4D<T> matrix(shape_);
            auto value_element([&](size_t i, size_t j, size_t k, size_t l) {
                for (size_t z = 0; z < size(2); ++z)
                    proto::value(matrix)(i, j, k, l) += proto::value(*this)(i, j, z, l) * proto::value(matr1)(z, k, l);
            });
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(3) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    const auto& out_l = out.dim(3);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                                    value_element(i, j, k, l);
                });

            return matrix;
        }
        if constexpr(Index == 'l') {
            BOOST_ASSERT_MSG(size(2) <= matr1.size(1),
                             "MATRICES 3D * 2D SIZE IS NOT EQUAL, Index == i");
            SizeMatrix3D_context const sizes(size(3), size(3), size(3));
            proto::eval(proto::as_expr<Matrix3D_domain>(matr1), sizes);
            Matrix4D<T> matrix(shape_);
            auto value_element([&](size_t i, size_t j, size_t k, size_t l) {
                for (size_t z = 0; z < size(3); ++z)
                    proto::value(matrix)(i, j, k, l) += proto::value(*this)(i, j, k, z) * proto::value(matr1)(z, l, l);
            });
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(3) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    const auto& out_l = out.dim(3);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                                    value_element(i, j, k, l);
                });

            return matrix;
        }
    }
    template< typename Expr >
    Matrix4D<T> operator * (Expr const& matr1) {
        BOOST_ASSERT_MSG(size(0) == matr1.size(0) && size(1) == matr1.size(1) && size(2) == matr1.size(2) && size(3) == matr1.size(3),
                         "MATRICES SIZE IS NOT EQUAL");
        SizeMatrix4D_context const sizes(size(0), size(1), size(2), size(3));
        proto::eval(proto::as_expr<Matrix4D_domain>(matr1), sizes);
        Matrix4D<T> matrix(shape_);
        for (size_t i = 0; i < size(0); ++i)
            for (size_t j = 0; j < size(1); ++j)
                for (size_t k = 0; k < size(2); ++k)
                for (size_t l = 0; l < size(3); ++l) {
                    proto::value(matrix)(i, j, k, l) =
                            tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(3)), proto::value(matrix)(i, j, k, l),
                                 [=](const tbb::blocked_range<size_t>& r, T tmp) {
                                     for (size_t z = r.begin(); z != r.end(); ++z) {
                                         tmp += proto::value(*this)(i, j, k, z) * matr1(z, j, k, l);
                                     } return tmp; }, std::plus<T>());
                    proto::value(matrix)(i, j, k, l) =
                            tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(0)), proto::value(matrix)(i, j, k, l),
                                 [=](const tbb::blocked_range<size_t>& r, T tmp) {
                                     for (size_t z = r.begin(); z != r.end(); ++z) {
                                         tmp += proto::value(*this)(i, j, z, l) * matr1(z, j, k, l);
                                     } return tmp; }, std::plus<T>());
                    proto::value(matrix)(i, j, k, l) =
                            tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), proto::value(matrix)(i, j, k, l),
                                 [=](const tbb::blocked_range<size_t>& r, T tmp) {
                                     for (size_t z = r.begin(); z != r.end(); ++z) {
                                         tmp += proto::value(*this)(i, z, k, l) * matr1(z, j, k, l);
                                     } return tmp; }, std::plus<T>());
                    proto::value(matrix)(i, j, k, l) =
                            tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), proto::value(matrix)(i, j, k, l),
                                 [=](const tbb::blocked_range<size_t>& r, T tmp) {
                                     for (size_t z = r.begin(); z != r.end(); ++z) {
                                         tmp += proto::value(*this)(i, j, k, z) * matr1(i, z, k, l);
                                     } return tmp; }, std::plus<T>());
                    proto::value(matrix)(i, j, k, l) =
                            tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), proto::value(matrix)(i, j, k, l),
                                 [=](const tbb::blocked_range<size_t>& r, T tmp) {
                                     for (size_t z = r.begin(); z != r.end(); ++z) {
                                         tmp += proto::value(*this)(i, j, z, l) * matr1(i, z, k, l);
                                     } return tmp; }, std::plus<T>());
                    proto::value(matrix)(i, j, k, l) =
                            tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), proto::value(matrix)(i, j, k, l),
                                 [=](const tbb::blocked_range<size_t>& r, T tmp) {
                                     for (size_t z = r.begin(); z != r.end(); ++z) {
                                         tmp += proto::value(*this)(i, j, k, z) * matr1(i, j, z, l);
                                     } return tmp; }, std::plus<T>());
                }

        return matrix;
    }
    Matrix4D<T> operator + (const T val) const {
        Matrix4D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
            const auto& out_i = out.dim(0);
            const auto& out_j = out.dim(1);
            const auto& out_k = out.dim(2);
            const auto& out_l = out.dim(3);
            for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                    for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                        for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                            proto::value(matrix)(i, j, k, l) = proto::value(*this)(i, j, k, l) + val;
        });
        return matrix;
    }
    Matrix4D<T> operator * (const T val) const {
        Matrix4D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
            const auto& out_i = out.dim(0);
            const auto& out_j = out.dim(1);
            const auto& out_k = out.dim(2);
            const auto& out_l = out.dim(3);
            for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                    for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                        for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                            proto::value(matrix)(i, j, k, l) = proto::value(*this)(i, j, k, l) * val;
        });
        return matrix;
    }
    Matrix4D<T> operator / (const T val) const {
        Matrix4D<T> matrix(shape_);
        tbb::parallel_for(range_tbb({0, size(0)}, {0, size(1)}, {0, size(2)}, {0, size(3)}),
        [&](const range_tbb& out){
            const auto& out_i = out.dim(0);
            const auto& out_j = out.dim(1);
            const auto& out_k = out.dim(2);
            const auto& out_l = out.dim(3);
            for(size_t i = out_i.begin(); i < out_i.end(); ++i)
                for(size_t j = out_j.begin(); j < out_j.end(); ++j)
                    for(size_t k = out_k.begin(); k < out_k.end(); ++k)
                        for(size_t l = out_l.begin(); l < out_l.end(); ++l)
                            proto::value(matrix)(i, j, k, l) = proto::value(*this)(i, j, k, l) / val;
        });
        return matrix;
    }
    template<char Index>
    Matrix3D<T> transversal_matrix(int N) {
        static_assert((Index == 'i') || (Index == 'j') || (Index == 'k') || (Index == 'l'), "Не совпадение индексов");
        typename array_type::index_gen indices;
        if constexpr(Index == 'i') {
            std::array<size_t, 3> shi{ {size(1), size(2), size(3)} };
            Matrix3D<T> tmp(shi);
            tmp = proto::value(*this)[indices[N][range(0, size(1))][range(0, size(2))][range(0, size(3))]];
            return tmp;
        } else if constexpr(Index == 'j') {
            std::array<size_t, 3> shj{ {size(0), size(2), size(3)} };
            Matrix3D<T> tmp(shj);
            tmp = proto::value(*this)[indices[range(0, size(0))][N][range(0, size(2))][range(0, size(3))]];
            return tmp;
        } else if constexpr(Index == 'k') {
            std::array<size_t, 3> shk{ {size(0), size(1), size(3)} };
            Matrix3D<T> tmp(shk);
            tmp = proto::value(*this)[indices[range(0, size(0))][range(0, size(1))][N][range(0, size(3))]];
            return tmp;
        } else if constexpr(Index == 'l') {
            std::array<size_t, 3> shl{ {size(0), size(1), size(2)} };
            Matrix3D<T> tmp(shl);
            tmp = proto::value(*this)[indices[range(0, size(0))][range(0, size(1))][range(0, size(2))][N]];
            return tmp;
        }
    }
    template<char Index>
    std::vector<T> transversal_vector(int indx)  {
        static_assert((Index == 'i') || (Index == 'j') || (Index == 'k') || (Index == 'l'), "Не совпадение индексов");
        std::vector<T> transversal_vector;
        transversal_vector.reserve(size(0)*size(1)*size(2)*size(3));
        std::function<size_t(size_t, size_t)> permutation_rec;
        permutation_rec = [&permutation_rec](size_t N, size_t K) {
            size_t tmp = 0;
            if( N <= 0 || K <= 0) {
                tmp = 0;
            } else {
                tmp = (N - 1)*permutation_rec(N - 1, K) + permutation_rec(N - 1, K - 1);
            }
            return tmp;
        };

        if constexpr(Index == 'i') {
            for (size_t j = 0; j != size(1); ++j) {
                for (size_t k = 0; k != size(2); ++k) {
                    for (size_t l = 0; l != size(3); ++l) {
                        if(permutation_rec(size(0), indx)%2 == 0){
                            transversal_vector.push_back(transversal_matrix<'i'>(indx)(j, k, l));
                        }
                    }
                }
            }
        } else if constexpr(Index == 'j') {
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t k = 0; k != size(2); ++k) {
                    for (size_t l = 0; l != size(3); ++l) {
                        if(permutation_rec(size(1), indx)%2 == 0){
                            transversal_vector.push_back(transversal_matrix<'j'>(indx)(i, k, l));
                        }
                    }
                }
            }
        } else if constexpr(Index == 'k') {
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t j = 0; j != size(1); ++j) {
                    for (size_t l = 0; l != size(3); ++l) {
                        if(permutation_rec(size(2), indx)%2 == 0){
                            transversal_vector.push_back(transversal_matrix<'k'>(indx)(i, j, l));
                        }
                    }
                }
            }
        } else if constexpr(Index == 'l') {
            for (size_t i = 0; i != size(0); ++i) {
                for (size_t j = 0; j != size(1); ++j) {
                    for (size_t k = 0; k != size(2); ++k) {
                        if(permutation_rec(size(3), indx)%2 == 0){
                            transversal_vector.push_back(transversal_matrix<'k'>(indx)(i, j, k));
                        }
                    }
                }
            }
        }
        return transversal_vector;
    }
    template<char Index>
    T DET_orient(int indx) {
        BOOST_STATIC_ASSERT_MSG((Index == 'i') || (Index == 'j') || (Index == 'k') || (Index == 'l'), "Не совпадение индексов");
        if constexpr(Index == 'i') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(0)), T(1),
                    [=](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector<'i'>(indx)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if constexpr(Index == 'j') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), T(1),
                    [=](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector<'j'>(indx)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if constexpr(Index == 'k') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(2)), T(1),
                    [=](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector<'k'>(indx)[i];
                        } return tmp; }, std::multiplies<T>());
        } else if constexpr(Index == 'l') {
            return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(3)), T(1),
                    [=](const tbb::blocked_range<size_t>& r, T tmp) {
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            tmp *= transversal_vector<'l'>(indx)[i];
                        } return tmp; }, std::multiplies<T>());
        }
    }
    T DET_FULL() {
        return tbb::parallel_reduce(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(3) }), T(0),
                [=](const range_tbb& out, T tmp) {
                const auto& out_i = out.dim(0);
                const auto& out_j = out.dim(1);
                const auto& out_k = out.dim(2);
                const auto& out_l = out.dim(3);
                for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                    for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                        for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                            for (size_t l = out_l.begin(); l < out_l.end(); ++l)
                                tmp += std::pow(-1, _my::invers_loop<4>({i, j, k, l})[i] + _my::invers_loop<4>({i, j, k, l})[j] + _my::invers_loop<4>({i, j, k, l})[k] + _my::invers_loop<4>({i, j, k, l})[l])
                                        *DET_orient<'i'>(i)*DET_orient<'j'>(j)*DET_orient<'k'>(k)*DET_orient<'l'>(l);
                return tmp; }, std::plus<T>() );
    }

    friend std::ostream& operator << (std::ostream& os, const Matrix4D<T>& A){
            for(const auto& x : proto::value(A)) {
                for(const auto& y : x) {
                    for(const auto& z : y) {
                        for(const auto& f : z) os << f << "\t";
                        os << std::endl;
                    } os << std::endl;
                } os << std::endl;
            } os << std::endl;
            return os;
        }
private:

};

}
