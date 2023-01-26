#pragma once
#include "../include/Matrix2D.hpp"

proto::terminal< std::ostream & >::type cout_ = {std::cout};

namespace _spatial {

    template<typename Expr>
    struct Matrix3D_expr;

    struct Matrix3DGrammar : proto::or_<
        proto::terminal<matrix3_<proto::_>>,
//        proto::plus<Matrix3DGrammar, Matrix3DGrammar>,
        proto::minus<Matrix3DGrammar, Matrix3DGrammar>,
        proto::negate< Matrix3DGrammar>,
        proto::less_equal< Matrix3DGrammar, Matrix3DGrammar>,
        proto::greater_equal< Matrix3DGrammar, Matrix3DGrammar>,
        proto::less< Matrix3DGrammar, Matrix3DGrammar>,
        proto::greater< Matrix3DGrammar, Matrix3DGrammar>,
        proto::not_equal_to< Matrix3DGrammar, Matrix3DGrammar>
//        proto::equal_to< Matrix3DGrammar, Matrix3DGrammar>
    > {};
    struct Matrix3D_domain : proto::domain<proto::generator<Matrix3D_expr>, Matrix3DGrammar>{};

    struct Matrix3D_context : proto::callable_context< Matrix3D_context const > {
        Matrix3D_context(size_t i, size_t j, size_t k) : i(i), j(j), k(k) {}

        template<typename Expr, typename Tag = typename Expr::proto_tag>
        struct eval : proto::default_eval<Expr, Matrix3D_context> {};

        template<typename Expr>
        struct eval<Expr, proto::tag::terminal> {
            using result_type = typename proto::result_of::value<Expr>::type::value_type;

            result_type operator()(const Expr& expr, Matrix3D_context& ctx) const {
                return proto::value(expr)(ctx.i, ctx.j, ctx.k);
            }

        };

    public:
        size_t i;
        size_t j;
        size_t k;
    };

    struct SizeMatrix3D_context {
        SizeMatrix3D_context(size_t Ni, size_t Nj, size_t Nk)
            : NI(Ni), NJ(Nj), NK(Nk) {}
        template<typename Expr, typename EnableIf = void>
        struct eval : proto::null_eval<Expr, SizeMatrix3D_context const> {};

        template<typename Expr>
        struct eval<Expr, typename boost::enable_if<
            proto::matches<Expr, proto::terminal<matrix3_<proto::_> > >
        >::type
        >
        {
            typedef void result_type;

            result_type operator ()(Expr& expr, SizeMatrix3D_context const& ctx) const
            {
                if (ctx.NI != proto::value(expr).shape()[0]) {
                    throw std::runtime_error("Матрицы не совпадают в размерности по индексу i");
                }
                else if (ctx.NJ != proto::value(expr).shape()[1]) {
                    throw std::runtime_error("Матрицы не совпадают в размерности по индексу j");
                }
                else if (ctx.NK != proto::value(expr).shape()[2]) {
                    throw std::runtime_error("Матрицы не совпадают в размерности по индексу k");
                }
            }
        };

        size_t NI;
        size_t NJ;
        size_t NK;
    };

    template<typename Expr>
    struct Matrix3D_expr : proto::extends<Expr, Matrix3D_expr<Expr>, Matrix3D_domain> {
        Matrix3D_expr(const Expr& expr = Expr()) : Matrix3D_expr::proto_extends(expr) {}

        typename proto::result_of::eval<Expr, Matrix3D_context>::type
            operator () (size_t i, size_t j, size_t k) const {
            Matrix3D_context ctx(i, j, k);
            return proto::eval(*this, ctx);
        }
    };
    template<typename T>
    struct Matrix3D : Matrix3D_expr<typename proto::terminal< matrix3_<T>>::type> {
        using expr_type = typename proto::terminal< matrix3_<T>>::type;
        using range_tbb = tbb::blocked_rangeNd<size_t, 3>;
        using range_tbb4 = tbb::blocked_rangeNd<size_t, 4>;

        using array_type = typename matrix3_<T>::array_type;
        using range = boost::multi_array_types::index_range;
        using sub_matrix_view2D = typename matrix3_<T>::sub_matrix_view2D;
        using sub_matrix_view3D = typename matrix3_<T>::sub_matrix_view3D;

        const std::array<size_t, 3>& shape_;
        static constexpr size_t dim = 3;
        constexpr Matrix3D(const std::array<size_t, 3>& shape) :
            Matrix3D_expr<expr_type>(expr_type::make(matrix3_<T>(shape))), shape_(shape) {

        }
        template<typename Expr, typename  = std::enable_if_t<std::is_convertible_v<Expr::expr_type, T>>>
        constexpr Matrix3D(const Expr& matrix) :
            Matrix3D_expr<expr_type>(expr_type::make(matrix3_<T>(matrix.shape_))), shape_(matrix.shape_) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(matrix), sizes);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(*this)(i, j, k) = matrix(i, j, k);
                });

        }
        size_t size(size_t i) const {
            BOOST_ASSERT_MSG((i < 3), "Error i >= 4");
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
                                    proto::value(*this)(i, j, k) = dist(gen);
                }
                if constexpr(!std::is_integral_v<T> && std::is_same_v<T, double>) {
                    boost::random::uniform_real_distribution<> dist{double(min), double(max)};
                    for(size_t i = 0; i < size(0); ++i)
                        for(size_t j = 0; j < size(1); ++j)
                            for(size_t k = 0; k < size(2); ++k)
                                    proto::value(*this)(i, j, k) = dist(gen);
                }
        }
        template<typename F, typename = std::enable_if_t<std::is_same_v<T, std::complex<typename _my::is_type<F>::type>> &&
                                                          std::is_same_v<F, typename _my::is_type<F>::type>>>
        void Random(F min, F max) {
            std::time_t now = std::time(0);
            boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
            boost::random::uniform_real_distribution<> dist{double(min), double(max)};
            for(size_t i = 0; i < size(0); ++i)
                for(size_t j = 0; j < size(1); ++j)
                    for(size_t k = 0; k < size(2); ++k)
                        proto::value(*this)(i, j, k) = std::complex(dist(gen), dist(gen));
        }
        template< typename Expr >
        Matrix3D<T>& operator = (Expr const& expr) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(expr), sizes);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(*this)(i, j, k) = expr(i, j, k);
                });
            return *this;
        }
        Matrix3D<T>& operator = (sub_matrix_view3D const& expr) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(expr), sizes);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(*this)(i, j, k) = expr[i][j][k];
                });
            return *this;
        }
        template< typename Expr >
        Matrix3D<T>& operator += (Expr const& expr) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(expr), sizes);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(*this)(i, j, k) += expr(i, j, k);
                });
            return *this;
        }
        template< typename Expr >
        bool operator == (Expr const& expr) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(expr), sizes);
            size_t tmp_ = tbb::parallel_reduce(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }), size_t(0),
                    [=](const range_tbb& out, size_t tmp) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                tmp += proto::value(*this)(i, j, k) != expr(i, j, k) ? 1 : 0;
                    return tmp; }, std::plus<size_t>() );
            size_t tmp2_ = tbb::parallel_reduce(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }), size_t(0),
                    [=](const range_tbb& out, size_t tmp) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                tmp += proto::value(*this)(i, j, k) == expr(i, j, k) ? 1 : 0;
                    return tmp; }, std::plus<size_t>() );
            std::cout << "EQUAL = " << tmp_ << "; " << tmp2_ << std::endl;
            if(tmp_ == 0) return true;
            else return false;
        }
        template< typename Expr >
        Matrix3D<T> operator + (Expr const& expr) {
            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(expr), sizes);
            Matrix3D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(matrix)(i, j, k) = proto::value(*this)(i, j, k) + expr(i, j, k);
                });
            return matrix;
        }
        Matrix3D<T> operator + (const T val) const {
            Matrix3D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(matrix)(i, j, k) = proto::value(*this)(i, j, k) + val;
                });
            return matrix;
        }
        template<char Index, typename F, typename =
                 std::enable_if_t<std::is_same_v<T, F> && (Index == 'i' || Index == 'j' || Index == 'k')> >
        Matrix3D<T> product(Matrix2D<F> const& matr1) {
            if constexpr(Index == 'i') {
                BOOST_ASSERT_MSG(size(0) <= matr1.size(1),
                                 "MATRICES 3D * 2D SIZE IS NOT EQUAL, Index == i");
                SizeMatrix2D_context const sizes(size(0), size(0));
                proto::eval(proto::as_expr<Matrix2D_domain>(matr1), sizes);
                Matrix3D<T> matrix(shape_);
                auto value_element([&](size_t i, size_t j, size_t k) {
                    for (size_t l = 0; l < size(0); ++l)
                        proto::value(matrix)(i, j, k) += proto::value(*this)(l, j, k) * proto::value(matr1)(l, i);
                });
                tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                    [&](const range_tbb& out) {
                        const auto& out_i = out.dim(0);
                        const auto& out_j = out.dim(1);
                        const auto& out_k = out.dim(2);
                        for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                            for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                                for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                        value_element(i, j, k);
                    });
                return matrix;
            }
            if constexpr(Index == 'j') {
                BOOST_ASSERT_MSG(size(0) <= matr1.size(1),
                                 "MATRICES 3D * 2D SIZE IS NOT EQUAL, Index == j");
                SizeMatrix2D_context const sizes(size(1), size(1));
                proto::eval(proto::as_expr<Matrix2D_domain>(matr1), sizes);
                Matrix3D<T> matrix(shape_);
                auto value_element([&](size_t i, size_t j, size_t k) {
                    for (size_t l = 0; l < size(1); ++l)
                        proto::value(matrix)(i, j, k) += proto::value(*this)(i, l, k) * proto::value(matr1)(l, j);
                });
                tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                    [&](const range_tbb& out) {
                        const auto& out_i = out.dim(0);
                        const auto& out_j = out.dim(1);
                        const auto& out_k = out.dim(2);
                        for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                            for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                                for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                        value_element(i, j, k);
                    });
                return matrix;
            }
            if constexpr(Index == 'k') {
                BOOST_ASSERT_MSG(size(2) <= matr1.size(1),
                                 "MATRICES 3D * 2D SIZE IS NOT EQUAL, Index == i");
                SizeMatrix2D_context const sizes(size(2), size(2));
                proto::eval(proto::as_expr<Matrix2D_domain>(matr1), sizes);
                Matrix3D<T> matrix(shape_);
                auto value_element([&](size_t i, size_t j, size_t k) {
                    for (size_t l = 0; l < size(2); ++l)
                        proto::value(matrix)(i, j, k) += proto::value(*this)(i, j, l) * proto::value(matr1)(l, k);
                });
                tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                    [&](const range_tbb& out) {
                        const auto& out_i = out.dim(0);
                        const auto& out_j = out.dim(1);
                        const auto& out_k = out.dim(2);
                        for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                            for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                                for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                        value_element(i, j, k);
                    });
                return matrix;
            }
        }
        template<char Index, typename F, typename =
                 std::enable_if_t<std::is_same_v<T, F> && (Index == 'i' || Index == 'j' || Index == 'k')> >
        Matrix4D<T> product(Matrix3D<F> const& matr1) {

            SizeMatrix3D_context const sizes(size(0), size(1), size(2));
            proto::eval(proto::as_expr<Matrix3D_domain>(matr1), sizes);
            if constexpr(Index == 'i') {
                BOOST_ASSERT_MSG(size(0) <= matr1.size(1),
                                 "MATRICES 3D * 2D SIZE IS NOT EQUAL, Index == i");
                std::array<size_t, 4> shape4D = { {size(0), size(0), size(1), size(2)} };
                Matrix4D<T> matrix(shape4D);
                auto value_element([&](size_t k1, size_t k2, size_t k3, size_t k4) {
                    for (size_t l = 0; l < size(2); ++l)
                        proto::value(matrix)(k1, k2, k3, k4) += proto::value(*this)(l, k2, k3) * proto::value(matr1)(l, k1, k4);
                });
                tbb::parallel_for(range_tbb4({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }, { 0, size(2) }),
                    [&](const range_tbb4& out) {
                        const auto& out_k1 = out.dim(0);
                        const auto& out_k2 = out.dim(1);
                        const auto& out_k3 = out.dim(2);
                        const auto& out_k4 = out.dim(3);
                        for (size_t k1 = out_k1.begin(); k1 < out_k1.end(); ++k1)
                            for (size_t k2 = out_k2.begin(); k2 < out_k2.end(); ++k2)
                                for (size_t k3 = out_k3.begin(); k3 < out_k3.end(); ++k3)
                                    for (size_t k4 = out_k4.begin(); k4 < out_k4.end(); ++k4)
                                        value_element(k1, k2, k3, k4);
                    });
                return matrix;
            }
            if constexpr(Index == 'j') {
                BOOST_ASSERT_MSG(size(0) <= matr1.size(1) && size(0) == size(1),
                                 "MATRICES 3D * 3D SIZE IS NOT EQUAL, Index == j");
                std::array<size_t, 4> shape4D = { {size(0), size(1), size(1), size(2)} };
                Matrix4D<T> matrix(shape4D);
                auto value_element([&](size_t k1, size_t k2, size_t k3, size_t k4) {
                    for (size_t l = 0; l < size(1); ++l)
                        proto::value(matrix)(k1, k2, k3, k4) += proto::value(*this)(k1, l, k3) * proto::value(matr1)(l, k2, k4);
                });
                tbb::parallel_for(range_tbb4({ 0, size(0) }, { 0, size(1) }, { 0, size(1) }, { 0, size(2) }),
                    [&](const range_tbb4& out) {
                        const auto& out_k1 = out.dim(0);
                        const auto& out_k2 = out.dim(1);
                        const auto& out_k3 = out.dim(2);
                        const auto& out_k4 = out.dim(3);
                        for (size_t k1 = out_k1.begin(); k1 < out_k1.end(); ++k1)
                            for (size_t k2 = out_k2.begin(); k2 < out_k2.end(); ++k2)
                                for (size_t k3 = out_k3.begin(); k3 < out_k3.end(); ++k3)
                                    for (size_t k4 = out_k4.begin(); k4 < out_k4.end(); ++k4)
                                        value_element(k1, k2, k3, k4);
                    });
                return matrix;
            }
            if constexpr(Index == 'k') {
                BOOST_ASSERT_MSG(size(0) == matr1.size(0) && size(0) == size(2),
                                 "MATRICES 3D * 3D SIZE IS NOT EQUAL, Index == i");
                std::array<size_t, 4> shape4D = { {size(0), size(1), size(2), size(2)} };
                Matrix4D<T> matrix(shape4D);
                auto value_element([&](size_t k1, size_t k2, size_t k3, size_t k4) {
                    for (size_t l = 0; l < size(2); ++l)
                        proto::value(matrix)(k1, k2, k3, k4) += proto::value(*this)(k1, k2, l) * proto::value(matr1)(l, k3, k4);
                });
                tbb::parallel_for(range_tbb4({ 0, size(0) }, { 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                    [&](const range_tbb4& out) {
                        const auto& out_k1 = out.dim(0);
                        const auto& out_k2 = out.dim(1);
                        const auto& out_k3 = out.dim(2);
                        const auto& out_k4 = out.dim(3);
                        for (size_t k1 = out_k1.begin(); k1 < out_k1.end(); ++k1)
                            for (size_t k2 = out_k2.begin(); k2 < out_k2.end(); ++k2)
                                for (size_t k3 = out_k3.begin(); k3 < out_k3.end(); ++k3)
                                    for (size_t k4 = out_k4.begin(); k4 < out_k4.end(); ++k4)
                                        value_element(k1, k2, k3, k4);
                    });
                return matrix;
            }
        }

        Matrix3D<T> operator * (const T val) const {
            Matrix3D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(matrix)(i, j, k) = proto::value(*this)(i, j, k) * val;
                });
            return matrix;
        }
        Matrix3D<T> operator / (const T val) const {
            Matrix3D<T> matrix(shape_);
            tbb::parallel_for(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }),
                [&](const range_tbb& out) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                proto::value(matrix)(i, j, k) = proto::value(*this)(i, j, k) / val;
                });
            return matrix;
        }

        Matrix2D<T> transversal_matrix(char index, int N) {
            BOOST_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') , "Не совпадение индексов");
            typename array_type::index_gen indices;
            if (index == 'i') {
                std::array<size_t, 2> shi{ {size(1), size(2)} };
                Matrix2D<T> tmp(shi);
                tmp = proto::value(*this)[indices[N][range(0, size(1))][range(0, size(2))]];
                return tmp;
            } else if (index == 'j') {
                std::array<size_t, 2> shj{ {size(0), size(2)} };
                Matrix2D<T> tmp(shj);
                tmp = proto::value(*this)[indices[range(0, size(0))][N][range(0, size(2))]];
                return tmp;
            } else if (index == 'k') {
                std::array<size_t, 2> shk{ {size(0), size(1)} };
                Matrix2D<T> tmp(shk);
                tmp = proto::value(*this)[indices[range(0, size(0))][range(0, size(1))][N]];
                return tmp;
            }
        }
        std::vector<T> transversal_vector(char index, int N)  {
            BOOST_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') , "Не совпадение индексов");
            typename array_type::index_gen indices;
            std::vector<T> transversal_vector;
            transversal_vector.reserve(size(0)*size(1)*size(2));
            struct index2D { int J; int K; };
            auto predicate([s_ = std::max({size(0), size(1), size(2)})](int IndexN, std::vector<index2D>& indexJK) ->bool {
                int j{};
                do {
                    j = IndexN == 0 ? IndexN : IndexN-1;
                    if(j >= 0 && j < indexJK.size() - 1)
                        while ((indexJK[j].J >= indexJK[j + 1].J) && (indexJK[j].K >= indexJK[j + 1].K)) j--;
                    else j++;
                } while (j > (int)s_ - 1);
                //std::cout << "INDEX J: " << j << std::endl;
                if (j < 0) return false;
                else return true;
            });
            if (index == 'i') {
                std::vector<index2D> tmp;
                tmp.reserve(size(0)*size(1)*size(2));
                for (size_t j = 0; j != size(1); ++j) {
                    for (size_t k = 0; k != size(2); ++k) {
                        tmp.push_back({int(j), int(k)});
                        if(predicate(int(N), tmp)){
                            transversal_vector.push_back(transversal_matrix('i', N)(j, k));
                        }
                    }
                }
            } else if (index == 'j') {
                std::vector<index2D> tmp;
                tmp.reserve(size(0)*size(1)*size(2));
                for (size_t i = 0; i != size(0); ++i) {
                    for (size_t k = 0; k != size(2); ++k) {
                        tmp.push_back({int(i), int(k)});
                        if(predicate(int(N), tmp)){
                            transversal_vector.push_back(transversal_matrix('j', N)(i, k));
                        }
                    }
                }
            } else if (index == 'k') {
                std::vector<index2D> tmp;
                tmp.reserve(size(0)*size(1)*size(2));
                for (size_t i = 0; i != size(0); ++i) {
                    for (size_t j = 0; j != size(1); ++j) {
                        tmp.push_back({int(i), int(j)});
                        if(predicate(int(N), tmp)){
                            transversal_vector.push_back(transversal_matrix('k', N)(i, j));
                        }
                    }
                }
            }
            return transversal_vector;
        }
        T DET_orient(char index, int N) {
            BOOST_ASSERT_MSG((index == 'i') || (index == 'j') || (index == 'k') , "Не совпадение индексов");
            if(index == 'i') {
                return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(0)), T(1),
                        [=](const tbb::blocked_range<size_t>& r, T tmp) {
                            for (size_t i = r.begin(); i != r.end(); ++i) {
                                tmp *= transversal_vector('i', N)[i];
                            } return tmp; }, std::multiplies<T>());
            } else if(index == 'j') {
                return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(1)), T(1),
                        [=](const tbb::blocked_range<size_t>& r, T tmp) {
                            for (size_t i = r.begin(); i != r.end(); ++i) {
                                tmp *= transversal_vector('j', N)[i];
                            } return tmp; }, std::multiplies<T>());
            } else if(index == 'k') {
                return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size(2)), T(1),
                        [=](const tbb::blocked_range<size_t>& r, T tmp) {
                            for (size_t i = r.begin(); i != r.end(); ++i) {
                                tmp *= transversal_vector('k', N)[i];
                            } return tmp; }, std::multiplies<T>());
            }
        }
        T DET_FULL() {
            return tbb::parallel_reduce(range_tbb({ 0, size(0) }, { 0, size(1) }, { 0, size(2) }), T(0),
                    [=](const range_tbb& out, T tmp) {
                    const auto& out_i = out.dim(0);
                    const auto& out_j = out.dim(1);
                    const auto& out_k = out.dim(2);
                    for (size_t i = out_i.begin(); i < out_i.end(); ++i)
                        for (size_t j = out_j.begin(); j < out_j.end(); ++j)
                            for (size_t k = out_k.begin(); k < out_k.end(); ++k)
                                tmp += std::pow(-1, _my::invers<3>(i, j, k))*DET_orient('i', i)*DET_orient('j', j)*DET_orient('k', k);
                    return tmp; }, std::plus<T>() );
        }
        template< typename Expr >
        void evaluate( Expr const & expr )
        {
            proto::default_context ctx;
            proto::eval(expr, ctx);
        }

        friend std::ostream& operator << (std::ostream& os, const Matrix3D<T>& A) {
            for (const auto& x : proto::value(A)) {
                for (const auto& y : x) {
                    for (const auto& f : y) os << f << "\t";
                    os << std::endl;
                } os << std::endl;
            } os << std::endl;
            return os;
        }
    private:

    };

}
