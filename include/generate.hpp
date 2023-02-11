#ifndef GENERATE_HPP
#define GENERATE_HPP

#include "Matrix_type.hpp"


struct PASSAGES_ITEM {
    int default_zone_v       = 0;
    int penalty           = 4;
    int penalty_fromRoute = -1;
    size_t Length            = 0;
    friend std::ostream& operator << (std::ostream& os, const PASSAGES_ITEM& A) {
        os << "{ default_zone_v = " << A.default_zone_v << "; penalty = " << A.penalty
           << "; penalty_fromRoute = " << A.penalty_fromRoute
           << "; Length = " << A.Length << " }" << std::endl;;
        return os;
    }
};
struct KNOT {
    int X = -1;
    int Y = -1;
    int Z = -1;
    friend std::ostream& operator << (std::ostream& os, const KNOT& A) {
        os << "{ X = " << A.X << "; Y = " << A.Y
           << "; Z = " << A.Z << " }" << std::endl;;
        return os;
    }
};
static constexpr int indexSimplex[3][6] = {
                                            {-1, 1, 0, 0, 0, 0},
                                            {0, 0, -1, 1, 0, 0},
                                            {0, 0, 0, 0, -1, 1}
                                          };


struct for_one {
private:
    const MATRIX<3, int>& MAZE;
    MATRIX<3, PASSAGES_ITEM>& matrix;
    const int& key;
    size_t x; size_t y; size_t z;
public:
    for_one(const MATRIX<3, int>& MAZE, MATRIX<3, PASSAGES_ITEM>& matrix, const int& key, size_t x, size_t y, size_t z) :
        MAZE(MAZE), matrix(matrix), key(key), x(x), y(y), z(z)  {}
    template <int J>
    void apply() {
        int index_X = x + indexSimplex[0][J];
        int index_Y = y + indexSimplex[1][J];
        int index_Z = z + indexSimplex[2][J];
        if(index_X < MAZE.size(0) && index_Y < MAZE.size(1) && index_Z < MAZE.size(2)) {
            proto::value(matrix)(x, y, z).penalty += MAZE(index_X, index_Y, index_Z) == 6 ?
                        10 : (key == 0 && MAZE(index_X, index_Y, index_Z) == 4 ? 1 : 0);
        }
    }
};

inline void for_(const MATRIX<3, int>& MAZE, MATRIX<3, PASSAGES_ITEM>& matrix, const int& key, size_t x, size_t y, size_t z) {
    for_one closure(MAZE, matrix, key, x, y, z);
    _my::meta_loop<6>(closure);
}

class generate
{
    using range_tbb = tbb::blocked_rangeNd<size_t, 3>;
    static constexpr std::array<size_t, 4> penaltyWeight = { {0, 100000, 1000000, 100000000} };
    static constexpr size_t penaltyNearZero = (penaltyWeight[1] + penaltyWeight[2])/2;
    const MATRIX<3, int>& MAZE;
public:
    generate(const MATRIX<3, int>& MAZE_) : MAZE(MAZE_) {}
private:
    MATRIX<3, PASSAGES_ITEM> create_simplex() const {
        const std::array<size_t, 3> shape3D = { {MAZE.size(0), MAZE.size(1), MAZE.size(2)} };

        MATRIX<3, PASSAGES_ITEM> matrix(shape3D);
        auto lambda_bind([this, &matrix] (size_t x, size_t y, size_t z) {
            const int& key = MAZE(x, y, z)/2;
            proto::value(matrix)(x, y, z).default_zone_v = key;
            proto::value(matrix)(x, y, z).penalty = key >= 1 && key <= 6 ? penaltyWeight[key] : 0;
            for_(MAZE, matrix, key, x, y, z);
            if(z < MAZE.size(2) - 1) {
                if(MAZE(x, y, z + 1) == 0) {
                    proto::value(matrix)(x, y, z).penalty += penaltyNearZero;
                }
            }
        });
        tbb::parallel_for(range_tbb({ 0, MAZE.size(0) }, { 0, MAZE.size(1) }, { 0, MAZE.size(2) }),
            [&](const range_tbb& out) {
                const auto& out_i = out.dim(0);
                const auto& out_j = out.dim(1);
                const auto& out_k = out.dim(2);
                for (size_t x = out_i.begin(); x < out_i.end(); ++x)
                    for (size_t y = out_j.begin(); y < out_j.end(); ++y)
                        for (size_t z = out_k.begin(); z < out_k.end(); ++z)
                            lambda_bind(x, y, z);
            });

        return matrix;
    }
public:
    int path_generate(const KNOT& start, const KNOT& finish) const {
        auto tmp_arr = create_simplex();
        const std::array<int, 8> arr_dl = {{0, 100, 1000, 1004, 1000, 1004, 1414, 1417}};
        int _hole_Xmin, _hole_Ymin, _hole_Zmin;
        int _hole_Xmax, _hole_Ymax, _hole_Zmax;
        double x_interval, y_interval, z_interval;
        size_t penalty;
        bool modify = false;
        const std::array<size_t, 2> shape2D = { {8, 3} };

        MATRIX<2, int*> matrix(shape2D);
        for(size_t i = 0; i < matrix.size(0); ++i) {
            if(i < 4) {
                proto::value(matrix)(i, 0) = &_hole_Xmin;
                if(i < 2) {
                    proto::value(matrix)(i, 1) = &_hole_Ymin;
                } else {
                    proto::value(matrix)(i, 1) = &_hole_Ymax;
                }
                if(i%2 == 0) proto::value(matrix)(i, 2) = &_hole_Zmin;
                else proto::value(matrix)(i, 2) = &_hole_Zmax;
            } else {
                proto::value(matrix)(i, 0) = &_hole_Xmax;
                if(i < 6) {
                    proto::value(matrix)(i, 1) = &_hole_Ymin;
                } else {
                    proto::value(matrix)(i, 1) = &_hole_Ymax;
                }
                if(i%2 == 0) proto::value(matrix)(i, 2) = &_hole_Zmin;
                else proto::value(matrix)(i, 2) = &_hole_Zmax;
            }
        }
        std::cout << tmp_arr << std::endl;
        if(tmp_arr(start.X, start.Y, start.Z).default_zone_v == 0 || tmp_arr(finish.X, finish.Y, finish.Z).default_zone_v == 0) {
            return -1;
        }
        proto::value(tmp_arr)(start.X, start.Y, start.Z).penalty_fromRoute = proto::value(tmp_arr)(start.X, start.Y, start.Z).penalty;
        proto::value(tmp_arr)(start.X, start.Y, start.Z).Length = 0;
        MATRIX<3, PASSAGES_ITEM> matrixCopy(tmp_arr);
        MATRIX<3, PASSAGES_ITEM> matrixDest(tmp_arr);
        MATRIX<3, PASSAGES_ITEM> matrixSrc(tmp_arr);
        int Xmin = std::max(0, start.X);
        int Ymin = std::max(0, start.Y);
        int Zmin = std::max(0, start.Z);

        int Xmax = std::min(int(tmp_arr.size(0) - 1), finish.X);
        int Ymax = std::min(int(tmp_arr.size(1) - 1), finish.Y);
        int Zmax = std::min(int(tmp_arr.size(2) - 1), finish.Z);

        x_interval = (Xmin + Xmax)/2.;
        y_interval = (Ymin + Ymax)/2.;
        z_interval = (Zmin + Zmax)/2.;

        int hole_Xmin = Xmax; int hole_Xmax = Xmin;
        int hole_Ymin = Ymax; int hole_Ymax = Ymin;
        int hole_Zmin = Zmax; int hole_Zmax = Zmin;

        int PMax = tmp_arr.size(0)*(tmp_arr.size(1) - 1)*(tmp_arr.size(2)/2 + 1)/2 + 1;

        for(int p = 1; p < PMax; ++p) {
            modify = false;
            int _Xmin = tmp_arr.size(0); int _Ymin = tmp_arr.size(1);int _Zmin = tmp_arr.size(2);
            int _hole_Xmin = Xmin; int _hole_Xmax = Xmax;
            int _hole_Ymin = Ymin; int _hole_Ymax = Ymax;
            int _hole_Zmin = Zmin; int _hole_Zmax = Zmax;
            for(int z = Zmin; z <= Zmax; ++z) {
                int Z_env_min = std::max(0, z - 1);
                int Z_env_max = std::min(int(tmp_arr.size(2)), z + 1);
                int Z_YN = z*tmp_arr.size(1);
                for(int y = Ymin; y <= Ymax; ++y) {
                    int Y_env_min = std::max(0, y - 1);
                    int Y_env_max = std::min(int(tmp_arr.size(1)), y + 1);
                    int Z_YN_Y_XN = (Z_YN + y)*tmp_arr.size(0);
                    for(int x = Xmin; x <= Xmax; ++x) {
                        if(z >= hole_Zmin && z <= hole_Zmax &&
                           y >= hole_Ymin && y <= hole_Ymax &&
                           x >= hole_Xmin && x <= hole_Xmax) continue;
                        int X_env_min = std::max(0, x - 1);
                        int X_env_max = std::min(int(tmp_arr.size(0)), x + 1);
                        if(matrixDest(x, y, z).default_zone_v < 0) continue;
                        matrixDest(x, y, z) = matrixSrc(x, y, z);
                        int indexHole = (z > z_interval)*4 + (y > y_interval)*2 + (x > x_interval);
                        for(int z_env = Z_env_min; z_env <= Z_env_max; ++z_env) {

                        }
                    }
                }
            }
        }
    }
};

#endif // GENERATE_HPP
