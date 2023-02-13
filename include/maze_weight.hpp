#pragma once

#include <boost/mpl/if.hpp>
#include <boost/unordered_map.hpp>
#include <array>
#include <complex>

namespace mpl = boost::mpl;


static constexpr uint8_t FLOOR = 1;   //00 00 00 01;
static constexpr uint8_t EAST  = 2;   //00 00 00 10;
static constexpr uint8_t NORTH = 4;   //00 00 01 00;
static constexpr uint8_t WEST  = 8;   //00 00 10 00;
static constexpr uint8_t SOUTH = 16;  //00 01 00 00;
static constexpr uint8_t CEIL  = 32;  //00 10 00 00;
static constexpr uint8_t FORW  = 64;  //01 00 00 00;
static constexpr uint8_t BACKW = 128; //10 00 00 00;

template<typename T>
struct Cmpare_complex {
  Cmpare_complex(bool flip) : flip(flip) {}

  bool operator()(const T& lhs, const T& rhs) const {
      if constexpr(std::is_same_v<T, int>) {
         return flip ? lhs < rhs : rhs < lhs;
      }
      if constexpr(std::is_same_v<T, std::complex<int>>) {
         return flip ? std::norm(lhs) < std::norm(rhs) : std::norm(rhs) < std::norm(lhs);
      }
  }

  bool flip;
};

template<typename T, template<typename> class Cont>
using Base_Cont = std::enable_if_t<(Cont<T>::dim == 3 || Cont<T>::dim == 4),
                                    typename mpl::if_c<std::is_same_v<T, int>,
                                                            Cont<int>,
                                    typename mpl::if_c<std::is_same_v<T, std::complex<int> >,
                                                            Cont<std::complex<int> >,
                                                            mpl::false_
                                             >::type
                                             >::type
                                             >;

template<typename T, template<typename> class Cont>
class maze_weight : public Base_Cont<T, Cont>


{
    using vector_pos_type = typename mpl::if_c<(Cont<T>::dim == 3),
                                                std::vector<std::tuple<int, int, int>>,
                                                std::vector<std::tuple<int, int, int, int>>
                                                >::type;
    using passages_type = typename mpl::if_c<(Cont<T>::dim == 3),
                                              std::vector<std::tuple<int, int, int, int, int, int>>,
                                              std::vector<std::tuple<int, int, int, int, int, int, int, int>>
                                              >::type;
    using base_type = Base_Cont<T, Cont>;
    passages_type passages;
private:
    void calc_cell_values() {
        if constexpr(base_type::dim == 3) {

            int i_from, i_to;
            int j_from, j_to;
            int k_from, k_to;
            if constexpr(std::is_same_v<T, int>) {
                for(size_t i = 0; i < this->size(0); ++i)
                    for(size_t j = 0; j < this->size(1); ++j)
                        for(size_t k = 0; k < this->size(2); ++k)
                            proto::value(*this)(i, j, k) = 6;
                for(auto& pass : passages) {
                    std::tie(i_from, j_from, k_from, i_to, j_to, k_to) = pass;
                    uint8_t cell_num_from = proto::value(*this)(i_from, j_from, k_from);
                    uint8_t cell_num_to   = proto::value(*this)(i_to, j_to, k_to);
                    if(i_from != i_to) {
                        proto::value(*this)(i_from, j_from, k_from) = cell_num_from & ~(SOUTH);
                        proto::value(*this)(i_to, j_to, k_to)       = cell_num_to   & ~(NORTH);
                    } else if(j_from != j_to) {
                        proto::value(*this)(i_from, j_from, k_from) = cell_num_from & ~(EAST);
                        proto::value(*this)(i_to, j_to, k_to)       = cell_num_to   & ~(WEST);
                    } else if(k_from != k_to) {
                        proto::value(*this)(i_from, j_from, k_from) = cell_num_from & ~(CEIL);
                        proto::value(*this)(i_to, j_to, k_to)       = cell_num_to   & ~(FLOOR);
                    }
                }
            } else {
                for(size_t i = 0; i < this->size(0); ++i)
                    for(size_t j = 0; j < this->size(1); ++j)
                        for(size_t k = 0; k < this->size(2); ++k)
                            proto::value(*this)(i, j, k) = std::complex<int >(6, 12);
                for(auto& pass : passages) {
                    std::tie(i_from, j_from, k_from, i_to, j_to, k_to) = pass;
                    std::complex<uint8_t > cell_num_from = proto::value(*this)(i_from, j_from, k_from);
                    std::complex<uint8_t > cell_num_to   = proto::value(*this)(i_to, j_to, k_to);
                    if(i_from != i_to) {
                        proto::value(*this)(i_from, j_from, k_from) = std::complex<uint8_t >(cell_num_from.real() & ~(SOUTH), cell_num_from.imag() & ~(SOUTH));
                        proto::value(*this)(i_to, j_to, k_to)       = std::complex<uint8_t >(cell_num_to.real()   & ~(NORTH), cell_num_to.imag()   & ~(NORTH));
                    } else if(j_from != j_to) {
                        proto::value(*this)(i_from, j_from, k_from) = std::complex<uint8_t >(cell_num_from.real() & ~(EAST), cell_num_from.imag() & ~(EAST));
                        proto::value(*this)(i_to, j_to, k_to)       = std::complex<uint8_t >(cell_num_to.real()   & ~(WEST), cell_num_to.imag()   & ~(WEST));
                    } else if(k_from != k_to) {
                        proto::value(*this)(i_from, j_from, k_from) = std::complex<uint8_t >(cell_num_from.real() & ~(CEIL), cell_num_from.imag() & ~(CEIL));
                        proto::value(*this)(i_to, j_to, k_to)       = std::complex<uint8_t >(cell_num_to.real()   & ~(FLOOR), cell_num_to.imag()   & ~(FLOOR));
                    }
                }
            }
        } else {
            int i_from, i_to;
            int j_from, j_to;
            int k_from, k_to;
            int l_from, l_to;
            if constexpr(std::is_same_v<T, int>) {
                for(size_t i = 0; i < this->size(0); ++i)
                    for(size_t j = 0; j < this->size(1); ++j)
                        for(size_t k = 0; k < this->size(2); ++k)
                            for(size_t l = 0; l < this->size(3); ++l)
                                proto::value(*this)(i, j, k, l) = 6;
                for(auto& pass : passages) {
                    std::tie(i_from, j_from, k_from, l_from, i_to, j_to, k_to, l_to) = pass;
                    uint8_t cell_num_from = proto::value(*this)(i_from, j_from, k_from, l_from);
                    uint8_t cell_num_to   = proto::value(*this)(i_to, j_to, k_to, l_to);
                    if(i_from != i_to) {
                        proto::value(*this)(i_from, j_from, k_from, l_from) = cell_num_from & ~(SOUTH);
                        proto::value(*this)(i_to, j_to, k_to, l_to)         = cell_num_to   & ~(NORTH);
                    } else if(j_from != j_to) {
                        proto::value(*this)(i_from, j_from, k_from, l_from) = cell_num_from & ~(EAST);
                        proto::value(*this)(i_to, j_to, k_to, l_to)         = cell_num_to   & ~(WEST);
                    } else if(k_from != k_to) {
                        proto::value(*this)(i_from, j_from, k_from, l_from) = cell_num_from & ~(CEIL);
                        proto::value(*this)(i_to, j_to, k_to, l_to)         = cell_num_to   & ~(FLOOR);
                    } else if(l_from != l_to) {
                        proto::value(*this)(i_from, j_from, k_from, l_from) = cell_num_from & ~(FORW);
                        proto::value(*this)(i_to, j_to, k_to, l_to)         = cell_num_to   & ~(BACKW);
                    }
                }
            } else {
                for(size_t i = 0; i < this->size(0); ++i)
                    for(size_t j = 0; j < this->size(1); ++j)
                        for(size_t k = 0; k < this->size(2); ++k)
                            for(size_t l = 0; l < this->size(3); ++l)
                                proto::value(*this)(i, j, k, l) = std::complex<int >(6, 12);
                for(auto& pass : passages) {
                    std::tie(i_from, j_from, k_from, l_from, i_to, j_to, k_to, l_to) = pass;
                    std::complex<uint8_t > cell_num_from = proto::value(*this)(i_from, j_from, k_from, l_from);
                    std::complex<uint8_t > cell_num_to   = proto::value(*this)(i_to, j_to, k_to, l_to);
                    if(i_from != i_to) {
                        proto::value(*this)(i_from, j_from, k_from, l_from) = std::complex<uint8_t >(cell_num_from.real() & ~(SOUTH), cell_num_from.imag() & ~(SOUTH));
                        proto::value(*this)(i_to, j_to, k_to, l_to)         = std::complex<uint8_t >(cell_num_to.real()   & ~(NORTH), cell_num_to.imag()   & ~(NORTH));
                    } else if(j_from != j_to) {
                        proto::value(*this)(i_from, j_from, k_from, l_from) = std::complex<uint8_t >(cell_num_from.real() & ~(EAST), cell_num_from.imag() & ~(EAST));
                        proto::value(*this)(i_to, j_to, k_to, l_to)         = std::complex<uint8_t >(cell_num_to.real()   & ~(WEST), cell_num_to.imag()   & ~(WEST));
                    } else if(k_from != k_to) {
                        proto::value(*this)(i_from, j_from, k_from, l_from) = std::complex<uint8_t >(cell_num_from.real() & ~(CEIL), cell_num_from.imag() & ~(CEIL));
                        proto::value(*this)(i_to, j_to, k_to, l_to)         = std::complex<uint8_t >(cell_num_to.real()   & ~(FLOOR), cell_num_to.imag()   & ~(FLOOR));
                    } else if(l_from != l_to) {
                        proto::value(*this)(i_from, j_from, k_from, l_from) = std::complex<uint8_t >(cell_num_from.real() & ~(FORW), cell_num_from.imag() & ~(FORW));
                        proto::value(*this)(i_to, j_to, k_to, l_to)         = std::complex<uint8_t >(cell_num_to.real()   & ~(BACKW), cell_num_to.imag()   & ~(BACKW));
                    }
                }
            }
        }
    }
public:
    maze_weight() = delete;
    template<typename F = T, template<typename> class Arr = Cont, typename std::enable_if_t<Arr<T>::dim == 3>* = nullptr>
    maze_weight(const std::array<size_t, 3>& shape) : base_type(shape) {
        if constexpr(std::is_same_v<F, int>) {
            for(size_t i = 0; i < this->size(0); ++i)
                for(size_t j = 0; j < this->size(1); ++j)
                    for(size_t k = 0; k < this->size(2); ++k)
                        proto::value(*this)(i, j, k) += k + this->size(2)*(j + this->size(1)*i);
        } else {
            for(size_t i = 0; i < this->size(0); ++i)
                for(size_t j = 0; j < this->size(1); ++j)
                    for(size_t k = 0; k < this->size(2); ++k)
                        proto::value(*this)(i, j, k) += std::complex<int>(k + this->size(2)*(j + this->size(1)*i), k + this->size(2)*(j + this->size(1)*i));
        }
    }
    template<typename F = T, template<typename> class Arr = Cont, typename std::enable_if_t<Arr<T>::dim == 4>* = nullptr>
    maze_weight(const std::array<size_t, 4>& shape) : base_type(shape) {
        if constexpr(std::is_same_v<F, int>) {
            for(size_t i = 0; i < this->size(0); ++i)
                for(size_t j = 0; j < this->size(1); ++j)
                    for(size_t k = 0; k < this->size(2); ++k)
                        for(size_t l = 0; l < this->size(3); ++l)
                            proto::value(*this)(i, j, k, l) += l + this->size(3)*(k + this->size(2)*(j + this->size(1)*i));
        } else {
            for(size_t i = 0; i < this->size(0); ++i)
                for(size_t j = 0; j < this->size(1); ++j)
                    for(size_t k = 0; k < this->size(2); ++k)
                        for(size_t l = 0; l < this->size(3); ++l)
                            proto::value(*this)(i, j, k, l) += std::complex<int>(l + this->size(3)*(k + this->size(2)*(j + this->size(1)*i)), l + this->size(3)*(k + this->size(2)*(j + this->size(1)*i)));
        }
    }

    void calc_maze(const double horisontal_bias, const double vertical_bias) {
        const double east_wall_threshold = horisontal_bias;
        const double south_wall_threshold = vertical_bias;
        if constexpr(base_type::dim == 3) {
            T room_value;
            T east_room_value;
            T south_room_value;
            std::map<T, vector_pos_type, Cmpare_complex<T>> room_set(Cmpare_complex<T>(true));
            std::map<T, T, Cmpare_complex<T>> merged_room_sets(Cmpare_complex<T>(true));
            srand(time(NULL));
            for(size_t k = 0; k < this->size(2) - 1; k++) {
                for(size_t i = 0; i < this->size(0); i++) {
                    for(size_t j = 0; j < this->size(1); j++) {
                        room_value = proto::value(*this)(i, j, k);
                        if(merged_room_sets.find(room_value) != merged_room_sets.end()) {
                            room_value = merged_room_sets[room_value];
                        } else {
                            merged_room_sets[room_value] = room_value;
                        }
                        if(j < this->size(1) - 1) {
                            east_room_value = (merged_room_sets.find(proto::value(*this)(i, j + 1, k)) != merged_room_sets.end()) ?
                                                merged_room_sets[proto::value(*this)(i, j + 1, k)] : proto::value(*this)(i, j + 1, k);
                        }
                        if(i < this->size(0) - 1) {
                            south_room_value = (merged_room_sets.find(proto::value(*this)(i + 1, j, k)) != merged_room_sets.end()) ?
                                                merged_room_sets[proto::value(*this)(i + 1, j, k)] : proto::value(*this)(i + 1, j, k);
                        }
                        room_set[room_value].push_back(std::make_tuple(i, j, k));
                        if(j == this->size(1) - 1) {}
                        else if((rand() < east_wall_threshold*RAND_MAX) && (room_value != east_room_value)) {
                            vector_pos_type next_room_position_set = room_set[proto::value(*this)(i, j + 1, k)];
                            room_set[room_value].insert(room_set[room_value].end(), next_room_position_set.begin(), next_room_position_set.end());
                            room_set.erase(proto::value(*this)(i, j + 1, k));
                            merged_room_sets[proto::value(*this)(i, j + 1, k)] = room_value;
                            proto::value(*this)(i, j + 1, k) = room_value;
                            passages.push_back(std::make_tuple(i, j, k, i, j + 1, k));
                        }
                        if(i == this->size(0) - 1) {}
                        else if((rand() < south_wall_threshold*RAND_MAX) && (room_value != south_room_value)) {
                            proto::value(*this)(i + 1, j, k) = room_value;
                            passages.push_back(std::make_tuple(i, j, k, i + 1, j, k));
                        }
                    }
                }
                merged_room_sets.clear();
                for(auto entry : room_set) {
                    auto group = entry.second;
                    int x, y, z;
                    std::tie(y, x, z) = group.at(rand()%group.size());
                    passages.push_back(std::make_tuple(y, x, z, y, x, z + 1));
                }
                room_set.clear();
            }
            int floor = this->size(2) - 1;
            std::set<T, Cmpare_complex<T>> can_go_south(Cmpare_complex<T>(true));
            for(size_t i = 0; i < this->size(0) - 1; i++) {
                for(size_t j = 0; j < this->size(1); j++) {
                    room_value = proto::value(*this)(i, j, floor);
                    if(merged_room_sets.find(room_value) != merged_room_sets.end()) {
                        room_value = merged_room_sets[room_value];
                    } else {
                        merged_room_sets[room_value] = room_value;
                    }
                    if(j < this->size(1) - 1) {
                        east_room_value = (merged_room_sets.find(proto::value(*this)(i, j + 1, floor)) != merged_room_sets.end()) ?
                                            merged_room_sets[proto::value(*this)(i, j + 1, floor)] : proto::value(*this)(i, j + 1, floor);
                    }
                    south_room_value = (merged_room_sets.find(proto::value(*this)(i + 1, j, floor)) != merged_room_sets.end()) ?
                                        merged_room_sets[proto::value(*this)(i + 1, j, floor)] : proto::value(*this)(i + 1, j, floor);
                    if(j == this->size(1) - 1) {}
                    else if((rand() < east_wall_threshold*RAND_MAX) && (room_value != east_room_value)) {
                        passages.push_back(std::make_tuple(i, j, floor, i, j + 1, floor));
                        merged_room_sets[east_room_value] = room_value;
                        proto::value(*this)(i, j + 1, floor) = room_value;
                    }
                    if(((rand() < south_wall_threshold*RAND_MAX) || (can_go_south.find(room_value) == can_go_south.end()))
                            && (room_value != south_room_value)) {
                        passages.push_back(std::make_tuple(i, j, floor, i + 1, j, floor));
                        merged_room_sets[south_room_value] = room_value;
                        proto::value(*this)(i + 1, j, floor) = room_value;
                        can_go_south.insert(room_value);
                    }
                }
                can_go_south.clear();
            }
            int row = this->size(0) - 1;
            for(size_t j = 0; j < this->size(1) - 1; j++) {
                room_value = proto::value(*this)(row, j, floor);
                if(merged_room_sets.find(room_value) != merged_room_sets.end()) {
                    room_value = merged_room_sets[room_value];
                } else {
                    merged_room_sets[room_value] = room_value;
                }
                east_room_value = (merged_room_sets.find(proto::value(*this)(row, j + 1, floor)) != merged_room_sets.end()) ?
                                    merged_room_sets[proto::value(*this)(row, j + 1, floor)] : proto::value(*this)(row, j + 1, floor);
                if(room_value != east_room_value) {
                    passages.push_back(std::make_tuple(row, j, floor, row, j + 1, floor));
                    merged_room_sets[east_room_value] = room_value;
                    proto::value(*this)(row, j + 1, floor) = room_value;
                }
            }
        } else {
            T room_value;
            T east_room_value;
            T south_room_value;
            std::map<T, vector_pos_type, Cmpare_complex<T>> room_set(Cmpare_complex<T>(true));
            std::map<T, T, Cmpare_complex<T>> merged_room_sets(Cmpare_complex<T>(true));
            srand(time(NULL));
            for(size_t l = 0; l < this->size(3) - 1; l++) {
                for(size_t i = 0; i < this->size(0); i++) {
                    for(size_t j = 0; j < this->size(1); j++)
                        for(size_t k = 0; k < this->size(2); k++) {
                        room_value = proto::value(*this)(i, j, k, l);
                        if(merged_room_sets.find(room_value) != merged_room_sets.end()) {
                            room_value = merged_room_sets[room_value];
                        } else {
                            merged_room_sets[room_value] = room_value;
                        }
                        if(k < this->size(2) - 1) {
                            east_room_value = (merged_room_sets.find(proto::value(*this)(i, j, k + 1, l)) != merged_room_sets.end()) ?
                                                merged_room_sets[proto::value(*this)(i, j, k + 1, l)] : proto::value(*this)(i, j, k + 1, l);
                        }
                        if(j < this->size(1) - 1) {
                            east_room_value = (merged_room_sets.find(proto::value(*this)(i, j + 1, k, l)) != merged_room_sets.end()) ?
                                                merged_room_sets[proto::value(*this)(i, j + 1, k, l)] : proto::value(*this)(i, j + 1, k, l);
                        }
                        if(i < this->size(0) - 1) {
                            south_room_value = (merged_room_sets.find(proto::value(*this)(i + 1, j, k, l)) != merged_room_sets.end()) ?
                                                merged_room_sets[proto::value(*this)(i + 1, j, k, l)] : proto::value(*this)(i + 1, j, k, l);
                        }
                        room_set[room_value].push_back(std::make_tuple(i, j, k, l));
                        if(k == this->size(2) - 1) {}
                        else if((rand() < east_wall_threshold*RAND_MAX) && (room_value != east_room_value)) {
                            vector_pos_type next_room_position_set = room_set[proto::value(*this)(i, j, k + 1, l)];
                            room_set[room_value].insert(room_set[room_value].end(), next_room_position_set.begin(), next_room_position_set.end());
                            room_set.erase(proto::value(*this)(i, j, k + 1, l));
                            merged_room_sets[proto::value(*this)(i, j, k + 1, l)] = room_value;
                            proto::value(*this)(i, j, k + 1, l) = room_value;
                            passages.push_back(std::make_tuple(i, j, k, l, i, j, k + 1, l));
                        }
                        if(j == this->size(1) - 1) {}
                        else if((rand() < east_wall_threshold*RAND_MAX) && (room_value != east_room_value)) {
                            vector_pos_type next_room_position_set = room_set[proto::value(*this)(i, j + 1, k, l)];
                            room_set[room_value].insert(room_set[room_value].end(), next_room_position_set.begin(), next_room_position_set.end());
                            room_set.erase(proto::value(*this)(i, j + 1, k, l));
                            merged_room_sets[proto::value(*this)(i, j + 1, k, l)] = room_value;
                            proto::value(*this)(i, j + 1, k, l) = room_value;
                            passages.push_back(std::make_tuple(i, j, k, l, i, j + 1, k, l));
                        }
                        if(i == this->size(0) - 1) {}
                        else if((rand() < south_wall_threshold*RAND_MAX) && (room_value != south_room_value)) {
                            proto::value(*this)(i + 1, j, k, l) = room_value;
                            passages.push_back(std::make_tuple(i, j, k, l, i + 1, j, k, l));
                        }
                    }
                }
                merged_room_sets.clear();
                for(auto entry : room_set) {
                    auto group = entry.second;
                    int x, y, z, t;
                    std::tie(y, x, z, t) = group.at(rand()%group.size());
                    passages.push_back(std::make_tuple(y, x, z, t, y, x, z, t + 1));
                }
                room_set.clear();
            }
            int floor = this->size(2) - 1;
            std::set<T, Cmpare_complex<T>> can_go_south(Cmpare_complex<T>(true));
            for(size_t i = 0; i < this->size(0) - 1; i++) {
                for(size_t j = 0; j < this->size(1); j++)
                    for(size_t l = 0; l < this->size(3); l++) {
                        room_value = proto::value(*this)(i, j, floor, l);
                        if(merged_room_sets.find(room_value) != merged_room_sets.end()) {
                            room_value = merged_room_sets[room_value];
                        } else {
                            merged_room_sets[room_value] = room_value;
                        }
                        if(l < this->size(3) - 1) {
                            east_room_value = (merged_room_sets.find(proto::value(*this)(i, j, floor, l + 1)) != merged_room_sets.end()) ?
                                                merged_room_sets[proto::value(*this)(i, j, floor, l + 1)] : proto::value(*this)(i, j, floor, l + 1);
                        }
                        south_room_value = (merged_room_sets.find(proto::value(*this)(i + 1, j, floor, l)) != merged_room_sets.end()) ?
                                            merged_room_sets[proto::value(*this)(i + 1, j, floor, l)] : proto::value(*this)(i + 1, j, floor, l);
                        if(j < this->size(1) - 1) {
                            east_room_value = (merged_room_sets.find(proto::value(*this)(i, j + 1, floor, l)) != merged_room_sets.end()) ?
                                                merged_room_sets[proto::value(*this)(i, j + 1, floor, l)] : proto::value(*this)(i, j + 1, floor, l);
                        }
                        south_room_value = (merged_room_sets.find(proto::value(*this)(i + 1, j, floor, l)) != merged_room_sets.end()) ?
                                            merged_room_sets[proto::value(*this)(i + 1, j, floor, l)] : proto::value(*this)(i + 1, j, floor, l);
                        if(l == this->size(3) - 1) {}
                        else if((rand() < east_wall_threshold*RAND_MAX) && (room_value != east_room_value)) {
                            passages.push_back(std::make_tuple(i, j, floor, l, i, j, floor, l + 1));
                            merged_room_sets[east_room_value] = room_value;
                            proto::value(*this)(i, j, floor, l + 1) = room_value;
                        }
                        if(j == this->size(1) - 1) {}
                        else if((rand() < east_wall_threshold*RAND_MAX) && (room_value != east_room_value)) {
                            passages.push_back(std::make_tuple(i, j, floor, l, i, j + 1, floor, l));
                            merged_room_sets[east_room_value] = room_value;
                            proto::value(*this)(i, j + 1, floor, l) = room_value;
                        }
                        if(((rand() < south_wall_threshold*RAND_MAX) || (can_go_south.find(room_value) == can_go_south.end()))
                                && (room_value != south_room_value)) {
                            passages.push_back(std::make_tuple(i, j, floor, l, i + 1, j, floor, l));
                            merged_room_sets[south_room_value] = room_value;
                            proto::value(*this)(i + 1, j, floor, l) = room_value;
                            can_go_south.insert(room_value);
                        }
                }
                can_go_south.clear();
            }
            int row = this->size(0) - 1;
            for(size_t j = 0; j < this->size(1) - 1; j++)
                for(size_t l = 0; l < this->size(3); l++) {
                    room_value = proto::value(*this)(row, j, floor, l);
                    if(merged_room_sets.find(room_value) != merged_room_sets.end()) {
                        room_value = merged_room_sets[room_value];
                    } else {
                        merged_room_sets[room_value] = room_value;
                    }
                    if(l < this->size(3) - 1) {
                        east_room_value = (merged_room_sets.find(proto::value(*this)(row, j, floor, l + 1)) != merged_room_sets.end()) ?
                        merged_room_sets[proto::value(*this)(row, j, floor, l + 1)] : proto::value(*this)(row, j, floor, l + 1);
                    }
                    east_room_value = (merged_room_sets.find(proto::value(*this)(row, j + 1, floor, l)) != merged_room_sets.end()) ?
                                        merged_room_sets[proto::value(*this)(row, j + 1, floor, l)] : proto::value(*this)(row, j + 1, floor, l);
                    if(room_value != east_room_value) {
                        passages.push_back(std::make_tuple(row, j, floor, l, row, j + 1, floor, l));
                        merged_room_sets[east_room_value] = room_value;
                        proto::value(*this)(row, j + 1, floor, l) = room_value;
                    }
            }
        }
        calc_cell_values();
    }
    base_type& get_maze() {
        return *this;
    }
};
