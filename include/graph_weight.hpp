#ifndef GRAPH_WEIGHT_HPP
#define GRAPH_WEIGHT_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/visitors.hpp>


template<typename T>
using graph_type = boost::adjacency_list<boost::vecS,
                                         boost::vecS,
                                         boost::bidirectionalS, T>;
static constexpr size_t vertex_count = 27;

namespace _spatial {

template<typename T>
struct Matrix3D;
}

#include <iosfwd>
namespace detail {
template<class GraphT>
class vertex_write
{
    const GraphT& graph;
public:
    explicit vertex_write(const GraphT& graph) : graph(graph) {}
    template<typename VertexDiscr>
    void operator ()(std::ostream& out, const VertexDiscr& dd) const {
        out << " [KNOT = \"" << boost::get(boost::vertex_bundle, graph)[dd] << "\"]";
    }
};
}
namespace _graph {

template<typename T>
struct local_baricentric_star
{
private:
    using descriptor_t = typename boost::graph_traits<graph_type<T>>::vertex_descriptor;
    using vertex_it_t  = typename boost::graph_traits<graph_type<T>>::vertex_iterator;
    using edge_it_t    = typename boost::graph_traits<graph_type<T>>::edge_iterator;
    using edge_desc_t  = typename boost::graph_traits<graph_type<T>>::edge_descriptor;

    _spatial::Matrix3D<T>& node;
    size_t X;
    size_t Y;
    size_t Z;
    graph_type<T> graph;
    descriptor_t source;

    descriptor_t gorizont_front;
    descriptor_t gorizont_back;
    descriptor_t gorizont_right;
    descriptor_t gorizont_left;

    descriptor_t gorizont_front_right;
    descriptor_t gorizont_back_right;
    descriptor_t gorizont_front_left;
    descriptor_t gorizont_back_left;
public:
    local_baricentric_star(_spatial::Matrix3D<T>& node, size_t X, size_t Y, size_t Z) : node(node), X(X), Y(Y), Z(Z) {
        source         = boost::add_vertex(node(X, Y, Z), graph);

        gorizont_front = X < node.size(0) - 1 ? boost::add_vertex(node(X + 1, Y, Z), graph) :
                                                boost::add_vertex(node(node.size(0) - 1, Y, Z), graph);
        gorizont_back  = X > 0 ? boost::add_vertex(node(X - 1, Y, Z), graph) :
                                 boost::add_vertex(node(0, Y, Z), graph);
        gorizont_right = Y < node.size(1) - 1 ? boost::add_vertex(node(X, Y + 1, Z), graph) :
                                                boost::add_vertex(node(X, node.size(1) - 1, Z), graph);
        gorizont_left  = Y > 0 ? boost::add_vertex(node(X, Y - 1, Z), graph) :
                                 boost::add_vertex(node(X, 0, Z), graph);

        gorizont_front_right = X < node.size(0) - 1 && Y < node.size(1) - 1 ?
                    boost::add_vertex(node(X + 1, Y + 1, Z), graph) : boost::add_vertex(node(node.size(0) - 1, node.size(1) - 1, Z), graph);
        gorizont_back_right  = X > 0 && Y < node.size(1) - 1 ?
                    boost::add_vertex(node(X - 1, Y + 1, Z), graph) : boost::add_vertex(node(0, node.size(1) - 1, Z), graph);
        gorizont_front_left  = X < node.size(0) - 1 && Y > 0 ?
                    boost::add_vertex(node(X + 1, Y - 1, Z), graph) : boost::add_vertex(node(node.size(0) - 1, 0, Z), graph);
        gorizont_back_left   = X > 0 && Y > 0 ?
                    boost::add_vertex(node(X - 1, Y - 1, Z), graph) : boost::add_vertex(node(0, 0, Z), graph);

        boost::add_edge(source, gorizont_front, graph);
        boost::add_edge(source, gorizont_back, graph);
        boost::add_edge(source, gorizont_right, graph);
        boost::add_edge(source, gorizont_left, graph);

        boost::add_edge(source, gorizont_front_right, graph);
        boost::add_edge(source, gorizont_back_right, graph);
        boost::add_edge(source, gorizont_front_left, graph);
        boost::add_edge(source, gorizont_back_left, graph);

//        std::array<int, 9> distance{{0}};

//        boost::breadth_first_search(graph, gorizont_back_left,
//                                    boost::visitor(
//                                        boost::make_bfs_visitor(
//                                            boost::record_distances(distance.begin(), boost::on_tree_edge{}))));

//        std::copy(distance.begin(), distance.end(), std::ostream_iterator<int>(std::cout, "\n"));
    }
    void find_print(size_t X, size_t Y, size_t Z) const {
        vertex_it_t it, end;
        boost::tie(it, end) = boost::vertices(graph);
        for(; it != end; ++it) {
            const descriptor_t desk = *it;
            const T& vertex = boost::get(boost::vertex_bundle, graph)[desk];
            if(vertex == node(X, Y, Z)) {
                break;
            }
        }
        BOOST_ASSERT_MSG(it != end, "Граф не рабочий.");
        std::cout << "Вершина с координатами: {"
                  << X << "; " << Y
                  << "; " << Z << "} = " << node(X, Y, Z) << std::endl;
    }
    friend std::ostream& operator << (std::ostream& out, const local_baricentric_star<T>& g) {
        detail::vertex_write<graph_type<T>> vw(g.graph);
        boost::write_graphviz(out, g.graph, vw);
        return out;
    }

};
}

#endif // GRAPH_WEIGHT_HPP
