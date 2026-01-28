#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include <omp.h>

#include <gbwtgraph/cached_gbwtgraph.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class GraphOperations : public ::testing::Test
{
public:
  gbwt::GBWT index;
  NaiveGraph source;
  GBWTGraph graph;
  CachedGBWTGraph cached_graph;
  std::vector<node_type> correct_nodes;
  std::set<gbwt_edge> correct_edges;

  GraphOperations()
  {
  }

  void SetUp() override
  {
    this->index = build_gbwt_index();
    this->source = build_naive_graph(false);
    this->graph = GBWTGraph(this->index, this->source);
    this->cached_graph = CachedGBWTGraph(this->graph);

    this->correct_nodes =
    {
      { nid_t(1), "G" },
      { nid_t(2), "A" },
      { nid_t(4), "GGG" },
      { nid_t(5), "T" },
      { nid_t(6), "A" },
      { nid_t(7), "C" },
      { nid_t(8), "A" },
      { nid_t(9), "A" }
    };

    this->correct_edges =
    {
      { gbwt::Node::encode(1, false), gbwt::Node::encode(2, false) },
      { gbwt::Node::encode(1, false), gbwt::Node::encode(4, false) },
      { gbwt::Node::encode(2, false), gbwt::Node::encode(4, false) },
      { gbwt::Node::encode(4, false), gbwt::Node::encode(5, false) },
      { gbwt::Node::encode(5, false), gbwt::Node::encode(6, false) },
      { gbwt::Node::encode(6, false), gbwt::Node::encode(7, false) },
      { gbwt::Node::encode(6, false), gbwt::Node::encode(8, false) },
      { gbwt::Node::encode(7, false), gbwt::Node::encode(9, false) },
      { gbwt::Node::encode(8, false), gbwt::Node::encode(9, false) }
    };
  }
};

TEST_F(GraphOperations, EmptyGraph)
{
  gbwt::GBWT empty_index;
  NaiveGraph empty_source;
  GBWTGraph empty_graph(empty_index, empty_source);
  CachedGBWTGraph empty_cache(empty_graph);

  std::vector<node_type> nodes;
  std::set<gbwt_edge> edges;
  handle_graph_handles(empty_cache, nodes, true);
  handle_graph_nodes(empty_cache, nodes);
  handle_graph_edges(empty_cache, nodes, edges);
}

// HandleGraph interface for the default test graph.
TEST_F(GraphOperations, HandleGraphTests)
{
  handle_graph_handles(this->cached_graph, this->correct_nodes, true);
  handle_graph_nodes(this->cached_graph, this->correct_nodes);
  handle_graph_edges(this->cached_graph, this->correct_nodes, this->correct_edges);
}

//------------------------------------------------------------------------------

} // namespace
