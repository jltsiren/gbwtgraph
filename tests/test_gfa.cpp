#include <gtest/gtest.h>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class GFAConstruction : public ::testing::Test
{
public:
  gbwt::GBWT index, gfa_index;
  SequenceSource source;
  GBWTGraph graph, gfa_graph;

  GFAConstruction()
  {
  }

  void SetUp() override
  {
    this->index = build_gbwt_index();
    build_source(this->source);
    this->graph = GBWTGraph(this->index, this->source);

    auto gfa_parse = gfa_to_gbwt("example.gfa");
    this->gfa_index = *(gfa_parse.first);
    this->gfa_graph = GBWTGraph(this->gfa_index, *(gfa_parse.second));
  }
};

TEST_F(GFAConstruction, GBWTComparison)
{
  ASSERT_EQ(this->index.size(), this->gfa_index.size()) << "Different data size";
  ASSERT_EQ(this->index.sequences(), this->gfa_index.sequences()) << "Number of sequences";
  ASSERT_EQ(this->index.sigma(), this->gfa_index.sigma()) << "Different alphabet size";
  ASSERT_EQ(this->index.effective(), this->gfa_index.effective()) << "Different effective alphabet size";
  ASSERT_EQ(this->index.samples(), this->gfa_index.samples()) << "Different number of samples";
}

TEST_F(GFAConstruction, GraphComparison)
{
  ASSERT_EQ(this->graph.header, this->gfa_graph.header) << "Graph headers are not identical";
  this->graph.for_each_handle([&](const handle_t& handle)
  {
    nid_t node_id = this->graph.get_id(handle);
    ASSERT_TRUE(this->gfa_graph.has_node(this->graph.get_id(handle))) << "GFA graph does not have node " << node_id;
    EXPECT_EQ(this->graph.get_sequence(this->graph.get_handle(node_id, false)),
              this->gfa_graph.get_sequence(this->gfa_graph.get_handle(node_id, false))) << "Wrong forward sequence for node " << node_id;
    EXPECT_EQ(this->graph.get_sequence(this->graph.get_handle(node_id, true)),
              this->gfa_graph.get_sequence(this->gfa_graph.get_handle(node_id, true))) << "Wrong reverse sequence for node " << node_id;
  });
  this->graph.for_each_edge([&](const edge_t& edge)
  {
    nid_t id_from = this->graph.get_id(edge.first), id_to = this->graph.get_id(edge.second);
    bool rev_from = this->graph.get_is_reverse(edge.first), rev_to = this->graph.get_is_reverse(edge.second);
    handle_t from = this->gfa_graph.get_handle(id_from, rev_from);
    handle_t to = this->gfa_graph.get_handle(id_to, rev_to);
    EXPECT_TRUE(this->gfa_graph.has_edge(from, to)) << "GFA graph does not have the edge ((" << id_from << ", " << rev_from <<"), (" << id_to << ", " << rev_to << "))";
  });
}

//------------------------------------------------------------------------------

} // namespace
