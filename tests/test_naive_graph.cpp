#include <gtest/gtest.h>

#include <gbwtgraph/naive_graph.h>

#include "shared.h"

#include <set>

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class NaiveGraphTest : public ::testing::Test
{
public:
  std::vector<node_type> correct_nodes;
  std::vector<translation_type> correct_translation;
  std::set<gbwt_edge> correct_edges;

  void SetUp() override
  {
    this->correct_nodes =
    {
      { nid_t(1), "G" },
      { nid_t(2), "A" },
      { nid_t(3), "T" },
      { nid_t(4), "GGG" },
      { nid_t(5), "T" },
      { nid_t(6), "A" },
      { nid_t(7), "C" },
      { nid_t(8), "A" },
      { nid_t(9), "A" }
    };

    this->correct_translation =
    {
      { "s1", { 1, 2 } },
      { "s2", { 2, 3 } },
      { "s3", { 3, 4 } },
      { "s4", { 4, 6 } },
      { "s5", { 6, 7 } },
      { "s6", { 7, 8 } },
      { "s7", { 8, 9 } },
      { "s8", { 9, 10 } },
    };

    this->correct_edges =
    {
      { gbwt::Node::encode(1, false), gbwt::Node::encode(2, false) },
      { gbwt::Node::encode(2, false), gbwt::Node::encode(4, false) },
      { gbwt::Node::encode(4, false), gbwt::Node::encode(5, false) },
      { gbwt::Node::encode(5, false), gbwt::Node::encode(6, false) },
      { gbwt::Node::encode(6, false), gbwt::Node::encode(7, false) },
      { gbwt::Node::encode(6, false), gbwt::Node::encode(8, false) },
      { gbwt::Node::encode(7, false), gbwt::Node::encode(9, false) },
      { gbwt::Node::encode(8, false), gbwt::Node::encode(9, false) }
    };
  }

  void check_sequence_views(const NaiveGraph& graph, const std::vector<node_type>& nodes) const
  {
    for(const node_type& node : nodes)
    {
      std::string_view view = graph.get_sequence_view(node.first);
      EXPECT_EQ(std::string(view), node.second) << "Invalid sequence view for node " << node.first;
    }
  }

  void check_translation(const NaiveGraph& graph, bool expect_translation) const
  {
    ASSERT_EQ(graph.uses_translation(), expect_translation) << "Segment translation is not used as expected";
    if(graph.uses_translation())
    {
      ASSERT_EQ(graph.get_segment_count(), correct_translation.size()) << "Invalid number of segments";
      for(const translation_type& translation : correct_translation)
      {
        EXPECT_TRUE(graph.has_segment(translation.first)) << "Segment " << translation.first << " is missing from the sequence source";
        EXPECT_TRUE(graph.has_node_or_segment(translation.first)) << "Node/segment " << translation.first << " is missing from the sequence source";
        EXPECT_EQ(graph.translate(translation.first), translation.second) << "Invalid translation for " << translation.first;
        // Here we assume that segment names are not numeric in the tests.
        std::string node_str = std::to_string(translation.second.first);
        EXPECT_FALSE(graph.has_segment(node_str)) << "Node " << node_str << " of segment " << translation.first << " is identified as a segment";
        EXPECT_FALSE(graph.has_node_or_segment(node_str)) << "Node " << node_str << " of segment " << translation.first << " is identified as a node/segment";
      }
      {
        std::string missing = "missing_segment";
        EXPECT_FALSE(graph.has_segment(missing)) << "Missing segment is reported present";
        EXPECT_FALSE(graph.has_node_or_segment(missing)) << "Missing node/segment is reported present";
        EXPECT_EQ(graph.translate(missing), NaiveGraph::no_translation()) << "A translation is returned for a missing segment";
      }
    }
    else
    {
      for(nid_t id = graph.min_node_id(); id <= graph.max_node_id(); id++)
      {
        std::string segment = std::to_string(id);
        std::pair<nid_t, nid_t> translation(id, id + 1);
        EXPECT_EQ(graph.translate(segment), translation) << "Invalid translation for " << segment;
        EXPECT_FALSE(graph.has_segment(segment)) << "Node " << segment << " is identified as a segment";
        EXPECT_TRUE(graph.has_node_or_segment(segment)) << "Node " << segment << " is not identified as a node/segment";
      }
      {
        std::string missing = "0";
        EXPECT_EQ(graph.translate(missing), NaiveGraph::no_translation()) << "A translation is returned for a missing segment";
      }
    }
  }

  void check_inverse_translation(const NaiveGraph& graph, const std::function<bool(std::pair<nid_t, nid_t>)>& is_present) const
  {
    auto result = graph.invert_translation(is_present);
    ASSERT_EQ(result.first.size(), graph.get_segment_count()) << "Invalid number of segments in inverse translation";
    ASSERT_EQ(result.second.size(), size_t(graph.max_node_id() + 1)) << "Invalid number of nodes in inverse translation";
    ASSERT_EQ(result.first.size(), result.second.ones()) << "Inconsistent number of segments in inverse translation";
    auto iter = result.second.one_begin();
    nid_t start = iter->second;
    for(size_t i = 0; i < result.first.size(); i++)
    {
      ++iter;
      nid_t limit = iter->second;
      std::string segment = result.first.str(i);
      if(is_present(std::make_pair(start, limit)))
      {
        std::pair<nid_t, nid_t> translation = graph.translate(segment);
        EXPECT_EQ(start, translation.first) << "Invalid start for segment " << segment;
        EXPECT_EQ(limit, translation.second) << "Invalid limit for segment " << segment;
      }
      else
      {
        EXPECT_TRUE(segment.empty()) << "Got a name for an unused segment from " << start << " to " << limit;
      }
      start = limit;
    }
  }
};

TEST_F(NaiveGraphTest, EmptyGraph)
{
  NaiveGraph graph;
  std::vector<node_type> nodes;
  std::set<gbwt_edge> edges;
  handle_graph_handles(graph, nodes, false);
  handle_graph_nodes(graph, nodes);
  handle_graph_edges(graph, nodes, edges);
  this->check_sequence_views(graph, nodes);
  this->check_translation(graph, false);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });
}

TEST_F(NaiveGraphTest, AddNodes)
{
  NaiveGraph graph = build_naive_graph(false);
  handle_graph_handles(graph, this->correct_nodes, false);
  handle_graph_nodes(graph, this->correct_nodes);
  handle_graph_edges(graph, this->correct_nodes, this->correct_edges);
  this->check_sequence_views(graph, this->correct_nodes);
  this->check_translation(graph, false);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });
}

TEST_F(NaiveGraphTest, TranslateSegments)
{
  NaiveGraph graph = build_naive_graph(true);
  handle_graph_handles(graph, this->correct_nodes, false);
  handle_graph_nodes(graph, this->correct_nodes);
  handle_graph_edges(graph, this->correct_nodes, this->correct_edges);
  this->check_sequence_views(graph, this->correct_nodes);
  this->check_translation(graph, true);  
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });

  // Check an inverse translation without the names for s4 and s5.
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t> nodes) -> bool
  {
    return (nodes.first < 4 || nodes.first > 6);
  });
}

//------------------------------------------------------------------------------

} // namespace
