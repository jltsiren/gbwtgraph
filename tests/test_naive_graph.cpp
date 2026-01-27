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
  using gbwt_edge = std::pair<gbwt::node_type, gbwt::node_type>;
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

  void check_handles(const NaiveGraph& graph, bool expect_nodes) const
  {
    std::vector<node_type> empty;
    const std::vector<node_type>& correct_nodes = (expect_nodes ? this->correct_nodes : empty);

    size_t count = 0;
    std::set<nid_t> visited_ids;
    graph.for_each_handle([&](const handle_t& handle)
    {
      nid_t id = graph.get_id(handle);
      EXPECT_FALSE(graph.get_is_reverse(handle)) << "Handle at index " << count << " is in reverse orientation for node " << id;
      handle_t flipped = graph.flip(handle);
      EXPECT_EQ(graph.get_id(flipped), id) << "Flipped handle at index " << count << " has incorrect id for node " << id;
      EXPECT_TRUE(graph.get_is_reverse(flipped)) << "Flipped handle at index " << count << " is in forward orientation for node " << id;
      count++; visited_ids.insert(id);
    });
    EXPECT_EQ(count, correct_nodes.size()) << "Incorrect number of handles visited";

    for(const node_type& node : correct_nodes)
    {
      EXPECT_TRUE(visited_ids.find(node.first) != visited_ids.end()) << "Node id " << node.first << " was not visited";
    }
  }

  void check_nodes(const NaiveGraph& graph, bool expect_nodes) const
  {
    std::vector<node_type> empty;
    const std::vector<node_type>& correct_nodes = (expect_nodes ? this->correct_nodes : empty);
    nid_t min_id = (expect_nodes ? correct_nodes.front().first : nid_t(std::numeric_limits<nid_t>::max()));
    nid_t max_id = (expect_nodes ? correct_nodes.back().first : 0);

    ASSERT_EQ(graph.get_node_count(), correct_nodes.size()) << "Incorrect number of nodes";
    EXPECT_EQ(graph.min_node_id(), min_id) << "Incorrect minimum node id";
    EXPECT_EQ(graph.max_node_id(), max_id) << "Incorrect maximum node id";

    for(const node_type& node : correct_nodes)
    {
      ASSERT_TRUE(graph.has_node(node.first)) << "Node id " << node.first << " is missing from the sequence source";
      handle_t handle = graph.get_handle(node.first, false);
      EXPECT_EQ(graph.get_length(handle), node.second.length()) << "Invalid sequence length for node " << node.first;
      EXPECT_EQ(graph.get_sequence(handle), node.second) << "Invalid sequence for node " << node.first;
      std::string reverse = reverse_complement(node.second);
      handle_t rev_handle = graph.get_handle(node.first, true);
      EXPECT_EQ(graph.get_length(rev_handle), node.second.length()) << "Invalid sequence length for reverse node " << node.first;
      EXPECT_EQ(graph.get_sequence(rev_handle), reverse) << "Invalid sequence for reverse node " << node.first;

      std::string_view view = graph.get_sequence_view(node.first);
      EXPECT_EQ(std::string(view), node.second) << "Invalid sequence view for node " << node.first;

      for(size_t i = 0; i < node.second.length(); i++)
      {
        EXPECT_EQ(graph.get_base(handle, i), node.second[i]) << "Invalid base at index " << i << " for node " << node.first;
        EXPECT_EQ(graph.get_base(rev_handle, i), reverse[i]) << "Invalid base at index " << i << " for reverse node " << node.first;
        EXPECT_EQ(graph.get_subsequence(handle, i, 2), node.second.substr(i, 2)) << "Invalid subsequence at index " << i << " for node " << node.first;
        EXPECT_EQ(graph.get_subsequence(rev_handle, i, 2), reverse.substr(i, 2)) << "Invalid subsequence at index " << i << " for reverse node " << node.first;
      }
    }

    nid_t missing_id = graph.max_node_id() + 100;
    EXPECT_FALSE(graph.has_node(missing_id)) << "Missing node id " << missing_id << " is reported present";
    std::string missing_id_str = std::to_string(missing_id);
    EXPECT_FALSE(graph.has_segment(missing_id_str)) << "Missing segment " << missing_id_str << " is reported present";
  }

  void check_edges(const NaiveGraph& graph, bool expect_edges) const
  {
    std::set<gbwt_edge> empty;
    const std::set<gbwt_edge>& correct_edges = (expect_edges ? this->correct_edges : empty);
    std::set<gbwt_edge> reverse_edges;
    for(const auto& edge : correct_edges)
    {
      reverse_edges.insert({ gbwt::Node::reverse(edge.second), gbwt::Node::reverse(edge.first) });
    }
    ASSERT_EQ(graph.get_edge_count(), correct_edges.size()) << "Incorrect number of edges";

    // follow_edges and get_degree
    std::set<gbwt_edge> fw_succ, fw_pred, rev_succ, rev_pred;
    for(const auto& node : correct_nodes)
    {
      handle_t forward_handle = NaiveGraph::node_to_handle(gbwt::Node::encode(node.first, false));
      handle_t reverse_handle = NaiveGraph::node_to_handle(gbwt::Node::encode(node.first, true));
      size_t fw_out = 0, fw_in = 0, rev_out = 0, rev_in = 0;
      graph.follow_edges(forward_handle, false, [&](const handle_t& handle)
      {
        fw_succ.insert(gbwt_edge(GBWTGraph::handle_to_node(forward_handle), GBWTGraph::handle_to_node(handle)));
        fw_out++;
      });
      graph.follow_edges(forward_handle, true, [&](const handle_t& handle)
      {
        fw_pred.insert(gbwt_edge(GBWTGraph::handle_to_node(handle), GBWTGraph::handle_to_node(forward_handle)));
        fw_in++;
      });
      graph.follow_edges(reverse_handle, false, [&](const handle_t& handle)
      {
        rev_succ.insert(gbwt_edge(GBWTGraph::handle_to_node(reverse_handle), GBWTGraph::handle_to_node(handle)));
        rev_out++;
      });
      graph.follow_edges(reverse_handle, true, [&](const handle_t& handle)
      {
        rev_pred.insert(gbwt_edge(GBWTGraph::handle_to_node(handle), GBWTGraph::handle_to_node(reverse_handle)));
        rev_in++;
      });
      EXPECT_EQ(graph.get_degree(forward_handle, false), fw_out) << "Wrong outdegree for forward handle " << node.first;
      EXPECT_EQ(graph.get_degree(forward_handle, true), fw_in) << "Wrong indegree for forward handle " << node.first;
      EXPECT_EQ(graph.get_degree(reverse_handle, false), rev_out) << "Wrong outdegree for reverse handle " << node.first;
      EXPECT_EQ(graph.get_degree(reverse_handle, true), rev_in) << "Wrong indegree for reverse handle " << node.first;
    }
    EXPECT_EQ(fw_succ, correct_edges) << "Wrong forward successors";
    EXPECT_EQ(fw_pred, correct_edges) << "Wrong forward predecessors";
    EXPECT_EQ(rev_succ, reverse_edges) << "Wrong reverse successors";
    EXPECT_EQ(rev_pred, reverse_edges) << "Wrong reverse predecessors";

    // Presence / absence of edges
    for(nid_t from = graph.min_node_id(); from <= graph.max_node_id(); from++)
    {
      for(nid_t to = graph.min_node_id(); to <= graph.max_node_id(); to++)
      {
        for(bool from_rev : { false, true })
        {
          for(bool to_rev : { false, true })
          {
            handle_t from_handle = graph.get_handle(from, from_rev);
            handle_t to_handle = graph.get_handle(to, to_rev);
            gbwt_edge edge(GBWTGraph::handle_to_node(from_handle), GBWTGraph::handle_to_node(to_handle));
            bool should_have = (correct_edges.find(edge) != correct_edges.end());
            should_have |= (reverse_edges.find(edge) != reverse_edges.end());
            EXPECT_EQ(graph.has_edge(from_handle, to_handle), should_have) <<
              "has_edge() failed with (" << from << ", " << from_rev << ") to (" << to << ", " << to_rev <<")";
          }
        }
      }
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
        EXPECT_EQ(graph.translate(translation.first), translation.second) << "Invalid translation for " << translation.first;
        // Here we assume that segment names are not numeric in the tests.
        std::string node_str = std::to_string(translation.second.first);
        EXPECT_FALSE(graph.has_segment(node_str)) << "Node " << node_str << " of segment " << translation.first << " is identified as a segment";
      }
      {
        std::string missing = "missing_segment";
        EXPECT_FALSE(graph.has_segment(missing)) << "Missing segment is reported present";
        EXPECT_EQ(graph.translate(missing), SequenceSource::invalid_translation()) << "A translation is returned for a missing segment";
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
  this->check_handles(graph, false);
  this->check_nodes(graph, false);
  this->check_edges(graph, false);
  this->check_translation(graph, false);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });
}

TEST_F(NaiveGraphTest, AddNodes)
{
  NaiveGraph graph = build_naive_graph(false);
  this->check_handles(graph, true);
  this->check_nodes(graph, true);
  this->check_edges(graph, true);
  this->check_translation(graph, false);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });
}

TEST_F(NaiveGraphTest, TranslateSegments)
{
  NaiveGraph graph = build_naive_graph(true);
  this->check_handles(graph, true);
  this->check_nodes(graph, true);
  this->check_edges(graph, true);
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
