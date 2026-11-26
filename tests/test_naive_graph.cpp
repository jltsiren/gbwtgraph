#include <gtest/gtest.h>

#include <gbwtgraph/naive_graph.h>

#include "shared.h"

#include <cctype>
#include <limits>
#include <map>
#include <set>
#include <tuple>

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

  // Nodes created with `create_node()` are not stored in superstrings.
  constexpr static size_t NO_SUPERSTRING = std::numeric_limits<size_t>::max();

  // Where the sequence of a node is stored in a graph built with `build_dedup_graph()`.
  struct DedupSource
  {
    size_t superstring;
    size_t offset;
  };

  std::vector<std::string> correct_superstrings { "GATGGG", "TACAA" };
  std::map<nid_t, DedupSource> correct_dedup_sources;

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

    this->correct_dedup_sources =
    {
      { nid_t(1), { 0, 0 } },
      { nid_t(2), { 0, 1 } },
      { nid_t(3), { NO_SUPERSTRING, 0 } },
      { nid_t(4), { 0, 3 } },
      { nid_t(5), { 1, 0 } },
      { nid_t(6), { NO_SUPERSTRING, 0 } },
      { nid_t(7), { 1, 2 } },
      { nid_t(8), { 1, 1 } },
      { nid_t(9), { 1, 1 } }, // Same substring as node 8.
    };
  }

  static std::vector<char> to_vector(const std::string& str, bool lower_case)
  {
    std::vector<char> result(str.begin(), str.end());
    if(lower_case)
    {
      for(char& c : result) { c = std::tolower(static_cast<unsigned char>(c)); }
    }
    return result;
  }

  static std::string to_string(const std::string& str, bool lower_case)
  {
    std::vector<char> result = to_vector(str, lower_case);
    return std::string(result.begin(), result.end());
  }

  // Builds the same graph as `build_naive_graph(false, lower_case)`, but with
  // most node sequences stored in superstrings according to `correct_dedup_sources`.
  NaiveGraph build_dedup_graph(bool lower_case) const
  {
    NaiveGraph graph;

    graph.add_superstring(to_vector(this->correct_superstrings[0], lower_case));
    graph.create_dedup_node(1, 0, 0, 1);
    graph.create_dedup_node(2, 0, 1, 1);
    graph.create_node(3, to_string("T", lower_case));
    graph.create_dedup_node(4, 0, 3, 3);

    // Add the second superstring after some deduplicated nodes have already been created.
    graph.add_superstring(to_vector(this->correct_superstrings[1], lower_case));
    graph.create_dedup_node(5, 1, 0, 1);
    graph.create_node(6, to_string("A", lower_case));
    graph.create_dedup_node(7, 1, 2, 1);
    graph.create_dedup_node(8, 1, 1, 1);
    graph.create_dedup_node(9, 1, 1, 1);

    for(const gbwt_edge& edge : this->correct_edges)
    {
      graph.create_edge(edge.first, edge.second);
    }
    graph.remove_duplicate_edges();

    return graph;
  }

  void check_sequence_views(const NaiveGraph& graph, const std::vector<node_type>& nodes) const
  {
    for(const node_type& node : nodes)
    {
      std::string_view view = graph.get_sequence_view(node.first);
      EXPECT_EQ(std::string(view), node.second) << "Invalid sequence view for node " << node.first;
    }
  }

  // Checks the sequence sources for all nodes, assuming that they were stored
  // in superstrings according to `correct_dedup_sources` if `dedup` is set.
  void check_sequence_sources(const NaiveGraph& graph, const std::vector<node_type>& nodes, bool dedup) const
  {
    const std::vector<std::vector<char>>& superstrings = graph.get_superstrings();
    ASSERT_EQ(superstrings.size(), graph.get_num_superstrings()) << "Inconsistent number of superstrings";
    for(const node_type& node : nodes)
    {
      NaiveGraph::StringSource source = graph.get_sequence_source(node.first);
      EXPECT_EQ(std::string(source.view), node.second) << "Invalid sequence source view for node " << node.first;
      EXPECT_EQ(source.length(), node.second.length()) << "Invalid sequence source length for node " << node.first;
      EXPECT_TRUE(source.valid(superstrings)) << "Invalid sequence source for node " << node.first;

      auto iter = this->correct_dedup_sources.find(node.first);
      if(dedup && iter != this->correct_dedup_sources.end() && iter->second.superstring != NO_SUPERSTRING)
      {
        EXPECT_EQ(source.index, iter->second.superstring) << "Invalid superstring for node " << node.first;
        EXPECT_FALSE(source.linearized(superstrings)) << "Deduplicated node " << node.first << " would be linearized";
        EXPECT_EQ(source.offset(superstrings), iter->second.offset) << "Invalid superstring offset for node " << node.first;
      }
      else
      {
        EXPECT_EQ(source.index, graph.get_num_superstrings()) << "Node " << node.first << " is not in the concatenated sequences";
        EXPECT_TRUE(source.linearized(superstrings)) << "Node " << node.first << " would not be linearized";
      }
    }

    // Missing node.
    {
      nid_t missing = graph.max_node_id() + 1;
      NaiveGraph::StringSource source = graph.get_sequence_source(missing);
      EXPECT_TRUE(source.view.empty()) << "Non-empty sequence source for a missing node";
      EXPECT_EQ(source.index, graph.get_num_superstrings()) << "Invalid superstring index for a missing node";
      EXPECT_TRUE(source.valid(superstrings)) << "Invalid sequence source for a missing node";
      EXPECT_TRUE(source.linearized(superstrings)) << "A missing node would not be linearized";
    }
  }

  // Checks the superstrings in a graph built with `build_dedup_graph()`.
  // Superstrings should be normalized to upper case, even outside node sequences.
  void check_superstrings(const NaiveGraph& graph) const
  {
    ASSERT_EQ(graph.get_num_superstrings(), this->correct_superstrings.size()) << "Invalid number of superstrings";
    for(size_t i = 0; i < graph.get_num_superstrings(); i++)
    {
      EXPECT_EQ(std::string(as_view(graph.get_superstring(i))), this->correct_superstrings[i]) << "Invalid superstring " << i;
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
  this->check_sequence_sources(graph, nodes, false);
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
  this->check_sequence_sources(graph, this->correct_nodes, false);
  this->check_translation(graph, false);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });
}

TEST_F(NaiveGraphTest, AddNodesLowerCase)
{
  NaiveGraph graph = build_naive_graph(false, true);
  handle_graph_handles(graph, this->correct_nodes, false);
  handle_graph_nodes(graph, this->correct_nodes);
  handle_graph_edges(graph, this->correct_nodes, this->correct_edges);
  this->check_sequence_views(graph, this->correct_nodes);
  this->check_sequence_sources(graph, this->correct_nodes, false);
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
  this->check_sequence_sources(graph, this->correct_nodes, false);
  this->check_translation(graph, true);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });

  // Check an inverse translation without the names for s4 and s5.
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t> nodes) -> bool
  {
    return (nodes.first < 4 || nodes.first > 6);
  });
}

TEST_F(NaiveGraphTest, TranslateSegmentsLowerCase)
{
  NaiveGraph graph = build_naive_graph(true, true);
  handle_graph_handles(graph, this->correct_nodes, false);
  handle_graph_nodes(graph, this->correct_nodes);
  handle_graph_edges(graph, this->correct_nodes, this->correct_edges);
  this->check_sequence_views(graph, this->correct_nodes);
  this->check_sequence_sources(graph, this->correct_nodes, false);
  this->check_translation(graph, true);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });

  // Check an inverse translation without the names for s4 and s5.
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t> nodes) -> bool
  {
    return (nodes.first < 4 || nodes.first > 6);
  });
}

//------------------------------------------------------------------------------

TEST_F(NaiveGraphTest, DedupNodes)
{
  NaiveGraph graph = this->build_dedup_graph(false);
  handle_graph_handles(graph, this->correct_nodes, false);
  handle_graph_nodes(graph, this->correct_nodes);
  handle_graph_edges(graph, this->correct_nodes, this->correct_edges);
  this->check_sequence_views(graph, this->correct_nodes);
  this->check_sequence_sources(graph, this->correct_nodes, true);
  this->check_superstrings(graph);
  this->check_translation(graph, false);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });
}

TEST_F(NaiveGraphTest, DedupNodesLowerCase)
{
  NaiveGraph graph = this->build_dedup_graph(true);
  handle_graph_handles(graph, this->correct_nodes, false);
  handle_graph_nodes(graph, this->correct_nodes);
  handle_graph_edges(graph, this->correct_nodes, this->correct_edges);
  this->check_sequence_views(graph, this->correct_nodes);
  this->check_sequence_sources(graph, this->correct_nodes, true);
  this->check_superstrings(graph);
  this->check_translation(graph, false);
  this->check_inverse_translation(graph, [](std::pair<nid_t, nid_t>) -> bool { return true; });
}

TEST_F(NaiveGraphTest, OnlyDedupNodes)
{
  std::string superstring = "GATTACA";
  std::vector<node_type> nodes =
  {
    { nid_t(2), "GAT" },
    { nid_t(3), "TTA" }, // Overlaps node 2.
    { nid_t(5), "GATTACA" }, // The entire superstring.
  };

  NaiveGraph graph;
  graph.add_superstring(std::vector<char>(superstring.begin(), superstring.end()));
  graph.create_dedup_node(2, 0, 0, 3);
  graph.create_dedup_node(3, 0, 2, 3);
  graph.create_dedup_node(5, 0, 0, 7);
  ASSERT_EQ(graph.get_node_count(), nodes.size()) << "Invalid number of nodes";
  EXPECT_EQ(graph.min_node_id(), 2) << "Invalid minimum node id";
  EXPECT_EQ(graph.max_node_id(), 5) << "Invalid maximum node id";
  handle_graph_nodes(graph, nodes);
  this->check_sequence_views(graph, nodes);

  const std::vector<std::vector<char>>& superstrings = graph.get_superstrings();
  std::vector<size_t> offsets = { 0, 2, 0 };
  for(size_t i = 0; i < nodes.size(); i++)
  {
    NaiveGraph::StringSource source = graph.get_sequence_source(nodes[i].first);
    EXPECT_EQ(source.index, size_t(0)) << "Invalid superstring for node " << nodes[i].first;
    EXPECT_EQ(source.offset(superstrings), offsets[i]) << "Invalid superstring offset for node " << nodes[i].first;
    EXPECT_FALSE(source.linearized(superstrings)) << "Deduplicated node " << nodes[i].first << " would be linearized";
  }
}

TEST_F(NaiveGraphTest, InvalidDedupNodes)
{
  NaiveGraph graph;

  // No superstrings yet.
  EXPECT_THROW(graph.create_dedup_node(1, 0, 0, 1), std::runtime_error) << "Created a node in a missing superstring";
  EXPECT_EQ(graph.get_node_count(), size_t(0)) << "A node was created in a missing superstring";

  std::string superstring = "GATTACA";
  graph.add_superstring(std::vector<char>(superstring.begin(), superstring.end()));
  graph.add_superstring(std::vector<char>()); // An empty superstring.
  graph.create_node(1, "C");
  graph.create_dedup_node(2, 0, 0, 3);
  std::vector<node_type> nodes = { { nid_t(1), "C" }, { nid_t(2), "GAT" } };

  // (node id, superstring index, offset, length, test case)
  std::vector<std::tuple<nid_t, size_t, size_t, size_t, std::string>> invalid =
  {
    { nid_t(3), 0, 0, 0, "empty sequence" },
    { nid_t(3), 2, 0, 1, "missing superstring" },
    { nid_t(3), 1, 0, 1, "empty superstring" },
    { nid_t(3), 0, 7, 1, "offset at the end" },
    { nid_t(3), 0, 8, 1, "offset past the end" },
    { nid_t(3), 0, 5, 3, "sequence past the end" },
    { nid_t(3), 0, 0, 8, "sequence longer than superstring" },
    { nid_t(3), 0, 1, std::numeric_limits<size_t>::max(), "maximum length" },
    { nid_t(1), 0, 1, 1, "duplicate of a concatenated node" },
    { nid_t(2), 0, 4, 3, "duplicate of a deduplicated node" },
  };
  for(const auto& test : invalid)
  {
    const std::string& test_case = std::get<4>(test);
    EXPECT_THROW(graph.create_dedup_node(std::get<0>(test), std::get<1>(test), std::get<2>(test), std::get<3>(test)), std::runtime_error) << "No exception for " << test_case;
    EXPECT_EQ(graph.get_node_count(), nodes.size()) << "Invalid number of nodes after " << test_case;
    EXPECT_FALSE(graph.has_node(3)) << "Node 3 was created after " << test_case;
    this->check_sequence_views(graph, nodes);
  }

  // A concatenated node cannot replace a deduplicated node.
  EXPECT_THROW(graph.create_node(2, "A"), std::runtime_error) << "Replaced a deduplicated node with a concatenated node";
  EXPECT_EQ(graph.get_node_count(), nodes.size()) << "Invalid number of nodes after replacing a deduplicated node";
  this->check_sequence_views(graph, nodes);

  // The superstrings are unchanged.
  ASSERT_EQ(graph.get_num_superstrings(), size_t(2)) << "Invalid number of superstrings";
  EXPECT_EQ(std::string(as_view(graph.get_superstring(0))), superstring) << "The first superstring was modified";
  EXPECT_TRUE(graph.get_superstring(1).empty()) << "The empty superstring was modified";
}

TEST_F(NaiveGraphTest, CopyAndMoveDedupNodes)
{
  NaiveGraph original = this->build_dedup_graph(false);

  // Sequence sources in the copy must point to the superstrings of the copy.
  NaiveGraph copy(original);
  this->check_sequence_views(copy, this->correct_nodes);
  this->check_sequence_sources(copy, this->correct_nodes, true);
  for(const node_type& node : this->correct_nodes)
  {
    EXPECT_NE(copy.get_sequence_view(node.first).data(), original.get_sequence_view(node.first).data()) << "The copy shares the sequence of node " << node.first;
  }

  NaiveGraph assigned;
  assigned = original;
  this->check_sequence_views(assigned, this->correct_nodes);
  this->check_sequence_sources(assigned, this->correct_nodes, true);

  NaiveGraph moved(std::move(copy));
  this->check_sequence_views(moved, this->correct_nodes);
  this->check_sequence_sources(moved, this->correct_nodes, true);

  NaiveGraph move_assigned;
  move_assigned = std::move(assigned);
  this->check_sequence_views(move_assigned, this->correct_nodes);
  this->check_sequence_sources(move_assigned, this->correct_nodes, true);
}

//------------------------------------------------------------------------------

/*
  Build a `DedupStringArray` from a `NaiveGraph` in the same way as `GBWTGraph`
  would build its sequence array: one string per node id in a contiguous range
  (with empty strings for missing ids), optionally with reverse complements.
*/

class NaiveGraphToDedupTest : public NaiveGraphTest
{
public:
  // Node ids from 0 to max_node_id(), with id 0 missing.
  DedupStringArray build_array(const NaiveGraph& graph, bool reverse_complements) const
  {
    size_t n = graph.max_node_id() + 1;
    return DedupStringArray
    (
      n, reverse_complements, graph.get_superstrings(),
      [&](size_t i) -> NaiveGraph::StringSource { return graph.get_sequence_source(nid_t(i)); }
    );
  }

  void check_array(const NaiveGraph& graph, const DedupStringArray& array, bool reverse_complements, bool dedup, const std::string& test_case) const
  {
    size_t n = graph.max_node_id() + 1;
    size_t mult = (reverse_complements ? 2 : 1);
    ASSERT_EQ(array.size(), n * mult) << test_case << ": Invalid number of strings";
    EXPECT_EQ(array.contains_reverse_complements(), reverse_complements) << test_case << ": Invalid reverse complement flag";

    std::vector<std::string> truth(n);
    for(const node_type& node : this->correct_nodes) { truth[node.first] = node.second; }

    std::vector<std::string> all_strings;
    size_t linearized_length = 0, total_length = 0;
    for(size_t i = 0; i < n; i++)
    {
      bool in_superstring = false;
      auto iter = this->correct_dedup_sources.find(nid_t(i));
      if(dedup && iter != this->correct_dedup_sources.end() && iter->second.superstring != NO_SUPERSTRING)
      {
        in_superstring = true;
      }
      else { linearized_length += truth[i].length(); }
      total_length += truth[i].length();

      EXPECT_EQ(array.str(i * mult), truth[i]) << test_case << ": Invalid sequence for node " << i;
      EXPECT_EQ(array.linearized(i * mult), !in_superstring) << test_case << ": Invalid linearization status for node " << i;
      all_strings.push_back(truth[i]);
      if(reverse_complements)
      {
        std::string rc = reverse_complement(truth[i]);
        EXPECT_EQ(array.str(i * mult + 1), rc) << test_case << ": Invalid reverse complement for node " << i;
        EXPECT_EQ(array.linearized(i * mult + 1), !in_superstring) << test_case << ": Invalid linearization status for the reverse complement of node " << i;
        all_strings.push_back(rc);
      }
    }

    size_t superstring_length = 0;
    for(const std::vector<char>& superstring : graph.get_superstrings()) { superstring_length += superstring.size(); }
    EXPECT_EQ(array.length(), total_length * mult) << test_case << ": Invalid total length";
    EXPECT_EQ(array.linearized_length, linearized_length * mult) << test_case << ": Invalid linearized length";
    EXPECT_EQ(array.data_array_length, (linearized_length + superstring_length) * mult) << test_case << ": Invalid data array length";

    // Compare to a deduplicated array built from a naive string array.
    DedupStringArray from_naive(gbwt::StringArray(all_strings), reverse_complements);
    EXPECT_TRUE(array.equivalent(from_naive)) << test_case << ": Not equivalent to an array built from a naive string array";

    // Linearization should not change the strings.
    DedupStringArray linearized = array;
    linearized.linearize();
    EXPECT_TRUE(linearized.linear()) << test_case << ": Linearization failed";
    EXPECT_TRUE(linearized.equivalent(array)) << test_case << ": Linearization changed the strings";
  }
};

TEST_F(NaiveGraphToDedupTest, WithoutSuperstrings)
{
  for(bool rc : { false, true })
  {
    for(bool translation : { false, true })
    {
      NaiveGraph graph = build_naive_graph(translation);
      DedupStringArray array = this->build_array(graph, rc);
      std::string test_case = std::string("Translation ") + (translation ? "on" : "off") + ", reverse complements " + (rc ? "on" : "off");
      this->check_array(graph, array, rc, false, test_case);
      EXPECT_TRUE(array.linear()) << test_case << ": The array is not linear";
    }
  }
}

TEST_F(NaiveGraphToDedupTest, WithSuperstrings)
{
  for(bool rc : { false, true })
  {
    for(bool lower_case : { false, true })
    {
      NaiveGraph graph = this->build_dedup_graph(lower_case);
      DedupStringArray array = this->build_array(graph, rc);
      std::string test_case = std::string("Lower case ") + (lower_case ? "on" : "off") + ", reverse complements " + (rc ? "on" : "off");
      this->check_array(graph, array, rc, true, test_case);
      EXPECT_FALSE(array.linear()) << test_case << ": The array is linear";
    }
  }
}

//------------------------------------------------------------------------------

} // namespace
