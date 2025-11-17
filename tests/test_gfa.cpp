#include <gtest/gtest.h>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>

#include <handlegraph/algorithms/canonical_gfa.hpp>

#include <sstream>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class GFAConstruction : public ::testing::Test
{
public:
  gbwt::GBWT index;
  SequenceSource source;
  GBWTGraph graph;

  GFAConstruction()
  {
  }

  void SetUp() override
  {
    this->index = build_gbwt_index();
    build_source(this->source);
    this->graph = GBWTGraph(this->index, this->source);
  }

  void check_gbwt(const gbwt::GBWT& gfa_index, const gbwt::GBWT* truth) const
  {
    ASSERT_EQ(gfa_index.size(), truth->size()) << "Different data size";
    ASSERT_EQ(gfa_index.sequences(), truth->sequences()) << "Number of sequences";
    ASSERT_EQ(gfa_index.sigma(), truth->sigma()) << "Different alphabet size";
    ASSERT_EQ(gfa_index.effective(), truth->effective()) << "Different effective alphabet size";
    ASSERT_EQ(gfa_index.samples(), truth->samples()) << "Different number of samples";
  }

  void check_graph(const GBWTGraph& gfa_graph, const GBWTGraph* truth) const
  {
    // Only check node counts. Headers may be different due to flags that are out of
    // scope for this test.
    ASSERT_EQ(gfa_graph.get_node_count(), truth->get_node_count()) << "Node counts are not identical";

    truth->for_each_handle([&](const handle_t& handle)
    {
      nid_t node_id = truth->get_id(handle);
      ASSERT_TRUE(gfa_graph.has_node(gfa_graph.get_id(handle))) << "GFA graph does not have node " << node_id;
      EXPECT_EQ(gfa_graph.get_sequence(gfa_graph.get_handle(node_id, false)),
                truth->get_sequence(truth->get_handle(node_id, false))) << "Wrong forward sequence for node " << node_id;
      EXPECT_EQ(gfa_graph.get_sequence(gfa_graph.get_handle(node_id, true)),
                truth->get_sequence(truth->get_handle(node_id, true))) << "Wrong reverse sequence for node " << node_id;
    });

    truth->for_each_edge([&](const edge_t& edge)
    {
      nid_t id_from = truth->get_id(edge.first), id_to = truth->get_id(edge.second);
      bool rev_from = truth->get_is_reverse(edge.first), rev_to = truth->get_is_reverse(edge.second);
      handle_t from = gfa_graph.get_handle(id_from, rev_from);
      handle_t to = gfa_graph.get_handle(id_to, rev_to);
      EXPECT_TRUE(gfa_graph.has_edge(from, to)) << "GFA graph does not have the edge ((" << id_from << ", " << rev_from <<"), (" << id_to << ", " << rev_to << "))";
    });
  }

  void check_paths(const GBWTGraph& gfa_graph, const GBWTGraph* truth) const
  {
    std::unordered_map<std::string, handlegraph::PathSense> truth_paths;
    truth->for_each_path_matching(nullptr, nullptr, nullptr, [&](const path_handle_t truth_path)
    {
      truth_paths[truth->get_path_name(truth_path)] = truth->get_sense(truth_path);
    });

    std::unordered_map<std::string, handlegraph::PathSense> observed_paths;
    gfa_graph.for_each_path_matching(nullptr, nullptr, nullptr, [&](const path_handle_t gfa_path)
    {
      observed_paths[gfa_graph.get_path_name(gfa_path)] = gfa_graph.get_sense(gfa_path);
    });
    
    for (auto& kv : observed_paths) {
      EXPECT_TRUE(truth_paths.count(kv.first)) << "GFA has extraneous path " << kv.first;
    }
    for (auto& kv : truth_paths) {
      ASSERT_TRUE(observed_paths.count(kv.first)) << "GFA is missing truth path " << kv.first;
      EXPECT_EQ(observed_paths[kv.first], kv.second) << "GFA has wrong sense for path " << kv.first;
    }
  }

  void check_translation(const SequenceSource& source, const std::vector<translation_type>& truth) const
  {
    ASSERT_TRUE(source.uses_translation()) << "Segment translation not in use";
    ASSERT_GE(source.segment_translation.size(), truth.size()) << "Too few segments in the translation";
    for(const translation_type& translation : truth)
    {
      EXPECT_EQ(source.get_translation(translation.first), translation.second) << "Invalid translation for " << translation.first;
    }
  }

  void check_no_translation(const GBWTGraph& graph) const
  {
    ASSERT_FALSE(graph.has_segment_names()) << "The graph has segment names";
    graph.for_each_handle([&](const handle_t& handle)
    {
      nid_t id = graph.get_id(handle);
      std::pair<std::string, size_t> correct(std::to_string(id), 0);
      EXPECT_EQ(graph.get_segment_name_and_offset(handle), correct) << "Invalid null translation for node " << id;
      EXPECT_EQ(graph.get_segment_name(handle), correct.first) << "Invalid null segment name for node " << id;
      EXPECT_EQ(graph.get_segment_offset(handle), correct.second) << "Invalid null segment offset for node " << id;
      handle_t reverse = graph.flip(handle);
      EXPECT_EQ(graph.get_segment_name_and_offset(reverse), correct) << "Invalid null translation for reverse node " << id;
      EXPECT_EQ(graph.get_segment_name(reverse), correct.first) << "Invalid null segment name for reverse node " << id;
      EXPECT_EQ(graph.get_segment_offset(reverse), correct.second) << "Invalid null segment offset for reverse node " << id;
      auto segment = graph.get_segment(handle);
      std::pair<nid_t, nid_t> range(id, id + 1);
      EXPECT_EQ(segment.first, correct.first) << "Invalid null segment containing node " << id;
      EXPECT_EQ(segment.second, range) << "Invalid node range for null segment containing node " << id;
    });
  }

  void check_translation(const GBWTGraph& graph, const std::vector<translation_type>& truth) const
  {
    ASSERT_TRUE(graph.has_segment_names()) << "The graph has no segment names";

    // Forward orientation.
    for(const translation_type& translation : truth)
    {
      size_t offset = 0;
      for(nid_t id = translation.second.first; id < translation.second.second; id++)
      {
        handle_t handle = graph.get_handle(id, false);
        size_t node_length = graph.get_length(handle);
        ASSERT_TRUE(node_length > 0) << "Node appears empty";

        // Check the SegmentHandleGraph API
        std::pair<std::string, size_t> correct(translation.first, offset);
        EXPECT_EQ(graph.get_segment_name_and_offset(handle), correct) << "Invalid translation for node " << id;
        EXPECT_EQ(graph.get_segment_name(handle), correct.first) << "Invalid segment name for node " << id;
        EXPECT_EQ(graph.get_segment_offset(handle), correct.second) << "Invalid segment offset for node " << id;

        // Check the NamedNodeBackTranslation API
        for (size_t start = 0; start + 1 < node_length; start++) {
          for (size_t length = 1; start + length <= node_length; length++) {
            // For each possible range on this node, make a range.
            oriented_node_range_t graph_range {id, false, start, length};
            // And translate it back to GFA coordinates
            std::vector<oriented_node_range_t> back_translated = graph.translate_back(graph_range);

            EXPECT_EQ(back_translated.size(), 1ul) << "GFA segments appear merged in graph";
            EXPECT_EQ(graph.get_back_graph_node_name(std::get<0>(back_translated[0])), translation.first) << "Graph node belongs to wrong GFA segment";
            ASSERT_FALSE(std::get<1>(back_translated[0])) << "Graph node is backward in GFA segment";
            EXPECT_EQ(std::get<2>(back_translated[0]), offset + start) << "Graph node range begins at wrong offset on GFA segment";
            EXPECT_EQ(std::get<3>(back_translated[0]), length) << "Graph node range translated to GFA range of wrong length";
          }
        }

        offset += node_length;
      }
    }

    // Reverse orientation.
    for(const translation_type& translation : truth)
    {
      size_t offset = 0;
      for(nid_t i = translation.second.second; i > translation.second.first; i--)
      {
        nid_t id = i - 1;
        handle_t handle = graph.get_handle(id, true);
        size_t node_length = graph.get_length(handle);
        ASSERT_TRUE(node_length > 0) << "Node appears empty";

        // Check the SegmentHandleGraph API
        std::pair<std::string, size_t> correct(translation.first, offset);
        EXPECT_EQ(graph.get_segment_name_and_offset(handle), correct) << "Invalid translation for reverse node " << id;
        EXPECT_EQ(graph.get_segment_name(handle), correct.first) << "Invalid segment name for reverse node " << id;
        EXPECT_EQ(graph.get_segment_offset(handle), correct.second) << "Invalid segment offset for reverse node " << id;

        // Check the NamedNodeBackTranslation API
        for (size_t start = 0; start + 1 < node_length; start++) {
          for (size_t length = 1; start + length <= node_length; length++) {
            // For each possible range on this node, make a range.
            oriented_node_range_t graph_range {id, true, start, length};
            // And translate it back to GFA coordinates
            std::vector<oriented_node_range_t> back_translated = graph.translate_back(graph_range);

            EXPECT_EQ(back_translated.size(), 1ul) << "GFA segments appear merged in graph";
            EXPECT_EQ(graph.get_back_graph_node_name(std::get<0>(back_translated[0])), translation.first) << "Graph node belongs to wrong GFA segment";
            ASSERT_TRUE(std::get<1>(back_translated[0])) << "Graph node is backward in GFA segment";
            EXPECT_EQ(std::get<2>(back_translated[0]), offset + start) << "Graph node range begins at wrong offset on GFA segment";
            EXPECT_EQ(std::get<3>(back_translated[0]), length) << "Graph node range translated to GFA range of wrong length";
          }
        }

        offset += node_length;
      }
    }

    // Entire segments.
    for(const translation_type& translation : truth)
    {
      for(nid_t id = translation.second.first; id < translation.second.second; id++)
      {
        EXPECT_EQ(graph.get_segment(graph.get_handle(id, false)), translation) << "Invalid segment for node " << id;
      }
    }

    // For each segment.
    bool ok = true;
    auto iter = truth.begin();
    graph.for_each_segment([&](const std::string& name, std::pair<nid_t, nid_t> nodes) -> bool
    {
      if(iter == truth.end() || name != iter->first || nodes != iter->second)
      {
        ok = false; return false;
      }
      ++iter;
      return true;
    });
    ASSERT_TRUE(ok) << "for_each_segment() did not find the right translations";
    EXPECT_EQ(iter, truth.end()) << "for_each_segment() did not find all translations";
  }

  void check_links(const GBWTGraph& graph, const std::vector<edge_t>& edges) const
  {
    bool ok = true;
    auto iter = edges.begin();
    graph.for_each_link([&](const edge_t& edge, const std::string& from, const std::string& to) -> bool
    {
      std::string from_segment = graph.get_segment_name(edge.first);
      std::string to_segment = graph.get_segment_name(edge.second);
      if(iter == edges.end() || edge != *iter || from != from_segment || to != to_segment)
      {
        ok = false; return false;
      }
      ++iter;
      return true;
    });
    ASSERT_TRUE(ok) << "for_each_link() did not find the right links";
    EXPECT_EQ(iter, edges.end()) << "for_each_links() did not find all links";
  }
};

TEST_F(GFAConstruction, NormalGraph)
{
  auto gfa_parse = gfa_to_gbwt("gfas/example.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  ASSERT_FALSE(gfa_parse.second->uses_translation()) << "Unnecessary segment translation";

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_no_translation(graph);
}

TEST_F(GFAConstruction, WithZeroSegment)
{
  auto gfa_parse = gfa_to_gbwt("gfas/example_0-based.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  std::vector<translation_type> translation =
  {
    { "0", { 1, 2 } },
    { "1", { 2, 3 } },
    { "3", { 4, 5 } },
    { "4", { 5, 6 } },
    { "5", { 6, 7 } },
    { "6", { 7, 8 } },
    { "7", { 8, 9 } },
    { "8", { 9, 10 } },
  };
  this->check_translation(*(gfa_parse.second), translation);

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_translation(graph, translation);

  std::vector<edge_t> links =
  {
    edge_t(graph.get_handle(1, false), graph.get_handle(2, false)),
    edge_t(graph.get_handle(1, false), graph.get_handle(4, false)),
    edge_t(graph.get_handle(2, false), graph.get_handle(4, false)),
    edge_t(graph.get_handle(4, false), graph.get_handle(5, false)),
    edge_t(graph.get_handle(5, false), graph.get_handle(6, false)),
    edge_t(graph.get_handle(6, false), graph.get_handle(7, false)),
    edge_t(graph.get_handle(6, false), graph.get_handle(8, false)),
    edge_t(graph.get_handle(7, false), graph.get_handle(9, false)),
    edge_t(graph.get_handle(8, false), graph.get_handle(9, false)),
  };
  this->check_links(graph, links);
}

TEST_F(GFAConstruction, StringSegmentNames)
{
  auto gfa_parse = gfa_to_gbwt("gfas/example_str-names.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  std::vector<translation_type> translation =
  {
    { "s1", { 1, 2 } },
    { "s2", { 2, 3 } },
    { "s4", { 4, 5 } },
    { "s5", { 5, 6 } },
    { "s6", { 6, 7 } },
    { "s7", { 7, 8 } },
    { "s8", { 8, 9 } },
    { "s9", { 9, 10 } },
  };
  this->check_translation(*(gfa_parse.second), translation);

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_translation(graph, translation);

  std::vector<edge_t> links =
  {
    edge_t(graph.get_handle(1, false), graph.get_handle(2, false)),
    edge_t(graph.get_handle(1, false), graph.get_handle(4, false)),
    edge_t(graph.get_handle(2, false), graph.get_handle(4, false)),
    edge_t(graph.get_handle(4, false), graph.get_handle(5, false)),
    edge_t(graph.get_handle(5, false), graph.get_handle(6, false)),
    edge_t(graph.get_handle(6, false), graph.get_handle(7, false)),
    edge_t(graph.get_handle(6, false), graph.get_handle(8, false)),
    edge_t(graph.get_handle(7, false), graph.get_handle(9, false)),
    edge_t(graph.get_handle(8, false), graph.get_handle(9, false)),
  };
  this->check_links(graph, links);
}

TEST_F(GFAConstruction, SegmentChopping)
{
  GFAParsingParameters parameters;
  parameters.max_node_length = 3;
  auto gfa_parse = gfa_to_gbwt("gfas/example_chopping.gfa", parameters);
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  std::vector<translation_type> translation =
  {
    { "1", { 1, 2 } },
    { "2", { 2, 3 } },
    { "4", { 4, 6 } },
    { "6", { 6, 7 } },
    { "7", { 7, 8 } },
    { "8", { 8, 9 } },
    { "9", { 9, 10 } },
  };
  this->check_translation(*(gfa_parse.second), translation);

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_translation(graph, translation);

  std::vector<edge_t> links =
  {
    edge_t(graph.get_handle(1, false), graph.get_handle(2, false)),
    edge_t(graph.get_handle(1, false), graph.get_handle(4, false)),
    edge_t(graph.get_handle(2, false), graph.get_handle(4, false)),
    edge_t(graph.get_handle(5, false), graph.get_handle(6, false)),
    edge_t(graph.get_handle(6, false), graph.get_handle(7, false)),
    edge_t(graph.get_handle(6, false), graph.get_handle(8, false)),
    edge_t(graph.get_handle(7, false), graph.get_handle(9, false)),
    edge_t(graph.get_handle(8, false), graph.get_handle(9, false)),
  };
  this->check_links(graph, links);
}

class GFAConstructionReversal : public GFAConstruction
{
public:

  void SetUp() override
  {
    // We need to parse the P line in our truth GFA as a haplotype, so we can
    // match the W line in the test GFA.
    GFAParsingParameters parameters;
    parameters.path_name_formats.clear();
    parameters.path_name_formats.emplace_back(
      GFAParsingParameters::PAN_SN_REGEX,
      GFAParsingParameters::PAN_SN_FIELDS,
      GFAParsingParameters::PAN_SN_SENSE
    );
    auto truth = gfa_to_gbwt("gfas/reversal.gfa", parameters);
    this->index = *truth.first;
    this->source = *truth.second;
    this->graph = GBWTGraph(this->index, this->source);
  }
};

TEST_F(GFAConstructionReversal, ChoppingWithReversal)
{
  GFAParsingParameters parameters;
  parameters.max_node_length = 3;
  auto gfa_parse = gfa_to_gbwt("gfas/reversal_chopping.gfa", parameters);
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  std::vector<translation_type> translation =
  {
    { "1", { 1, 2 } },
    { "2", { 2, 3 } },
    { "3", { 3, 4 } },
    { "4", { 4, 6 } },
  };
  this->check_translation(*(gfa_parse.second), translation);

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_translation(graph, translation);

  std::vector<edge_t> links =
  {
    edge_t(graph.get_handle(1, false), graph.get_handle(2, false)),
    edge_t(graph.get_handle(1, false), graph.get_handle(3, false)),
    edge_t(graph.get_handle(2, false), graph.get_handle(5, true)),
    edge_t(graph.get_handle(3, false), graph.get_handle(4, false)),
  };
  this->check_links(graph, links);
}

TEST_F(GFAConstructionReversal, WalksWithReversal)
{
  auto gfa_parse = gfa_to_gbwt("gfas/reversal_walks.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  ASSERT_FALSE(gfa_parse.second->uses_translation()) << "Unnecessary segment translation";

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_paths(graph, &(this->graph));
  this->check_no_translation(graph);
}

class GFAConstructionWalks : public GFAConstruction
{
public:
  void SetUp() override
  {
    this->index = build_gbwt_example_walks();
    build_source(this->source);
    this->graph = GBWTGraph(this->index, this->source);
  }
};

TEST_F(GFAConstructionWalks, WalksAndPaths)
{
  auto gfa_parse = gfa_to_gbwt("gfas/example_walks.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  ASSERT_FALSE(gfa_parse.second->uses_translation()) << "Unnecessary segment translation";

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_paths(graph, &(this->graph));
  this->check_no_translation(graph);
}

TEST_F(GFAConstructionWalks, WalksOnly)
{
  auto gfa_parse = gfa_to_gbwt("gfas/example_walks-only.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  ASSERT_FALSE(gfa_parse.second->uses_translation()) << "Unnecessary segment translation";

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_paths(graph, &(this->graph));
  this->check_no_translation(graph);
}

class GFAConstructionReference : public GFAConstruction
{
public:
  void SetUp() override
  {
    this->index = build_gbwt_example_reference();
    build_source(this->source);
    this->graph = GBWTGraph(this->index, this->source);
  }
};

TEST_F(GFAConstructionReference, ReferencePaths)
{
  GFAParsingParameters parameters;
  // Parse panSN paths if possible.
  parameters.path_name_formats.emplace_front(GFAParsingParameters::PAN_SN_REGEX, GFAParsingParameters::PAN_SN_FIELDS, GFAParsingParameters::PAN_SN_SENSE);
  auto gfa_parse = gfa_to_gbwt("gfas/example_reference.gfa", parameters);
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  ASSERT_FALSE(gfa_parse.second->uses_translation()) << "Unnecessary segment translation";

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
  this->check_paths(graph, &(this->graph));
  this->check_no_translation(graph);
}

//------------------------------------------------------------------------------

class GBWTSubgraph : public ::testing::Test
{
public:
  gbwt::GBWT select_paths(const gbwt::GBWT& source, gbwt::size_type skip) const
  {
    gbwt::GBWTBuilder builder(sdsl::bits::length(source.sigma() - 1), source.size());
    for(gbwt::size_type path_id = 0; path_id < source.metadata.paths(); path_id += 1 + skip)
    {
      gbwt::vector_type path = source.extract(gbwt::Path::encode(path_id, false));
      builder.insert(path, true);
    }
    builder.finish();
    return builder.index;
  }

  void check_subgraph(const GBWTGraph& graph, const GBWTGraph& subgraph) const
  {
    bool nodes_ok = true;
    bool sequences_ok = true;
    bool reverse_complements_ok = true;
    subgraph.for_each_handle([&](const handle_t& handle)
    {
      if(!(graph.has_node(subgraph.get_id(handle)))) { nodes_ok = false; return; }
      if(subgraph.get_sequence(handle) != graph.get_sequence(handle)) { sequences_ok = false; }
      handle_t flipped = subgraph.flip(handle);
      if(subgraph.get_sequence(flipped) != graph.get_sequence(flipped)) { reverse_complements_ok = false; }
    });
    ASSERT_TRUE(nodes_ok) << "Some nodes were missing from the supergraph";
    ASSERT_TRUE(sequences_ok) << "Some sequences were not identical";
    ASSERT_TRUE(reverse_complements_ok) << "Some reverse complement sequences were not identical";

    bool edges_ok = true;
    subgraph.for_each_edge([&](const edge_t& edge)
    {
      if(!(graph.has_edge(edge.first, edge.second))) { edges_ok = false; }
    });
    ASSERT_TRUE(edges_ok) << "Some edges were missing from the supergraph";
  }

  void check_translation(const GBWTGraph& graph, const GBWTGraph& subgraph) const
  {
    ASSERT_EQ(subgraph.has_segment_names(), graph.has_segment_names()) << "Node-to-segment translation mismatch";
    if(!(subgraph.has_segment_names())) { return; }

    bool translation_ok = true;
    subgraph.for_each_handle([&](const handle_t& handle)
    {
      auto subgraph_translation = subgraph.get_segment(handle);
      auto graph_translation = graph.get_segment(handle);
      if(subgraph_translation != graph_translation) { translation_ok = false; }
    });
    ASSERT_TRUE(translation_ok) << "Some translations were not identical";
  }
};

TEST_F(GBWTSubgraph, WithoutTranslation)
{
  auto gfa_parse = gfa_to_gbwt("gfas/for_subgraph.gfa");
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));
  gbwt::GBWT selected = this->select_paths(*(gfa_parse.first), 1);
  GBWTGraph subgraph = graph.subgraph(selected);

  ASSERT_NO_THROW(subgraph.sanity_checks()) << "The subgraph failed sanity checks";
  this->check_subgraph(graph, subgraph);
  this->check_translation(graph, subgraph);
}

TEST_F(GBWTSubgraph, WithTranslation)
{
  GFAParsingParameters parameters;
  parameters.max_node_length = 3;
  auto gfa_parse = gfa_to_gbwt("gfas/for_subgraph.gfa", parameters);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));
  gbwt::GBWT selected = this->select_paths(*(gfa_parse.first), 1);
  GBWTGraph subgraph = graph.subgraph(selected);

  ASSERT_NO_THROW(subgraph.sanity_checks()) << "The subgraph failed sanity checks";
  this->check_subgraph(graph, subgraph);
  this->check_translation(graph, subgraph);
}

//------------------------------------------------------------------------------

class GFAExtraction : public ::testing::Test
{
public:
  void extract_gfa(const GBWTGraph& graph, const GraphName* graph_name, const std::string& filename, const GFAExtractionParameters& parameters) const
  {
    std::ofstream out(filename, std::ios_base::binary);
    gbwt_to_gfa(graph, graph_name, out, parameters);
    out.close();
  }

  void compare_gfas(const std::string& test_file, const std::string& truth_file, const std::string& name) const
  {
    std::vector<std::string> test_rows, truth_rows;
    gbwt::readRows(test_file, test_rows, false);
    gbwt::readRows(truth_file, truth_rows, false);

    ASSERT_EQ(test_rows.size(), truth_rows.size()) << name << ": Invalid number of rows";
    size_t line = 0;
    for(auto iter = test_rows.begin(), truth = truth_rows.begin(); iter != test_rows.end(); ++iter, ++truth)
    {
      EXPECT_EQ(*iter, *truth) << name << ": Invalid line " << line;
      line++;
    }
  }
};

TEST_F(GFAExtraction, Components)
{
  std::string input = "gfas/components_walks.gfa";
  auto gfa_parse = gfa_to_gbwt(input);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  std::string output = gbwt::TempFile::getName("gfa-extraction");
  GFAExtractionParameters parameters;
  this->extract_gfa(graph, nullptr, output, parameters);

  this->compare_gfas(output, input, "Components");
  gbwt::TempFile::remove(output);
}

TEST_F(GFAExtraction, PathsAndWalks)
{
  std::string input = "gfas/example_walks.gfa";
  auto gfa_parse = gfa_to_gbwt(input);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  std::string output = gbwt::TempFile::getName("gfa-extraction");
  GFAExtractionParameters parameters;
  this->extract_gfa(graph, nullptr, output, parameters);

  this->compare_gfas(output, input, "Paths and walks");
  gbwt::TempFile::remove(output);
}

TEST_F(GFAExtraction, CacheRecords)
{
  std::string input = "gfas/components_walks.gfa";
  auto gfa_parse = gfa_to_gbwt(input);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  std::vector<size_t> cache_limits = { 0, 1, 2, 4, 8, 16 };
  for(size_t limit : cache_limits)
  {
    std::string output = gbwt::TempFile::getName("gfa-extraction");
    GFAExtractionParameters parameters;
    parameters.large_record_bytes = limit;
    this->extract_gfa(graph, nullptr, output, parameters);
    std::string name = "Cache records " + std::to_string(limit);
    this->compare_gfas(output, input, name);
    gbwt::TempFile::remove(output);
  }
}

TEST_F(GFAExtraction, PathModes)
{
  std::string input = "gfas/default.gfa";
  auto gfa_parse = gfa_to_gbwt(input);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  // Default mode.
  {
    std::string truth = "gfas/default.gfa";
    std::string output = gbwt::TempFile::getName("gfa-modes");
    GFAExtractionParameters parameters; parameters.mode = GFAExtractionParameters::mode_default;
    this->extract_gfa(graph, nullptr, output, parameters);
    this->compare_gfas(output, truth, "Default");
    gbwt::TempFile::remove(output);
  }

  // PanSN mode.
  {
    std::string truth = "gfas/pan-sn.gfa";
    std::string output = gbwt::TempFile::getName("gfa-modes");
    GFAExtractionParameters parameters; parameters.mode = GFAExtractionParameters::mode_pan_sn;
    this->extract_gfa(graph, nullptr, output, parameters);
    this->compare_gfas(output, truth, "Default");
    gbwt::TempFile::remove(output);
  }

  // Reference-only mode.
  {
    std::string truth = "gfas/ref-only.gfa";
    std::string output = gbwt::TempFile::getName("gfa-modes");
    GFAExtractionParameters parameters; parameters.mode = GFAExtractionParameters::mode_ref_only;
    this->extract_gfa(graph, nullptr, output, parameters);
    this->compare_gfas(output, truth, "Default");
    gbwt::TempFile::remove(output);
  }
}

TEST_F(GFAExtraction, Translation)
{
  GFAParsingParameters parameters; parameters.max_node_length = 3;
  std::string input = "gfas/example_chopping.gfa";
  auto gfa_parse = gfa_to_gbwt(input, parameters);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  // Use translation.
  {
    std::string truth = "gfas/example_from_chopping.gfa";
    std::string output = gbwt::TempFile::getName("gfa-translation");
    GFAExtractionParameters parameters; parameters.use_translation = true;
    this->extract_gfa(graph, nullptr, output, parameters);
    this->compare_gfas(output, truth, "With translation");
    gbwt::TempFile::remove(output);
  }

  // No translation.
  {
    std::string truth = "gfas/example_chopped.gfa";
    std::string output = gbwt::TempFile::getName("gfa-translation");
    GFAExtractionParameters parameters; parameters.use_translation = false;
    this->extract_gfa(graph, nullptr, output, parameters);
    this->compare_gfas(output, truth, "Without translation");
    gbwt::TempFile::remove(output);
  }
}

TEST_F(GFAExtraction, CanonicalGFA)
{
  std::vector<std::string> inputs = { "gfas/components_walks.gfa", "gfas/example_walks.gfa", "gfas/example_chopping.gfa" };

  for(const std::string& input : inputs)
  {
    auto gfa_parse = gfa_to_gbwt(input);
    GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

    std::stringstream truth_stream;
    handlegraph::algorithms::canonical_gfa(graph, truth_stream, true);
    std::string truth = truth_stream.str();

    std::stringstream output_stream;
    gbwt_to_canonical_gfa(graph, output_stream);
    std::string output = output_stream.str();

    ASSERT_EQ(output, truth) << "Canonical GFA mismatch for " << input;
  }
}

//------------------------------------------------------------------------------

class GBWTMetadata : public ::testing::Test
{
public:
  std::string path_name_regex;
  std::string samples_and_haplotypes;
  std::string contigs_and_fragments;

  std::vector<std::string> names;
  std::vector<gbwt::PathName::path_name_type> ids;

  GBWTMetadata()
  {
  }

  void SetUp() override
  {
    this->path_name_regex = "(.)(.)";
    this->samples_and_haplotypes = "-SH";
    this->contigs_and_fragments = "-CF";
    this->names = { "A", "B" };
    this->ids =
    {
      static_cast<gbwt::PathName::path_name_type>(1),
      static_cast<gbwt::PathName::path_name_type>(2)
    };
  }

  void check_metadata(const gbwt::Metadata& metadata, const gbwt::Metadata& expected)
  {
    // Samples.
    ASSERT_EQ(metadata.samples(), expected.samples()) << "Wrong number of samples";
    if(expected.hasSampleNames())
    {
      ASSERT_TRUE(metadata.hasSampleNames()) << "No sample names";
      bool names_ok = true;
      for(gbwt::size_type i = 0; i < metadata.samples(); i++)
      {
        if(metadata.sample(i) != expected.sample(i)) { names_ok = false; }
      }
      EXPECT_TRUE(names_ok) << "Invalid sample names";
    }
    else
    {
      EXPECT_FALSE(metadata.hasSampleNames()) << "Sample names were created";
    }

    // Haplotypes.
    EXPECT_EQ(metadata.haplotypes(), expected.haplotypes()) << "Wrong number of haplotypes";

    // Contigs.
    ASSERT_EQ(metadata.contigs(), expected.contigs()) << "Wrong number of contigs";
    if(expected.hasContigNames())
    {
      ASSERT_TRUE(metadata.hasContigNames()) << "No contig names";
      bool names_ok = true;
      for(gbwt::size_type i = 0; i < metadata.contigs(); i++)
      {
        if(metadata.contig(i) != expected.contig(i)) { names_ok = false; }
      }
      EXPECT_TRUE(names_ok) << "Invalid contig names";
    }
    else
    {
      EXPECT_FALSE(metadata.hasContigNames()) << "Contig names were created";
    }

    // Paths.
    ASSERT_EQ(metadata.paths(), expected.paths()) << "Wrong number of paths";
    if(expected.hasPathNames())
    {
      ASSERT_TRUE(metadata.hasPathNames()) << "No path names";
      bool names_ok = true;
      for(gbwt::size_type i = 0; i < metadata.paths(); i++)
      {
        if(metadata.path(i) != expected.path(i)) { names_ok = false; }
      }
      EXPECT_TRUE(names_ok) << "Invalid path names";
    }
    else
    {
      EXPECT_FALSE(metadata.hasPathNames()) << "Path names were created";
    }
  }
};

TEST_F(GBWTMetadata, SamplesAndHaplotypes)
{
  GFAParsingParameters parameters;
  parameters.path_name_formats.clear();
  parameters.path_name_formats.emplace_back(this->path_name_regex, this->samples_and_haplotypes, PathSense::HAPLOTYPE);
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa", parameters);
  const gbwt::GBWT& index = *(gfa_parse.first);

  gbwt::Metadata expected_metadata;
  expected_metadata.setSamples(this->names);
  expected_metadata.setHaplotypes(4);
  expected_metadata.setContigs(1);
  expected_metadata.addPath(0, 0, 1, 0);
  expected_metadata.addPath(0, 0, 2, 0);
  expected_metadata.addPath(1, 0, 1, 0);
  expected_metadata.addPath(1, 0, 2, 0);

  ASSERT_TRUE(index.hasMetadata()) << "No GBWT metadata was created";
  this->check_metadata(index.metadata, expected_metadata);
}

TEST_F(GBWTMetadata, ContigsAndFragments)
{
  GFAParsingParameters parameters;
  parameters.path_name_formats.clear();
  parameters.path_name_formats.emplace_back(this->path_name_regex, this->contigs_and_fragments, PathSense::HAPLOTYPE);
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa", parameters);
  const gbwt::GBWT& index = *(gfa_parse.first);

  gbwt::Metadata expected_metadata;
  expected_metadata.setSamples(1);
  expected_metadata.setHaplotypes(1);
  expected_metadata.setContigs(this->names);
  expected_metadata.addPath(0, 0, GBWTGraph::NO_PHASE, 1);
  expected_metadata.addPath(0, 0, GBWTGraph::NO_PHASE, 2);
  expected_metadata.addPath(0, 1, GBWTGraph::NO_PHASE, 1);
  expected_metadata.addPath(0, 1, GBWTGraph::NO_PHASE, 2);

  ASSERT_TRUE(index.hasMetadata()) << "No GBWT metadata was created";
  this->check_metadata(index.metadata, expected_metadata);
}

TEST_F(GBWTMetadata, Walks)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components_walks.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);

  gbwt::Metadata expected_metadata;
  std::vector<std::string> samples = { "sample" };
  expected_metadata.setSamples(samples);
  expected_metadata.setHaplotypes(2);
  expected_metadata.setContigs(this->names);
  expected_metadata.addPath(0, 0, 1, 0);
  expected_metadata.addPath(0, 0, 2, 0);
  expected_metadata.addPath(0, 1, 1, 0);
  expected_metadata.addPath(0, 1, 2, 0);

  ASSERT_TRUE(index.hasMetadata()) << "No GBWT metadata was created";
  this->check_metadata(index.metadata, expected_metadata);
}

TEST_F(GBWTMetadata, WalksNoInterval)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components_walks_no_interval.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);

  gbwt::Metadata expected_metadata;
  std::vector<std::string> samples = { "sample" };
  expected_metadata.setSamples(samples);
  expected_metadata.setHaplotypes(2);
  expected_metadata.setContigs(this->names);
  expected_metadata.addPath(0, 0, 1, 0);
  expected_metadata.addPath(0, 0, 2, 0);
  expected_metadata.addPath(0, 1, 1, 0);
  expected_metadata.addPath(0, 1, 2, 0);

  ASSERT_TRUE(index.hasMetadata()) << "No GBWT metadata was created";
  this->check_metadata(index.metadata, expected_metadata);
}

TEST_F(GBWTMetadata, WalksAndPaths)
{
  auto gfa_parse = gfa_to_gbwt("gfas/example_walks.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);

  gbwt::Metadata expected_metadata;
  std::vector<std::string> samples = { REFERENCE_PATH_SAMPLE_NAME, "short", "alt" };
  expected_metadata.setSamples(samples);
  expected_metadata.setHaplotypes(4);
  std::vector<std::string> contigs = { "short", "alt1", "alt2", "chr" };
  expected_metadata.setContigs(contigs);
  expected_metadata.addPath(0, 0, GBWTGraph::NO_PHASE, 0);
  expected_metadata.addPath(0, 1, GBWTGraph::NO_PHASE, 0);
  expected_metadata.addPath(0, 2, GBWTGraph::NO_PHASE, 0);
  expected_metadata.addPath(1, 3, 1, 0);
  expected_metadata.addPath(2, 3, 0, 0);
  expected_metadata.addPath(1, 3, 2, 0);

  ASSERT_TRUE(index.hasMetadata()) << "No GBWT metadata was created";
  this->check_metadata(index.metadata, expected_metadata);
}

//------------------------------------------------------------------------------

} // namespace
