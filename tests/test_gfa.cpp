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
    ASSERT_EQ(gfa_graph.header, truth->header) << "Graph headers are not identical";
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

  void check_translation(const SequenceSource& source, const std::vector<translation_type>& truth) const
  {
    ASSERT_TRUE(source.uses_translation()) << "Segment translation not in use";
    ASSERT_EQ(source.segment_translation.size(), truth.size()) << "Invalid number of segments";
    for(const translation_type& translation : truth)
    {
      EXPECT_EQ(source.get_translation(translation.first), translation.second) << "Invalid translation for " << translation.first;
    }
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
}

TEST_F(GFAConstruction, WalksAndPaths)
{
  auto gfa_parse = gfa_to_gbwt("gfas/example_walks.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  ASSERT_FALSE(gfa_parse.second->uses_translation()) << "Unnecessary segment translation";

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
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
    { "2", { 3, 4 } },
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
    { "s3", { 3, 4 } },
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
    { "3", { 3, 4 } },
    { "4", { 4, 6 } },
    { "6", { 6, 7 } },
    { "7", { 7, 8 } },
    { "8", { 8, 9 } },
    { "9", { 9, 10 } },
  };
  this->check_translation(*(gfa_parse.second), translation);

  this->check_gbwt(index, &(this->index));
  this->check_graph(graph, &(this->graph));
}

TEST_F(GFAConstruction, ChoppingWithReversal)
{
  auto truth = gfa_to_gbwt("gfas/reversal.gfa");
  GBWTGraph truth_graph(*(truth.first), *(truth.second));

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

  this->check_gbwt(index, truth.first.get());
  this->check_graph(graph, &truth_graph);
}

TEST_F(GFAConstruction, WalksWithReversal)
{
  auto truth = gfa_to_gbwt("gfas/reversal.gfa");
  GBWTGraph truth_graph(*(truth.first), *(truth.second));

  auto gfa_parse = gfa_to_gbwt("gfas/reversal_walks.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(*(gfa_parse.first), *(gfa_parse.second));

  ASSERT_FALSE(gfa_parse.second->uses_translation()) << "Unnecessary segment translation";

  this->check_gbwt(index, truth.first.get());
  this->check_graph(graph, &truth_graph);
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
  parameters.path_name_regex = this->path_name_regex;
  parameters.path_name_fields = this->samples_and_haplotypes;
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
  parameters.path_name_regex = this->path_name_regex;
  parameters.path_name_fields = this->contigs_and_fragments;
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa", parameters);
  const gbwt::GBWT& index = *(gfa_parse.first);

  gbwt::Metadata expected_metadata;
  expected_metadata.setSamples(1);
  expected_metadata.setHaplotypes(1);
  expected_metadata.setContigs(this->names);
  expected_metadata.addPath(0, 0, 0, 1);
  expected_metadata.addPath(0, 0, 0, 2);
  expected_metadata.addPath(0, 1, 0, 1);
  expected_metadata.addPath(0, 1, 0, 2);

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

//------------------------------------------------------------------------------

} // namespace
