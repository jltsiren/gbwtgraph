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
  auto gfa_parse = gfa_to_gbwt("components.gfa", parameters);
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
  auto gfa_parse = gfa_to_gbwt("components.gfa", parameters);
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

//------------------------------------------------------------------------------

} // namespace
