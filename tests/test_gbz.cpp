#include <gtest/gtest.h>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class GBZSerialization : public ::testing::Test
{
public:
  std::unique_ptr<GBZ> create_gbz()
  {
    SequenceSource source; build_source(source);
    return std::make_unique<GBZ>(build_gbwt_index(), source);
  }

  void check_gbz(const GBZ& gbz, const GBZ& truth, bool check_tags = true) const
  {
    // GBZ
    ASSERT_EQ(gbz.header, truth.header) << "GBZ: Invalid header";
    if(check_tags)
    {
      ASSERT_EQ(gbz.tags, truth.tags) << "GBZ: Invalid tags";
    }

    // GBWT
    ASSERT_EQ(gbz.index.size(), truth.index.size()) << "GBWT: Invalid size";
    ASSERT_EQ(gbz.index.sequences(), truth.index.sequences()) << "GBWT: Invalid number of sequences";
    ASSERT_EQ(gbz.index.sigma(), truth.index.sigma()) << "GBWT: Invalid alphabet size";
    ASSERT_EQ(gbz.index.effective(), truth.index.effective()) << "GBWT: Invalid effective alphabet size";
    ASSERT_EQ(gbz.index.samples(), truth.index.samples()) << "GBWT: Invalid number of samples";

    // Graph
    ASSERT_EQ(gbz.graph.header, truth.graph.header) << "Graph: Invalid header";
    ASSERT_EQ(gbz.graph.sequences, truth.graph.sequences) << "Graph: Invalid sequences";
    ASSERT_EQ(gbz.graph.real_nodes, truth.graph.real_nodes) << "Graph: Invalid real nodes";
    ASSERT_EQ(gbz.graph.segments, truth.graph.segments) << "Graph: Invalid segments";
    ASSERT_EQ(gbz.graph.node_to_segment, truth.graph.node_to_segment) << "Graph: Invalid node-to-segment mapping";
  }
};

TEST_F(GBZSerialization, Empty)
{
  GBZ empty;
  size_t expected_size = empty.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
  std::string filename = gbwt::TempFile::getName("gbz");
  sdsl::simple_sds::serialize_to(empty, filename);

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  size_t bytes = gbwt::fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid file size";
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, empty);

  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, NonEmpty)
{
  std::unique_ptr<GBZ> original = this->create_gbz();
  size_t expected_size = original->simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
  std::string filename = gbwt::TempFile::getName("gbz");
  sdsl::simple_sds::serialize_to(*original, filename);

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  size_t bytes = gbwt::fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid file size";
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, *original);

  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, ExternalObjects)
{
  // Serialization into separate GBWT / GBWTGraph files does not preserve tags.
  // We therefore do not check file sizes.
  std::unique_ptr<GBZ> original = this->create_gbz();
  std::string filename = gbwt::TempFile::getName("gbz");
  std::ofstream out(filename, std::ios_base::binary);
  GBZ::simple_sds_serialize(original->index, original->graph, out);
  out.close();

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, *original, false);

  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, CopyAndSerialize)
{
  std::string filename = gbwt::TempFile::getName("gbz");
  {
    std::unique_ptr<GBZ> original = this->create_gbz();
    GBZ copied = *original; original.reset();
    sdsl::simple_sds::serialize_to(copied, filename);
  }
  {
    std::unique_ptr<GBZ> truth = this->create_gbz();
    GBZ loaded; sdsl::simple_sds::load_from(loaded, filename);
    this->check_gbz(loaded, *truth);
  }
  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, MoveAndSerialize)
{
  std::string filename = gbwt::TempFile::getName("gbz");
  {
    std::unique_ptr<GBZ> original = this->create_gbz();
    GBZ moved = std::move(*original); original.reset();
    sdsl::simple_sds::serialize_to(moved, filename);
  }
  {
    std::unique_ptr<GBZ> truth = this->create_gbz();
    GBZ loaded; sdsl::simple_sds::load_from(loaded, filename);
    this->check_gbz(loaded, *truth);
  }
  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, SwapAndSerialize)
{
  std::string filename = gbwt::TempFile::getName("gbz");
  {
    std::unique_ptr<GBZ> original = this->create_gbz();
    GBZ swapped;
    swapped.swap(*original); original.reset();
    sdsl::simple_sds::serialize_to(swapped, filename);
  }
  {
    std::unique_ptr<GBZ> truth = this->create_gbz();
    GBZ loaded; sdsl::simple_sds::load_from(loaded, filename);
    this->check_gbz(loaded, *truth);
  }
  gbwt::TempFile::remove(filename);
}

//------------------------------------------------------------------------------

class GBZFunctionality : public ::testing::Test
{
public:
  GBZ build_gbz(const std::string& graph_name)
  {
    auto parse = gfa_to_gbwt(graph_name);
    return GBZ(parse.first, parse.second);
  }

  void set_reference_samples(GBZ& gbz, const std::unordered_set<std::string>& samples, size_t expected, const std::string& test)
  {
    size_t present = gbz.set_reference_samples(samples);
    ASSERT_EQ(present, expected) << test << ": Unexpected number of sample names present in the graph";
  }

  void check_named_paths(const GBZ& gbz, const std::unordered_set<std::string>& true_samples, size_t expected_paths, const std::string& test)
  {
    ASSERT_EQ(gbz.named_paths(), expected_paths) << test << ": Invalid number of named paths";

    const std::unordered_set<std::string>& samples = gbz.get_reference_samples();
    ASSERT_EQ(samples.size(), true_samples.size()) << test << ": Invalid number of reference samples";
    for(const std::string& sample : true_samples)
    {
      ASSERT_TRUE(samples.find(sample) != samples.end()) << test << ": Missing reference sample " << sample;
    }
  }

  // The test GFA files have their graph names computed with the reference pggname implementation.
  // Tests in test_gfa.cpp ensure that our implementation computes the same names.
  void check_graph_name(GBZ& gbz, const GraphName* parent, bool missing_name, const std::string& test_case)
  {
    EXPECT_EQ(gbz.pggname().empty(), missing_name) << "Unexpected graph name presence for " << test_case;
    EXPECT_EQ(gbz.translation_target(), "") << "Translation target should not be set for " << test_case;
    GraphName old = gbz.graph_name();
    EXPECT_EQ(gbz.pggname(), old.name()) << "Graph name mismatch for " << test_case;

    gbz.compute_pggname(parent);
    GraphName recomputed = gbz.graph_name();
    EXPECT_TRUE(recomputed.has_name()) << "Recomputed graph name missing for " << test_case;
    if(!missing_name)
    {
      EXPECT_EQ(recomputed, old) << "Graph name changed after recomputation for " << test_case;
    }
  }
};

TEST_F(GBZFunctionality, ReferenceSamples)
{
  GBZ gbz = this->build_gbz("gfas/components_ref.gfa");
  std::unordered_set<std::string> samples { "ref" };
  std::unordered_set<std::string> true_samples = samples;
  this->check_named_paths(gbz, true_samples, 2, "Initial graph");

  samples.erase("ref");
  samples.insert("sample");
  true_samples = samples;
  this->set_reference_samples(gbz, samples, 1, "New sample");
  this->check_named_paths(gbz, true_samples, 4, "New sample");

  samples.insert("ref");
  true_samples = samples;
  this->set_reference_samples(gbz, samples, 2, "Both samples");
  this->check_named_paths(gbz, true_samples, 6, "Both samples");

  samples.erase("sample");
  true_samples = samples;
  samples.insert("missing");
  this->set_reference_samples(gbz, samples, 1, "Invalid sample");
  this->check_named_paths(gbz, true_samples, 2, "Invalid sample");
}

TEST_F(GBZFunctionality, GraphNames)
{
  // Additional tests with parent graph relationships are in test_gfa.cpp.

  // Constructor from gfa_to_gbwt() output.
  {
    std::unique_ptr<gbwt::GBWT> index = std::make_unique<gbwt::GBWT>(build_gbwt_index());
    std::unique_ptr<SequenceSource> source = std::make_unique<SequenceSource>();
    build_source(*source);
    GBZ gbz(index, source);
    this->check_graph_name(gbz, nullptr, false, "gfa_to_gbwt() output");
  }

  // Constructor from GBWT and SequenceSource.
  {
    gbwt::GBWT index = build_gbwt_index();
    SequenceSource source; build_source(source);
    GBZ gbz(index, source);
    this->check_graph_name(gbz, nullptr, false, "GBWT and SequenceSource");
  }

  // Subgraph construction.
  {
    gbwt::GBWT index = build_gbwt_index();
    SequenceSource source; build_source(source);
    GBZ supergraph(index, source);
    GraphName parent = supergraph.graph_name();

    std::vector<gbwt::vector_type> paths { alt_path };
    gbwt::GBWT sub_index = build_gbwt(paths);
    GBZ gbz(std::move(sub_index), supergraph);
    this->check_graph_name(gbz, &parent, false, "subgraph construction");
  }

  // Constructor from GBWT and HandleGraph; graph name must be set manually.
  {
    gbwt::GBWT parent_index = build_gbwt_index();
    SequenceSource source; build_source(source);
    GBWTGraph parent_graph(parent_index, source);

    gbwt::GBWT index = parent_index;
    GBZ gbz(std::move(index), parent_graph);
    this->check_graph_name(gbz, nullptr, true, "GBWT and HandleGraph");
  }
}

//------------------------------------------------------------------------------

} // namespace
