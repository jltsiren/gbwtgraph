#include <gtest/gtest.h>

#include <gbwtgraph/gbz.h>

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

  void check_gbz(const GBZ& gbz, const GBZ& truth) const
  {
    // GBZ
    ASSERT_EQ(gbz.header, truth.header) << "GBZ: Invalid header";
    ASSERT_EQ(gbz.tags, truth.tags) << "GBZ: Invalid tags";

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
  std::unique_ptr<GBZ> original = this->create_gbz();
  size_t expected_size = original->simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
  std::string filename = gbwt::TempFile::getName("gbz");
  std::ofstream out(filename, std::ios_base::binary);
  GBZ::simple_sds_serialize(original->index, original->graph, out);
  out.close();

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  size_t bytes = gbwt::fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid file size";
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, *original);

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

} // namespace
