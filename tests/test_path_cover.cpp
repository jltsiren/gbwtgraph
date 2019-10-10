#include <gtest/gtest.h>

#include <set>
#include <vector>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>
#include <gbwtgraph/path_cover.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class PathCoverTest : public ::testing::Test
{
public:
  gbwt::GBWT index;
  GBWTGraph graph;
  size_t components;
  std::vector<std::set<gbwt::vector_type>> correct_paths;

  PathCoverTest()
  {
  }

  void SetUp() override
  {
    auto gfa_parse = gfa_to_gbwt("components.gfa");
    this->index = *(gfa_parse.first);
    this->graph = GBWTGraph(this->index, *(gfa_parse.second));
    this->components = 2;
  }
};

TEST_F(PathCoverTest, CorrectPaths)
{
  size_t paths_per_component = 4;
  size_t context_length = 3;
  gbwt::size_type expected_sequences = this->components * paths_per_component * 2;
  std::vector<std::set<gbwt::vector_type>> correct_paths =
  {
    {
      {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(12, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(15, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(17, false))
      },
      {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(12, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(16, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(17, false))
      },
      {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(13, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(15, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(17, false))
      },
      {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(13, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(16, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(17, false))
      }
    },
    {
      {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(22, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(24, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(25, false))
      },
      {
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(22, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(24, false)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(23, true)),
        static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, true))
      }
    }
  };

  gbwt::GBWT cover = path_cover_gbwt(this->graph, paths_per_component, context_length);
  ASSERT_EQ(cover.sequences(), expected_sequences) << "Wrong number of sequences in the path cover GBWT";

  // We insert the smaller of a path and its reverse complement to handle paths
  // that flip the orientation.
  std::vector<std::set<gbwt::vector_type>> result(this->components);
  for(size_t i = 0; i < this->components; i++)
  {
    for(size_t j = 0; j < paths_per_component; j++)
    {
      size_t seq_id = 2 * (i * paths_per_component + j);
      gbwt::vector_type forward = cover.extract(seq_id), reverse;
      gbwt::reversePath(forward, reverse);
      result[i].insert(std::min(forward, reverse));
    }
  }
  for(size_t i = 0; i < this->components; i++)
  {
    ASSERT_EQ(result[i].size(), correct_paths[i].size()) << "Wrong number of distinct paths for component " << i;
    auto result_iter = result[i].begin();
    auto correct_iter = correct_paths[i].begin();
    while(result_iter != result[i].end())
    {
      EXPECT_EQ(*result_iter, *correct_iter) << "Wrong path in component " << i;
      ++result_iter; ++correct_iter;
    }
  }
}

TEST_F(PathCoverTest, Metadata)
{
  size_t paths_per_component = 4;
  size_t context_length = 3;
  size_t expected_paths = paths_per_component * this->components;

  gbwt::GBWT cover = path_cover_gbwt(this->graph, paths_per_component, context_length);
  ASSERT_TRUE(cover.hasMetadata()) << "Path cover GBWT contains no metadata";
  EXPECT_EQ(cover.metadata.samples(), paths_per_component) << "Wrong number of samples in the metadata";
  EXPECT_EQ(cover.metadata.contigs(), this->components) << "Wrong number of contigs in the metadata";
  EXPECT_EQ(cover.metadata.haplotypes(), expected_paths) << "Wrong number of haplotypes in the metadata";
  EXPECT_TRUE(cover.metadata.hasPathNames()) << "No path names in the metadata";
}

//------------------------------------------------------------------------------

} // namespace
