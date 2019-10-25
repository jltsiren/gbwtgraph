#include <gtest/gtest.h>

#include <set>
#include <stack>
#include <vector>

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
  EXPECT_EQ(cover.metadata.haplotypes(), paths_per_component) << "Wrong number of haplotypes in the metadata";
  EXPECT_TRUE(cover.metadata.hasPathNames()) << "No path names in the metadata";
  EXPECT_EQ(cover.metadata.paths(), expected_paths) << "Wrong number of path names in the metadata";
}

//------------------------------------------------------------------------------

class LocalHaplotypesTest : public ::testing::Test
{
public:
  gbwt::GBWT index;
  GBWTGraph graph;
  size_t components;

  LocalHaplotypesTest()
  {
  }

  void SetUp() override
  {
    auto gfa_parse = gfa_to_gbwt("components.gfa");
    this->index = *(gfa_parse.first);
    this->graph = GBWTGraph(this->index, *(gfa_parse.second));
    this->components = 2;
  }

  struct SearchState
  {
    gbwt::SearchState left, right;
    size_t depth;

    SearchState successor(const gbwt::GBWT& left_index, const gbwt::GBWT& right_index, gbwt::node_type to)
    {
      return { left_index.extend(this->left, to), right_index.extend(this->right, to), static_cast<size_t>(this->depth + 1) };
    }

    bool ok() const { return (this->left.empty() == this->right.empty()); }

    gbwt::node_type node() const { return this->left.node; }
  };

  // Does the second index have the same depth-k extensions as the first?
  bool same_extensions(const gbwt::GBWT& baseline, const gbwt::GBWT& candidate, size_t k)
  {
    for(gbwt::node_type node = baseline.firstNode(); node < baseline.sigma(); node++)
    {
      std::stack<SearchState> states;
      states.push({ baseline.find(node), candidate.find(node), static_cast<size_t>(1)});
      while(!(states.empty()))
      {
        SearchState state = states.top(); states.pop();
        if(!(state.ok())) { return false; }
        if(state.depth >= k) { continue; }
        for(auto edge : baseline.edges(state.node()))
        {
          if(edge.first == gbwt::ENDMARKER) { continue; }
          states.push(state.successor(baseline, candidate, edge.first));
        }
      }
    }
    return true;
  }
};

TEST_F(LocalHaplotypesTest, CorrectSubPaths)
{
  size_t paths_per_component = 4;
  size_t context_length = 3;
  gbwt::size_type expected_sequences = this->components * paths_per_component * 2;

  gbwt::GBWT cover = local_haplotypes(this->graph, paths_per_component, context_length);
  ASSERT_EQ(cover.sequences(), expected_sequences) << "Wrong number of sequences in the local haplotype GBWT";
  ASSERT_EQ(cover.sigma(), this->index.sigma()) << "Wrong alphabet size in the local haplotype GBWT";
  ASSERT_EQ(cover.effective(), this->index.effective()) << "Wrong effective alphabet size in the local haplotype GBWT";

  bool all_correct_subpaths = this->same_extensions(this->index, cover, context_length);
  EXPECT_TRUE(all_correct_subpaths) << "Missing " << context_length << "-subpaths in the local haplotype GBWT";
  bool no_extra_subpaths = this->same_extensions(cover, this->index, context_length);
  EXPECT_TRUE(no_extra_subpaths) << "Additional " << context_length << "-subpaths in the local haplotype GBWT";
}

TEST_F(LocalHaplotypesTest, Metadata)
{
  size_t paths_per_component = 4;
  size_t context_length = 3;
  size_t expected_paths = paths_per_component * this->components;

  gbwt::GBWT cover = local_haplotypes(this->graph, paths_per_component, context_length);
  ASSERT_TRUE(cover.hasMetadata()) << "Local haplotype GBWT contains no metadata";
  EXPECT_EQ(cover.metadata.samples(), paths_per_component) << "Wrong number of samples in the metadata";
  EXPECT_EQ(cover.metadata.contigs(), this->components) << "Wrong number of contigs in the metadata";
  EXPECT_EQ(cover.metadata.haplotypes(), paths_per_component) << "Wrong number of haplotypes in the metadata";
  EXPECT_TRUE(cover.metadata.hasPathNames()) << "No path names in the metadata";
  EXPECT_EQ(cover.metadata.paths(), expected_paths) << "Wrong number of path names in the metadata";
}

TEST_F(LocalHaplotypesTest, Frequencies)
{
  size_t paths_per_component = 4;
  size_t context_length = 3;
  std::vector<gbwt::node_type> frequent_path
  {
    gbwt::Node::encode(21, false),
    gbwt::Node::encode(22, false),
    gbwt::Node::encode(24, false)
  };
  std::vector<gbwt::node_type> rare_path
  {
    gbwt::Node::encode(22, false),
    gbwt::Node::encode(24, false),
    gbwt::Node::encode(25, false)
  };

  gbwt::GBWT cover = local_haplotypes(this->graph, paths_per_component, context_length);
  gbwt::SearchState frequent_state = cover.find(frequent_path.begin(), frequent_path.end());
  gbwt::SearchState rare_state = cover.find(rare_path.begin(), rare_path.end());
  EXPECT_GE(frequent_state.size(), rare_state.size()) << "Local haplotype frequencies do not reflect true frequencies";
}

//------------------------------------------------------------------------------

} // namespace
