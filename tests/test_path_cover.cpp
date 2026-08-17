#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/gfa.h>
#include <gbwtgraph/path_cover.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

// Included reference paths for each component from default.gfa.
std::vector<std::vector<gbwt::vector_type>> generic_paths =
{
  {
    {
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(12, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(15, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(17, false))
    }
  },
  {
    {
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(22, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(24, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(25, false))
    }
  }
};

// Included reference paths for each component from components_ref.gfa.
std::vector<std::vector<gbwt::vector_type>> reference_paths =
{
  {
    {
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(13, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(15, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(17, false))
    }
  },
  {
    {
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(23, false)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(24, true)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(22, true)),
      static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, true))
    }
  }
};

// (n = 4, k = 3)-paths generated for components*.gfa / default.gfa.
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

std::string sample_name(size_t i)
{
  return "path_cover_" + std::to_string(i);
}

void check_metadata(
  const gbwt::GBWT& cover,
  const PathCoverParameters& params,
  const std::vector<std::string>& ref_samples,
  const std::vector<std::string>& contig_names,
  const std::vector<std::vector<gbwt::vector_type>>& included_paths)
{
  bool include_paths = !(included_paths.empty());
  ASSERT_TRUE(cover.hasMetadata()) << "Path cover GBWT contains no metadata";

  gbwt::size_type expected_samples = params.num_paths + ref_samples.size();
  EXPECT_EQ(cover.metadata.samples(), expected_samples) << "Wrong number of samples in the metadata";
  ASSERT_TRUE(cover.metadata.hasSampleNames()) << "No sample names in the metadata";
  for(size_t i = 0; i < ref_samples.size(); i++)
  {
    EXPECT_EQ(cover.metadata.sample(i), ref_samples[i]) << "Wrong reference sample name " << i << " in the metadata";
  }
  for(size_t i = 0; i < params.num_paths; i++)
  {
    EXPECT_EQ(cover.metadata.sample(i + ref_samples.size()), sample_name(i)) << "Wrong path cover sample name " << i << " in the metadata";
  }

  EXPECT_EQ(cover.metadata.contigs(), contig_names.size()) << "Wrong number of contigs in the metadata";
  ASSERT_TRUE(cover.metadata.hasContigNames()) << "No contig names in the metadata";
  for(size_t i = 0; i < contig_names.size(); i++)
  {
    EXPECT_EQ(cover.metadata.contig(i), contig_names[i]) << "Wrong contig name " << i << " in the metadata";
  }

  EXPECT_EQ(cover.metadata.haplotypes(), params.num_paths + ref_samples.size()) << "Wrong number of haplotypes in the metadata";

  gbwt::size_type expected_paths = contig_names.size() * params.num_paths;
  if(include_paths)
  {
    for(size_t i = 0; i < included_paths.size(); i++)
    {
      expected_paths += included_paths[i].size();
    }
  }
  EXPECT_EQ(cover.metadata.paths(), expected_paths) << "Wrong number of path names in the metadata";
  ASSERT_TRUE(cover.metadata.hasPathNames()) << "No path names in the metadata";
  size_t path_offset = 0;
  for(size_t i = 0; i < contig_names.size(); i++)
  {
    if(include_paths)
    {
      for(size_t j = 0; j < included_paths[i].size(); j++)
      {
        size_t path_id = path_offset + j;
        gbwt::PathName path_name = cover.metadata.path(path_id);
        size_t phase = (ref_samples[j] == GENERIC_PATH_SAMPLE_NAME ? GBWTGraph::NO_PHASE : 0);
        gbwt::PathName correct =
        {
          gbwt::PathName::path_name_type(j),
          gbwt::PathName::path_name_type(i),
          gbwt::PathName::path_name_type(phase),
          0
        };
        EXPECT_EQ(path_name, correct) << "Wrong included path name " << path_id << " in the metadata";
      }
      path_offset += included_paths[i].size();
    }
    for(size_t j = 0; j < params.num_paths; j++)
    {
      size_t path_id = path_offset + j;
      gbwt::PathName path_name = cover.metadata.path(path_id);
      gbwt::PathName correct =
      {
        gbwt::PathName::path_name_type(ref_samples.size() + j),
        gbwt::PathName::path_name_type(i),
        0, 0
      };
      EXPECT_EQ(path_name, correct) << "Wrong path cover path name " << path_id << " in the metadata";
    }
    path_offset += params.num_paths;
  }
}

//------------------------------------------------------------------------------

class PathStorageTest : public ::testing::Test
{
public:
  gbwt::GBWT index1;
  GBWTGraph graph1;
  
  gbwt::GBWT index2;
  GBWTGraph graph2;

  NaiveGraph source;
  size_t node_width;
  
  std::unordered_map<std::string, PathSense> paths1 {
    {"GRCh38#0#chr1", PathSense::REFERENCE},
    {"GRCh37#0#chr1", PathSense::REFERENCE},
    {"coolgene", PathSense::GENERIC},
    {"sample1#1#chr1#0", PathSense::HAPLOTYPE},
    {"sample1#2#chr1#0", PathSense::HAPLOTYPE},
    {"CHM13#0#chr1#0", PathSense::HAPLOTYPE}
  };
  
  std::unordered_map<std::string, PathSense> paths2 {
    {"GRCh38#1#chr1", PathSense::REFERENCE},
    {"GRCh37#2#chr1", PathSense::REFERENCE},
    {"coolergene", PathSense::GENERIC},
    {"sample1#1#chr2#0", PathSense::HAPLOTYPE},
    {"sample1#2#chr2#0", PathSense::HAPLOTYPE},
    {"CHM13v2#0#chr1#0", PathSense::HAPLOTYPE}
  };
  
  std::unordered_set<PathSense> all_senses {PathSense::GENERIC, PathSense::REFERENCE, PathSense::HAPLOTYPE};
  
  std::unordered_set<PathSense> named_senses {PathSense::GENERIC, PathSense::REFERENCE};

  PathStorageTest()
  {
  }

  void SetUp() override
  {
    // Need to parse PanSN
    GFAParsingParameters parameters;
    parameters.path_name_formats.emplace_front(
      GFAParsingParameters::PAN_SN_REGEX,
      GFAParsingParameters::PAN_SN_FIELDS,
      GFAParsingParameters::PAN_SN_SENSE
    );
    
    auto gfa_parse1 = gfa_to_gbwt("gfas/example_reference.gfa", parameters);
    this->index1 = *(gfa_parse1.first);
    this->source = *gfa_parse1.second;
    this->graph1 = GBWTGraph(this->index1, this->source);
    this->node_width = gbwt::bit_length(this->index1.sigma() - 1);
    
    // Grab another graph with different paths but (we assume) the same node ID space.
    auto gfa_parse2 = gfa_to_gbwt("gfas/example_more_reference.gfa", parameters);
    this->index2 = *(gfa_parse2.first);
    this->graph2 = GBWTGraph(this->index2, this->source);
  }
};

// Make sure the right senses of the right paths are present.
void 
check_stored_paths(
  const GBWTGraph& constructed, 
  const std::vector<const std::unordered_map<std::string, PathSense>*>& path_lists,
  const std::unordered_set<PathSense>& wanted_senses,
  const std::unordered_set<std::string>& unwanted_names
)
{
  for(auto path_list : path_lists)
  {
    for(auto& kv : *path_list)
    {
      if(wanted_senses.count(kv.second) && !unwanted_names.count(kv.first))
      {
        // Should have been copied over
        EXPECT_TRUE(constructed.has_path(kv.first)) << "Wanted path not copied: " << kv.first;
        if(constructed.has_path(kv.first))
        {
          EXPECT_EQ(constructed.get_sense(constructed.get_path_handle(kv.first)), kv.second) << "Wrong path sense: " << kv.first;
        }
      }
      else
      {
        // Should not have been copied over
        EXPECT_FALSE(constructed.has_path(kv.first)) << "Unwanted path copied: " << kv.first;
      }
    }
  }
}

TEST_F(PathStorageTest, StoreNamedPathsOneGraph)
{
  gbwt::GBWTBuilder builder(this->node_width);
  builder.index.addMetadata();
  store_named_paths(builder, this->graph1, nullptr);
  builder.finish();
  // Static-ify the GBWT so it doesn't happen in a temporary that GBWTGraph will keep a pointer to.
  gbwt::GBWT built(builder.index);
  GBWTGraph constructed(built, this->source);
  
  ASSERT_TRUE(constructed.index->hasMetadata()) << "Index missing metadata";
  ASSERT_TRUE(constructed.index->metadata.hasPathNames()) << "Index missing path names";
  ASSERT_TRUE(constructed.index->metadata.hasSampleNames()) << "Index missing sample names";
  ASSERT_TRUE(constructed.index->metadata.hasContigNames()) << "Index missing contig names";
  
  EXPECT_EQ(constructed.index->metadata.sample_names.size(), (gbwt::size_type) 3) << "Index has wrong number of samples";
  EXPECT_EQ(constructed.index->metadata.contig_names.size(), (gbwt::size_type) 2) << "Index has wrong number of contigs";
  EXPECT_LT(constructed.index->metadata.sample(GENERIC_PATH_SAMPLE_NAME), constructed.index->metadata.sample_names.size()) << "Index is missing generic path sample";
  EXPECT_LT(constructed.index->metadata.contig("chr1"), constructed.index->metadata.contig_names.size()) << "Index is missing chr1 contig";
  
  check_stored_paths(constructed, {&this->paths1}, this->named_senses, {});
}

TEST_F(PathStorageTest, StoreNamedPathsTwoGraphs)
{
  gbwt::GBWTBuilder builder(this->node_width);
  builder.index.addMetadata();
  store_named_paths(builder, this->graph1, nullptr);
  store_named_paths(builder, this->graph2, nullptr);
  builder.finish();
  // Static-ify the GBWT so it doesn't happen in a temporary that GBWTGraph will keep a pointer to.
  gbwt::GBWT built(builder.index);
  GBWTGraph constructed(built, this->source);
  
  EXPECT_EQ(constructed.index->metadata.sample_names.size(), (gbwt::size_type) 3) << "Index has wrong number of samples";
  EXPECT_EQ(constructed.index->metadata.contig_names.size(), (gbwt::size_type) 3) << "Index has wrong number of contigs";
  EXPECT_LT(constructed.index->metadata.contig("coolergene"), constructed.index->metadata.contig_names.size()) << "Index is missing coolergene contig";
  
  check_stored_paths(constructed, {&this->paths1, &this->paths2}, this->named_senses, {});
}

TEST_F(PathStorageTest, StoreAllPathsTwoGraphs)
{
  gbwt::GBWTBuilder builder(this->node_width);
  builder.index.addMetadata();
  store_paths(builder, this->graph1, this->all_senses, nullptr);
  store_paths(builder, this->graph2, this->all_senses, nullptr);
  builder.finish();
  // Static-ify the GBWT so it doesn't happen in a temporary that GBWTGraph will keep a pointer to.
  gbwt::GBWT built(builder.index);
  GBWTGraph constructed(built, this->source);
  
  EXPECT_EQ(constructed.index->metadata.sample_names.size(), (gbwt::size_type) 6) << "Index has wrong number of samples";
  EXPECT_EQ(constructed.index->metadata.contig_names.size(), (gbwt::size_type) 4) << "Index has wrong number of contigs";
  EXPECT_LT(constructed.index->metadata.sample("CHM13v2"), constructed.index->metadata.sample_names.size()) << "Index is missing CHM13v2 sample";
  EXPECT_LT(constructed.index->metadata.contig("chr2"), constructed.index->metadata.contig_names.size()) << "Index is missing chr2 contig";
  
  check_stored_paths(constructed, {&this->paths1, &this->paths2}, this->all_senses, {});
}

TEST_F(PathStorageTest, StoreAllPathsExceptTwoGraphs)
{
  std::unordered_set<std::string> unwanted_names {"coolgene", "GRCh38#0#chr1"};
  std::function<bool(const path_handle_t&)> filter1 = [&](const path_handle_t& path_handle)
  {
    return !unwanted_names.count(this->graph1.get_path_name(path_handle));
  };
  
  gbwt::GBWTBuilder builder(this->node_width);
  builder.index.addMetadata();
  store_paths(builder, this->graph1, this->all_senses, &filter1);
  store_paths(builder, this->graph2, this->all_senses, nullptr);
  builder.finish();
  // Static-ify the GBWT so it doesn't happen in a temporary that GBWTGraph will keep a pointer to.
  gbwt::GBWT built(builder.index);
  GBWTGraph constructed(built, this->source);
  
  check_stored_paths(constructed, {&this->paths1, &this->paths2}, this->all_senses, unwanted_names);
}

//------------------------------------------------------------------------------

class PathCoverTest : public ::testing::Test
{
public:
  size_t components;

  PathCoverTest()
  {
  }

  void SetUp() override
  {
    this->components = correct_paths.size();
  }

  void check_paths(
    const gbwt::GBWT& cover,
    const PathCoverParameters& params,
    const std::vector<std::vector<gbwt::vector_type>>& included_paths)
  {
    bool include_paths = !(included_paths.empty());

    gbwt::size_type expected_sequences = this->components * params.num_paths * 2;
    if(include_paths)
    {
      for(size_t i = 0; i < included_paths.size(); i++)
      {
        expected_sequences += 2 * included_paths[i].size();
      }
    }
    ASSERT_EQ(cover.sequences(), expected_sequences) << "Wrong number of sequences in the path cover GBWT";

    // Check that the included paths are there.
    size_t path_offset = 0;
    if(include_paths)
    {
      for(size_t i = 0; i < this->components; i++)
      {
        for(size_t j = 0; j < included_paths[i].size(); j++)
        {
          gbwt::size_type seq_id = gbwt::Path::encode(path_offset + j, false);
          gbwt::vector_type path = cover.extract(seq_id);
          EXPECT_EQ(path, included_paths[i][j]) << "Wrong included path " << j << " in component " << i;
        }
        path_offset += included_paths[i].size() + params.num_paths;
      }
    }

    // Check the paths in the path cover. We insert the smaller of a path and
    // its reverse complement to handle paths that flip the orientation.
    std::vector<std::set<gbwt::vector_type>> result(this->components);
    path_offset = 0;
    for(size_t i = 0; i < this->components; i++)
    {
      if(include_paths) { path_offset += included_paths[i].size(); }
      for(size_t j = 0; j < params.num_paths; j++)
      {
        size_t seq_id = gbwt::Path::encode(path_offset + j, false);
        gbwt::vector_type forward = cover.extract(seq_id), reverse;
        gbwt::reversePath(forward, reverse);
        result[i].insert(std::min(forward, reverse));
      }
      path_offset += params.num_paths;
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
};

TEST_F(PathCoverTest, SingleThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  gbwt::GBWT cover = path_cover_gbwt(graph, params);
  this->check_paths(cover, params, {});
  check_metadata(cover, params, {}, { "A1", "B1" }, {});
}

TEST_F(PathCoverTest, MultiThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  params.parallel_jobs = 2;
  gbwt::GBWT cover = path_cover_gbwt(graph, params);
  this->check_paths(cover, params, {});
  check_metadata(cover, params, {}, { "A1", "B1" }, {});
}

TEST_F(PathCoverTest, IncludeGenericSingleThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/default.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  gbwt::GBWT cover = path_cover_gbwt(graph, params, true);
  this->check_paths(cover, params, generic_paths);
  check_metadata(cover, params, { GENERIC_PATH_SAMPLE_NAME }, { "A", "B" }, generic_paths);
}

TEST_F(PathCoverTest, IncludeGenericMultiThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/default.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  params.parallel_jobs = 2;
  gbwt::GBWT cover = path_cover_gbwt(graph, params, true);
  this->check_paths(cover, params, generic_paths);
  check_metadata(cover, params, { GENERIC_PATH_SAMPLE_NAME }, { "A", "B" }, generic_paths);
}

TEST_F(PathCoverTest, IncludeReferenceSingleThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components_ref.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  gbwt::GBWT cover = path_cover_gbwt(graph, params, true);
  this->check_paths(cover, params, reference_paths);
  check_metadata(cover, params, { "ref" }, { "A", "B" }, reference_paths);
}

TEST_F(PathCoverTest, IncludeReferenceMultiThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components_ref.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  params.parallel_jobs = 2;
  gbwt::GBWT cover = path_cover_gbwt(graph, params, true);
  this->check_paths(cover, params, reference_paths);
  check_metadata(cover, params, { "ref" }, { "A", "B" }, reference_paths);
}

//------------------------------------------------------------------------------

class LocalHaplotypesTest : public ::testing::Test
{
public:
  size_t components;

  LocalHaplotypesTest()
  {
  }

  void SetUp() override
  {
    this->components = correct_paths.size();
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

  void check_cover(
    const gbwt::GBWT& index,
    const gbwt::GBWT& cover,
    const PathCoverParameters& params,
    const std::vector<std::vector<gbwt::vector_type>>& included_paths)
  {
    gbwt::size_type expected_sequences = this->components * params.num_paths * 2;
    for(size_t i = 0; i < included_paths.size(); i++)
    {
      expected_sequences += 2 * included_paths[i].size();
    }

    ASSERT_EQ(cover.sequences(), expected_sequences) << "Wrong number of sequences in the local haplotype GBWT";
    ASSERT_EQ(cover.sigma(), index.sigma()) << "Wrong alphabet size in the local haplotype GBWT";
    ASSERT_EQ(cover.effective(), index.effective()) << "Wrong effective alphabet size in the local haplotype GBWT";

    bool all_correct_subpaths = this->same_extensions(index, cover, params.context);
    EXPECT_TRUE(all_correct_subpaths) << "Missing " << params.context << "-subpaths in the local haplotype GBWT";
    bool no_extra_subpaths = this->same_extensions(cover, index, params.context);
    EXPECT_TRUE(no_extra_subpaths) << "Additional " << params.context << "-subpaths in the local haplotype GBWT";
  }
};

TEST_F(LocalHaplotypesTest, SingleThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  gbwt::GBWT cover = local_haplotypes(graph, index, params);
  this->check_cover(index, cover, params, {});
  check_metadata(cover, params, {}, { "A1", "B1" }, {});
}

TEST_F(LocalHaplotypesTest, MultiThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  params.parallel_jobs = 2;
  gbwt::GBWT cover = local_haplotypes(graph, index, params);
  this->check_cover(index, cover, params, {});
  check_metadata(cover, params, {}, { "A1", "B1" }, {});
}

TEST_F(LocalHaplotypesTest, IncludeGenericSingleThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/default.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  gbwt::GBWT cover = local_haplotypes(graph, index, params, true);
  this->check_cover(index, cover, params, generic_paths);
  check_metadata(cover, params, { GENERIC_PATH_SAMPLE_NAME }, { "A", "B" }, generic_paths);
}

TEST_F(LocalHaplotypesTest, IncludeGenericMultiThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/default.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  params.parallel_jobs = 2;
  gbwt::GBWT cover = local_haplotypes(graph, index, params, true);
  this->check_cover(index, cover, params, generic_paths);
  check_metadata(cover, params, { GENERIC_PATH_SAMPLE_NAME }, { "A", "B" }, generic_paths);
}

TEST_F(LocalHaplotypesTest, IncludeReferenceSingleThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components_ref.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  gbwt::GBWT cover = local_haplotypes(graph, index, params, true);
  this->check_cover(index, cover, params, reference_paths);
  check_metadata(cover, params, { "ref" }, { "A", "B" }, reference_paths);
}

TEST_F(LocalHaplotypesTest, IncludeReferenceMultiThreaded)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components_ref.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  params.parallel_jobs = 2;
  gbwt::GBWT cover = local_haplotypes(graph, index, params, true);
  this->check_cover(index, cover, params, reference_paths);
  check_metadata(cover, params, { "ref" }, { "A", "B" }, reference_paths);
}

TEST_F(LocalHaplotypesTest, Frequencies)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
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

  gbwt::GBWT cover = local_haplotypes(graph, index, params);
  gbwt::SearchState frequent_state = cover.find(frequent_path.begin(), frequent_path.end());
  gbwt::SearchState rare_state = cover.find(rare_path.begin(), rare_path.end());
  EXPECT_GE(frequent_state.size(), rare_state.size()) << "Local haplotype frequencies do not reflect true frequencies";
}

TEST_F(LocalHaplotypesTest, RevertToPathCover)
{
  auto gfa_parse = gfa_to_gbwt("gfas/components.gfa");
  const gbwt::GBWT& index = *(gfa_parse.first);
  GBWTGraph graph(index, *(gfa_parse.second));

  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  gbwt::size_type expected_sequences = this->components * params.num_paths * 2;

  gbwt::GBWT haplotype_cover = local_haplotypes(graph, index, params);
  ASSERT_EQ(haplotype_cover.sequences(), expected_sequences) << "Wrong number of sequences in the local haplotype GBWT";
  gbwt::GBWT path_cover = path_cover_gbwt(graph, params);
  ASSERT_EQ(path_cover.sequences(), expected_sequences) << "Wrong number of sequences in the path cover GBWT";

  auto first_only = gfa_to_gbwt("gfas/components_first.gfa");
  gbwt::GBWT mixed_cover = local_haplotypes(graph, *(first_only.first), params);
  ASSERT_EQ(mixed_cover.sequences(), expected_sequences) << "Wrong number of sequences in the mixed cover GBWT";

  // For the first component, we should have the same paths as with local haplotypes.
  for(size_t i = 0; i < params.num_paths; i++)
  {
    gbwt::vector_type path = mixed_cover.extract(gbwt::Path::encode(i, false));
    gbwt::vector_type correct_path = haplotype_cover.extract(gbwt::Path::encode(i, false));
    EXPECT_EQ(path, correct_path) << "Wrong path " << i << " in the first component";
  }

  // For the second component, we should have the same paths as with path cover.
  for(size_t i = 0; i < params.num_paths; i++)
  {
    gbwt::vector_type path = mixed_cover.extract(gbwt::Path::encode(params.num_paths + i, false));
    gbwt::vector_type correct_path = path_cover.extract(gbwt::Path::encode(params.num_paths + i, false));
    EXPECT_EQ(path, correct_path) << "Wrong path " << i << " in the second component";
  }
}

//------------------------------------------------------------------------------

class AugmentTest : public ::testing::Test
{
public:
  gbwt::GBWT index;
  GBWTGraph graph;
  std::vector<std::vector<nid_t>> components;
  size_t samples;

  AugmentTest()
  {
  }

  void SetUp() override
  {
    auto gfa_parse = gfa_to_gbwt("gfas/components.gfa");
    this->index = *(gfa_parse.first);
    this->graph = GBWTGraph(this->index, *(gfa_parse.second));
    this->components = weakly_connected_components(this->graph);
    this->samples = 2;
  }

  gbwt::DynamicGBWT create_gbwt(const std::set<size_t>& components_present) const
  {
    size_t node_width = gbwt::bit_length(this->index.sigma() - 1);
    size_t total_length = this->index.size();
    gbwt::GBWTBuilder builder(node_width, total_length);
    builder.index.addMetadata();

    std::vector<size_t> samples_per_component(this->components.size(), 0);
    std::vector<size_t> component_to_rank(this->components.size(), 0);
    for(size_t i = 0, rank = 0; i < this->components.size(); i++)
    {
      component_to_rank[i] = rank;
      if(components_present.find(i) != components_present.end()) { rank++; }
    }

    for(gbwt::size_type i = 0; i < this->index.sequences(); i += 2)
    {
      gbwt::vector_type sequence = this->index.extract(i);
      nid_t start_node = gbwt::Node::id(sequence.front());
      size_t component = 0;
      while(component < this->components.size())
      {
        const std::vector<nid_t>& curr = this->components[component];
        if(std::find(curr.begin(), curr.end(), start_node) != curr.end()) { break; }
        component++;
      }
      if(components_present.find(component) != components_present.end())
      {
        builder.insert(sequence, true);
        builder.index.metadata.addPath(samples_per_component[component], component_to_rank[component], 0, 0);
        samples_per_component[component]++;
      }
    }

    builder.finish();
    if(!(components_present.empty()))
    {
      std::vector<std::string> sample_names, contig_names;
      for(size_t i = 0; i < this->samples; i++)
      {
        sample_names.push_back("sample_" + std::to_string(i));
      }
      builder.index.metadata.setSamples(sample_names);
      for(size_t contig : components_present)
      {
        contig_names.push_back("contig_" + std::to_string(contig));
      }
      builder.index.metadata.setContigs(contig_names);
      builder.index.metadata.setHaplotypes(this->samples);
    }
    return builder.index;
  }

  std::string set_name(const std::set<size_t>& components) const
  {
    std::stringstream ss;
    ss << "{";
    for(size_t i : components) { ss << " " << i; }
    ss << " }";
    return ss.str();
  }
};

TEST_F(AugmentTest, CorrectPaths)
{
  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  std::vector<std::set<size_t>> component_sets
  {
    { }, { 0 }, { 1 }, { 0, 1 }
  };

  for(const std::set<size_t>& components_present : component_sets)
  {
    gbwt::DynamicGBWT original = this->create_gbwt(components_present);
    gbwt::DynamicGBWT augmented = original;
    std::string name = this->set_name(components_present);

    // Augment the GBWT.
    size_t expected_components = this->components.size() - components_present.size();
    size_t covered_components = augment_gbwt(this->graph, augmented, params);
    ASSERT_EQ(covered_components, expected_components) << "Wrong number of covered components for components " << name;
    size_t expected_sequences = original.sequences() + 2 * expected_components * params.num_paths;
    ASSERT_EQ(augmented.sequences(), expected_sequences) << "Wrong number of sequences for components " << name;

    // Check original paths.
    size_t path_id = 0;
    while(path_id < original.metadata.paths())
    {
      size_t i = gbwt::Path::encode(path_id, false);
      bool same_paths = (original.extract(i) == augmented.extract(i));
      EXPECT_TRUE(same_paths) << "Wrong original path " << path_id << " for components " << name;
      path_id++;
    }

    // Check generated paths.
    for(size_t component = 0; component < this->components.size(); component++)
    {
      if(components_present.find(component) != components_present.end()) { continue; }
      for(size_t i = 0; i < params.num_paths; i++)
      {
        gbwt::vector_type forward = augmented.extract(gbwt::Path::encode(path_id, false)), reverse;
        gbwt::reversePath(forward, reverse);
        gbwt::vector_type path = std::min(forward, reverse);
        bool correct_path = (correct_paths[component].find(path) != correct_paths[component].end());
        EXPECT_TRUE(correct_path) << "Wrong augmented path " << i << " in component " << component << " for components " << name;
        path_id++;
      }
    }
  }
}

TEST_F(AugmentTest, PathNames)
{
  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  std::vector<std::set<size_t>> component_sets
  {
    { }, { 0 }, { 1 }, { 0, 1 }
  };

  for(const std::set<size_t>& components_present : component_sets)
  {
    gbwt::DynamicGBWT original = this->create_gbwt(components_present);
    gbwt::DynamicGBWT augmented = original;
    std::string name = this->set_name(components_present);

    // Augment the GBWT.
    size_t expected_components = this->components.size() - components_present.size();
    augment_gbwt(this->graph, augmented, params);

    // Check metadata.
    ASSERT_TRUE(augmented.hasMetadata()) << "No metadata for components " << name;
    ASSERT_TRUE(augmented.metadata.hasPathNames()) << "No path names for components " << name;
    size_t expected_paths = original.metadata.paths() + expected_components * params.num_paths;
    ASSERT_EQ(augmented.metadata.paths(), expected_paths) << "Wrong number of path names for components " << name;

    // Check original paths.
    size_t path_id = 0;
    while(path_id < original.metadata.paths())
    {
      bool same_names = (original.metadata.path(path_id) == augmented.metadata.path(path_id));
      EXPECT_TRUE(same_names) << "Wrong original path name " << path_id << " for components " << name;
      path_id++;
    }

    // Check generated paths.
    size_t rank = 0;
    for(size_t component = 0; component < this->components.size(); component++)
    {
      if(components_present.find(component) != components_present.end()) { continue; }
      for(size_t i = 0; i < params.num_paths; i++)
      {
        gbwt::PathName path_name = augmented.metadata.path(path_id);
        EXPECT_EQ(path_name.sample, original.metadata.samples() + i) << "Wrong sample for augmented path " << i << " in component " << component << " for components " << name;
        EXPECT_EQ(path_name.contig, components_present.size() + rank) << "Wrong contig for augmented path " << i << " in component " << component << " for components " << name;
        path_id++;
      }
      rank++;
    }
  }
}

TEST_F(AugmentTest, SamplesAndContigs)
{
  PathCoverParameters params;
  params.num_paths = 4; params.context = 3;
  std::vector<std::set<size_t>> component_sets
  {
    { }, { 0 }, { 1 }, { 0, 1 }
  };

  for(const std::set<size_t>& components_present : component_sets)
  {
    gbwt::DynamicGBWT original = this->create_gbwt(components_present);
    gbwt::DynamicGBWT augmented = original;
    std::string name = this->set_name(components_present);

    // Augment the GBWT.
    size_t expected_components = this->components.size() - components_present.size();
    augment_gbwt(this->graph, augmented, params);

    // Check metadata.
    ASSERT_TRUE(augmented.hasMetadata()) << "No metadata for components " << name;
    ASSERT_TRUE(augmented.metadata.hasSampleNames()) << "No sample names for components " << name;
    size_t expected_samples = original.metadata.samples() + (expected_components == 0 ? 0 : params.num_paths);
    ASSERT_EQ(augmented.metadata.samples(), expected_samples) << "Wrong number of samples for components " << name;
    ASSERT_TRUE(augmented.metadata.hasContigNames()) << "No contig names for components " << name;
    ASSERT_EQ(augmented.metadata.contigs(), this->components.size()) << "Wrong number of contigs for components " << name;

    // Check sample names.
    size_t sample_id = 0;
    while(sample_id < original.metadata.samples())
    {
      std::string expected_name = original.metadata.sample(sample_id);
      EXPECT_EQ(augmented.metadata.sample(sample_id), expected_name) << "Wrong name for original sample " << sample_id << " for components " << name;
      sample_id++;
    }
    if(expected_components > 0)
    {
      for(size_t i = 0; i < params.num_paths; i++)
      {
        std::string expected_name = "path_cover_" + std::to_string(i);
        EXPECT_EQ(augmented.metadata.sample(sample_id), expected_name) << "Wrong name for augmented sample " << i << " for components " << name;
        sample_id++;
      }
    }

    // Check contig names.
    size_t contig_id = 0;
    for(size_t i = 0; i < this->components.size(); i++)
    {
      if(components_present.find(i) != components_present.end())
      {
        std::string expected_name = "contig_" + std::to_string(i);
        EXPECT_EQ(augmented.metadata.contig(contig_id), expected_name) << "Wrong name for original contig " << i << " for components " << name;
        contig_id++;
      }
    }
    if(expected_components > 0)
    {
      for(size_t i = 0; i < this->components.size(); i++)
      {
        if(components_present.find(i) == components_present.end())
        {
          std::string expected_name = "component_" + std::to_string(i);
          EXPECT_EQ(augmented.metadata.contig(contig_id), expected_name) << "Wrong name for augmented contig " << i << " for components " << name;
          contig_id++;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

} // namespace
