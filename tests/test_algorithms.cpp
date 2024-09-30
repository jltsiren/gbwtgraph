#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <unordered_set>
#include <vector>

#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/gfa.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class ComponentTest : public ::testing::Test
{
public:
  gbwt::GBWT index;
  GBWTGraph graph;
  size_t components;
  std::vector<std::set<gbwt::vector_type>> correct_paths;

  ComponentTest()
  {
  }

  void SetUp() override
  {
    auto gfa_parse = gfa_to_gbwt("gfas/components.gfa");
    this->index = *(gfa_parse.first);
    this->graph = GBWTGraph(this->index, *(gfa_parse.second));
    this->components = 2;
  }

  void check_path_assignments(
    const GBWTGraph& graph,
    const std::vector<std::vector<path_handle_t>>& result,
    const std::vector<std::vector<size_t>>& correct,
    const std::string& filename) const
  {
    ASSERT_EQ(result.size(), correct.size()) << "Wrong number of components in " << filename;
    for(size_t i = 0; i < result.size(); i++)
    {
      const std::vector<path_handle_t>& paths = result[i];
      const std::vector<size_t>& correct_paths = correct[i];
      ASSERT_EQ(paths.size(), correct_paths.size()) << "Wrong number of paths in component " << i << " in " << filename;
      for(size_t j = 0; j < paths.size(); j++)
      {
        size_t path_id = graph.handle_to_path(paths[j]);
        EXPECT_EQ(path_id, correct_paths[j]) << "Incorrect path " << j << " in component " << i << " in " << filename;
      }
    }
  }
};

TEST_F(ComponentTest, Components)
{
  std::vector<std::vector<nid_t>> correct_components =
  {
    { 11, 12, 13, 14, 15, 16, 17 },
    { 21, 22, 23, 24, 25 }
  };
  std::vector<std::vector<nid_t>> result = weakly_connected_components(this->graph);
  ASSERT_EQ(result.size(), correct_components.size()) << "Wrong number of components";

  for(size_t i = 0; i < result.size(); i++)
  {
    ASSERT_EQ(result[i].size(), correct_components[i].size()) << "Wrong number of nodes in component " << i;
    auto result_iter = result[i].begin();
    auto correct_iter = correct_components[i].begin();
    while(result_iter != result[i].end())
    {
      EXPECT_EQ(*result_iter, *correct_iter) << "Incorrect node in component " << i;
      ++result_iter; ++correct_iter;
    }
  }
}

TEST_F(ComponentTest, HeadNodes)
{
  std::vector<std::vector<nid_t>> correct_heads =
  {
    { 11 },
    { }
  };
  std::vector<std::vector<nid_t>> components = weakly_connected_components(this->graph);
  ASSERT_EQ(components.size(), correct_heads.size()) << "Wrong number of components";

  // For both components, add a version with a nonexistent node id.
  components.emplace_back(components[0]);
  components.back().push_back(42);
  correct_heads.emplace_back(correct_heads[0]);
  components.emplace_back(components[1]);
  components.back().push_back(42);
  correct_heads.emplace_back(correct_heads[1]);

  for(size_t i = 0; i < components.size(); i++)
  {
    std::vector<nid_t> heads = is_nice_and_acyclic(this->graph, components[i]);
    ASSERT_EQ(heads.size(), correct_heads[i].size()) << "Wrong number of head nodes in component " << i;
    auto result_iter = heads.begin();
    auto correct_iter = correct_heads[i].begin();
    while(result_iter != heads.end())
    {
      EXPECT_EQ(*result_iter, *correct_iter) << "Incorrect head node in component " << i;
      ++result_iter; ++correct_iter;
    }
  }
}

TEST_F(ComponentTest, ConstructionJobs)
{
  std::vector<std::vector<nid_t>> components = weakly_connected_components(this->graph);
  std::vector<std::pair<size_t, size_t>> bounds_and_jobs { { 0, 2 }, { 11, 2 }, { 12, 1 } };
  for(auto params : bounds_and_jobs)
  {
    ConstructionJobs jobs = gbwt_construction_jobs(this->graph, params.first);
    ASSERT_EQ(jobs.size(), params.second) << "Invalid number of jobs with size bound " << params.first;
    if(params.second == 1)
    {
      ASSERT_EQ(jobs.job_size(0), this->graph.get_node_count()) << "Invalid size for job 0 with size bound " << params.first;
    }
    ASSERT_EQ(jobs.components(), components.size()) << "Invalid number of components with size bound " << params.first;
    ASSERT_EQ(jobs.weakly_connected_components, components) << "Invalid components with size bound " << params.first;

    // For each component, check that the mappings for each node in the component are correct.
    for(size_t i = 0; i < components.size(); i++)
    {
      for(nid_t id : components[i])
      {
        EXPECT_EQ(jobs.component(id), i) << "Invalid component for node " << id << " with size bound " << params.first;
      }
      if(params.second == components.size())
      {
        ASSERT_EQ(jobs.job_size(i), components[i].size()) << "Invalid size for job " << i << " with size bound " << params.first;
      }
      size_t expected_job = (params.second == 1 ? 0 : i);
      for(nid_t id : components[i])
      {
        EXPECT_EQ(jobs.job(id), expected_job) << "Invalid job for node " << id << " with size bound " << params.first;
      }
    }

    // For each job, check that the mappings for each component in the job are correct.
    std::vector<std::vector<size_t>> components_per_job = jobs.components_per_job();
    for(size_t job_id = 0; job_id < components_per_job.size(); job_id++)
    {
      for(size_t component_id : components_per_job[job_id])
      {
        EXPECT_EQ(jobs.job_for_component(component_id), job_id) << "Invalid job for component " << component_id << " with size bound " << params.first;
      }
    }
  }
}

TEST_F(ComponentTest, ContigNames)
{
  std::vector<std::pair<std::string, std::vector<std::string>>> files_and_components =
  {
    { "gfas/components.gfa", { "A1", "B1" } },
    { "gfas/components_walks.gfa", { "component_0", "component_1" } },
    { "gfas/components_ref.gfa", { "A", "B" } },
    { "gfas/default.gfa", { "A", "B" } }
  };

  for(auto params : files_and_components)
  {
    auto gfa_parse = gfa_to_gbwt(params.first);
    const gbwt::GBWT& index = *(gfa_parse.first);
    GBWTGraph graph(index, *(gfa_parse.second));

    ConstructionJobs jobs = gbwt_construction_jobs(graph, 0);
    std::vector<std::string> names = jobs.contig_names(graph);
    ASSERT_EQ(names.size(), params.second.size()) << "Wrong number of contig names in " << params.first;
    for(size_t i = 0; i < names.size(); i++)
    {
      EXPECT_EQ(names[i], params.second[i]) << "Incorrect contig name " << i << " in " << params.first;
    }
  }
}
TEST_F(ComponentTest, ContigNamesWithFilter)
{
  std::vector<std::pair<std::string, std::vector<std::string>>> files_and_components =
  {
    { "gfas/components.gfa", { "A1", "B1" } },
    { "gfas/components_walks.gfa", { "component_0", "component_1" } },
    { "gfas/components_ref.gfa", { "component_0", "component_1" } },
    { "gfas/default.gfa", { "A", "B" } }
  };

  for(auto params : files_and_components)
  {
    auto gfa_parse = gfa_to_gbwt(params.first);
    const gbwt::GBWT& index = *(gfa_parse.first);
    GBWTGraph graph(index, *(gfa_parse.second));

    std::function<bool(const path_handle_t&)> generic_filter = [&](const path_handle_t& path) -> bool
    {
      return (graph.get_sense(path) == PathSense::GENERIC);
    };

    ConstructionJobs jobs = gbwt_construction_jobs(graph, 0);
    std::vector<std::string> names = jobs.contig_names(graph, generic_filter);
    ASSERT_EQ(names.size(), params.second.size()) << "Wrong number of contig names in " << params.first;
    for(size_t i = 0; i < names.size(); i++)
    {
      EXPECT_EQ(names[i], params.second[i]) << "Incorrect contig name " << i << " in " << params.first;
    }
  }
}

TEST_F(ComponentTest, AssignPaths)
{
  std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>> files_and_paths =
  {
    { "gfas/components.gfa", { { 0, 1 }, { 2, 3 } } },
    { "gfas/components_walks.gfa", { { }, { } } },
    { "gfas/components_ref.gfa", { { 0 }, { 3 } } },
    { "gfas/default.gfa", { { 0 }, { 3 } } }
  };

  for(auto params : files_and_paths)
  {
    auto gfa_parse = gfa_to_gbwt(params.first);
    const gbwt::GBWT& index = *(gfa_parse.first);
    GBWTGraph graph(index, *(gfa_parse.second));

    ConstructionJobs jobs = gbwt_construction_jobs(graph, 0);
    auto result = assign_paths(graph, jobs, nullptr, nullptr);
    this->check_path_assignments(graph, result, params.second, params.first);
  }
}

TEST_F(ComponentTest, AssignPathsWithFilter)
{
  std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>> files_and_paths =
  {
    { "gfas/components.gfa", { { 0, 1 }, { 2, 3 } } },
    { "gfas/components_walks.gfa", { { }, { } } },
    { "gfas/components_ref.gfa", { { }, { } } },
    { "gfas/default.gfa", { { 0 }, { 3 } } }
  };

  for(auto params : files_and_paths)
  {
    auto gfa_parse = gfa_to_gbwt(params.first);
    const gbwt::GBWT& index = *(gfa_parse.first);
    GBWTGraph graph(index, *(gfa_parse.second));

    std::function<bool(const path_handle_t&)> generic_filter = [&](const path_handle_t& path) -> bool
    {
      return (graph.get_sense(path) == PathSense::GENERIC);
    };

    ConstructionJobs jobs = gbwt_construction_jobs(graph, 0);
    auto result = assign_paths(graph, jobs, nullptr, &generic_filter);
    this->check_path_assignments(graph, result, params.second, params.first);
  }
}

TEST_F(ComponentTest, InsertPaths)
{
  std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>> files_and_paths =
  {
    { "gfas/components.gfa", { { 0, 1 }, { 2, 3 } } },
    { "gfas/components_walks.gfa", { { }, { } } },
    { "gfas/components_ref.gfa", { { 0 }, { 3 } } },
    { "gfas/default.gfa", { { 0 }, { 3 } } }
  };

  for(auto params : files_and_paths)
  {
    auto gfa_parse = gfa_to_gbwt(params.first);
    const gbwt::GBWT& index = *(gfa_parse.first);
    GBWTGraph graph(index, *(gfa_parse.second));
    MetadataBuilder metadata_builder;

    // Assign paths to jobs.
    ConstructionJobs jobs = gbwt_construction_jobs(graph, 0);
    auto assigned = assign_paths(graph, jobs, &metadata_builder, nullptr);

    // Build a GBWT with the assigned paths.
    gbwt::size_type node_width = sdsl::bits::length(gbwt::Node::encode(graph.max_node_id(), true));
    gbwt::GBWTBuilder builder(node_width, index.size());
    for(size_t job = 0; job < assigned.size(); job++)
    {
      insert_paths(graph, assigned[job], builder, job, false);
    }
    builder.finish();
    gbwt::GBWT built = builder.index;

    // Add metadata.
    built.addMetadata();
    built.metadata = metadata_builder.get_metadata();
    if(index.tags.contains(REFERENCE_SAMPLE_LIST_GBWT_TAG))
    {
      built.tags.set(REFERENCE_SAMPLE_LIST_GBWT_TAG, index.tags.get(REFERENCE_SAMPLE_LIST_GBWT_TAG));
    }
    GBWTGraph built_graph(built, graph);

    // Check that the paths are correct.
    gbwt::size_type path_id = 0;
    for(size_t job = 0; job < assigned.size(); job++)
    {
      for(size_t i = 0; i < assigned[job].size(); i++)
      {
        path_handle_t handle = built_graph.path_to_handle(path_id);
        path_handle_t old_handle = assigned[job][i];
        EXPECT_EQ(built_graph.get_path_name(handle), graph.get_path_name(old_handle)) << "Incorrect path name for path " << path_id << " in job " << job;
        EXPECT_EQ(built_graph.get_sense(handle), graph.get_sense(old_handle)) << "Incorrect path sense for path " << path_id << " in job " << job;
        gbwt::vector_type built_path = built.extract(gbwt::Path::encode(path_id, false));
        gbwt::vector_type correct_path = index.extract(gbwt::Path::encode(graph.handle_to_path(old_handle), false));
        EXPECT_EQ(built_path, correct_path) << "Incorrect path sequence for path " << path_id << " in job " << job;
        path_id++;
      }
    }
  }
}

//------------------------------------------------------------------------------

class TopologicalOrderTest : public ::testing::Test
{
public:
  gbwt::GBWT index;
  GBWTGraph graph;

  TopologicalOrderTest()
  {
  }

  void SetUp() override
  {
    auto gfa_parse = gfa_to_gbwt("gfas/cyclic.gfa");
    this->index = *(gfa_parse.first);
    this->graph = GBWTGraph(this->index, *(gfa_parse.second));
  }

  void check_subgraph(const std::unordered_set<nid_t>& subgraph, bool acyclic) const
  {
    std::vector<handle_t> order = topological_order(this->graph, subgraph);
    if(!acyclic)
    {
      ASSERT_TRUE(order.empty()) << "Non-empty order for a subgraph containing cycles";
      return;
    }

    // Determine the node ids that do not exist in the graph.
    size_t missing_nodes = 0;
    for(nid_t node : subgraph)
    {
      if(!(this->graph.has_node(node))) { missing_nodes++; }
    }

    ASSERT_EQ(order.size(), 2 * (subgraph.size() - missing_nodes)) << "Wrong number of handles in the order";
    for(nid_t node : subgraph)
    {
      if(!(this->graph.has_node(node))) { continue; }
      for(bool orientation : { false, true })
      {
        handle_t from = this->graph.get_handle(node, orientation);
        auto from_iter = std::find(order.begin(), order.end(), from);
        ASSERT_NE(from_iter, order.end()) << "Node " << node << ", orientation " << orientation << " not found in the order";
        bool ok = this->graph.follow_edges(from, false, [&](const handle_t& to) -> bool
        {
          if(subgraph.find(this->graph.get_id(to)) == subgraph.end()) { return true; }
          auto to_iter = std::find(order.begin(), order.end(), to);
          if(to_iter == order.end()) { return false; }
          return (from_iter < to_iter);
        });
        EXPECT_TRUE(ok) << "Constraints not satisfied for node " << node << ", orientation " << orientation;
      }
    }
  }
};

TEST_F(TopologicalOrderTest, SingleComponent)
{
  std::unordered_set<nid_t> subgraph =
  {
    static_cast<nid_t>(1),
    static_cast<nid_t>(2),
    static_cast<nid_t>(4),
    static_cast<nid_t>(5),
    static_cast<nid_t>(6)
  };
  this->check_subgraph(subgraph, true);
}

TEST_F(TopologicalOrderTest, TwoComponents)
{
  std::unordered_set<nid_t> subgraph =
  {
    static_cast<nid_t>(1),
    static_cast<nid_t>(2),
    static_cast<nid_t>(4),
    static_cast<nid_t>(6),
    static_cast<nid_t>(7),
    static_cast<nid_t>(8),
    static_cast<nid_t>(9)
  };
  this->check_subgraph(subgraph, true);
}

TEST_F(TopologicalOrderTest, CyclicComponent)
{
  std::unordered_set<nid_t> subgraph =
  {
    static_cast<nid_t>(2),
    static_cast<nid_t>(4),
    static_cast<nid_t>(5),
    static_cast<nid_t>(6),
    static_cast<nid_t>(8)
  };
  this->check_subgraph(subgraph, false);
}

TEST_F(TopologicalOrderTest, MissingNodes)
{
  std::unordered_set<nid_t> subgraph =
  {
    static_cast<nid_t>(1),
    static_cast<nid_t>(2),
    static_cast<nid_t>(4),
    static_cast<nid_t>(5),
    static_cast<nid_t>(6),
    static_cast<nid_t>(42)
  };
  this->check_subgraph(subgraph, true);
}

//------------------------------------------------------------------------------

class LCSTest : public ::testing::Test
{
public:
  gbwt::GBWT index;
  GBWTGraph graph;
  std::mt19937_64 rng;
  std::vector<nid_t> node_ids;

  LCSTest() :
    rng(0xACDCABBABABA ^ std::random_device()())
  {
  }

  void SetUp() override
  {
    // This graph contains nodes 1 to 8 with variable length sequences.
    auto gfa_parse = gfa_to_gbwt("gfas/for_subgraph.gfa");
    this->index = *(gfa_parse.first);
    this->graph = GBWTGraph(this->index, *(gfa_parse.second));
    this->graph.for_each_handle([&](const handle_t& handle)
    {
      this->node_ids.push_back(this->graph.get_id(handle));
    });
  }

  void compare_paths(const gbwt::vector_type& left, const gbwt::vector_type& right, const std::string& name)
  {
    std::vector<std::pair<size_t, size_t>> result = path_lcs(this->graph, left, right);
    for(size_t i = 0; i < result.size(); i++)
    {
      ASSERT_LT(result[i].first, left.size()) << "Invalid left LCS position " << i << " for " << name;
      ASSERT_LT(result[i].second, right.size()) << "Invalid right LCS position " << i << " for " << name;
      if(i > 0)
      {
        ASSERT_LT(result[i - 1].first, result[i].first) << "Non-increasing left LCS position " << i << " for " << name;
        ASSERT_LT(result[i - 1].second, result[i].second) << "Non-increasing right LCS position " << i << " for " << name;
      }
      ASSERT_EQ(left[result[i].first], right[result[i].second]) << "Incorrect LCS element " << i << " for " << name;
    }
  }

  gbwt::vector_type forward_path(const std::vector<nid_t>& nodes)
  {
    gbwt::vector_type result;
    for(nid_t node : nodes)
    {
      result.push_back(gbwt::Node::encode(node, false));
    }
    return result;
  }

  gbwt::vector_type random_path(size_t length)
  {
    gbwt::vector_type result;
    std::uniform_int_distribution<size_t> dist(0, this->node_ids.size() - 1);
    for(size_t i = 0; i < length; i++)
    {
      nid_t node = this->node_ids[dist(this->rng)];
      result.push_back(gbwt::Node::encode(node, false));
    }
    return result;
  }

  gbwt::vector_type mutate_path(const gbwt::vector_type& path, double rate)
  {
    gbwt::vector_type result = path;
    std::uniform_real_distribution<double> mutate(0.0, 1.0);
    std::uniform_int_distribution<size_t> node_dist(0, this->node_ids.size() - 1);
    for(size_t i = 0; i < result.size(); i++)
    {
      if(mutate(this->rng) < rate)
      {
        nid_t node = this->node_ids[node_dist(this->rng)];
        result[i] = gbwt::Node::encode(node, false);
      }
    }
    return result;
  }
};

TEST_F(LCSTest, EmptyLCS)
{
  gbwt::vector_type empty;
  gbwt::vector_type path = this->index.extract(gbwt::Path::encode(0, false));

  ASSERT_TRUE(path_lcs(this->graph, empty, empty).empty()) << "Nonempty LCS for two empty paths";
  ASSERT_TRUE(path_lcs(this->graph, empty, path).empty()) << "Nonempty LCS for an empty first path";
  ASSERT_TRUE(path_lcs(this->graph, path, empty).empty()) << "Nonempty LCS for an empty second path";

  gbwt::vector_type reverse = this->index.extract(gbwt::Path::encode(0, true));
  ASSERT_TRUE(path_lcs(this->graph, path, reverse).empty()) << "Nonempty LCS for non-overlapping paths";
}

TEST_F(LCSTest, IndenticalLCS)
{
  for(gbwt::size_type i = 0; i < this->index.sequences(); i++)
  {
    gbwt::vector_type path = this->index.extract(i);
    std::vector<std::pair<size_t, size_t>> truth;
    for(size_t j = 0; j < path.size(); j++)
    {
      truth.push_back(std::make_pair(j, j));
    }
    std::vector<std::pair<size_t, size_t>> result = path_lcs(this->graph, path, path);
    gbwt::size_type id = gbwt::Path::id(i);
    bool orientation = gbwt::Path::is_reverse(i);
    ASSERT_EQ(result, truth) << "Incorrect LCS for path " << id << ", orientation " << orientation << " with itself";
  }
}

TEST_F(LCSTest, FirstLastLCS)
{
  gbwt::vector_type path = this->forward_path({ 1, 2, 3, 4, 5 });

  gbwt::vector_type no_no = this->forward_path({ 6, 2, 4, 7 });
  std::vector<std::pair<size_t, size_t>> no_no_left = { { 1, 1 }, { 2, 3 } };
  ASSERT_EQ(path_lcs(this->graph, no_no, path), no_no_left) << "Incorrect LCS for no first, no last, left";
  std::vector<std::pair<size_t, size_t>> no_no_right = { { 1, 1 }, { 3, 2 } };
  ASSERT_EQ(path_lcs(this->graph, path, no_no), no_no_right) << "Incorrect LCS for no first, no last, right";

  gbwt::vector_type no_yes = this->forward_path({ 6, 2, 4, 5 });
  std::vector<std::pair<size_t, size_t>> no_yes_left = { { 1, 1 }, { 2, 3 }, { 3, 4 } };
  ASSERT_EQ(path_lcs(this->graph, no_yes, path), no_yes_left) << "Incorrect LCS for no first, yes last, left";
  std::vector<std::pair<size_t, size_t>> no_yes_right = { { 1, 1 }, { 3, 2 }, { 4, 3 } };
  ASSERT_EQ(path_lcs(this->graph, path, no_yes), no_yes_right) << "Incorrect LCS for no first, yes last, right";

  gbwt::vector_type yes_no = this->forward_path({ 1, 2, 4, 7 });
  std::vector<std::pair<size_t, size_t>> yes_no_left = { { 0, 0 }, { 1, 1 }, { 2, 3 } };
  ASSERT_EQ(path_lcs(this->graph, yes_no, path), yes_no_left) << "Incorrect LCS for yes first, no last, left";
  std::vector<std::pair<size_t, size_t>> yes_no_right = { { 0, 0 }, { 1, 1 }, { 3, 2 } };
  ASSERT_EQ(path_lcs(this->graph, path, yes_no), yes_no_right) << "Incorrect LCS for yes first, no last, right";

  gbwt::vector_type yes_yes = this->forward_path({ 1, 2, 4, 5 });
  std::vector<std::pair<size_t, size_t>> yes_yes_left = { { 0, 0 }, { 1, 1 }, { 2, 3 }, { 3, 4 } };
  ASSERT_EQ(path_lcs(this->graph, yes_yes, path), yes_yes_left) << "Incorrect LCS for yes first, yes last, left";
  std::vector<std::pair<size_t, size_t>> yes_yes_right = { { 0, 0 }, { 1, 1 }, { 3, 2 }, { 4, 3 } };
  ASSERT_EQ(path_lcs(this->graph, path, yes_yes), yes_yes_right) << "Incorrect LCS for yes first, yes last, right";
}

TEST_F(LCSTest, MinimalMaximalLCS)
{
  gbwt::vector_type three_odd = this->forward_path({ 1, 3, 5 });
  gbwt::vector_type four_odd = this->forward_path({ 1, 3, 5, 7 });
  gbwt::vector_type three_even = this->forward_path({ 4, 6, 8 });
  gbwt::vector_type four_even = this->forward_path({ 2, 4, 6, 8 });

  ASSERT_TRUE(path_lcs(this->graph, four_odd, four_even).empty()) << "Nonempty LCS for non-overlapping paths";

  std::vector<std::pair<size_t, size_t>> four_four { { 0, 0 }, { 1, 1 }, { 2, 2 }, { 3, 3 } };
  ASSERT_EQ(path_lcs(this->graph, four_odd, four_odd), four_four) << "Incorrect LCS for identical paths";

  std::vector<std::pair<size_t, size_t>> prefix { { 0, 0 }, { 1, 1 }, { 2, 2 } };
  ASSERT_EQ(path_lcs(this->graph, three_odd, four_odd), prefix) << "Incorrect LCS for prefix / full";
  ASSERT_EQ(path_lcs(this->graph, four_odd, three_odd), prefix) << "Incorrect LCS for full / prefix";

  std::vector<std::pair<size_t, size_t>> left_suffix { { 0, 1 }, { 1, 2 }, { 2, 3 } };
  ASSERT_EQ(path_lcs(this->graph, three_even, four_even), left_suffix) << "Incorrect LCS for suffix / full";
  std::vector<std::pair<size_t, size_t>> right_suffix { { 1, 0 }, { 2, 1 }, { 3, 2 } };
  ASSERT_EQ(path_lcs(this->graph, four_even, three_even), right_suffix) << "Incorrect LCS for full / suffix";
}

TEST_F(LCSTest, RandomLCS)
{
  for(size_t length : { 10, 30, 100, 300 })
  {
    gbwt::vector_type left = this->random_path(length);
    gbwt::vector_type right = this->random_path(length);
    std::string name = "(length " + std::to_string(length) + ")";
    this->compare_paths(left, right, name);
  }
}

TEST_F(LCSTest, RealLCS)
{
  gbwt::size_type paths = this->index.metadata.paths();
  for(gbwt::size_type left_id = 0; left_id < paths; left_id++)
  {
    for(gbwt::size_type right_id = 0; right_id < paths; right_id++)
    {
      for(bool orientation : { false, true })
      {
        gbwt::vector_type left = this->index.extract(gbwt::Path::encode(left_id, orientation));
        gbwt::vector_type right = this->index.extract(gbwt::Path::encode(right_id, orientation));
        std::string name = "(paths " + std::to_string(left_id) + " and " + std::to_string(right_id) + " " + (orientation ? "reverse" : "forward") + ")";
        this->compare_paths(left, right, name);
      }
    }
  }
}

TEST_F(LCSTest, LCSWeights)
{
  // This is an artificial scenario where node lengths matter.
  // Nodes 3 and 6 are long nodes, while nodes 5, 7, and 8 are short nodes.
  gbwt::vector_type left = this->forward_path({ 2, 3, 6, 5, 7, 8 });
  gbwt::vector_type right = this->forward_path({ 2, 5, 7, 8, 3, 6 });
  std::vector<std::pair<size_t, size_t>> truth = { { 0, 0 }, { 1, 4 }, { 2, 5 } };
  std::vector<std::pair<size_t, size_t>> result = path_lcs(this->graph, left, right);
  ASSERT_EQ(result, truth) << "Incorrect LCS weighted by node lengths";
}

TEST_F(LCSTest, SimilarLCS)
{
  std::vector<std::pair<size_t, double>> params = { { 100, 0.02 }, { 1000, 0.01 }, { 10000, 0.005 }, { 30000, 0.005 } };
  for(std::pair<size_t, double> test : params)
  {
    gbwt::vector_type left = this->random_path(test.first);
    gbwt::vector_type right = this->mutate_path(left, test.second);
    std::string name = "(length " + std::to_string(test.first) + ", mutation rate " + std::to_string(test.second) + ")";
    this->compare_paths(left, right, name);
  }
}

//------------------------------------------------------------------------------

} // namespace
