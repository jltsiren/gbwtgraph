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

} // namespace
