#include <gtest/gtest.h>

#include <set>
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
    auto gfa_parse = gfa_to_gbwt("components.gfa");
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

//------------------------------------------------------------------------------

class TopologicalOrderTest : public ::testing::Test
{
public:
  gbwt::GBWT index;
  GBWTGraph graph;
  size_t components;
  std::vector<std::set<gbwt::vector_type>> correct_paths;

  TopologicalOrderTest()
  {
  }

  void SetUp() override
  {
    auto gfa_parse = gfa_to_gbwt("components.gfa"); // FIXME a new graph with a single component containing a cycle
    this->index = *(gfa_parse.first);
    this->graph = GBWTGraph(this->index, *(gfa_parse.second));
    this->components = 2;
  }
};

TEST_F(TopologicalOrderTest, SingleComponent)
{
  // Define subset
  // Define constraints
  // Get topological order
  // Check that the size is correct
  // Check that all constraints hold
}

TEST_F(TopologicalOrderTest, TwoComponents)
{
  
}

TEST_F(TopologicalOrderTest, CyclicComponent)
{
  
}

//------------------------------------------------------------------------------

} // namespace
