#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include <omp.h>

#include <gbwtgraph/cached_gbwtgraph.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class GraphOperations : public ::testing::Test
{
public:
  typedef std::pair<gbwt::node_type, gbwt::node_type> gbwt_edge;

  gbwt::GBWT index;
  SequenceSource source;
  GBWTGraph graph;
  CachedGBWTGraph cached_graph;

  GraphOperations()
  {
  }

  void SetUp() override
  {
    this->index = build_gbwt_index();
    build_source(this->source);
    this->graph = GBWTGraph(this->index, this->source);
    this->cached_graph = CachedGBWTGraph(this->graph);
  }
};

TEST_F(GraphOperations, EmptyGraph)
{
  gbwt::GBWT empty_index;
  SequenceSource empty_source;
  GBWTGraph empty_graph(empty_index, empty_source);
  CachedGBWTGraph empty_cache(empty_graph);
  EXPECT_EQ(empty_cache.get_node_count(), static_cast<size_t>(0)) << "Empty graph contains nodes";
}

TEST_F(GraphOperations, CorrectNodes)
{
  ASSERT_EQ(this->cached_graph.get_node_count(), this->graph.get_node_count()) << "Wrong number of nodes";
  EXPECT_EQ(this->cached_graph.min_node_id(), this->graph.min_node_id()) << "Wrong minimum node id";
  EXPECT_EQ(this->cached_graph.max_node_id(), this->graph.max_node_id()) << "Wrong maximum node id";
  for(nid_t id = this->cached_graph.min_node_id(); id <= this->cached_graph.max_node_id(); id++)
  {
    EXPECT_EQ(this->cached_graph.has_node(id), this->graph.has_node(id)) << "Node set incorrect at " << id;
  }
}

TEST_F(GraphOperations, Handles)
{
  for(nid_t id = this->cached_graph.min_node_id(); id <= this->cached_graph.max_node_id(); id++)
  {
    if(!(this->cached_graph.has_node(id))) { continue; }
    handle_t forward_handle = this->cached_graph.get_handle(id, false);
    handle_t reverse_handle = this->cached_graph.get_handle(id, true);
    EXPECT_EQ(this->cached_graph.get_id(forward_handle), id) << "Wrong node id for forward handle " << id;
    EXPECT_FALSE(this->cached_graph.get_is_reverse(forward_handle)) << "Forward handle " << id << " is reverse";
    EXPECT_EQ(this->cached_graph.get_id(reverse_handle), id) << "Wrong node id for reverse handle " << id;
    EXPECT_TRUE(this->cached_graph.get_is_reverse(reverse_handle)) << "Reverse handle " << id << " is not reverse";
    handle_t flipped_fw = this->cached_graph.flip(forward_handle);
    handle_t flipped_rev = this->cached_graph.flip(reverse_handle);
    EXPECT_NE(as_integer(forward_handle), as_integer(reverse_handle)) << "Forward and reverse handles are identical";
    EXPECT_EQ(as_integer(flipped_fw), as_integer(reverse_handle)) << "Flipped forward handle is not identical to the reverse handle";
    EXPECT_EQ(as_integer(flipped_rev), as_integer(forward_handle)) << "Flipped reverse handle is not identical to the forward handle";
  }
}

TEST_F(GraphOperations, Sequences)
{
  for(nid_t id = this->cached_graph.min_node_id(); id <= this->cached_graph.max_node_id(); id++)
  {
    if(!(this->cached_graph.has_node(id))) { continue; }
    handle_t gbwt_fw = this->cached_graph.get_handle(id, false);
    handle_t gbwt_rev = this->cached_graph.get_handle(id, true);
    handle_t source_fw = this->source.get_handle(id, false);
    std::string source_str = this->source.get_sequence(source_fw);
    EXPECT_EQ(this->cached_graph.get_length(gbwt_fw), source_str.length()) << "Wrong forward length at node " << id;
    EXPECT_EQ(this->cached_graph.get_sequence(gbwt_fw), source_str) << "Wrong forward sequence at node " << id;
    std::string source_rev = reverse_complement(source_str);
    EXPECT_EQ(this->cached_graph.get_length(gbwt_rev), source_rev.length()) << "Wrong reverse length at node " << id;
    EXPECT_EQ(this->cached_graph.get_sequence(gbwt_rev), source_rev) << "Wrong reverse sequence at node " << id;
  }
}

TEST_F(GraphOperations, Substrings)
{
  for(nid_t id = this->cached_graph.min_node_id(); id <= this->cached_graph.max_node_id(); id++)
  {
    if(!(this->cached_graph.has_node(id))) { continue; }
    handle_t fw = this->cached_graph.get_handle(id, false);
    handle_t rev = this->cached_graph.get_handle(id, true);
    std::string fw_str = this->cached_graph.get_sequence(fw);
    std::string rev_str = this->cached_graph.get_sequence(rev);
    ASSERT_EQ(fw_str.length(), rev_str.length()) << "Forward and reverse sequences have different lengths at node " << id;
    for(size_t i = 0; i < fw_str.length(); i++)
    {
      EXPECT_EQ(this->cached_graph.get_base(fw, i), fw_str[i]) << "Wrong forward base " << i << " at node " << id;
      EXPECT_EQ(this->cached_graph.get_base(rev, i), rev_str[i]) << "Wrong reverse base " << i << " at node " << id;
      EXPECT_EQ(this->cached_graph.get_subsequence(fw, i, 2), fw_str.substr(i, 2)) << "Wrong forward substring " << i << " at node " << id;
      EXPECT_EQ(this->cached_graph.get_subsequence(rev, i, 2), rev_str.substr(i, 2)) << "Wrong reverse substring " << i << " at node " << id;
    }
  }
}

TEST_F(GraphOperations, Edges)
{
  std::set<gbwt_edge> correct_edges, reverse_edges;
  for(nid_t id = this->cached_graph.min_node_id(); id <= this->cached_graph.max_node_id(); id++)
  {
    if(!(this->graph.has_node(id))) { continue; }
    handle_t forward_handle = this->graph.get_handle(id, false);
    handle_t reverse_handle = this->graph.get_handle(id, true);
    this->graph.follow_edges(forward_handle, false, [&](const handle_t& handle) {
      correct_edges.insert(gbwt_edge(GBWTGraph::handle_to_node(forward_handle), GBWTGraph::handle_to_node(handle)));
    });
    this->graph.follow_edges(reverse_handle, false, [&](const handle_t& handle) {
      reverse_edges.insert(gbwt_edge(GBWTGraph::handle_to_node(reverse_handle), GBWTGraph::handle_to_node(handle)));
    });
  }

  std::set<gbwt_edge> fw_succ, fw_pred, rev_succ, rev_pred;
  for(nid_t id = this->cached_graph.min_node_id(); id <= this->cached_graph.max_node_id(); id++)
  {
    if(!(this->cached_graph.has_node(id))) { continue; }
    handle_t forward_handle = this->cached_graph.get_handle(id, false);
    handle_t reverse_handle = this->cached_graph.get_handle(id, true);
    size_t fw_out = 0, fw_in = 0, rev_out = 0, rev_in = 0;
    this->cached_graph.follow_edges(forward_handle, false, [&](const handle_t& handle) {
      fw_succ.insert(gbwt_edge(CachedGBWTGraph::handle_to_node(forward_handle), CachedGBWTGraph::handle_to_node(handle)));
      fw_out++;
    });
    this->cached_graph.follow_edges(forward_handle, true, [&](const handle_t& handle) {
      fw_pred.insert(gbwt_edge(CachedGBWTGraph::handle_to_node(handle), CachedGBWTGraph::handle_to_node(forward_handle)));
      fw_in++;
    });
    this->cached_graph.follow_edges(reverse_handle, false, [&](const handle_t& handle) {
      rev_succ.insert(gbwt_edge(CachedGBWTGraph::handle_to_node(reverse_handle), CachedGBWTGraph::handle_to_node(handle)));
      rev_out++;
    });
    this->cached_graph.follow_edges(reverse_handle, true, [&](const handle_t& handle) {
      rev_pred.insert(gbwt_edge(CachedGBWTGraph::handle_to_node(handle), CachedGBWTGraph::handle_to_node(reverse_handle)));
      rev_in++;
    });
    EXPECT_EQ(this->cached_graph.get_degree(forward_handle, false), fw_out) << "Wrong outdegree for forward handle " << id;
    EXPECT_EQ(this->cached_graph.get_degree(forward_handle, true), fw_in) << "Wrong indegree for forward handle " << id;
    EXPECT_EQ(this->cached_graph.get_degree(reverse_handle, false), rev_out) << "Wrong outdegree for reverse handle " << id;
    EXPECT_EQ(this->cached_graph.get_degree(reverse_handle, true), rev_in) << "Wrong indegree for reverse handle " << id;
  }
  EXPECT_EQ(fw_succ, correct_edges) << "Wrong forward successors";
  EXPECT_EQ(fw_pred, correct_edges) << "Wrong forward predecessors";
  EXPECT_EQ(rev_succ, reverse_edges) << "Wrong reverse successors";
  EXPECT_EQ(rev_pred, reverse_edges) << "Wrong reverse predecessors";

  for(nid_t from = this->cached_graph.min_node_id(); from <= this->cached_graph.max_node_id(); from++)
  {
    for(nid_t to = this->cached_graph.min_node_id(); to <= this->cached_graph.max_node_id(); to++)
    {
      for(bool from_rev : { false, true })
      {
        for(bool to_rev : { false, true })
        {
          handle_t from_handle = this->cached_graph.get_handle(from, from_rev);
          handle_t to_handle = this->cached_graph.get_handle(to, to_rev);
          gbwt_edge edge(CachedGBWTGraph::handle_to_node(from_handle), CachedGBWTGraph::handle_to_node(to_handle));
          bool should_have = (correct_edges.find(edge) != correct_edges.end());
          should_have |= (reverse_edges.find(edge) != reverse_edges.end());
          EXPECT_EQ(this->cached_graph.has_edge(from_handle, to_handle), should_have) <<
            "has_edge() failed with (" << from << ", " << from_rev << ") to (" << to << ", " << to_rev <<")";
        }
      }
    }
  }
}

TEST_F(GraphOperations, ForEachHandle)
{
  std::vector<handle_t> correct_handles, found_handles;
  this->graph.for_each_handle([&](const handle_t& handle)
  {
    correct_handles.push_back(handle);
  }, false);
  this->cached_graph.for_each_handle([&](const handle_t& handle)
  {
    found_handles.push_back(handle);
  }, false);
  ASSERT_EQ(found_handles, correct_handles) << "Sequential: Wrong handles in the graph";

  found_handles.clear();
  int old_thread_count = omp_get_max_threads();
  omp_set_num_threads(2);
  this->cached_graph.for_each_handle([&](const handle_t& handle)
  {
    #pragma omp critical
    {
      found_handles.push_back(handle);
    }
  }, false);
  omp_set_num_threads(old_thread_count);
  std::sort(correct_handles.begin(), correct_handles.end());
  std::sort(found_handles.begin(), found_handles.end());
  ASSERT_EQ(found_handles, correct_handles) << "Parallel: Wrong handles in the graph";
}

//------------------------------------------------------------------------------

} // namespace
