#include <gtest/gtest.h>

#include <fstream>
#include <set>
#include <utility>
#include <vector>

#include <arpa/inet.h>
#include <omp.h>

#include <gbwtgraph/gbwtgraph.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

/*
  Test most GBWTGraph operations. SegmentHandleGraph / node-to-segment translation
  is GFA-specific functionality that is currently tested in test_gfa.cpp.
*/

class GraphOperations : public ::testing::Test
{
public:
  typedef std::pair<gbwt::node_type, gbwt::node_type> gbwt_edge;

  gbwt::GBWT index;
  NaiveGraph source;
  GBWTGraph graph;
  std::vector<node_type> correct_nodes;
  std::set<gbwt_edge> correct_edges;
  std::set<gbwt::vector_type> correct_paths;
  std::map<std::string, gbwt::vector_type> correct_named_paths;
  std::map<std::string, gbwt::vector_type> correct_reference_paths;
  std::map<std::string, gbwt::vector_type> correct_haplotype_paths;
  std::map<std::string, std::string> correct_sample_name;
  std::map<std::string, std::string> correct_locus_name;
  std::map<std::string, size_t> correct_haplotype_number;
  std::map<std::string, size_t> correct_phase_block_number;
  std::map<std::string, handlegraph::subrange_t> correct_subrange;

  void SetUp() override
  {
    this->index = build_gbwt_index_with_named_paths();
    this->source = build_naive_graph(false);
    this->graph = GBWTGraph(this->index, this->source);

    this->correct_nodes =
    {
      { nid_t(1), "G" },
      { nid_t(2), "A" },
      { nid_t(4), "GGG" },
      { nid_t(5), "T" },
      { nid_t(6), "A" },
      { nid_t(7), "C" },
      { nid_t(8), "A" },
      { nid_t(9), "A" }
    };

    this->correct_edges =
    {
      { gbwt::Node::encode(1, false), gbwt::Node::encode(2, false) },
      { gbwt::Node::encode(1, false), gbwt::Node::encode(4, false) },
      { gbwt::Node::encode(2, false), gbwt::Node::encode(4, false) },
      { gbwt::Node::encode(4, false), gbwt::Node::encode(5, false) },
      { gbwt::Node::encode(5, false), gbwt::Node::encode(6, false) },
      { gbwt::Node::encode(6, false), gbwt::Node::encode(7, false) },
      { gbwt::Node::encode(6, false), gbwt::Node::encode(8, false) },
      { gbwt::Node::encode(7, false), gbwt::Node::encode(9, false) },
      { gbwt::Node::encode(8, false), gbwt::Node::encode(9, false) }
    };

    this->correct_paths = { alt_path, short_path };

    // Path order is short, alt, short, empty, empty as ref 1, ref 2, sample, empty 1, empty 2
    this->correct_named_paths = {{"chr1", short_path}, {"chr2", alt_path}, {"empty1", empty_path}, {"GRCh38#empty2", empty_path}};

    this->correct_reference_paths = {{"GRCh38#empty2", empty_path}};
    
    this->correct_haplotype_paths = {{"Jouni Sirén#0#chr1#0", short_path}};
    
    this->correct_sample_name = {{"chr1", handlegraph::PathMetadata::NO_SAMPLE_NAME},
                                 {"chr2", handlegraph::PathMetadata::NO_SAMPLE_NAME},
                                 {"empty1", handlegraph::PathMetadata::NO_SAMPLE_NAME},
                                 {"GRCh38#empty2", "GRCh38"},
                                 {"Jouni Sirén#0#chr1#0", "Jouni Sirén"}};
                                 
    this->correct_locus_name = {{"chr1", "chr1"},
                                {"chr2", "chr2"},
                                {"empty1", "empty1"},
                                {"GRCh38#empty2", "empty2"},
                                {"Jouni Sirén#0#chr1#0", "chr1"}};
    
    this->correct_haplotype_number = {{"chr1", handlegraph::PathMetadata::NO_HAPLOTYPE},
                                      {"chr2", handlegraph::PathMetadata::NO_HAPLOTYPE},
                                      {"empty1", handlegraph::PathMetadata::NO_HAPLOTYPE},
                                      {"GRCh38#empty2", handlegraph::PathMetadata::NO_HAPLOTYPE},
                                      {"Jouni Sirén#0#chr1#0", 0}};
    
    this->correct_phase_block_number = {{"chr1", handlegraph::PathMetadata::NO_PHASE_BLOCK},
                                        {"chr2", handlegraph::PathMetadata::NO_PHASE_BLOCK},
                                        {"empty1", handlegraph::PathMetadata::NO_PHASE_BLOCK},
                                        {"GRCh38#empty2", handlegraph::PathMetadata::NO_PHASE_BLOCK},
                                        {"Jouni Sirén#0#chr1#0", 0}};
    
    this->correct_subrange = {{"chr1", handlegraph::PathMetadata::NO_SUBRANGE},
                              {"chr2", handlegraph::PathMetadata::NO_SUBRANGE},
                              {"empty1", handlegraph::PathMetadata::NO_SUBRANGE},
                              {"GRCh38#empty2", handlegraph::PathMetadata::NO_SUBRANGE},
                              {"Jouni Sirén#0#chr1#0", handlegraph::PathMetadata::NO_SUBRANGE}};
    
  }
};

TEST_F(GraphOperations, EmptyGraph)
{
  gbwt::GBWT empty_index;
  NaiveGraph empty_source;
  GBWTGraph empty_graph(empty_index, empty_source);
  EXPECT_EQ(empty_graph.get_node_count(), static_cast<size_t>(0)) << "Empty graph contains nodes";
  EXPECT_FALSE(empty_graph.has_segment_names()) << "Empty graph has segment names";

  std::vector<node_type> nodes;
  std::set<gbwt_edge> edges;
  handle_graph_handles(empty_graph, nodes, true);
  handle_graph_nodes(empty_graph, nodes);
  handle_graph_edges(empty_graph, nodes, edges);
}

// HandleGraph interface for the default test graph.
TEST_F(GraphOperations, HandleGraphTests)
{
  handle_graph_handles(this->graph, this->correct_nodes, true);
  handle_graph_nodes(this->graph, this->correct_nodes);
  handle_graph_edges(this->graph, this->correct_nodes, this->correct_edges);
}

// Custom interface for the default test graph.
TEST_F(GraphOperations, CustomInterface)
{
  EXPECT_FALSE(this->graph.has_segment_names()) << "Graph has segment names";

  // Check that handles are as expected.
  for(const node_type& node : this->correct_nodes)
  {
    handle_t fw_handle = this->graph.get_handle(node.first, false);
    handle_t rev_handle = this->graph.get_handle(node.first, true);
    gbwt::node_type fw_node = gbwt::Node::encode(node.first, false);
    gbwt::node_type rev_node = gbwt::Node::encode(node.first, true);
    handle_t fw_from_node = GBWTGraph::node_to_handle(fw_node);
    handle_t rev_from_node = GBWTGraph::node_to_handle(rev_node);
    gbwt::node_type fw_from_handle = GBWTGraph::handle_to_node(fw_handle);
    gbwt::node_type rev_from_handle = GBWTGraph::handle_to_node(rev_handle);
    EXPECT_EQ(as_integer(fw_handle), as_integer(fw_from_node)) << "Forward handle incorrect for node " << node.first;
    EXPECT_EQ(as_integer(rev_handle), as_integer(rev_from_node)) << "Reverse handle incorrect for node " << node.first;
    EXPECT_EQ(fw_node, fw_from_handle) << "Incorrect GBWT node id from forward handle for node " << node.first;
    EXPECT_EQ(rev_node, rev_from_handle) << "Incorrect GBWT node id from reverse handle for node " << node.first;
  }

  // Sequence views; starts/ends with.
  for(const node_type& node : this->correct_nodes)
  {
    handle_t fw_handle = this->graph.get_handle(node.first, false);
    EXPECT_EQ(this->graph.get_sequence_view(fw_handle), std::string_view(node.second)) << "Invalid sequence view for forward node " << node.first;
    EXPECT_TRUE(this->graph.starts_with(fw_handle, node.second.front())) << "Wrong first character for node " << node.first;
    EXPECT_FALSE(this->graph.starts_with(fw_handle, 'x')) << "First character is 'x' for node " << node.first;
    EXPECT_TRUE(this->graph.ends_with(fw_handle, node.second.back())) << "Wrong last character for node " << node.first;
    EXPECT_FALSE(this->graph.ends_with(fw_handle, 'x')) << "Last character is 'x' for node " << node.first;

    handle_t rev_handle = this->graph.get_handle(node.first, true);
    std::string reverse = node.second;
    reverse_complement_in_place(reverse);
    EXPECT_EQ(this->graph.get_sequence_view(rev_handle), std::string_view(reverse)) << "Invalid sequence view for reverse node " << node.first;
    EXPECT_TRUE(this->graph.starts_with(rev_handle, reverse.front())) << "Wrong first character for reverse node " << node.first;
    EXPECT_FALSE(this->graph.starts_with(rev_handle, 'x')) << "First character is 'x' for reverse node " << node.first;
    EXPECT_TRUE(this->graph.ends_with(rev_handle, reverse.back())) << "Wrong last character for reverse node " << node.first;
    EXPECT_FALSE(this->graph.ends_with(rev_handle, 'x')) << "Last character is 'x' for reverse node " << node.first;
  }
}

TEST_F(GraphOperations, FromHandleGraph)
{
  GBWTGraph copy(this->index, this->graph, nullptr);
  EXPECT_EQ(copy.header, this->graph.header) << "Invalid header";
  EXPECT_EQ(copy.sequences, this->graph.sequences) << "Invalid sequences";
  EXPECT_EQ(copy.real_nodes, this->graph.real_nodes) << "Invalid real_nodes";
  EXPECT_FALSE(copy.has_segment_names()) << "Got segment names from HandleGraph";
}

TEST_F(GraphOperations, Paths)
{
  ASSERT_EQ(this->graph.get_path_count(), this->correct_named_paths.size() + this->correct_haplotype_paths.size()) << "Wrong number of paths";
}

TEST_F(GraphOperations, NamedPaths)
{
  EXPECT_FALSE(this->graph.has_path("SirNotAppearingInThisGraph")) << "Named path that shouldn't exist appears to exist";

  // We need to count the steps we see on each node.
  std::map<handlegraph::nid_t, size_t> steps_on_node;

  for(auto& kv : this->correct_named_paths)
  {
    ASSERT_TRUE(this->graph.has_path(kv.first))
      << "Named path " << kv.first << " that should exist appears not to";
    // Grab the path
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    // Make sure we can recover the name
    ASSERT_EQ(this->graph.get_path_name(path_handle), kv.first)
      << "Named path " << kv.first << " does not appear to have its own name";
    // Trace along the truth path
    handlegraph::step_handle_t step_handle = this->graph.path_begin(path_handle);
    // Remember the previous step
    handlegraph::step_handle_t expected_previous;
    // And count our step index in the truth
    size_t index_in_path = 0;
    while (index_in_path < kv.second.size())
    {
      // Make sure we are still on the correct path
      ASSERT_EQ(this->graph.get_path_handle_of_step(step_handle), path_handle)
        << "Step " << index_in_path << " of path " << kv.first
        << " does not appear to belong to the path";

      // Check step node and orientation
      handlegraph::handle_t visited = this->graph.get_handle_of_step(step_handle);
      ASSERT_EQ(this->graph.get_id(visited), static_cast<nid_t>(gbwt::Node::id(kv.second[index_in_path])))
        << "Step " << index_in_path << " of path " << kv.first << " visits the wrong node";
      ASSERT_EQ(this->graph.get_is_reverse(visited), gbwt::Node::is_reverse(kv.second[index_in_path]))
        << "Step " << index_in_path << " of path " << kv.first << " visits the right node backward";
      steps_on_node[this->graph.get_id(visited)]++;

      // Get the path steps on the node and make sure one of them is us.
      bool found_self = false;
      auto visit_step = [&](const handlegraph::step_handle_t other_step) -> bool {
        EXPECT_FALSE(found_self) << "Step iteration did not stop early";
        // Node could be forward or backward in the path, but it must be this node.
        if (other_step == step_handle)
        {
          // Stop as soon as we see ourselves.
          found_self = true;
          return false;
        }
        return true;
      };
      // Try forward
      this->graph.for_each_step_on_handle(visited, visit_step);
      EXPECT_TRUE(found_self) << "Step iteration on handle did not find current step";
      // And try reverse
      found_self = false;
      this->graph.for_each_step_on_handle(this->graph.flip(visited), visit_step);
      EXPECT_TRUE(found_self) << "Step iteration on flipped handle did not find current step";

      if (index_in_path == 0)
      {
        // This is the very first entry in the path
        EXPECT_FALSE(this->graph.has_previous_step(step_handle))
          << "Initial step of path " << kv.first << " appears to have a predecessor";

        // We can still look left and we expect to see path_front_end()
        handlegraph::step_handle_t observed_previous = this->graph.get_previous_step(step_handle);
        handlegraph::step_handle_t front_end_step_handle = this->graph.path_front_end(path_handle);
        ASSERT_EQ(observed_previous, front_end_step_handle)
          << "Running off start of path " << kv.first
          << " did not yield path_front_end for the path";
      } else {
        // This isn't the first entry in the path, so we should have a previous entry.
        EXPECT_TRUE(this->graph.has_previous_step(step_handle))
          << "Step " << index_in_path << " of path " << kv.first
          << " appears to be the first";
        // Look left
        handlegraph::step_handle_t observed_previous = this->graph.get_previous_step(step_handle);
        ASSERT_EQ(observed_previous, expected_previous)
          << "Step " << index_in_path << " of path " << kv.first
          << " has wrong previous step";
      }

      // Save the previous visited handle
      expected_previous = step_handle;

      if (index_in_path + 1 == kv.second.size())
      {
        // This is the last entry in the path
        EXPECT_FALSE(this->graph.has_next_step(step_handle))
          << "Final step " << index_in_path << " of path " << kv.first
          << " appears to have a successor";

        // While here, we expect to match path_back()
        handlegraph::step_handle_t back_step_handle = this->graph.path_back(path_handle);
        ASSERT_EQ(step_handle, back_step_handle)
          << "Passing last step in path " << kv.first
          << " did not yield path_back for the path";

        // We can still advance and we expect to see path_end()
        step_handle = this->graph.get_next_step(step_handle);
        handlegraph::step_handle_t end_step_handle = this->graph.path_end(path_handle);
        ASSERT_EQ(step_handle, end_step_handle)
          << "Running off end of path " << kv.first
          << " did not yield path_end for the path";
      } else {
        // Advance along the path
        EXPECT_TRUE(this->graph.has_next_step(step_handle))
          << "Step " << index_in_path << " of path " << kv.first
          << " appears to be the last";
        step_handle = this->graph.get_next_step(step_handle);
      }
      // And advance along the true path
      index_in_path++;
    }
    EXPECT_EQ(this->graph.get_step_count(path_handle), kv.second.size())
      << "Path " << kv.first << " appears to have the wrong number of steps";
  }
  this->graph.for_each_handle([&](const handlegraph::handle_t& here)
  {
    // Make sure the paths we traced agree with the per-node step counts.
    EXPECT_EQ(steps_on_node[this->graph.get_id(here)], this->graph.get_step_count(here))
      << "Node " << this->graph.get_id(here) << " has the wrong number of steps on it";
  });
}

TEST_F(GraphOperations, PathMetadata)
{
  
  // Make sure iteration of haplotype paths works.
  std::set<std::string> haplotype_paths_seen;
  this->graph.for_each_path_of_sense(handlegraph::PathSense::HAPLOTYPE, [&](const handlegraph::path_handle_t& path_handle)
  {
    auto name = this->graph.get_path_name(path_handle);
    haplotype_paths_seen.insert(name);
    EXPECT_TRUE(this->correct_haplotype_paths.count(name)) << "Unexpected haplotype path " << name;
  });
  EXPECT_EQ(haplotype_paths_seen.size(), this->correct_haplotype_paths.size())
    << "Found wrong number of haplotype paths";
    
  // Make sure iteration of reference paths works (even if we can't store any).
  std::set<std::string> reference_paths_seen;
  this->graph.for_each_path_of_sense(handlegraph::PathSense::REFERENCE, [&](const handlegraph::path_handle_t& path_handle)
  {
    auto name = this->graph.get_path_name(path_handle);
    reference_paths_seen.insert(name);
    EXPECT_TRUE(this->correct_reference_paths.count(name)) << "Unexpected reference path " << name;
  });
  EXPECT_EQ(reference_paths_seen.size(), this->correct_reference_paths.size())
    << "Found wrong number of reference paths";
    
  // Make sure iteration of generic paths works.
  std::set<std::string> generic_paths_seen;
  this->graph.for_each_path_of_sense(handlegraph::PathSense::GENERIC, [&](const handlegraph::path_handle_t& path_handle)
  {
    auto name = this->graph.get_path_name(path_handle);
    generic_paths_seen.insert(name);
    EXPECT_TRUE(this->correct_named_paths.count(name)) << "Unexpected haplotype path " << name;
  });
  EXPECT_EQ(generic_paths_seen.size(), this->correct_named_paths.size() - this->correct_reference_paths.size())
    << "Found wrong number of generic paths";
  
  for(auto& kv : this->correct_haplotype_paths)
  {
    ASSERT_TRUE(this->graph.has_path(kv.first))
      << "Haplotype path " << kv.first << " that should exist appears not to";
    // Grab the path
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    // Make sure we can recover the name
    ASSERT_EQ(this->graph.get_path_name(path_handle), kv.first)
      << "Haplotype path " << kv.first << " does not appear to have its own name";
      
    // Make sure we can get the right start and end.
    // We assume tracing works the same as for named paths.
    handlegraph::step_handle_t front_handle = this->graph.path_begin(path_handle);
    auto front_visit = this->graph.get_handle_of_step(front_handle);
    ASSERT_EQ(this->graph.get_id(front_visit), static_cast<nid_t>(gbwt::Node::id(kv.second.front())))
        << "Path " << kv.first << " starts at the wrong node";
    ASSERT_EQ(this->graph.get_is_reverse(front_visit), gbwt::Node::is_reverse(kv.second.front()))
      << "Path " << kv.first << " starts at the right node backward";
    handlegraph::step_handle_t back_handle = this->graph.path_back(path_handle);
    auto back_visit = this->graph.get_handle_of_step(back_handle);
    ASSERT_EQ(this->graph.get_id(back_visit), static_cast<nid_t>(gbwt::Node::id(kv.second.back())))
        << "Path " << kv.first << " ends at the wrong node";
    ASSERT_EQ(this->graph.get_is_reverse(back_visit), gbwt::Node::is_reverse(kv.second.back()))
      << "Path " << kv.first << " ends at the right node backward";
    
    // Make sure we can see steps
    bool found_front_step = false;
    this->graph.for_each_step_of_sense(front_visit, handlegraph::PathSense::HAPLOTYPE, [&](const handlegraph::step_handle_t& step)
    {
      if(step == front_handle)
      {
        found_front_step = true;
        return false;
      }
      return true;
    });
    EXPECT_TRUE(found_front_step) << "Front step of " << kv.first << " not visible on node";
    bool found_back_step = false;
    this->graph.for_each_step_of_sense(back_visit, handlegraph::PathSense::HAPLOTYPE, [&](const handlegraph::step_handle_t& step)
    {
      if(step == back_handle)
      {
        found_back_step = true;
        return false;
      }
      return true;
    });
    EXPECT_TRUE(found_back_step) << "Back step of " << kv.first << " not visible on node";
    
    // And check size.
    EXPECT_EQ(this->graph.get_step_count(path_handle), kv.second.size())
      << "Path " << kv.first << " appears to have the wrong number of steps";
      
    // Check sense
    EXPECT_EQ(this->graph.get_sense(path_handle), handlegraph::PathSense::HAPLOTYPE)
      << "Haplotype has wrong sense";
  }
  
  for(auto& kv : this->correct_reference_paths)
  {
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    // Check sense of reference paths
    EXPECT_EQ(this->graph.get_sense(path_handle), handlegraph::PathSense::REFERENCE)
      << "Named path has wrong sense";
  }
  
  for(auto& kv : this->correct_named_paths)
  {
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    // Check sense of named paths
    auto true_sense = this->correct_reference_paths.count(kv.first) ? handlegraph::PathSense::REFERENCE : handlegraph::PathSense::GENERIC;
    EXPECT_EQ(this->graph.get_sense(path_handle), true_sense)
      << "Named path has wrong sense";
  }
  
  // Now check metadata across path types.
  
  for(auto& kv : this->correct_sample_name)
  {
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    EXPECT_EQ(this->graph.get_sample_name(path_handle), kv.second)
      << "Path " << kv.first << " has wrong sample name";
  }
  
  for(auto& kv : this->correct_locus_name)
  {
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    EXPECT_EQ(this->graph.get_locus_name(path_handle), kv.second)
      << "Path " << kv.first << " has wrong sample name";
  }
  
  for(auto& kv : this->correct_haplotype_number)
  {
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    EXPECT_EQ(this->graph.get_haplotype(path_handle), kv.second)
      << "Path " << kv.first << " has wrong haplotype number";
  }
  
  for(auto& kv : this->correct_phase_block_number)
  {
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    EXPECT_EQ(this->graph.get_phase_block(path_handle), kv.second)
      << "Path " << kv.first << " has wrong phase block number";
  }
  
  for(auto& kv : this->correct_subrange)
  {
    handlegraph::path_handle_t path_handle = this->graph.get_path_handle(kv.first);
    EXPECT_EQ(this->graph.get_subrange(path_handle), kv.second)
      << "Path " << kv.first << " has wrong subrange";
  }
}

TEST_F(GraphOperations, ForwardTraversal)
{
  typedef std::pair<gbwt::SearchState, gbwt::vector_type> state_type;
  std::vector<gbwt::vector_type> found_paths;
  std::stack<state_type> states;

  // Extract all paths starting from node 1 in forward orientation.
  handle_t first_node = this->graph.get_handle(1, false);
  states.push({ this->graph.get_state(first_node),
                { static_cast<gbwt::vector_type::value_type>(GBWTGraph::handle_to_node(first_node)) } });
  while(!states.empty())
  {
    state_type curr = states.top();
    states.pop();
    bool extend_success = false;
    this->graph.follow_paths(curr.first, [&](const gbwt::SearchState& next_search) -> bool
    {
      if(!next_search.empty())
      {
        extend_success = true;
        gbwt::vector_type next_path = curr.second;
        next_path.push_back(next_search.node);
        states.push({ next_search, next_path });
      }
      return true;
    });
    if(!extend_success) { found_paths.push_back(curr.second); }
  }

  ASSERT_EQ(found_paths.size(), correct_paths.size()) << "Found a wrong number of paths";
  for(size_t i = 0; i < found_paths.size(); i++)
  {
    const gbwt::vector_type& path = found_paths[i];
    EXPECT_TRUE(correct_paths.find(path) != correct_paths.end()) << "Path " << i << " is incorrect";
  }
}

TEST_F(GraphOperations, BidirectionalTraversal)
{
  typedef std::pair<gbwt::BidirectionalState, gbwt::vector_type> state_type;
  std::vector<gbwt::vector_type> found_paths;
  std::stack<state_type> fw_states, rev_states;

  // Extract all paths visiting node 4 in forward orientation.
  handle_t first_node = this->graph.get_handle(4, false);
  fw_states.push({ this->graph.get_bd_state(first_node),
                 { static_cast<gbwt::vector_type::value_type>(GBWTGraph::handle_to_node(first_node)) } });
  while(!fw_states.empty())
  {
    state_type curr = fw_states.top();
    fw_states.pop();
    bool extend_success = false;
    this->graph.follow_paths(curr.first, false, [&](const gbwt::BidirectionalState& next_search) -> bool
    {
      if(!next_search.empty())
      {
        extend_success = true;
        gbwt::vector_type next_path = curr.second;
        next_path.push_back(next_search.forward.node);
        fw_states.push({ next_search, next_path });
      }
      return true;
    });
    if (!extend_success) { rev_states.push(curr); }
  }
  while (!rev_states.empty())
  {
    state_type curr = rev_states.top();
    rev_states.pop();
    bool extend_success = false;
    this->graph.follow_paths(curr.first, true, [&](const gbwt::BidirectionalState& next_search) -> bool
    {
      if(!next_search.empty())
      {
        extend_success = true;
        gbwt::vector_type next_path
        {
          static_cast<gbwt::vector_type::value_type>(gbwt::Node::reverse(next_search.backward.node))
        };
        next_path.insert(next_path.end(), curr.second.begin(), curr.second.end());
        rev_states.push({ next_search, next_path });
      }
      return true;
    });
    if (!extend_success) { found_paths.push_back(curr.second); }
  }

  ASSERT_EQ(found_paths.size(), correct_paths.size()) << "Found a wrong number of paths";
  for(size_t i = 0; i < found_paths.size(); i++)
  {
    const gbwt::vector_type& path = found_paths[i];
    EXPECT_TRUE(correct_paths.find(path) != correct_paths.end()) << "Path " << i << " is incorrect";
  }
}

//------------------------------------------------------------------------------

class GraphSerialization : public ::testing::Test
{
public:
  gbwt::GBWT index;
  NaiveGraph source;
  GBWTGraph graph;

  GraphSerialization()
  {
  }

  void SetUp() override
  {
    this->index = build_gbwt_index();
    this->source = build_naive_graph(false);
    this->graph = GBWTGraph(this->index, this->source);
  }

  void check_graph(const GBWTGraph& graph, const GBWTGraph& truth) const
  {
    ASSERT_EQ(graph.header, truth.header) << "Serialization did not preserve the header";
    ASSERT_EQ(graph.sequences, truth.sequences) << "Serialization did not preserve the sequences";
    ASSERT_EQ(graph.real_nodes, truth.real_nodes) << "Serialization did not preserve the real nodes";
    ASSERT_EQ(graph.segments, truth.segments) << "Serialization did not preserve the segments";
    ASSERT_EQ(graph.node_to_segment, truth.node_to_segment) << "Serialization did not preserve the node-to-segment mapping";
  }
};

TEST_F(GraphSerialization, SerializeEmpty)
{
  GBWTGraph empty_graph;
  std::string filename = gbwt::TempFile::getName("gbwtgraph");
  empty_graph.serialize(filename);

  GBWTGraph duplicate_graph;
  duplicate_graph.deserialize(filename);
  this->check_graph(duplicate_graph, empty_graph);

  gbwt::TempFile::remove(filename);
}

TEST_F(GraphSerialization, CompressEmpty)
{
  gbwt::GBWT empty_gbwt;
  GBWTGraph empty_graph;
  empty_graph.set_gbwt(empty_gbwt);
  size_t expected_size = empty_graph.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
  std::string filename = gbwt::TempFile::getName("gbwtgraph");
  sdsl::simple_sds::serialize_to(empty_graph, filename);

  GBWTGraph duplicate_graph;
  std::ifstream in(filename, std::ios_base::binary);
  size_t bytes = gbwt::fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid file size";
  duplicate_graph.simple_sds_load(in, *(empty_graph.index));
  in.close();
  this->check_graph(duplicate_graph, empty_graph);

  gbwt::TempFile::remove(filename);
}

TEST_F(GraphSerialization, SerializeNonEmpty)
{
  std::string filename = gbwt::TempFile::getName("gbwtgraph");
  this->graph.serialize(filename);

  GBWTGraph duplicate_graph;
  duplicate_graph.deserialize(filename);
  duplicate_graph.set_gbwt(this->index);
  this->check_graph(duplicate_graph, this->graph);

  gbwt::TempFile::remove(filename);
}

TEST_F(GraphSerialization, CompressNonEmpty)
{
  size_t expected_size = this->graph.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
  std::string filename = gbwt::TempFile::getName("gbwtgraph");
  sdsl::simple_sds::serialize_to(this->graph, filename);

  GBWTGraph duplicate_graph;
  std::ifstream in(filename, std::ios_base::binary);
  size_t bytes = gbwt::fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid file size";
  duplicate_graph.simple_sds_load(in, this->index);
  in.close();
  this->check_graph(duplicate_graph, this->graph);

  gbwt::TempFile::remove(filename);
}

TEST_F(GraphSerialization, SerializeTranslation)
{
  NaiveGraph source = build_naive_graph(true);
  GBWTGraph graph(this->index, source);

  std::string filename = gbwt::TempFile::getName("gbwtgraph");
  graph.serialize(filename);

  GBWTGraph duplicate_graph;
  duplicate_graph.deserialize(filename);
  duplicate_graph.set_gbwt(this->index);
  this->check_graph(duplicate_graph, graph);

  gbwt::TempFile::remove(filename);
}

TEST_F(GraphSerialization, CompressTranslation)
{
  NaiveGraph source = build_naive_graph(true);
  GBWTGraph graph(this->index, source);
  size_t expected_size = graph.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
  std::string filename = gbwt::TempFile::getName("gbwtgraph");
  sdsl::simple_sds::serialize_to(graph, filename);

  GBWTGraph duplicate_graph;
  std::ifstream in(filename, std::ios_base::binary);
  size_t bytes = gbwt::fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid file size";
  duplicate_graph.simple_sds_load(in, this->index);
  in.close();
  this->check_graph(duplicate_graph, graph);

  gbwt::TempFile::remove(filename);
}

TEST_F(GraphSerialization, DecompressSerialized)
{
  NaiveGraph source = build_naive_graph(true);
  GBWTGraph graph(this->index, source);

  std::string filename = gbwt::TempFile::getName("gbwtgraph");
  graph.serialize(filename);

  GBWTGraph duplicate_graph;
  std::ifstream in(filename, std::ios_base::binary);
  std::uint32_t magic = 0;
  in.read(reinterpret_cast<char*>(&magic), sizeof(std::uint32_t));
  EXPECT_EQ(ntohl(magic), graph.get_magic_number()) << "Magic number missing from serialized graph";
  duplicate_graph.simple_sds_load(in, this->index);
  in.close();
  this->check_graph(duplicate_graph, graph);

  gbwt::TempFile::remove(filename);
}

//------------------------------------------------------------------------------

class ForEachWindow : public ::testing::Test
{
public:
  typedef std::pair<std::vector<handle_t>, std::string> kmer_type;

  gbwt::GBWT index;
  NaiveGraph source;
  GBWTGraph graph;
  std::set<kmer_type> correct_kmers;
  std::set<kmer_type> correct_nonredundant;

  ForEachWindow()
  {
  }

  void SetUp() override
  {
    this->index = build_gbwt_index();
    this->source = build_naive_graph(false);
    this->graph = GBWTGraph(this->index, this->source);

    this->correct_kmers =
    {
      // alt_path forward
      { { this->graph.get_handle(1, false), this->graph.get_handle(2, false), this->graph.get_handle(4, false) }, "GAG" },
      { { this->graph.get_handle(2, false), this->graph.get_handle(4, false) }, "AGG" },
      { { this->graph.get_handle(4, false), this->graph.get_handle(5, false), this->graph.get_handle(6, false) }, "GGGTA" },
      { { this->graph.get_handle(5, false), this->graph.get_handle(6, false), this->graph.get_handle(8, false) }, "TAA" },
      { { this->graph.get_handle(6, false), this->graph.get_handle(8, false), this->graph.get_handle(9, false) }, "AAA" },

      // alt_path reverse
      { { this->graph.get_handle(9, true), this->graph.get_handle(8, true), this->graph.get_handle(6, true) }, "TTT" },
      { { this->graph.get_handle(8, true), this->graph.get_handle(6, true), this->graph.get_handle(5, true) }, "TTA" },
      { { this->graph.get_handle(6, true), this->graph.get_handle(5, true), this->graph.get_handle(4, true) }, "TAC" },
      { { this->graph.get_handle(5, true), this->graph.get_handle(4, true) }, "ACC" },
      { { this->graph.get_handle(4, true), this->graph.get_handle(2, true), this->graph.get_handle(1, true) }, "CCCTC" },

      // short_path forward
      { { this->graph.get_handle(1, false), this->graph.get_handle(4, false) }, "GGG" },
      { { this->graph.get_handle(4, false), this->graph.get_handle(5, false), this->graph.get_handle(6, false) }, "GGGTA" },
      { { this->graph.get_handle(5, false), this->graph.get_handle(6, false), this->graph.get_handle(7, false) }, "TAC" },
      { { this->graph.get_handle(6, false), this->graph.get_handle(7, false), this->graph.get_handle(9, false) }, "ACA" },

      // short_path reverse
      { { this->graph.get_handle(9, true), this->graph.get_handle(7, true), this->graph.get_handle(6, true) }, "TGT" },
      { { this->graph.get_handle(7, true), this->graph.get_handle(6, true), this->graph.get_handle(5, true) }, "GTA" },
      { { this->graph.get_handle(6, true), this->graph.get_handle(5, true), this->graph.get_handle(4, true) }, "TAC" },
      { { this->graph.get_handle(5, true), this->graph.get_handle(4, true) }, "ACC" },
      { { this->graph.get_handle(4, true), this->graph.get_handle(1, true) }, "CCCC" }
    };

    this->correct_nonredundant =
    {
      // alt_path forward
      { { this->graph.get_handle(1, false), this->graph.get_handle(2, false), this->graph.get_handle(4, false) }, "GAG" },
      { { this->graph.get_handle(2, false), this->graph.get_handle(4, false) }, "AGG" },
      { { this->graph.get_handle(4, false) }, "GGG" },
      { { this->graph.get_handle(4, false), this->graph.get_handle(5, false), this->graph.get_handle(6, false) }, "GGTA" },
      { { this->graph.get_handle(5, false), this->graph.get_handle(6, false), this->graph.get_handle(8, false) }, "TAA" },
      { { this->graph.get_handle(6, false), this->graph.get_handle(8, false), this->graph.get_handle(9, false) }, "AAA" },

      // alt_path reverse
      { { this->graph.get_handle(9, true), this->graph.get_handle(8, true), this->graph.get_handle(6, true) }, "TTT" },
      { { this->graph.get_handle(8, true), this->graph.get_handle(6, true), this->graph.get_handle(5, true) }, "TTA" },
      { { this->graph.get_handle(6, true), this->graph.get_handle(5, true), this->graph.get_handle(4, true) }, "TAC" },
      { { this->graph.get_handle(5, true), this->graph.get_handle(4, true) }, "ACC" },
      { { this->graph.get_handle(4, true) }, "CCC" },
      { { this->graph.get_handle(4, true), this->graph.get_handle(2, true), this->graph.get_handle(1, true) }, "CCTC" },

      // short_path forward
      { { this->graph.get_handle(1, false), this->graph.get_handle(4, false) }, "GGG" },
      { { this->graph.get_handle(4, false) }, "GGG" },
      { { this->graph.get_handle(4, false), this->graph.get_handle(5, false), this->graph.get_handle(6, false) }, "GGTA" },
      { { this->graph.get_handle(5, false), this->graph.get_handle(6, false), this->graph.get_handle(7, false) }, "TAC" },
      { { this->graph.get_handle(6, false), this->graph.get_handle(7, false), this->graph.get_handle(9, false) }, "ACA" },

      // short_path reverse
      { { this->graph.get_handle(9, true), this->graph.get_handle(7, true), this->graph.get_handle(6, true) }, "TGT" },
      { { this->graph.get_handle(7, true), this->graph.get_handle(6, true), this->graph.get_handle(5, true) }, "GTA" },
      { { this->graph.get_handle(6, true), this->graph.get_handle(5, true), this->graph.get_handle(4, true) }, "TAC" },
      { { this->graph.get_handle(5, true), this->graph.get_handle(4, true) }, "ACC" },
      { { this->graph.get_handle(4, true) }, "CCC" },
      { { this->graph.get_handle(4, true), this->graph.get_handle(1, true) }, "CCC" }
    };
  }
};

TEST_F(ForEachWindow, KmerExtraction)
{
  // Extract all haplotype-consistent windows of length 3.
  std::set<kmer_type> found_kmers;
  for_each_haplotype_window(this->graph, 3, [&found_kmers](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    found_kmers.insert(kmer_type(traversal, seq));
  }, false);

  ASSERT_EQ(found_kmers.size(), this->correct_kmers.size()) << "Found a wrong number of kmers";
  for(const kmer_type& kmer : found_kmers)
  {
    EXPECT_TRUE(correct_kmers.find(kmer) != correct_kmers.end()) << "Kmer " << kmer.second << " is incorrect";
  }
}

TEST_F(ForEachWindow, NonRedundantKmers)
{
  // Extract all haplotype-consistent non-redundant windows of length 3.
  std::set<kmer_type> found_kmers;
  for_each_nonredundant_window(this->graph, 3, [&found_kmers](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    found_kmers.insert(kmer_type(traversal, seq));
  }, false);

  ASSERT_EQ(found_kmers.size(), this->correct_nonredundant.size()) << "Found a wrong number of kmers";
  for(const kmer_type& kmer : found_kmers)
  {
    EXPECT_TRUE(correct_nonredundant.find(kmer) != correct_kmers.end()) << "Kmer " << kmer.second << " is incorrect";
  }
}

//------------------------------------------------------------------------------

} // namespace
