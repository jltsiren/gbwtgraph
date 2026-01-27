#ifndef GBWTGRAPH_TESTS_SHARED_H
#define GBWTGRAPH_TESTS_SHARED_H

#include <array>
#include <limits>
#include <set>
#include <vector>

#include <gbwt/dynamic_gbwt.h>

#include <gbwtgraph/minimizer.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/naive_graph.h>

#include <gtest/gtest.h>

/*
  shared.h: Utility functions and data definitions shared between the tests.
*/

namespace
{

//------------------------------------------------------------------------------

using gbwtgraph::nid_t;
using gbwtgraph::handle_t;
using gbwtgraph::pos_t;

//------------------------------------------------------------------------------

// Common HandleGraph tests.

typedef std::pair<nid_t, std::string> node_type; // (id, sequence)
typedef std::pair<gbwt::node_type, gbwt::node_type> gbwt_edge;
typedef std::pair<std::string, std::pair<nid_t, nid_t>> translation_type;

// Set in_order to expect the handles in the given order.
void handle_graph_handles(const handlegraph::HandleGraph& graph, const std::vector<node_type>& nodes, bool in_order)
{
  std::set<nid_t> visited;
  size_t index = 0;
  graph.for_each_handle([&](const handle_t& handle)
  {
    nid_t id = graph.get_id(handle);
    if(in_order)
    {
      EXPECT_EQ(id, nodes[index].first) << "Incorrect handle at index " << index;
    }
    EXPECT_FALSE(graph.get_is_reverse(handle)) << "Handle at index " << index << " is in reverse orientation for node " << id;

    handle_t flipped = graph.flip(handle);
    EXPECT_EQ(graph.get_id(flipped), id) << "Flipped handle at index " << index << " has incorrect id for node " << id;
    EXPECT_TRUE(graph.get_is_reverse(flipped)) << "Flipped handle at index " << index << " is in forward orientation for node " << id;

    handle_t flipped_from_graph = graph.get_handle(id, true);
    EXPECT_EQ(flipped, flipped_from_graph) << "Flipped handle at index " << index << " does not match handle from get_handle for node " << id;

    handle_t double_flipped = graph.flip(flipped);
    EXPECT_EQ(double_flipped, handle) << "Double flipped handle at index " << index << " does not match original handle for node " << id;

    index++; visited.insert(id);
  }, false);
  EXPECT_EQ(index, nodes.size()) << "Incorrect number of handles visited";

  for(const node_type& node : nodes)
  {
    EXPECT_TRUE(visited.find(node.first) != visited.end()) << "Node id " << node.first << " was not visited";
  }

  // Also test parallel for_each_handle.
  int old_thread_count = omp_get_max_threads();
  omp_set_num_threads(2);
  std::array<std::vector<gbwt::node_type>, 2> thread_visits;
  graph.for_each_handle([&](const handle_t& handle) -> bool
  {
    int thread_id = omp_get_thread_num();
    nid_t id = graph.get_id(handle);
    bool is_reverse = graph.get_is_reverse(handle);
    thread_visits[thread_id].push_back(gbwt::Node::encode(id, is_reverse));
    return true;
  }, true);
  omp_set_num_threads(old_thread_count);
  std::vector<gbwt::node_type> combined_visits = thread_visits[0];
  combined_visits.insert(combined_visits.end(), thread_visits[1].begin(), thread_visits[1].end());
  gbwt::removeDuplicates(combined_visits, false);

  ASSERT_EQ(combined_visits.size(), visited.size()) << "Parallel for_each_handle visited wrong number of handles";
  for(size_t i = 0; i < combined_visits.size(); i++)
  {
    nid_t id = gbwt::Node::id(combined_visits[i]);
    EXPECT_TRUE(visited.find(id) != visited.end()) << "Wrong handle at index " << i;
  }
}

void handle_graph_nodes(const handlegraph::HandleGraph& graph, const std::vector<node_type>& nodes)
{
  nid_t min_id = std::numeric_limits<nid_t>::max(), max_id = 0;
  for(const node_type& node : nodes)
  {
    if(node.first < min_id) { min_id = node.first; }
    if(node.first > max_id) { max_id = node.first; }
  }

  ASSERT_EQ(graph.get_node_count(), nodes.size()) << "Incorrect number of nodes";
  EXPECT_EQ(graph.min_node_id(), min_id) << "Incorrect minimum node id";
  EXPECT_EQ(graph.max_node_id(), max_id) << "Incorrect maximum node id";

  for(const node_type& node : nodes)
  {
    ASSERT_TRUE(graph.has_node(node.first)) << "Node id " << node.first << " is missing";
    handle_t handle = graph.get_handle(node.first, false);
    EXPECT_EQ(graph.get_length(handle), node.second.length()) << "Invalid sequence length for node " << node.first;
    EXPECT_EQ(graph.get_sequence(handle), node.second) << "Invalid sequence for node " << node.first;
    std::string reverse = gbwtgraph::reverse_complement(node.second);
    handle_t rev_handle = graph.get_handle(node.first, true);
    EXPECT_EQ(graph.get_length(rev_handle), node.second.length()) << "Invalid sequence length for reverse node " << node.first;
    EXPECT_EQ(graph.get_sequence(rev_handle), reverse) << "Invalid sequence for reverse node " << node.first;

    for(size_t i = 0; i < node.second.length(); i++)
    {
      EXPECT_EQ(graph.get_base(handle, i), node.second[i]) << "Invalid base at index " << i << " for node " << node.first;
      EXPECT_EQ(graph.get_base(rev_handle, i), reverse[i]) << "Invalid base at index " << i << " for reverse node " << node.first;
      EXPECT_EQ(graph.get_subsequence(handle, i, 2), node.second.substr(i, 2)) << "Invalid subsequence at index " << i << " for node " << node.first;
      EXPECT_EQ(graph.get_subsequence(rev_handle, i, 2), reverse.substr(i, 2)) << "Invalid subsequence at index " << i << " for reverse node " << node.first;
    }
  }

  nid_t missing_id = graph.max_node_id() + 100;
  EXPECT_FALSE(graph.has_node(missing_id)) << "Missing node id " << missing_id << " is reported present";
}

void handle_graph_edges(const handlegraph::HandleGraph& graph, const std::vector<node_type>& nodes, const std::set<gbwt_edge>& edges)
{
  ASSERT_EQ(graph.get_edge_count(), edges.size()) << "Incorrect number of edges";
  std::set<gbwt_edge> reverse_edges;
  for(const auto& edge : edges)
  {
    reverse_edges.insert({ gbwt::Node::reverse(edge.second), gbwt::Node::reverse(edge.first) });
  }

  // follow_edges and get_degree
  std::set<gbwt_edge> fw_succ, fw_pred, rev_succ, rev_pred;
  for(const auto& node : nodes)
  {
    handle_t forward_handle = graph.get_handle(node.first, false);
    handle_t reverse_handle = graph.get_handle(node.first, true);
    size_t fw_out = 0, fw_in = 0, rev_out = 0, rev_in = 0;
    graph.follow_edges(forward_handle, false, [&](const handle_t& handle)
    {
      gbwt::node_type from = gbwt::Node::encode(node.first, false);
      gbwt::node_type to = gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle));
      fw_succ.insert(gbwt_edge(from, to));
      fw_out++;
    });
    graph.follow_edges(forward_handle, true, [&](const handle_t& handle)
    {
      gbwt::node_type from = gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle));
      gbwt::node_type to = gbwt::Node::encode(node.first, false);
      fw_pred.insert(gbwt_edge(from, to));
      fw_in++;
    });
    graph.follow_edges(reverse_handle, false, [&](const handle_t& handle)
    {
      gbwt::node_type from = gbwt::Node::encode(node.first, true);
      gbwt::node_type to = gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle));
      rev_succ.insert(gbwt_edge(from, to));
      rev_out++;
    });
    graph.follow_edges(reverse_handle, true, [&](const handle_t& handle)
    {
      gbwt::node_type from = gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle));
      gbwt::node_type to = gbwt::Node::encode(node.first, true);
      rev_pred.insert(gbwt_edge(from, to));
      rev_in++;
    });
    EXPECT_EQ(graph.get_degree(forward_handle, false), fw_out) << "Wrong outdegree for forward handle " << node.first;
    EXPECT_EQ(graph.get_degree(forward_handle, true), fw_in) << "Wrong indegree for forward handle " << node.first;
    EXPECT_EQ(graph.get_degree(reverse_handle, false), rev_out) << "Wrong outdegree for reverse handle " << node.first;
    EXPECT_EQ(graph.get_degree(reverse_handle, true), rev_in) << "Wrong indegree for reverse handle " << node.first;
  }
  EXPECT_EQ(fw_succ, edges) << "Wrong forward successors";
  EXPECT_EQ(fw_pred, edges) << "Wrong forward predecessors";
  EXPECT_EQ(rev_succ, reverse_edges) << "Wrong reverse successors";
  EXPECT_EQ(rev_pred, reverse_edges) << "Wrong reverse predecessors";

  // Presence / absence of edges
  for(nid_t from = graph.min_node_id(); from <= graph.max_node_id(); from++)
  {
    for(nid_t to = graph.min_node_id(); to <= graph.max_node_id(); to++)
    {
      for(bool from_rev : { false, true })
      {
        for(bool to_rev : { false, true })
        {
          handle_t from_handle = graph.get_handle(from, from_rev);
          handle_t to_handle = graph.get_handle(to, to_rev);
          gbwt_edge edge(gbwt::Node::encode(from, from_rev), gbwt::Node::encode(to, to_rev));
          bool should_have = (edges.find(edge) != edges.end());
          should_have |= (reverse_edges.find(edge) != reverse_edges.end());
          EXPECT_EQ(graph.has_edge(from_handle, to_handle), should_have) <<
            "has_edge() failed with (" << from << ", " << from_rev << ") to (" << to << ", " << to_rev <<")";
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

// Construction of the basic test case.

gbwt::vector_type alt_path
{
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(2, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(8, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
};

gbwt::vector_type short_path
{
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(7, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
};

gbwt::vector_type empty_path
{
};

// Build a GBWT index for the given paths without metadata.
inline gbwt::GBWT
build_gbwt(const std::vector<gbwt::vector_type>& paths)
{
  gbwt::size_type node_width = 1, total_length = 0;
  for(auto& path : paths)
  {
    for(auto node : path)
    {
      node_width = std::max(node_width, gbwt::size_type(sdsl::bits::length(gbwt::Node::encode(node, true))));
    }
    total_length += 2 * (path.size() + 1);
  }

  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::GBWTBuilder builder(node_width, total_length);
  for(auto& path : paths) { builder.insert(path, true); }
  builder.finish();

  return gbwt::GBWT(builder.index);
}

// Build a GBWT with three paths including a duplicate, plus n empty paths.
inline gbwt::GBWT
build_gbwt_index(size_t empty_paths = 0)
{
  std::vector<gbwt::vector_type> paths
  {
    short_path, alt_path, short_path
  };
  for(size_t i = 0; i < empty_paths; i++)
  {
    paths.push_back(empty_path);
  }
  return build_gbwt(paths);
}

// Build a GBWT with 6 paths, 3 duplicates of each.
// TODO: what's the "ref" here?
// This should be identical to gfas/example_walks.gfa
inline gbwt::GBWT
build_gbwt_index_with_ref()
{
  std::vector<gbwt::vector_type> paths
  {
    short_path, alt_path, alt_path,
    short_path, alt_path, short_path
  };
  return build_gbwt(paths);
}

// This should be identical to gfas/example_walks.gfa, modulo name ordering.
inline gbwt::GBWT
build_gbwt_example_walks()
{
  gbwt::GBWT built = build_gbwt_index_with_ref();
  
  built.addMetadata();

  // Name the set of samples, including a special ref one for generic paths
  built.metadata.setSamples({gbwtgraph::GENERIC_PATH_SAMPLE_NAME, "short", "alt"});

  // Name the set of contigs we are over.
  built.metadata.setContigs({"short", "alt1", "alt2", "chr"});

  gbwt::PathName p_short;
  p_short.sample = 0;
  p_short.contig = 0;
  p_short.phase = 0;
  p_short.count = 0;
  built.metadata.addPath(p_short);
  gbwt::PathName p_alt1;
  p_alt1.sample = 0;
  p_alt1.contig = 1;
  p_alt1.phase = 0;
  p_alt1.count = 0;
  built.metadata.addPath(p_alt1);
  gbwt::PathName p_alt2;
  p_alt2.sample = 0;
  p_alt2.contig = 2;
  p_alt2.phase = 0;
  p_alt2.count = 0;
  built.metadata.addPath(p_alt2);
  gbwt::PathName w_short1;
  w_short1.sample = 1;
  w_short1.contig = 3;
  w_short1.phase = 1;
  w_short1.count = 0;
  built.metadata.addPath(w_short1);
  gbwt::PathName w_alt;
  w_alt.sample = 2;
  w_alt.contig = 3;
  w_alt.phase = 0;
  w_alt.count = 0;
  built.metadata.addPath(w_alt);
  gbwt::PathName w_short2;
  w_short2.sample = 1;
  w_short2.contig = 3;
  w_short2.phase = 2;
  w_short2.count = 0;
  built.metadata.addPath(w_short2);

  // Record that we have 6 total haplotypes in the GBWT.
  built.metadata.setHaplotypes(6);

  return built;
}

// This should be identical to gfas/example_reference.gfa, modulo name ordering.
inline gbwt::GBWT
build_gbwt_example_reference()
{
  gbwt::GBWT built = build_gbwt_index_with_ref();
  
  built.addMetadata();

  // Name the set of samples, including a special one for generic paths
  built.metadata.setSamples({gbwtgraph::GENERIC_PATH_SAMPLE_NAME, "GRCh38", "GRCh37", "sample1", "CHM13"});
  
  // Set which are references
  built.tags.set(gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG, "GRCh37" + std::string(1, gbwtgraph::REFERENCE_SAMPLE_LIST_SEPARATOR) + "GRCh38");

  // Name the set of contigs we are over.
  built.metadata.setContigs({"chr1", "coolgene"});

  gbwt::PathName p_grch38_chr1;
  p_grch38_chr1.sample = 1;
  p_grch38_chr1.contig = 0;
  p_grch38_chr1.phase = 0;
  p_grch38_chr1.count = 0;
  built.metadata.addPath(p_grch38_chr1);
  gbwt::PathName p_grch37_chr1;
  p_grch37_chr1.sample = 2;
  p_grch37_chr1.contig = 0;
  p_grch37_chr1.phase = 0;
  p_grch37_chr1.count = 0;
  built.metadata.addPath(p_grch37_chr1);
  gbwt::PathName p_coolgene;
  p_coolgene.sample = 0;
  p_coolgene.contig = 1;
  p_coolgene.phase = 0;
  p_coolgene.count = 0;
  built.metadata.addPath(p_coolgene);
  gbwt::PathName w_sample1_1;
  w_sample1_1.sample = 3;
  w_sample1_1.contig = 0;
  w_sample1_1.phase = 1;
  w_sample1_1.count = 0;
  built.metadata.addPath(w_sample1_1);
  gbwt::PathName w_chm13;
  w_chm13.sample = 4;
  w_chm13.contig = 0;
  w_chm13.phase = 0;
  w_chm13.count = 0;
  built.metadata.addPath(w_chm13);
  gbwt::PathName w_sample1_2;
  w_sample1_2.sample = 3;
  w_sample1_2.contig = 0;
  w_sample1_2.phase = 2;
  w_sample1_2.count = 0;
  built.metadata.addPath(w_sample1_2);

  // Record that we have 6 total haplotypes in the GBWT.
  built.metadata.setHaplotypes(6);

  return built;
}

// Build a GBWT that contains 3 paths, two of which are named, reference paths in the metadata.
inline gbwt::GBWT
build_gbwt_index_with_named_paths()
{
  // Start with the no-metadata GBWT we usually use, but with some empty paths.
  gbwt::GBWT built = build_gbwt_index(2);

  built.addMetadata();

  // Name the set of samples, including a special one for generic paths
  built.metadata.setSamples({gbwtgraph::GENERIC_PATH_SAMPLE_NAME, "Jouni Sirén", "GRCh38"});

  // Name the set of contigs we are over.
  built.metadata.setContigs({"chr1", "chr2", "empty1", "empty2"});

  // We have a ref path on the first contig, and a ref path on the second
  // contig, and a sample path on the first contig.
  gbwt::PathName ref1;
  ref1.sample = 0;
  ref1.contig = 0;
  ref1.phase = 0;
  ref1.count = 0;
  built.metadata.addPath(ref1);
  gbwt::PathName ref2;
  ref2.sample = 0;
  ref2.contig = 1;
  ref2.phase = 0;
  ref2.count = 0;
  built.metadata.addPath(ref2);
  gbwt::PathName sample1;
  sample1.sample = 1;
  sample1.contig = 0;
  sample1.phase = 0;
  sample1.count = 0;
  built.metadata.addPath(sample1);
  gbwt::PathName empty1;
  empty1.sample = 0;
  empty1.contig = 2;
  empty1.phase = 0;
  empty1.count = 0;
  built.metadata.addPath(empty1);
  gbwt::PathName empty2;
  empty2.sample = 2;
  empty2.contig = 3;
  empty2.phase = gbwtgraph::GBWTGraph::NO_PHASE;
  empty2.count = 0;
  built.metadata.addPath(empty2);


  // Record that we have 5 total haplotypes in the GBWT.
  built.metadata.setHaplotypes(5);

  // Make GRCh38 a reference
  built.tags.set(gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG, "GRCh38");

  return built;
}

// Dump a GBWT and its metadata.
inline void
dump_gbwt(const gbwt::GBWT& built)
{
  // Check to make sure threads and metadata match.
  if (built.hasMetadata())
  {
    for(size_t path_num = 0; path_num < built.metadata.paths(); path_num++)
    {
      auto extracted = built.extract(gbwt::Path::encode(path_num, false));
      std::cerr << "GBWT stored forward path " << path_num << std::endl;
      auto path_name = built.metadata.path(path_num);
      std::cerr << "\tSample " << path_name.sample
        << " Contig " << path_name.contig
        << " Phase " << path_name.phase
        << " Count " << path_name.count << std::endl;
      for(auto& gbwt_node : extracted)
      {
        std::cerr << "\t" << gbwt::Node::id(gbwt_node) << " " << gbwt::Node::is_reverse(gbwt_node) << std::endl;
      }
    }
  }
}

// Builds a NaiveGraph that could contain `alt_path` and `short_path`.
inline gbwtgraph::NaiveGraph
build_naive_graph(bool with_translation)
{
  gbwtgraph::NaiveGraph graph;

  if(with_translation)
  {
    std::string seq = "GATGGGTACAA";
    graph.translate_segment("s1", std::string_view(seq.data() + 0, 1), 3);
    graph.translate_segment("s2", std::string_view(seq.data() + 1, 1), 3);
    graph.translate_segment("s3", std::string_view(seq.data() + 2, 1), 3);
    graph.translate_segment("s4", std::string_view(seq.data() + 3, 4), 3);
    graph.translate_segment("s5", std::string_view(seq.data() + 7, 1), 3);
    graph.translate_segment("s6", std::string_view(seq.data() + 8, 1), 3);
    graph.translate_segment("s7", std::string_view(seq.data() + 9, 1), 3);
    graph.translate_segment("s8", std::string_view(seq.data() + 10, 1), 3);
  }
  else
  {
    graph.create_node(1, "G");
    graph.create_node(2, "A");
    graph.create_node(3, "T");
    graph.create_node(4, "GGG");
    graph.create_node(5, "T");
    graph.create_node(6, "A");
    graph.create_node(7, "C");
    graph.create_node(8, "A");
    graph.create_node(9, "A");
  }

  // The paths are 1, 2, 4, 5, 6, 8, 9 and 1, 4, 5, 6, 7, 9.
  graph.create_edge(gbwt::Node::encode(1, false), gbwt::Node::encode(2, false));
  graph.create_edge(gbwt::Node::encode(2, false), gbwt::Node::encode(4, false));
  graph.create_edge(gbwt::Node::encode(4, false), gbwt::Node::encode(5, false));
  graph.create_edge(gbwt::Node::encode(5, false), gbwt::Node::encode(6, false));
  graph.create_edge(gbwt::Node::encode(6, false), gbwt::Node::encode(7, false));
  graph.create_edge(gbwt::Node::encode(6, false), gbwt::Node::encode(8, false));
  graph.create_edge(gbwt::Node::encode(7, false), gbwt::Node::encode(9, false));
  graph.create_edge(gbwt::Node::encode(8, false), gbwt::Node::encode(9, false));
  graph.remove_duplicate_edges();

  return graph;
}

//------------------------------------------------------------------------------

// Helpers for k-mer / minimizer indexes.

typedef std::pair<gbwtgraph::Position, std::vector<std::uint64_t>> owned_value_type;
typedef std::vector<std::uint64_t> owned_multi_value_type;

inline owned_value_type
create_value(gbwtgraph::KmerEncoding::value_type value, size_t payload_size)
{
  return std::make_pair(value.first, std::vector<std::uint64_t>(value.second, value.second + payload_size));
}

inline owned_value_type
create_value(pos_t pos, size_t payload_size, std::uint64_t payload)
{
  return std::make_pair(gbwtgraph::Position(pos), std::vector<std::uint64_t>(payload_size, payload));
}

inline void
append_value(owned_multi_value_type& values, const owned_value_type& value, size_t payload_size)
{
  constexpr size_t POS_SIZE = sizeof(gbwtgraph::Position) / sizeof(std::uint64_t);
  size_t offset = values.size();
  values.resize(offset + POS_SIZE + payload_size);
  value.first.write(values.data() + offset);
  offset += POS_SIZE;
  for(size_t j = 0; j < payload_size; j++) { values[offset + j] = value.second[j]; }
}

// If with_paths is true, the last word of payload encodes the set of paths that contain the hit.
// This may be a superset of the set encoded in the truth value.
// If paths A and B both contain subpath X, the kmer corresponding to X may or may not be a
// minimizer, depending on the context.
inline bool
same_value(gbwtgraph::KmerEncoding::value_type value, const owned_value_type& truth, size_t payload_size, bool with_paths)
{
  if(value.first != truth.first) { std::cerr << "Wrong pos" << std::endl; return false; }
  std::vector<std::uint64_t> payload(value.second, value.second + payload_size);

  if(with_paths && payload_size > 0)
  {
    if((payload.back() & truth.second.back()) != truth.second.back())
    {
      std::cerr << "Wrong paths in payload" << std::endl;
      std::cerr << "  Got: " << payload.back() << std::endl;
      std::cerr << "  Expected at least: " << truth.second.back() << std::endl;
      return false;
    }
    payload.back() = truth.second.back();
  }

  if(payload != truth.second)
  {
    std::cerr << "Wrong payload" << std::endl;
    std::cerr << "  Got:";
    for(auto p : payload) { std::cerr << " " << p; }
    std::cerr << std::endl;
    std::cerr << "  Expected:";
    for(auto p : truth.second) { std::cerr << " " << p; }
    std::cerr << std::endl;
    return false;
  }
  return true;
}

// See above for with_paths.
inline bool
same_values(gbwtgraph::KmerEncoding::multi_value_type values, const std::set<owned_value_type>& truth, size_t payload_size, bool with_paths)
{
  constexpr size_t POS_SIZE = sizeof(gbwtgraph::Position) / sizeof(std::uint64_t);
  if(values.second != truth.size())
  {
    std::cerr << "Wrong size: " << values.second << " != " << truth.size() << std::endl;
    return false;
  }

  size_t value_offset = 0;
  for(auto& correct : truth)
  {
    gbwtgraph::Position pos(values.first[value_offset]);
    if(pos != correct.first) { std::cerr << "Wrong pos" << std::endl; return false; }
    value_offset += POS_SIZE;
    std::vector<std::uint64_t> payload(values.first + value_offset, values.first + value_offset + payload_size);

    if(with_paths && payload_size > 0)
    {
      if((payload.back() & correct.second.back()) != correct.second.back())
      {
        std::cerr << "Wrong paths in payload" << std::endl;
        std::cerr << "  Got: " << payload.back() << std::endl;
        std::cerr << "  Expected at least: " << correct.second.back() << std::endl;
        return false;
      }
      payload.back() = correct.second.back();
    }

    if(payload != correct.second)
    {
      std::cerr << "Wrong payload" << std::endl;
      std::cerr << "  Got:";
      for(auto p : payload) { std::cerr << " " << p; }
      std::cerr << std::endl;
      std::cerr << "  Expected:";
      for(auto p : correct.second) { std::cerr << " " << p; }
      std::cerr << std::endl;
      return false;
    }
    value_offset += payload_size;
  }

  return true;
}

template<class KeyType>
void insert_value(gbwtgraph::KmerIndex<KeyType>& index, KeyType key, const owned_value_type& value)
{
  index.insert(key, std::make_pair(value.first, value.second.data()));
}

template<class KeyType>
void insert_value
(
  gbwtgraph::MinimizerIndex<KeyType>& index,
  typename gbwtgraph::MinimizerIndex<KeyType>::minimizer_type key, const owned_value_type& value
)
{
  index.insert(key, std::make_pair(value.first, value.second.data()));
}

template<class KeyType>
gbwtgraph::Kmer<KeyType>
get_minimizer(KeyType key, gbwtgraph::offset_type offset = 0, bool orientation = false)
{
  return { key, key.hash(), offset, orientation };
}

template<class KeyType>
gbwtgraph::Kmer<KeyType>
get_minimizer(std::string key, gbwtgraph::offset_type offset = 0, bool orientation = false)
{
  return get_minimizer(KeyType::encode(key), offset, orientation);
}

//------------------------------------------------------------------------------

inline std::string
path_name_to_string(const gbwt::PathName& path_name)
{
  std::string result = "(" + std::to_string(path_name.sample)
    + ", " + std::to_string(path_name.contig)
    + ", " + std::to_string(path_name.phase)
    + ", " + std::to_string(path_name.count) + ")";
  return result;
}

inline std::string
path_to_string(const gbwtgraph::GBWTGraph& graph, const gbwt::vector_type& path)
{
  std::string str;
  for(gbwt::node_type node : path)
  {
    std::string_view view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(node));
    str.append(view.data(), view.size());
  }
  return str;
}

//------------------------------------------------------------------------------

} // anonymous namespace

#endif // GBWTGRAPH_TESTS_SHARED_H
