#include <gtest/gtest.h>

#include <gbwtgraph/subgraph.h>
#include <gbwtgraph/gfa.h>
#include <gbwtgraph/internal.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

GBZ
build_gbz(const std::string& graph_name)
{
  auto parse = gfa_to_gbwt(graph_name);
  return GBZ(parse.first, parse.second);
}

//------------------------------------------------------------------------------

class PathIndexTest : public ::testing::Test
{
public:
  std::vector<std::string> graphs;
  std::vector<std::vector<size_t>> indexed_path_lengths;

  PathIndexTest()
  {
    this->graphs = { "gfas/default.gfa", "gfas/components_ref.gfa" };
    this->indexed_path_lengths =
    {
      { 5, 4 },
      { 5, 5 }
    };
  }
};

TEST_F(PathIndexTest, BuildIndex)
{
  for(size_t i = 0; i < this->graphs.size(); i++)
  {
    GBZ gbz = build_gbz(this->graphs[i]);
    PathIndex index(gbz);
    ASSERT_EQ(index.paths(), this->indexed_path_lengths[i].size()) << "Invalid number of indexed paths for graph " << this->graphs[i];
    for(size_t j = 0; j < this->indexed_path_lengths[i].size(); j++)
    {
      path_handle_t handle = handlegraph::as_path_handle(j);
      EXPECT_EQ(index.path_length(handle), this->indexed_path_lengths[i][j]) << "Invalid length of path " << j << " for graph " << this->graphs[i];
    }
  }
}

TEST_F(PathIndexTest, IndexedPaths)
{
  for(auto& graph_name : this->graphs)
  {
    GBZ gbz = build_gbz(graph_name);
    for(size_t interval = 0; interval < 10; interval++)
    {
      PathIndex index(gbz, interval);
      for(size_t path_id = 0; path_id < index.paths(); path_id++)
      {
        size_t length = 0;
        path_handle_t handle = handlegraph::as_path_handle(path_id);
        std::vector<std::pair<size_t, gbwt::edge_type>> reference_positions = sample_path_positions(gbz, handle, interval, &length);
        size_t pos_offset = 0;
        for(size_t seq_offset = 0; seq_offset <= length; seq_offset++)
        {
          auto position = index.sampled_position(handle, seq_offset);
          if(pos_offset + 1 < reference_positions.size() && reference_positions[pos_offset + 1].first <= seq_offset) { pos_offset++; }
          ASSERT_EQ(position, reference_positions[pos_offset]) << "Invalid sampled position for path " << path_id << ", offset " << seq_offset << " in graph " << graph_name;
        }
      }
    }
  }
}

TEST_F(PathIndexTest, NonIndexedPaths)
{
  for(auto& graph_name : this->graphs)
  {
    GBZ gbz = build_gbz(graph_name);
    for(size_t interval = 0; interval < 10; interval++)
    {
      PathIndex index(gbz, interval);
      for(size_t path_id = index.paths(); path_id < gbz.index.metadata.paths(); path_id++)
      {
        path_handle_t handle = handlegraph::as_path_handle(path_id);
        for(size_t seq_offset = 0; seq_offset <= 10; seq_offset++)
        {
          auto position = index.sampled_position(handle, seq_offset);
          ASSERT_EQ(position.second, gbwt::invalid_edge()) << "Found a sampled position for non-indexed path " << path_id << ", offset " << seq_offset << " in graph " << graph_name;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

class AlignDivergingTest : public ::testing::Test
{
public:
  gbwt::GBWT index;
  GBWTGraph graph;

  void SetUp() override
  {
    std::vector<gbwt::vector_type> paths
    {
      forward_path({ 1, 2, 3, 4, 5, 6 })
    };
    this->index = build_gbwt(paths);

    std::vector<std::pair<nid_t, std::string>> nodes
    {
      { 1, "A" },
      { 2, "C" },
      { 3, "AA" },
      { 4, "CC" },
      { 5, "ACA" },
      { 6, "CAC" }
    };
    SequenceSource source;
    for(const auto& node : nodes)
    {
      source.add_node(node.first, node.second);
    }
    this->graph = GBWTGraph(this->index, source);
  }

  static gbwt::vector_type forward_path(const std::vector<nid_t>& nodes)
  {
    gbwt::vector_type result;
    for(nid_t node : nodes)
    {
      result.push_back(gbwt::Node::encode(node, false));
    }
    return result;
  }

  void check_edits(
    const gbwt::vector_type& path, const gbwt::vector_type& reference,
    const std::vector<std::pair<char, size_t>>& truth,
    const std::string& name
  ) const
  {
    std::vector<std::pair<char, size_t>> edits;
    align_diverging(
      this->graph,
      get_subpath(path, 0, path.size()),
      get_subpath(reference, 0, reference.size()),
      edits
    );
    ASSERT_EQ(edits.size(), truth.size()) << "Invalid number of edits for " << name;
    for(size_t i = 0; i < edits.size(); i++)
    {
      EXPECT_EQ(edits[i], truth[i]) << "Invalid edit " << i << " for " << name;
    }
  }
};

TEST_F(AlignDivergingTest, SpecialCases)
{
  gbwt::vector_type empty;
  gbwt::vector_type non_empty = forward_path({ 5, 6 });

  // (empty, empty)
  {
    std::vector<std::pair<char, size_t>> truth;
    this->check_edits(empty, empty, truth, "empty paths");
  }

  // (empty, non-empty)
  {
    std::vector<std::pair<char, size_t>> truth { { 'D', 6 } };
    this->check_edits(empty, non_empty, truth, "empty vs. non-empty paths");
  }

  // (non-empty, empty)
  {
    std::vector<std::pair<char, size_t>> truth { { 'I', 6 } };
    this->check_edits(non_empty, empty, truth, "non-empty vs. empty paths");
  }

  // (non-empty, non-empty)
  {
    std::vector<std::pair<char, size_t>> truth { { 'M', 6 } };
    this->check_edits(non_empty, non_empty, truth, "identical paths");
  }

  // Identical bases, different paths.
  {
    gbwt::vector_type path = forward_path({ 1, 5, 3 }); // AACAAA
    gbwt::vector_type reference = forward_path({ 3, 2, 1, 3 }); // AACAAA
    std::vector<std::pair<char, size_t>> truth { { 'M', 6 } };
    this->check_edits(path, reference, truth, "identical bases, different paths");
  }

  // Prefix + suffix length exceeds the length of the shorter path.
  {
    gbwt::vector_type short_path = forward_path({ 5, 5 }); // ACAACA
    gbwt::vector_type long_path = forward_path({ 1, 2, 3, 5, 3, 2, 1 }); // ACAAACAAACA
    std::vector<std::pair<char, size_t>> path_is_short { { 'M', 4 }, { 'D', 5 }, { 'M', 2 } };
    std::vector<std::pair<char, size_t>> path_is_long { { 'M', 4 }, { 'I', 5 }, { 'M', 2 } };
    this->check_edits(short_path, long_path, path_is_short, "prefix + suffix length exceeds path length");
    this->check_edits(long_path, short_path, path_is_long, "prefix + suffix length exceeds reference length");
  }
}

TEST_F(AlignDivergingTest, NoPrefixNoSuffix)
{
  // Insertion and deletion are in SpecialCases.

  // Mismatch + insertion.
  {
    gbwt::vector_type path = forward_path({ 1, 2, 5 }); // ACACA
    gbwt::vector_type reference = forward_path({ 4 }); // CC
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'I', 3 } };
    this->check_edits(path, reference, truth, "mismatch + insertion");
  }

  // Mismatch + deletion.
  {
    gbwt::vector_type path = forward_path({ 3 }); // AA
    gbwt::vector_type reference = forward_path({ 2, 1, 6 }); // CACAC
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'D', 3 } };
    this->check_edits(path, reference, truth, "mismatch + deletion");
  }

  // Insertion + deletion.
  {
    gbwt::vector_type path = forward_path({ 1, 2, 5, 3 }); // ACACAAA
    gbwt::vector_type reference = forward_path({ 4, 2, 1, 6 }); // CCCACAC
    std::vector<std::pair<char, size_t>> truth { { 'I', 7 }, { 'D', 7 } };
    this->check_edits(path, reference, truth, "insertion + deletion");
  }
}

TEST_F(AlignDivergingTest, NoPrefixWithSuffix)
{
  // Insertion.
  {
    gbwt::vector_type path = forward_path({ 1, 2, 5 }); // ACACA
    gbwt::vector_type reference = forward_path({ 2, 1 }); // CA
    std::vector<std::pair<char, size_t>> truth { { 'I', 3 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "insertion");
  }

  // Deletion.
  {
    gbwt::vector_type path = forward_path({ 1, 2 }); // AC
    gbwt::vector_type reference = forward_path({ 2, 1, 6 }); // CACAC
    std::vector<std::pair<char, size_t>> truth { { 'D', 3 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "deletion");
  }

  // Mismatch + insertion.
  {
    gbwt::vector_type path = forward_path({ 5, 2, 5 }); // ACACACA
    gbwt::vector_type reference = forward_path({ 4, 2, 1 }); // CCCA
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'I', 3 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "mismatch + insertion");
  }

  // Mismatch + deletion.
  {
    gbwt::vector_type path = forward_path({ 3, 1, 2 }); // AAAC
    gbwt::vector_type reference = forward_path({ 6, 1, 6 }); // CACACAC
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'D', 3 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "mismatch + deletion");
  }

  // Insertion + deletion.
  {
    gbwt::vector_type path = forward_path({ 3, 2, 5, 5 }); // AACACAACA
    gbwt::vector_type reference = forward_path({ 4, 2, 1, 6, 1 }); // CCCACACA
    std::vector<std::pair<char, size_t>> truth { { 'I', 6 }, { 'D', 5 }, { 'M', 3 } };
    this->check_edits(path, reference, truth, "insertion + deletion");
  }
}

TEST_F(AlignDivergingTest, WithPrefixNoSuffix)
{
  // Insertion.
  {
    gbwt::vector_type path = forward_path({ 5, 2, 1 }); // ACACA
    gbwt::vector_type reference = forward_path({ 1, 2 }); // AC
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'I', 3 } };
    this->check_edits(path, reference, truth, "insertion");
  }

  // Deletion.
  {
    gbwt::vector_type path = forward_path({ 2, 1 }); // CA
    gbwt::vector_type reference = forward_path({ 6, 1, 2 }); // CACAC
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'D', 3 } };
    this->check_edits(path, reference, truth, "deletion");
  }

  // Mismatch + insertion.
  {
    gbwt::vector_type path = forward_path({ 5, 2, 5 }); // ACACACA
    gbwt::vector_type reference = forward_path({ 1, 2, 4 }); // ACCC
    std::vector<std::pair<char, size_t>> truth { { 'M', 4 }, { 'I', 3 } };
    this->check_edits(path, reference, truth, "mismatch + insertion");
  }

  // Mismatch + deletion.
  {
    gbwt::vector_type path = forward_path({ 2, 1, 3 }); // CAAA
    gbwt::vector_type reference = forward_path({ 6, 1, 6 }); // CACACAC
    std::vector<std::pair<char, size_t>> truth { { 'M', 4 }, { 'D', 3 } };
    this->check_edits(path, reference, truth, "mismatch + deletion");
  }

  // Insertion + deletion.
  {
    gbwt::vector_type path = forward_path({ 3, 2, 5, 5 }); // AACACAACA
    gbwt::vector_type reference = forward_path({ 1, 1, 4, 2, 1, 6, 2 }); // AACCCACACCC
    std::vector<std::pair<char, size_t>> truth { { 'M', 3 }, { 'I', 6 }, { 'D', 7 } };
    this->check_edits(path, reference, truth, "insertion + deletion");
  }
}

TEST_F(AlignDivergingTest, WithPrefixWithSuffix)
{
  // Insertion.
  {
    gbwt::vector_type path = forward_path({ 5, 2, 5 }); // ACACACA
    gbwt::vector_type reference = forward_path({ 1, 4, 1 }); // ACCA
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'I', 3 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "insertion");
  }

  // Deletion.
  {
    gbwt::vector_type path = forward_path({ 2, 3, 2 }); // CAAC
    gbwt::vector_type reference = forward_path({ 6, 1, 6 }); // CACACAC
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'D', 3 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "deletion");
  }

  // Mismatch + insertion.
  {
    gbwt::vector_type path = forward_path({ 1, 2, 4, 6, 2, 1 }); // ACCCCACCA
    gbwt::vector_type reference = forward_path({ 5, 5 }); // ACAACA
    std::vector<std::pair<char, size_t>> truth { { 'M', 4 }, { 'I', 3 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "mismatch + insertion");
  }

  // Mismatch + deletion.
  {
    gbwt::vector_type path = forward_path({ 2, 3, 3, 2 }); // CAAAAC
    gbwt::vector_type reference = forward_path({ 2, 5, 3, 6 }); // CACAAACAC
    std::vector<std::pair<char, size_t>> truth { { 'M', 4 }, { 'D', 3 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "mismatch + deletion");
  }

  // Insertion + deletion.
  {
    gbwt::vector_type path = forward_path({ 1, 5, 4, 3, 4 }); // AACACCAACC
    gbwt::vector_type reference = forward_path({ 3, 1, 5, 1, 4, 2 }); // AAAACAACCC
    std::vector<std::pair<char, size_t>> truth { { 'M', 2 }, { 'I', 6 }, { 'D', 6 }, { 'M', 2 } };
    this->check_edits(path, reference, truth, "insertion + deletion");
  }
}

//------------------------------------------------------------------------------

class SubgraphQueryTest : public ::testing::Test
{
public:
  std::vector<std::string> graphs;
  std::vector<std::string> reference_samples;
  std::vector<size_t> reference_haplotypes;

  SubgraphQueryTest()
  {
    this->graphs = { "gfas/default.gfa", "gfas/components_ref.gfa" };
    this->reference_samples = { REFERENCE_PATH_SAMPLE_NAME, "ref" };
    this->reference_haplotypes = { GBWTGraph::NO_PHASE, 0 };
  }

  void find_path(const GBZ& gbz, size_t graph_number, const std::string& contig_name, path_handle_t& out) const
  {
    const gbwt::Metadata& metadata = gbz.index.metadata;
    std::vector<gbwt::size_type> path_ids = metadata.findPaths(metadata.sample(this->reference_samples[graph_number]), metadata.contig(contig_name));
    ASSERT_EQ(path_ids.size(), size_t(1)) << "Found " << path_ids.size() << " reference paths for contig " << contig_name << " in graph " << this->graphs[graph_number];
    out = gbz.graph.path_to_handle(path_ids.front());
  }

  Subgraph find_subgraph(const GBZ& gbz, const SubgraphQuery& query) const
  {
    PathIndex index(gbz);
    return Subgraph(gbz, &index, query);
  }

  void check_subgraph(size_t graph_number, const Subgraph& subgraph, const SubgraphQuery& query, const std::set<nid_t>& nodes, size_t path_count) const
  {
    ASSERT_EQ(subgraph.nodes.size(), nodes.size()) << "Invalid number of nodes for query " << query.to_string() << " in graph " << this->graphs[graph_number];
    for(nid_t node : nodes)
    {
      EXPECT_TRUE(subgraph.nodes.find(node) != subgraph.nodes.end()) << "Node " << node << " not found by query " << query.to_string() << " in graph " << this->graphs[graph_number];
    }
    ASSERT_EQ(subgraph.paths.size(), path_count) << "Invalid number of paths for query " << query.to_string() << " in graph " << this->graphs[graph_number];
  }

  void check_gfa(const GBZ& gbz, size_t graph_number, const Subgraph& subgraph, const SubgraphQuery& query, const std::vector<std::string>& gfa) const
  {
    std::ostringstream out;
    subgraph.to_gfa(gbz, out);
    std::string result = out.str();
    std::istringstream stream(result);
    std::vector<std::string> lines;
    std::string line;
    while(std::getline(stream, line)) { lines.push_back(line); }

    ASSERT_EQ(lines.size(), gfa.size()) << "Invalid number of lines in GFA output for query " << query.to_string() << " in graph " << this->graphs[graph_number];
    for(size_t i = 0; i < lines.size(); i++)
    {
      EXPECT_EQ(lines[i], gfa[i]) << "Invalid GFA line " << i << " for query " << query.to_string() << " in graph " << this->graphs[graph_number];
    }
  }
};

TEST_F(SubgraphQueryTest, AllHaplotypes)
{
  // All haplotypes within 1 bp of contig A, offset 2.
  std::string contig_name = "A";
  size_t offset = 2, context = 1;
  std::set<nid_t> nodes { 12, 13, 14, 15, 16 };
  size_t path_count = 3;

  std::vector<std::string> gfa
  {
    "H\tVN:Z:1.1\tRS:Z:", // Append reference sample to line 0.
    "S\t12\tA",
    "S\t13\tT",
    "S\t14\tT",
    "S\t15\tA",
    "S\t16\tC",
    "L\t12\t+\t14\t+\t0M",
    "L\t13\t+\t14\t+\t0M",
    "L\t14\t+\t15\t+\t0M",
    "L\t14\t+\t16\t+\t0M",
    "", // Set reference path manually on line 10.
    "W\tunknown\t1\tA\t0\t3\t>12>14>15\tCG:Z:3M", 
    "W\tunknown\t2\tA\t0\t3\t>13>14>16\tCG:Z:3M", 
  };
  std::vector<std::string> reference_paths
  {
    "W\t" + this->reference_samples[0] + "\t0\tA\t1\t4\t>12>14>15",
    "W\t" + this->reference_samples[1] + "\t0\tA\t1\t4\t>13>14>15",
  };

  for(size_t i = 0; i < this->graphs.size(); i++)
  {
    GBZ gbz = build_gbz(this->graphs[i]);
    gbwt::FullPathName path_name { this->reference_samples[i], contig_name, this->reference_haplotypes[i], 0 };
    SubgraphQuery query = SubgraphQuery::path_offset(path_name, offset, context, SubgraphQuery::all_haplotypes);
    Subgraph subgraph = this->find_subgraph(gbz, query);
    this->check_subgraph(i, subgraph, query, nodes, path_count);

    std::vector<std::string> true_gfa = gfa;
    true_gfa[0] += this->reference_samples[i];
    true_gfa[10] = reference_paths[i];
    this->check_gfa(gbz, i, subgraph, query, true_gfa);
  }
}

TEST_F(SubgraphQueryTest, DistinctHaplotypes)
{
  // Distinct haplotypes within 1 bp of contig A, offset 2.
  std::string contig_name = "A";
  size_t offset = 2, context = 1;
  std::set<nid_t> nodes { 12, 13, 14, 15, 16 };
  std::vector<size_t> path_counts { 2, 3 };

  std::vector<std::string> gfa
  {
    "H\tVN:Z:1.1\tRS:Z:", // Append reference sample to line 0.
    "S\t12\tA",
    "S\t13\tT",
    "S\t14\tT",
    "S\t15\tA",
    "S\t16\tC",
    "L\t12\t+\t14\t+\t0M",
    "L\t13\t+\t14\t+\t0M",
    "L\t14\t+\t15\t+\t0M",
    "L\t14\t+\t16\t+\t0M", // Append paths manually.
  };
  std::vector<std::vector<std::string>> gfa_paths
  {
    {
      "W\t" + this->reference_samples[0] + "\t0\tA\t1\t4\t>12>14>15\tWT:i:2",
      "W\tunknown\t1\tA\t0\t3\t>13>14>16\tWT:i:1\tCG:Z:3M", 
    },
    {
      "W\t" + this->reference_samples[1] + "\t0\tA\t1\t4\t>13>14>15\tWT:i:1",
      "W\tunknown\t1\tA\t0\t3\t>12>14>15\tWT:i:1\tCG:Z:3M", 
      "W\tunknown\t2\tA\t0\t3\t>13>14>16\tWT:i:1\tCG:Z:3M", 
    }
  };

  for(size_t i = 0; i < this->graphs.size(); i++)
  {
    GBZ gbz = build_gbz(this->graphs[i]);
    gbwt::FullPathName path_name { this->reference_samples[i], contig_name, this->reference_haplotypes[i], 0 };
    SubgraphQuery query = SubgraphQuery::path_offset(path_name, offset, context, SubgraphQuery::distinct_haplotypes);
    Subgraph subgraph = this->find_subgraph(gbz, query);
    this->check_subgraph(i, subgraph, query, nodes, path_counts[i]);

    std::vector<std::string> true_gfa = gfa;
    true_gfa[0] += this->reference_samples[i];
    true_gfa.insert(true_gfa.end(), gfa_paths[i].begin(), gfa_paths[i].end());
    this->check_gfa(gbz, i, subgraph, query, true_gfa);
  }
}

TEST_F(SubgraphQueryTest, NodeBased)
{
  // Distinct haplotypes within 2 bp of node 13.
  nid_t node_id = 13;
  size_t context = 2;
  std::set<nid_t> nodes { 11, 12, 13, 14, 15, 16 };
  std::vector<size_t> path_counts { 2, 3 };

  std::vector<std::string> gfa
  {
    "H\tVN:Z:1.1",
    "S\t11\tG",
    "S\t12\tA",
    "S\t13\tT",
    "S\t14\tT",
    "S\t15\tA",
    "S\t16\tC",
    "L\t11\t+\t12\t+\t0M",
    "L\t11\t+\t13\t+\t0M",
    "L\t12\t+\t14\t+\t0M",
    "L\t13\t+\t14\t+\t0M",
    "L\t14\t+\t15\t+\t0M",
    "L\t14\t+\t16\t+\t0M", // Append paths manually.
  };
  std::vector<std::vector<std::string>> gfa_paths
  {
    {
      "W\tunknown\t1\tunknown\t0\t4\t>11>12>14>15\tWT:i:2",
      "W\tunknown\t2\tunknown\t0\t4\t>11>13>14>16\tWT:i:1",
    },
    {
      "W\tunknown\t1\tunknown\t0\t4\t>11>12>14>15\tWT:i:1", 
      "W\tunknown\t2\tunknown\t0\t4\t>11>13>14>15\tWT:i:1",
      "W\tunknown\t3\tunknown\t0\t4\t>11>13>14>16\tWT:i:1",
    }
  };

  for(size_t i = 0; i < this->graphs.size(); i++)
  {
    GBZ gbz = build_gbz(this->graphs[i]);
    SubgraphQuery query = SubgraphQuery::node(node_id, context, SubgraphQuery::distinct_haplotypes);
    Subgraph subgraph = this->find_subgraph(gbz, query);
    this->check_subgraph(i, subgraph, query, nodes, path_counts[i]);

    std::vector<std::string> true_gfa = gfa;
    true_gfa.insert(true_gfa.end(), gfa_paths[i].begin(), gfa_paths[i].end());
    this->check_gfa(gbz, i, subgraph, query, true_gfa);
  }
}

TEST_F(SubgraphQueryTest, ReferenceOnly)
{
  // Reference haplotype within 1 bp of contig A, offset 2.
  std::string contig_name = "A";
  size_t offset = 2, context = 1;
  std::set<nid_t> nodes { 12, 13, 14, 15, 16 };
  size_t path_count = 1;

  std::vector<std::string> gfa
  {
    "H\tVN:Z:1.1\tRS:Z:", // Append reference sample to line 0.
    "S\t12\tA",
    "S\t13\tT",
    "S\t14\tT",
    "S\t15\tA",
    "S\t16\tC",
    "L\t12\t+\t14\t+\t0M",
    "L\t13\t+\t14\t+\t0M",
    "L\t14\t+\t15\t+\t0M",
    "L\t14\t+\t16\t+\t0M", // Append reference path manually.
  };
  std::vector<std::string> reference_paths
  {
    "W\t" + this->reference_samples[0] + "\t0\tA\t1\t4\t>12>14>15",
    "W\t" + this->reference_samples[1] + "\t0\tA\t1\t4\t>13>14>15",
  };

  for(size_t i = 0; i < this->graphs.size(); i++)
  {
    GBZ gbz = build_gbz(this->graphs[i]);
    gbwt::FullPathName path_name { this->reference_samples[i], contig_name, this->reference_haplotypes[i], 0 };
    SubgraphQuery query = SubgraphQuery::path_offset(path_name, offset, context, SubgraphQuery::reference_only);
    Subgraph subgraph = this->find_subgraph(gbz, query);
    this->check_subgraph(i, subgraph, query, nodes, path_count);

    std::vector<std::string> true_gfa = gfa;
    true_gfa[0] += this->reference_samples[i];
    true_gfa.push_back(reference_paths[i]);
    this->check_gfa(gbz, i, subgraph, query, true_gfa);
  }
}

TEST_F(SubgraphQueryTest, ContigB)
{
  // Distinct haplotypes within 1 bp of contig B, offset 2.
  std::string contig_name = "B";
  size_t offset = 2, context = 1;
  std::set<nid_t> nodes { 22, 23, 24, 25 };
  std::vector<size_t> path_counts { 2, 3 };

  std::vector<std::string> gfa
  {
    "H\tVN:Z:1.1\tRS:Z:", // Append reference sample to line 0.
    "S\t22\tA",
    "S\t23\tT",
    "S\t24\tT",
    "S\t25\tA",
    "L\t22\t+\t24\t+\t0M",
    "L\t23\t+\t24\t-\t0M",
    "L\t24\t+\t25\t+\t0M", // Append paths manually.
  };
  std::vector<std::vector<std::string>> gfa_paths
  {
    {
      "W\t" + this->reference_samples[0] + "\t0\tB\t1\t4\t>22>24>25\tWT:i:2",
      "W\tunknown\t1\tB\t0\t3\t>22>24<23\tWT:i:1\tCG:Z:3M",
    },
    {
      "W\t" + this->reference_samples[1] + "\t0\tB\t1\t4\t>23<24<22\tWT:i:1",
      "W\tunknown\t1\tB\t0\t3\t>22>24<23\tWT:i:2\tCG:Z:3M",
      "W\tunknown\t2\tB\t0\t3\t>22>24>25\tWT:i:1\tCG:Z:3M",
    }
  };

  for(size_t i = 0; i < this->graphs.size(); i++)
  {
    GBZ gbz = build_gbz(this->graphs[i]);
    gbwt::FullPathName path_name { this->reference_samples[i], contig_name, this->reference_haplotypes[i], 0 };
    SubgraphQuery query = SubgraphQuery::path_offset(path_name, offset, context, SubgraphQuery::distinct_haplotypes);
    Subgraph subgraph = this->find_subgraph(gbz, query);
    this->check_subgraph(i, subgraph, query, nodes, path_counts[i]);

    std::vector<std::string> true_gfa = gfa;
    true_gfa[0] += this->reference_samples[i];
    true_gfa.insert(true_gfa.end(), gfa_paths[i].begin(), gfa_paths[i].end());
    this->check_gfa(gbz, i, subgraph, query, true_gfa);
  }
}

TEST_F(SubgraphQueryTest, Fragmented)
{
  GBZ gbz = build_gbz("gfas/fragmented.gfa");
  gbwt::FullPathName path_name { "ref", "contig", 0, 0 };
  size_t context = 3;

  std::vector<size_t> offsets { 4, 11 };
  std::vector<std::set<nid_t>> nodes
  {
    { 1, 2, 3, 4, 5 },
    { 6, 7, 8, 9 }
  };
  std::vector<size_t> path_counts { 3, 2 };
  std::vector<std::vector<std::string>> gfas
  {
    {
      "H\tVN:Z:1.1\tRS:Z:ref",
      "S\t1\tGAT",
      "S\t2\tTA",
      "S\t3\tCA",
      "S\t4\tC",
      "S\t5\tG",
      "L\t1\t+\t2\t+\t0M",
      "L\t2\t+\t3\t+\t0M",
      "L\t3\t+\t4\t+\t0M",
      "L\t3\t+\t5\t+\t0M",
      "W\tref\t0\tcontig\t0\t7\t>1>2>3\tWT:i:1",
      "W\tunknown\t1\tcontig\t0\t8\t>1>2>3>4\tWT:i:1\tCG:Z:7M1I",
      "W\tunknown\t2\tcontig\t0\t8\t>1>2>3>5\tWT:i:1\tCG:Z:7M1I"
    },
    {
      "H\tVN:Z:1.1\tRS:Z:ref",
      "S\t6\tATTA",
      "S\t7\tC",
      "S\t8\tG",
      "S\t9\tAT",
      "L\t6\t+\t7\t+\t0M",
      "L\t6\t+\t8\t+\t0M",
      "L\t7\t+\t9\t+\t0M",
      "L\t8\t+\t9\t+\t0M",
      "W\tref\t0\tcontig\t8\t15\t>6>7>9\tWT:i:2",
      "W\tunknown\t1\tcontig\t0\t7\t>6>8>9\tWT:i:1\tCG:Z:7M"
    }
  };

  for(size_t i = 0; i < offsets.size(); i++)
  {
    SubgraphQuery query = SubgraphQuery::path_offset(path_name, offsets[i], context, SubgraphQuery::distinct_haplotypes);
    Subgraph subgraph = this->find_subgraph(gbz, query);
    this->check_subgraph(0, subgraph, query, nodes[i], path_counts[i]);
    this->check_gfa(gbz, 0, subgraph, query, gfas[i]);
  }
}

//------------------------------------------------------------------------------

} // namespace
