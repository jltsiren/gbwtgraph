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

  void check_subgraph(const GBZ& gbz, size_t graph_number, const Subgraph& subgraph, const SubgraphQuery& query, const std::set<nid_t>& nodes, size_t path_count) const
  {
    ASSERT_EQ(subgraph.nodes.size(), nodes.size()) << "Invalid number of nodes for query " << query.to_string(gbz) << " in graph " << this->graphs[graph_number];
    for(nid_t node : nodes)
    {
      EXPECT_TRUE(subgraph.nodes.find(node) != subgraph.nodes.end()) << "Node " << node << " not found by query " << query.to_string(gbz) << " in graph " << this->graphs[graph_number];
    }
    ASSERT_EQ(subgraph.paths.size(), path_count) << "Invalid number of paths for query " << query.to_string(gbz) << " in graph " << this->graphs[graph_number];
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

    ASSERT_EQ(lines.size(), gfa.size()) << "Invalid number of lines in GFA output for query " << query.to_string(gbz) << " in graph " << this->graphs[graph_number];
    for(size_t i = 0; i < lines.size(); i++)
    {
      EXPECT_EQ(lines[i], gfa[i]) << "Invalid GFA line " << i << " for query " << query.to_string(gbz) << " in graph " << this->graphs[graph_number];
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
    this->check_subgraph(gbz, i, subgraph, query, nodes, path_count);

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
    this->check_subgraph(gbz, i, subgraph, query, nodes, path_counts[i]);

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
    this->check_subgraph(gbz, i, subgraph, query, nodes, path_counts[i]);

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
    this->check_subgraph(gbz, i, subgraph, query, nodes, path_count);

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
    this->check_subgraph(gbz, i, subgraph, query, nodes, path_counts[i]);

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
    this->check_subgraph(gbz, 0, subgraph, query, nodes[i], path_counts[i]);
    this->check_gfa(gbz, 0, subgraph, query, gfas[i]);
  }
}

//------------------------------------------------------------------------------

} // namespace
