#ifndef GBWTGRAPH_SUBGRAPH_H
#define GBWTGRAPH_SUBGRAPH_H

#include <gbwtgraph/gbz.h>

#include <sdsl/sd_vector.hpp>

#include <iostream>
#include <set>
#include <vector>

/*
  subgraph.h: Subgraph extraction from a GBZ graph as in GBZ-base.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// FIXME tests
/*
  A data structure for indexing reference and generic paths in a GBZ graph for
  random access by sequence offset. For each path, the index stores GBWT
  positions corresponding to node starts once every `sample_interval` bp
  (default 1024 bp).
*/
class PathIndex
{
public:
  // Sample node starts once every 1024 bp by default.
  constexpr static size_t DEFAULT_SAMPLE_INTERVAL = 1024;

  // Build the index for the paths in the given GBZ graph.
  // The work is parallelized over the paths using OpenMP threads.
  explicit PathIndex(const GBZ& gbz, size_t sample_interval = DEFAULT_SAMPLE_INTERVAL);

  PathIndex(const PathIndex& source) = default;
  PathIndex(PathIndex&& source) = default;
  PathIndex& operator=(const PathIndex& source) = default;
  PathIndex& operator=(PathIndex&& source) = default;

  // For each `NamedPath` in the graph, mark the sampled sequence positions.
  std::vector<sdsl::sd_vector<>> sequence_positions;

  // Sampled GBWT positions for each `NamedPath` in the graph.
  std::vector<std::vector<gbwt::edge_type>> gbwt_positions;

  size_t paths() const { return this->sequence_positions.size(); }
  size_t path_length(size_t path_id) const { return this->sequence_positions[path_id].size(); }

  // Returns the last sampled position at or before `offset` in the path,
  // as a pair (sequence offset, GBWT position). If there is no such path,
  // returns (0, gbwt::invalid_edge()).
  std::pair<size_t, gbwt::edge_type> sampled_position(size_t path_id, size_t offset) const;
};

//------------------------------------------------------------------------------

// FIXME document
struct SubgraphQuery
{
  enum QueryType { path_offset_query, node_query };

  enum HaplotypeOutput { all_haplotypes, distinct_haplotypes, reference_only };

  // Type of the query.
  QueryType type;

  // Which haplotypes to output.
  HaplotypeOutput haplotype_output;

  // Path handle for a reference or generic path.
  path_handle_t path;

  // Query offset on the path.
  size_t offset;

  // Node id for the node query.
  nid_t node_id;

  // Context length around the query position (in bp).
  size_t context;

  // Creates a subgraph query based on a path offset.
  static SubgraphQuery path_offset(path_handle_t path, size_t offset, size_t context, HaplotypeOutput output);

  // Creates a subgraph query based on a path interval [`from`, `to`).
  // Throws `std::runtime_error` if the inteval is empty.
  static SubgraphQuery path_interval(path_handle_t path, size_t from, size_t to, size_t context, HaplotypeOutput output);

  // Creates a subgraph query based on a node id.
  static SubgraphQuery node(nid_t node, size_t context, HaplotypeOutput output);

  // Returns a string representation of the query.
  std::string to_string(const GBZ& gbz) const;
};

//------------------------------------------------------------------------------

// class Subgraph
// FIXME implement, document, tests
class Subgraph
{
public:
  // Build a subgraph from the given query.
  // If the query is based on a path, a path index must be provided.
  // Throws `std::runtime_error` on error.
  Subgraph(const GBZ& gbz, const PathIndex* path_index, const SubgraphQuery& query);

  Subgraph(const Subgraph& source) = default;
  Subgraph(Subgraph&& source) = default;
  Subgraph& operator=(const Subgraph& source) = default;
  Subgraph& operator=(Subgraph&& source) = default;

  // Node identifiers in the subgraph.
  std::set<nid_t> nodes;

  // Paths in the subgraph (in canonical orientation).
  std::vector<gbwt::vector_type> paths;

  // Length of each path (in bp).
  std::vector<size_t> path_lengths;

  // Weight (the number of duplicates) of each path.
  std::vector<size_t> path_weights;

  // CIGAR string for each path relative to the reference path.
  std::vector<std::string> path_cigars;

  // Offset of the reference path in `paths`.
  size_t reference_path;

  // Handle for the reference path.
  path_handle_t reference_handle;

  // Starting offset of the reference path fragment in this subgraph.
  // This is relative to the full path, ignoring any offset that may be stored
  // in path metadata.
  size_t reference_start;

  // Convert the subgraph to GFA.
  void to_gfa(const GBZ& gbz, std::ostream& out) const;

private:
  // Extract the paths within the subgraph and determine reference path information.
  void extract_paths(const GBZ& gbz, size_t query_offset, const std::pair<pos_t, gbwt::edge_type>& ref_pos);

  // Update the paths according to query type.
  void update_paths(const SubgraphQuery& query);

  gbwt::FullPathName reference_path_name(const GBZ& gbz) const;

  const size_t* weight(size_t path_id) const;
  const std::string* cigar(size_t path_id) const;
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_SUBGRAPH_H
