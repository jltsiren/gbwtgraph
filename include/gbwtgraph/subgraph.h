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

/*
  A data structure for indexing reference and generic paths in a GBZ graph for
  random access by sequence offset. For each path, the index stores GBWT
  positions corresponding to node starts once every `sample_interval` bp
  (default 1024 bp). If the path is a haplotype fragment, the positions are
  given relative to the start of the fragment, not the full haplotype.

  The index relies on the fact that path handles for generic / reference paths
  in a GBWTGraph are integers in the range [0, n), where n is the number of such
  paths in the graph.
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
  size_t path_length(path_handle_t handle) const { return this->sequence_positions[handlegraph::as_integer(handle)].size(); }

  // Returns the last sampled position at or before `offset` in the path,
  // as a pair (sequence offset, GBWT position). If there is no such path,
  // returns (0, gbwt::invalid_edge()).
  std::pair<size_t, gbwt::edge_type> sampled_position(path_handle_t handle, size_t offset) const;
};

//------------------------------------------------------------------------------

/*
  A query for extracting a subgraph from a GBZ graph. The subgraph contains all
  nodes and edges within the given context (in bp) around a path offset, path
  interval, or a node.

  The query can also extract all path fragments within the subgraph. If the
  query is based on a path, the corresponding fragment will have true metadata.
  All other paths are reported as unknown paths. It is possible to output all
  paths, only distinct paths (with the number of copies as the weight), or only
  the reference path.
*/
struct SubgraphQuery
{
  enum QueryType { path_offset_query, path_interval_query, node_query, invalid_query };

  enum HaplotypeOutput { all_haplotypes, distinct_haplotypes, reference_only };

  // Type of the query.
  QueryType type;

  // Which haplotypes to output.
  HaplotypeOutput haplotype_output;

  // A path offset/interval query as (sample, contig, haplotype, offset).
  gbwt::FullPathName path_query;

  // Node id for the node query.
  nid_t node_id;

  // Context length around the query position (in bp).
  size_t context;

  // Creates a subgraph query based on a haplotype offset.
  // The offset stored in `path_name` will be ignored.
  static SubgraphQuery path_offset(const gbwt::FullPathName& path_name, size_t offset, size_t context, HaplotypeOutput output);

  // Creates a subgraph query based on a haplotype interval [`from`, `to`).
  // The offset stored in `path_name` will be ignored.
  // Throws `std::runtime_error` if the interval is empty.
  static SubgraphQuery path_interval(const gbwt::FullPathName& path_name, size_t from, size_t to, size_t context, HaplotypeOutput output);

  // Creates a subgraph query based on a node id.
  static SubgraphQuery node(nid_t node, size_t context, HaplotypeOutput output);

  // Returns a string representation of the query.
  std::string to_string() const;
};

//------------------------------------------------------------------------------

/*
  A subgraph extracted from a GBZ graph using a `SubgraphQuery`. The subgraph
  is defined as the induced subgraph based on a stored set of node identifiers.
  There may be a defined reference path. There may be a weight (the number of
  duplicates) for each path and a CIGAR string relative to the reference for
  each non-reference path.

  At the moment, the only supported functionality is converting the subgraph
  to GFA format.
*/
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

  // Metadata for the reference path.
  gbwt::FullPathName reference_name;

  // Starting offset of the reference path fragment in this subgraph.
  // This is relative to the path fragment corresponding to `reference_handle`, not the full haplotype.
  // The starting offset of the path fragment in the full haplotype is `reference_name.offset`.
  size_t reference_start;

  // Convert the subgraph to GFA.
  void to_gfa(const GBZ& gbz, std::ostream& out) const;

private:
  // Extract the paths within the subgraph and determine reference path information.
  void extract_paths(const GBZ& gbz, const SubgraphQuery& query, size_t query_offset, const std::pair<pos_t, gbwt::edge_type>& ref_pos);

  // Update the paths according to query type.
  void update_paths(const SubgraphQuery& query);

  const size_t* weight(size_t path_id) const;
  const std::string* cigar(size_t path_id) const;
};

//------------------------------------------------------------------------------

/*
  An algorithm for aligning a subpath relative to a reference subpath. The two
  subpaths are assumed to be diverging. The alignment is done under the
  following assumptions:

  * There may be matching bases at the beginning and the end.
  * Any matching bases in the middle are considered accidental and treated as
    mismatches.
  * Scoring parameters are the same as in vg: match = 1, mismatch = -4,
    gap open = -6, gap extend = -1.
  * If multiple operations are used for the middle, the order is mismatch,
    insertion, deletion.

  The alignment is appended to the given vector of edit operation as pairs
  ('M', length), ('I', length), or ('D', length), as in CIGAR strings.
  Successive edits of the same type are merged.
*/
void align_diverging(
  const GBWTGraph& graph,
  subpath_type path, subpath_type reference,
  std::vector<std::pair<char, size_t>>& edits
);

/*
  Converts a sequence of edit operations to a CIGAR string.

  The edits are given as pairs (type, length).
*/
std::string to_cigar(const std::vector<std::pair<char, size_t>>& edits);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_SUBGRAPH_H
