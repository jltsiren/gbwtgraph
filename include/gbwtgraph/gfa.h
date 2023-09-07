#ifndef GBWTGRAPH_GFA_H
#define GBWTGRAPH_GFA_H

#include <memory>
#include <list>

#include <gbwt/dynamic_gbwt.h>

#include <gbwtgraph/gbwtgraph.h>

/*
  gfa.h: Tools for building GBWTGraph from GFA.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// TODO: Add sanity checks.
struct GFAParsingParameters
{
  // GBWT construction parameters. `node_width` is no longer used, as it can be determined
  // automatically. There are no sanity checks for `batch_size`.
  gbwt::size_type node_width = gbwt::WORD_BITS;
  gbwt::size_type batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;
  gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;

  // Chop segments longer than this into multiple nodes. Use 0 to disable chopping.
  size_t max_node_length = MAX_NODE_LENGTH;

  // To avoid creating too many jobs, combine small consecutive components into jobs
  // of at most `num_nodes / approximate_num_jobs` nodes. Value 0 is interpreted as 1.
  constexpr static size_t APPROXIMATE_NUM_JOBS = 32;
  size_t approximate_num_jobs = APPROXIMATE_NUM_JOBS;

  // Try to run this may construction jobs in parallel. Value 0 is interpreted as 1.
  size_t parallel_jobs = 1;

  // Determine GBWT batch size automatically. If the length of the longest path is `N`
  // segments, batch size will be the maximum of the default (100 million) and
  // `gbwt::DynamicGBWT::MIN_SEQUENCES_PER_BATCH * (N + 1)` but no more than GFA file
  // size in bytes. This should ensure that each batch consists of at least 10 paths
  // and their reverse complements. With heavy chopping, path length in nodes may be
  // much larger than `N`, and hence it may be useful to set the batch size manually.
  bool automatic_batch_size = true;

  bool show_progress = false;

  /*
    To parse path names, we use a string regex and a string listing field types.

    For the regex, each submatch (part of the regex separated by parentheses)
    is a field. The fields are numbered according to preorder traversal from
    left to right, with 0 corresponding to the entire path name.

    fields[i] maps field i to a GBWT path name component. Possible values are:

      S  sample name
      C  contig name
      H  haplotype identifier
      F  fragment identifier

    The values are case-insensitive. Any other character indicates that the field should
    not be used. If the string is too short, subsequent fields are not used. Each
    component may occur only once in the string.
  */
  const static std::string DEFAULT_REGEX;  // .*
  const static std::string DEFAULT_FIELDS; // C
  const static PathSense DEFAULT_SENSE;    // Defaults to generic
  const static std::string PAN_SN_REGEX;   // (.*)#(.*)#(.*)
  const static std::string PAN_SN_FIELDS;  // XSHC
  const static PathSense PAN_SN_SENSE;     // panSN names should default to haplotype

  struct PathNameParsingParameters
  {
    std::string regex;
    std::string fields;
    PathSense sense;

    PathNameParsingParameters() = default;
    PathNameParsingParameters(const PathNameParsingParameters& other) = default;
    PathNameParsingParameters(PathNameParsingParameters&& other) = default;
    PathNameParsingParameters& operator=(const PathNameParsingParameters& other) = default;
    PathNameParsingParameters& operator=(PathNameParsingParameters&& other) = default;

    PathNameParsingParameters(const std::string& path_name_regex, const std::string& path_name_fields, PathSense sense = PathSense::GENERIC);
  };

  // We consult these in order until one of the regular expressions matches.
  std::list<PathNameParsingParameters> path_name_formats {{DEFAULT_REGEX, DEFAULT_FIELDS, DEFAULT_SENSE}};

};

//------------------------------------------------------------------------------

struct GFAExtractionParameters
{
  // Use this many OpenMP threads for extracting paths and walks. Value 0 is interpreted
  // as 1.
  size_t num_threads = 1;
  size_t threads() const { return std::max(this->num_threads, size_t(1)); }

  // Cache GBWT records larger than this size to speed up decompression.
  constexpr static size_t LARGE_RECORD_BYTES = 1024;
  size_t large_record_bytes = LARGE_RECORD_BYTES;

  enum path_mode
  {
    mode_default,   // Named paths as P-lines, haplotype paths as W-lines.
    mode_pan_sn,    // All paths as P-lines with PanSN names.
    mode_ref_only,  // Named paths as P-lines.
  };
  path_mode mode = mode_default;

  static std::string mode_name(path_mode mode);
  static path_mode get_mode(const std::string& name);

  // Use the node-to-segment translation if it exists.
  bool use_translation = true;

  bool show_progress = false;
};

//------------------------------------------------------------------------------

/*
  Build GBWT from GFA P-lines and/or W-lines with the following assumptions:

    1. Links and paths have no overlaps between segments.
    2. There are no containments.

  If the construction fails, the function throws `std::runtime_error`. However,
  if there is an error during the multithreaded path/walk parsing stage,
  the construction exits with `std::exit()`.

  Before GBWT construction, the graph is partitioned into weakly connected
  components. The components are ordered by node ids, and contiguous ranges of
  components are assigned to jobs of roughly equal size. A separate GBWT index
  is built for each job, and the partial indexes are merged using the fast
  algorithm. Multiple jobs can be run in parallel.

  The construction is done in several passes over a memory-mapped GFA file. The
  function returns the GBWT index and a sequence source for GBWTGraph construction.

  If the GFA file contains both P-lines and W-lines, both will be used. In that
  case, P-lines will be interpreted as generic-sense paths and stored under a
  sample named only `REFERENCE_PATH_SAMPLE_NAME`, with the path name as contig
  name. If there are only P-lines, GBWT metadata will be parsed using the
  regular expression defined in the parameters.

  If there are segments longer than the maximum length specified in the parameters,
  such segments will be broken into nodes of that length. If segment identifiers are
  not positive integers, they will be translated into such identifiers. In both
  cases, the sequence source will contain a translation from segment names to
  ranges of node identifiers.
*/
std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<SequenceSource>>
gfa_to_gbwt(const std::string& gfa_filename, const GFAParsingParameters& parameters = GFAParsingParameters());

/*
  Writes the graph as GFA into the output stream in a normalized form. The lines are
  ordered in the following way:

  1. S-lines ordered by node ids.

  2. L-lines in canonical order. When using a single threads, the edges (from, to)
  are ordered by tuples (id(from), is_reverse(from), id(to), is_reverse(to)).
  All overlaps are `*`.

  3. P-lines for generic paths stored under the sample named only
  `REFERENCE_PATH_SAMPLE_NAME`. When using a single thread, the paths are ordered
  by GBWT path ids. All overlaps are `*`.

  4. W-lines for other paths. When using a single thread, the paths are ordered by
  GBWT path ids.

  If the GBWT does not contain path names, all GBWT paths will be written as P-lines.
*/
void gbwt_to_gfa(const GBWTGraph& graph, std::ostream& out, const GFAExtractionParameters& parameters = GFAExtractionParameters());

extern const std::string GFA_EXTENSION; // ".gfa"

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_GFA_H
