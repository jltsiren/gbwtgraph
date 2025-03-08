#ifndef GBWTGRAPH_UTILS_H
#define GBWTGRAPH_UTILS_H

#include <gbwt/gbwt.h>

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/serializable_handle_graph.hpp>
#include <handlegraph/named_node_back_translation.hpp>
#include <handlegraph/util.hpp>

#include <sdsl/sd_vector.hpp>

#include <functional>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

/*
  utils.h: Common utilities and definitions.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// We should be using std::hardware_destructive_interference_size, but Apple
// Clang does not support it and GCC really hates it. So we just use 128,
// which is a bit wasteful on Intel.
inline constexpr size_t CACHE_LINE_SIZE = 128;

//------------------------------------------------------------------------------

// Import type definitions from libhandlegraph.

using nid_t = handlegraph::nid_t;
using offset_t = handlegraph::offset_t;
using subrange_t = handlegraph::subrange_t;
using pos_t = handlegraph::pos_t;
using handle_t = handlegraph::handle_t;
using path_handle_t = handlegraph::path_handle_t;
using PathSense = handlegraph::PathSense;
using step_handle_t = handlegraph::step_handle_t;
using edge_t = handlegraph::edge_t;
using oriented_node_range_t = handlegraph::oriented_node_range_t;

using HandleGraph = handlegraph::HandleGraph;
using PathHandleGraph = handlegraph::PathHandleGraph;
using PathMetadata = handlegraph::PathMetadata;
using SerializableHandleGraph = handlegraph::SerializableHandleGraph;
using NamedNodeBackTranslation = handlegraph::NamedNodeBackTranslation;

//------------------------------------------------------------------------------

// In-place view of the sequence: (start, length).
// This is a quick replacement for std::string_view from C++17.
typedef std::pair<const char*, size_t> view_type;

inline view_type
str_to_view(const std::string& str)
{
  return view_type(str.data(), str.length());
}

// This is the equivalent of `view_type` or `std::string_view` for pats
// represented as `gbwt::vector_type`.
typedef std::pair<const gbwt::vector_type::value_type*, size_t> subpath_type;

// Returns a subpath corresponding the half-open interval [from, to) in the given path.
inline subpath_type
get_subpath(const gbwt::vector_type& path, size_t from, size_t to)
{
  return subpath_type(path.data() + from, to - from);
}

//------------------------------------------------------------------------------

// Custom exception that tells that something is wrong with the GBWT index.
class InvalidGBWT : public std::runtime_error
{
public:
  explicit InvalidGBWT(const std::string& message) : std::runtime_error(message) {}
  explicit InvalidGBWT(const char* message) : std::runtime_error(message) {}
};

//------------------------------------------------------------------------------

// Some tools may not work with nodes longer than this.
constexpr size_t MAX_NODE_LENGTH = 1024;

// Paths with this sample name (`_gbwt_ref`) are interpreted as generic named
// paths. All generic paths must have distinct contig names that can be used as
// path names.
extern const std::string REFERENCE_PATH_SAMPLE_NAME;

// Which paths are reference vs. haplotype paths is determined by storing the
// sample names in this GBWT tag
extern const std::string REFERENCE_SAMPLE_LIST_GBWT_TAG;

// Which paths are reference vs. haplotype paths is determined by storing the
// sample names in this GFA tag
extern const std::string REFERENCE_SAMPLE_LIST_GFA_TAG;

// What separator character is used in lists of reference sample names?
extern const char REFERENCE_SAMPLE_LIST_SEPARATOR;

// Cached information for a named path.
struct NamedPath
{
  gbwt::size_type id; // Original path id.
  gbwt::edge_type from, to; // First / last position on the path, or `gbwt::invalid_edge()` if empty.
  size_t length;
  PathSense sense;
};

//------------------------------------------------------------------------------

struct Version
{
  static std::string str(bool verbose = false);
  static void print(std::ostream& out, const std::string& tool_name, bool verbose = false, size_t new_lines = 2);

  constexpr static size_t MAJOR_VERSION     = 1;
  constexpr static size_t MINOR_VERSION     = 1;
  constexpr static size_t PATCH_VERSION     = 0;

  constexpr static size_t GBZ_VERSION       = 1;
  constexpr static size_t GRAPH_VERSION     = 3;
  constexpr static size_t MINIMIZER_VERSION = 10;

  const static std::string SOURCE_KEY; // source
  const static std::string SOURCE_VALUE; // jltsiren/gbwtgraph
};

//------------------------------------------------------------------------------

// Because we want to be able to work with path metadata with just the GBWT, we
// expose the utility functions for dealing with it. The packing format only
// has to touch these functions and MetadataBuilder, and some of the GBWTGraph
// search methods.

// Parse a path sense tag value to a collection of reference-sense sample names.
std::unordered_set<std::string> parse_reference_samples_tag(const char* cursor, const char* end);
// Parse a path sense tag value to a collection of reference-sense sample names.
std::unordered_set<std::string> parse_reference_samples_tag(const std::string& tag_value);
// Parse a path sense tag value to a collection of reference-sense sample names.
std::unordered_set<std::string> parse_reference_samples_tag(const view_type& tag_value);
// Parse the reference samples tag embedded in a GBWT index.
std::unordered_set<std::string> parse_reference_samples_tag(const gbwt::GBWT& index);

// Compose a reference sample name tag from a collection of reference-sense sample names.
std::string compose_reference_samples_tag(const std::unordered_set<std::string>& reference_samples);

// Determine the sense that a path ought to have, from stored metadata.
// Takes a reference path set from parse_reference_samples_tag, and handles identifying
// the string sample name and defaulting un-mentioned samples.
// Tolerates incomplete metadata.
PathSense get_path_sense(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, const std::unordered_set<std::string>& reference_samples);

// Determine the sense a path ought to have, given an index and the reference
// samples set parsed from it with parse_reference_samples_tag.
// Tolerates missing metadata.
PathSense get_path_sense(const gbwt::GBWT& index, gbwt::size_type path_number, const std::unordered_set<std::string>& reference_samples);

// Determine the sense that paths for a sample ought to have, given the sample number.
// Tolerates incomplete metadata.
PathSense get_sample_sense(const gbwt::Metadata& metadata, gbwt::size_type sample, const std::unordered_set<std::string>& reference_samples);

// Determine the sense that paths for a sample ought to have, given the string sample name.
// The metadata is assumed to be populated with sample names.
PathSense get_sample_sense(const std::string& sample_name, const std::unordered_set<std::string>& reference_samples);

// Determine the sample name that a path ought to present, from stored metadata.
// Tolerates incomplete metadata.
std::string get_path_sample_name(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense);

// Determine the sample name that a path ought to present, from stored metadata.
// Tolerates missing metadata.
std::string get_path_sample_name(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense);

// Determine the locus name that a path ought to present, from stored metadata.
// Tolerates incomplete metadata.
std::string get_path_locus_name(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense);

// Determine the locus name that a path ought to present, from stored metadata.
// Tolerates missing metadata.
std::string get_path_locus_name(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense);

// Determine the haplotype phase number that a path ought to present, from stored metadata.
// Tolerates incomplete metadata.
size_t get_path_haplotype(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense);

// Determine the haplotype phase number that a path ought to present, from stored metadata.
// Tolerates missing metadata.
size_t get_path_haplotype(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense);

// Determine the phase block number a path ought to present, from stored metadata.
// Tolerates incomplete metadata.
size_t get_path_phase_block(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense);

// Determine the phase block number a path ought to present, from stored metadata.
// Tolerates missing metadata.
size_t get_path_phase_block(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense);

// Determine the subrange that a path ought to present, from stored metadata.
// Tolerates incomplete metadata.
subrange_t get_path_subrange(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense);

// Determine the subrange that a path ought to present, from stored metadata.
// Tolerates missing metadata.
subrange_t get_path_subrange(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense);

// Get a libhandlegraph path name string from the given stored path metadata.
// Tolerates incomplete metadata.
std::string compose_path_name(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense);

// Get a libhandlegraph path name string from the given stored path metadata.
// Tolerates missing metadata.
std::string compose_path_name(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense);

// We also have write accessors.
// For write access to sample name, locus, haplotype, phase block, and
// subrange, see MetadataBuilder::add_path().

// Set the senses for samples' paths, into the given GBWT tag set.
void set_sample_path_senses(gbwt::Tags& tags, const std::unordered_map<std::string, PathSense>& senses);

//------------------------------------------------------------------------------

// Functions for pos_t manipulation.

inline pos_t
make_pos_t(nid_t id, bool is_rev, offset_t off)
{
  return std::make_tuple(id, is_rev, off);
}

inline nid_t
id(const pos_t& pos)
{
  return std::get<0>(pos);
}

inline bool
is_rev(const pos_t& pos)
{
  return std::get<1>(pos);
}

inline offset_t
offset(const pos_t& pos)
{
  return std::get<2>(pos);
}

inline nid_t&
get_id(pos_t& pos)
{
  return std::get<0>(pos);
}

inline bool&
get_is_rev(pos_t& pos)
{
  return std::get<1>(pos);
}

inline offset_t&
get_offset(pos_t& pos)
{
  return std::get<2>(pos);
}

inline bool
is_empty(const pos_t& pos)
{
  return (id(pos) == 0);
}

// Get a pos_t for the same base in the other orientation.
inline pos_t
reverse_base_pos(const pos_t& pos, size_t node_length)
{
  return make_pos_t(id(pos), !is_rev(pos), (node_length - 1) - offset(pos));
}

inline std::ostream&
operator<<(std::ostream& out, const pos_t& pos)
{
  return out << id(pos) << (is_rev(pos) ? "-" : "+") << offset(pos);
}

//------------------------------------------------------------------------------

/*
  Thomas Wang's integer hash function. In many implementations, std::hash
  is identity function for integers, which leads to performance issues.
*/

inline size_t
wang_hash_64(size_t key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

// Essentially boost::hash_combine.
inline size_t
hash(nid_t id, bool is_rev, size_t offset)
{
  size_t result = wang_hash_64(id);
  result ^= wang_hash_64(is_rev) + 0x9e3779b9 + (result << 6) + (result >> 2);
  result ^= wang_hash_64(offset) + 0x9e3779b9 + (result << 6) + (result >> 2);
  return result;
}

inline size_t
hash(const pos_t& pos)
{
  return hash(id(pos), is_rev(pos), offset(pos));
}

//------------------------------------------------------------------------------

// Utility functions.

std::string reverse_complement(const std::string& seq);
void reverse_complement_in_place(std::string& seq);

/*
  An edge is in canonical orientation, if:
  1. The destination node has a higher identifier than the source node.
  2. It is a self-loop with at least one end in forward orientation.
*/
bool edge_is_canonical(gbwt::node_type from, gbwt::node_type to);

/*
  A path is in canonical orientation, if:

  1. It is empty.
  2. Both the first and the last node are in forward orientation.
  3. The edge from the first node to the last node would be in canonical orientation.
*/
bool path_is_canonical(const gbwt::vector_type& path);

//------------------------------------------------------------------------------

/*
  An intermediate representation for building GBWTGraph from GFA. This class maps
  node ids to sequences and stores the translation from segment names to (ranges of)
  node ids.

  Nodes can be added with add_node(id, sequence) or translated from GFA segments
  with translate_segment(name, sequence, max_bp). These two approaches must not
  be mixed. In the latter case, the translation can be retrieved with
  get_translation(name).

  Nodes / segments with empty sequence are silently ignored.
*/
class SequenceSource
{
public:
  SequenceSource() : next_id(1) {}

  void swap(SequenceSource& another);

  void add_node(nid_t node_id, const std::string& sequence);
  void add_node(nid_t node_id, view_type sequence);

  bool has_node(nid_t id) const
  {
    return (this->nodes.find(id) != this->nodes.end());
  }

  size_t get_node_count() const { return this->nodes.size(); }

  size_t get_length(nid_t id) const
  {
    auto iter = this->nodes.find(id);
    if(iter == this->nodes.end()) { return 0; }
    return iter->second.second;
  }

  std::string get_sequence(nid_t id) const
  {
    auto iter = this->nodes.find(id);
    if(iter == this->nodes.end()) { return ""; }
    const char* ptr = this->sequences.data() + iter->second.first;
    return std::string(ptr, ptr + iter->second.second);
  }

  view_type get_sequence_view(nid_t id) const
  {
    auto iter = this->nodes.find(id);
    if(iter == this->nodes.end()) { return view_type(nullptr, 0); }
    const char* ptr = this->sequences.data() + iter->second.first;
    return view_type(ptr, iter->second.second);
  }

  // An empty translation or a translation that does not exist.
  constexpr static std::pair<nid_t, nid_t> invalid_translation()
  {
    return std::pair<nid_t, nid_t>(0, 0);
  }

  // Take a GFA segment (name, sequence). If the segment has not been translated
  // yet, break it into nodes of at most max_length bp each and assign them the
  // next unused node ids. Returns the node id range for the translated segment
  // or `invalid_translation()` on failure.
  std::pair<nid_t, nid_t> translate_segment(const std::string& name, view_type sequence, size_t max_length);

  bool uses_translation() const { return !(this->segment_translation.empty()); }

  // Returns a semiopen range of node ids, or `invalid_translation()` if there is
  // no such segment.
  std::pair<nid_t, nid_t> get_translation(const std::string& segment_name) const
  {
    auto iter = this->segment_translation.find(segment_name);
    if(iter == this->segment_translation.end()) { return invalid_translation(); }
    return iter->second;
  }

  // Translates the segment if translation is in use, or converts the segment
  // name into an integer `id` and returns `(id, id + 1)` otherwise.
  // Returns `invalid_translation()` on failure.
  std::pair<nid_t, nid_t> force_translate(const std::string& segment_name) const;

  // Returns `StringArray` of segment names and `sd_vector<>` mapping node ids to names.
  // If `is_present` returns false, the corresponding segment name will be empty.
  std::pair<gbwt::StringArray, sdsl::sd_vector<>> invert_translation(const std::function<bool(std::pair<nid_t, nid_t>)>& is_present) const;

  // (offset, length) for the node sequence.
  std::unordered_map<nid_t, std::pair<size_t, size_t>> nodes;
  std::vector<char> sequences;

  // If segment translation is enabled, this translates a segment identifier
  // into a semiopen range of node identifiers.
  std::unordered_map<std::string, std::pair<nid_t, nid_t>> segment_translation;
  nid_t next_id;

  const static std::string TRANSLATION_EXTENSION; // ".trans"
};

//------------------------------------------------------------------------------

/*
  A structure for building GBWT metadata in a way that is compatible with vg and
  the libhandlegraph path metadata model.

  All paths are assigned to a job, with 0 as the default job. When the GBWT
  metadata is extracted using get_metadata(), the paths are ordered by job.
  This corresponds to building partial GBWTs for each job and merging them in
  the same order.

  Note that this class is not thread-safe. When using multiple GBWT construction
  jobs, the best practice is to build the metadata first.

  Constructor and the methods for handling paths/walks throw `std::runtime_error`
  on failure.
*/
class MetadataBuilder
{
public:
  // Construct a MetadataBuilder with no path name formats.
  MetadataBuilder();

  // Construct a MetadataBuilder from pre-existing metadata and no path name formats.
  // Throws `std::runtime_error` if the metadata does not contain sample, contig,
  // and path names.
  MetadataBuilder(const gbwt::Metadata& source);

  // Construct a MetadataBuilder with one path name format.
  MetadataBuilder(
    const std::string& path_name_regex,
    const std::string& path_name_prefix,
    PathSense path_sense = PathSense::GENERIC);

  // An empty value.
  constexpr static size_t NO_FIELD = std::numeric_limits<size_t>::max();

  // A structure that parses path metadata from a string name.
  struct PathMetadataBuilder
  {
    std::regex parser;

    // Mapping from regex submatches to GBWT path name components.
    size_t sample_field, contig_field, haplotype_field, fragment_field;

    PathSense sense;

    PathMetadataBuilder(const std::string& path_name_regex, const std::string& path_name_prefix, PathSense path_sense);
  };

  // List of possible path name formats ordered by priority.
  std::vector<PathMetadataBuilder> path_name_formats;

  // GBWT metadata.
  std::map<std::string, size_t> sample_names, contig_names;
  std::set<std::pair<size_t, size_t>> haplotypes; // (sample id, phase id)
  std::vector<std::vector<gbwt::PathName>> path_names;
  std::map<gbwt::PathName, size_t> counts;

  // This can be set to true to enable warnings about the use of sample name
  // REFERENCE_PATH_SAMPLE_NAME in inappropriate situations.
  bool ref_path_sample_warning;

  // Register a format for parsing path names. Formats are tried in order until one matches.
  void add_path_name_format(const std::string& path_name_regex, const std::string& path_name_fields, PathSense path_sense);

  // Add a path defined by libhandlegraph metadata to the given job.
  // Doesn't create metadata for samples or contigs if the no-name sentinel is
  // passed for a sense that usually has them.
  void add_path(PathSense sense, const std::string& sample_name, const std::string& locus_name, size_t haplotype, size_t phase_block, const handlegraph::subrange_t& subrange, size_t job = 0);

  // Parse a path name using a regex to determine metadata, and assign it to the given job.
  void add_path(const std::string& name, size_t job = 0);

  // Add a path based on GFA walk metadata and assign it to the given job.
  void add_walk(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start, size_t job = 0);

  // Add a haplotype path and assign it to the given job.
  void add_haplotype(const std::string& sample, const std::string& contig, size_t haplotype, size_t fragment, size_t job = 0);

  // Add a generic path and assign it to the given job.
  void add_generic_path(const std::string& name, size_t job = 0);

  bool empty() const { return this->path_names.empty(); }

  // Build GBWT metadata from the current contents.
  // Paths come out ordered by job.
  gbwt::Metadata get_metadata() const;

  // Clears the existing metadata but leaves the path name formats intact.
  void clear();

private:
  // Add a path name to the given job, increasing the number of jobs if necessary.
  void add_path_name(const gbwt::PathName& path_name, size_t job);

  // Returns a vector of strings `result` such that `result[i] = "value"` if
  // `source["value"] = i`. The strings are ordered by the values. The
  // mapped values are assumed to be distinct and from `0` to `source.size() - 1`.
  static std::vector<std::string> map_to_vector(const std::map<std::string, size_t>& source);
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_UTILS_H
