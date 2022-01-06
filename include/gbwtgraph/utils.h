#ifndef GBWTGRAPH_UTILS_H
#define GBWTGRAPH_UTILS_H

#include <gbwt/support.h>

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/serializable_handle_graph.hpp>
#include <handlegraph/named_node_back_translation.hpp>
#include <handlegraph/util.hpp>

#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/simple_sds.hpp>

#include <functional>
#include <iostream>
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

// Import type definitions from libhandlegraph.

using nid_t = handlegraph::nid_t;
using off_t = handlegraph::off_t;
using pos_t = handlegraph::pos_t;
using handle_t = handlegraph::handle_t;
using path_handle_t = handlegraph::path_handle_t;
using step_handle_t = handlegraph::step_handle_t;
using edge_t = handlegraph::edge_t;
using oriented_node_range_t = handlegraph::oriented_node_range_t;

using HandleGraph = handlegraph::HandleGraph;
using PathHandleGraph = handlegraph::PathHandleGraph;
using SerializableHandleGraph = handlegraph::SerializableHandleGraph;
using NamedNodeBackTranslation = handlegraph::NamedNodeBackTranslation;

//------------------------------------------------------------------------------

// In-place view of the sequence: (start, length).
// This is a quick replacement to std::string_view from C++17.
typedef std::pair<const char*, size_t> view_type;

inline view_type
str_to_view(const std::string& str)
{
  return view_type(str.data(), str.length());
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

// Paths with this sample name (`_gbwt_ref`) are interpreted as reference paths.
// All such paths must have distinct contig names that can be used as path
// names.
extern const std::string REFERENCE_PATH_SAMPLE_NAME;

// Cached information for a reference path.
struct ReferencePath
{
  gbwt::size_type id; // Original path id.
  gbwt::edge_type from, to; // First / last position on the path, or `gbwt::invalid_edge()` if empty.
  size_t length;
};

//------------------------------------------------------------------------------

struct Version
{
  static std::string str(bool verbose = false);
  static void print(std::ostream& out, const std::string& tool_name, bool verbose = false, size_t new_lines = 2);

  constexpr static size_t MAJOR_VERSION     = 0;
  constexpr static size_t MINOR_VERSION     = 7;
  constexpr static size_t PATCH_VERSION     = 0;

  constexpr static size_t GBZ_VERSION       = 1;
  constexpr static size_t GRAPH_VERSION     = 3;
  constexpr static size_t MINIMIZER_VERSION = 7;

  const static std::string SOURCE_KEY; // source
  const static std::string SOURCE_VALUE; // jltsiren/gbwtgraph
};

//------------------------------------------------------------------------------

// Functions for pos_t manipulation.

inline pos_t
make_pos_t(nid_t id, bool is_rev, off_t off)
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

inline off_t
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

inline off_t&
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

//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------

  // Take a GFA segment (name, sequence). If the segment has not been translated
  // yet, break it into nodes of at most max_length bp each and assign them the
  // next unused node ids.
  void translate_segment(const std::string& name, view_type sequence, size_t max_length);

  bool uses_translation() const { return !(this->segment_translation.empty()); }

  // Returns a semiopen range of node ids, or (0, 0) if there is no such segment.
  std::pair<nid_t, nid_t> get_translation(const std::string& segment_name) const
  {
    auto iter = this->segment_translation.find(segment_name);
    if(iter == this->segment_translation.end()) { return std::pair<nid_t, nid_t>(0, 0); }
    return iter->second;
  }

  // Returns `StringArray` of segment names and `sd_vector<>` mapping node ids to names.
  // If `is_present` returns false, the corresponding segment name will be empty.
  std::pair<gbwt::StringArray, sdsl::sd_vector<>> invert_translation(const std::function<bool(std::pair<nid_t, nid_t>)>& is_present) const;

//------------------------------------------------------------------------------

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

} // namespace gbwtgraph

#endif // GBWTGRAPH_UTILS_H
