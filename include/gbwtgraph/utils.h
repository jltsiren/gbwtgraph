#ifndef GBWTGRAPH_UTILS_H
#define GBWTGRAPH_UTILS_H

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/serializable_handle_graph.hpp>
#include <handlegraph/util.hpp>

#include <iostream>
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
using edge_t = handlegraph::edge_t;

using HandleGraph = handlegraph::HandleGraph;
using SerializableHandleGraph = handlegraph::SerializableHandleGraph;

//------------------------------------------------------------------------------

// In-place view of the sequence: (start, length).
// This is a quick replacement to std::string_view from C++17.
typedef std::pair<const char*, size_t> view_type;

//------------------------------------------------------------------------------

// Some tools may not work with nodes longer than this.
constexpr size_t MAX_NODE_LENGTH = 1024;

//------------------------------------------------------------------------------

struct Version
{
  static std::string str(bool verbose = false);
  static void print(std::ostream& out, const std::string& tool_name, bool verbose = false, size_t new_lines = 2);

  constexpr static size_t MAJOR_VERSION     = 0;
  constexpr static size_t MINOR_VERSION     = 4;
  constexpr static size_t PATCH_VERSION     = 0;

  constexpr static size_t GRAPH_VERSION     = 1;
  constexpr static size_t MINIMIZER_VERSION = 7;
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
  A class that maps handles to strings. This can be used as a sequence source in
  GBWTGraph construction.

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

  void add_node(nid_t node_id, const std::string& sequence);
  void add_node(nid_t node_id, view_type sequence);

  // Take a GFA segment (name, sequence). If the segment has not been translated
  // yet, break it into nodes of at most max_length bp each and assign them the
  // next unused node ids.
  void translate_segment(const std::string& name, view_type sequence, size_t max_length);

  // Returns a semiopen range of node ids.
  std::pair<nid_t, nid_t> get_translation(const std::string& segment_name) const
  {
    auto iter = this->segment_translation.find(segment_name);
    if(iter == this->segment_translation.end()) { return std::pair<nid_t, nid_t>(0, 0); }
    return iter->second;
  }

  bool has_node(nid_t node_id) const
  {
    return (this->nodes.find(this->get_handle(node_id, false)) != this->nodes.end());
  }

  bool uses_translation() const { return !(this->segment_translation.empty()); }

  size_t get_node_count() const { return this->nodes.size(); }

  // This is compatible with GBWTGraph.
  handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const
  {
    return handlegraph::number_bool_packing::pack(node_id, is_reverse);
  }

  size_t get_length(const handle_t& handle) const
  {
    auto iter = this->nodes.find(handle);
    if(iter == this->nodes.end()) { return 0; }
    return iter->second.second;
  }

  std::string get_sequence(const handle_t& handle) const
  {
    auto iter = this->nodes.find(handle);
    if(iter == this->nodes.end()) { return ""; }
    const char* ptr = this->sequences.data() + iter->second.first;
    return std::string(ptr, ptr + iter->second.second);
  }

  view_type get_sequence_view(const handle_t& handle) const
  {
    auto iter = this->nodes.find(handle);
    if(iter == this->nodes.end()) { return view_type(nullptr, 0); }
    const char* ptr = this->sequences.data() + iter->second.first;
    return view_type(ptr, iter->second.second);
  }

  // Semiopen interval for the node sequence.
  std::unordered_map<handle_t, std::pair<size_t, size_t>> nodes;
  std::vector<char> sequences;

  // If segment translation is enabled, this translates a segment identifier
  // into a semiopen range of node identifiers.
  std::unordered_map<std::string, std::pair<nid_t, nid_t>> segment_translation;
  nid_t next_id;

  const static std::string TRANSLATION_EXTENSION; // ".trans"

private:
  SequenceSource(const SequenceSource&) = delete;
  SequenceSource& operator=(const SequenceSource&) = delete;
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_UTILS_H
