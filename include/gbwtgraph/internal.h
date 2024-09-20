#ifndef GBWTGRAPH_INTERNAL_H
#define GBWTGRAPH_INTERNAL_H

#include <gbwtgraph/gbz.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

/*
  internal.h: Internal support structures.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  A buffered TSV writer for single-threaded situations.
*/
struct TSVWriter
{
  explicit TSVWriter(std::ostream& out);
  ~TSVWriter();

  void put(char c)
  {
    this->buffer.push_back(c);
    if(this->buffer.size() >= BUFFER_SIZE) { this->flush(); }
  }
  void newline() { this->put('\n'); }
  void newfield() { this->put('\t'); }

  void write(view_type view);
  void write(const std::string& str) { this->write(view_type(str.data(), str.length())); }
  void write(size_t value)
  {
    std::string str = std::to_string(value);
    this->write(str);
  }

  void flush();

  // Buffer this many bytes;
  constexpr static size_t BUFFER_SIZE = 4 * 1048576;

  std::vector<char> buffer;
  std::ostream&     out;
};

//------------------------------------------------------------------------------

/*
  A buffered TSV writer where the buffer grows as necessary and must be flushed
  explicitly. Intended for multi-threaded situations.
*/
struct ManualTSVWriter
{
  explicit ManualTSVWriter(std::ostream& out);

  void put(char c) { this->buffer.push_back(c); }
  void newline() { this->put('\n'); }
  void newfield() { this->put('\t'); }

  void write(view_type view) { this->buffer.insert(this->buffer.end(), view.first, view.first + view.second); }
  void write(const std::string& str) { this->write(view_type(str.data(), str.length())); }
  void write(size_t value)
  {
    std::string str = std::to_string(value);
    this->write(str);
  }

  bool full() const { return this->buffer.size() >= BUFFER_FULL; }
  void flush();

  // Buffer this many bytes;
  constexpr static size_t BUFFER_SIZE = 4 * 1048576;

  // The buffer is (almost) full at this size.
  constexpr static size_t BUFFER_FULL = BUFFER_SIZE - 4096;

  std::vector<char> buffer;
  std::ostream&     out;
};

//------------------------------------------------------------------------------

/*
  In the following, the writer can be either TSVWriter or ManualTSVWriter. We
  only deal with nodes and edges (rather than segments and links). We further
  assume that the graph is small enough that caching segment names is not
  necessary.
*/

// Write a GFA header to the given writer.
// Set `reference_samples = nullptr` if there are no reference samples.
template<class Writer>
void
write_gfa_header(Writer& writer, const std::string* reference_samples)
{
  writer.put('H'); writer.newfield();
  writer.write(std::string("VN:Z:1.1"));
  if(reference_samples != nullptr)
  {
    writer.newfield();
    writer.write(REFERENCE_SAMPLE_LIST_GFA_TAG);
    writer.put(':');
    writer.put('Z');
    writer.put(':');
    writer.write(*reference_samples);
  }
  writer.newline();
}

// Write a node as a GFA segment to the given writer.
template<class Writer>
void
write_gfa_segment(Writer& writer, nid_t node_id, view_type sequence)
{
  writer.put('S'); writer.newfield();
  writer.write(node_id); writer.newfield();
  writer.write(sequence); writer.newline();
}

// Write an edge as a GFA link to the given writer.
template<class Writer>
void
write_gfa_link(Writer& writer, gbwt::node_type from, gbwt::node_type to)
{
  writer.put('L'); writer.newfield();
  writer.write(gbwt::Node::id(from)); writer.newfield();
  writer.put(gbwt::Node::is_reverse(from) ? '-' : '+'); writer.newfield();
  writer.write(gbwt::Node::id(to)); writer.newfield();
  writer.put(gbwt::Node::is_reverse(to) ? '-' : '+'); writer.newfield();
  writer.write("0M"); writer.newline();
}

// Write a path as a GFA walk to the given writer.
// Set `weight = nullptr` and/or `cigar = nullptr` if the information is not available.
template<class Writer>
void
write_gfa_walk(
  Writer& writer,
  const gbwt::vector_type& path, const gbwt::FullPathName& path_name,
  size_t length, const size_t* weight, const std::string* cigar
)
{
  writer.put('W'); writer.newfield();
  writer.write(path_name.sample_name); writer.newfield();
  // NOTE: It was a mistake to define NO_PHASE as unsigned -1.
  // GFA uses 0, which is much more convenient.
  size_t haplotype = (path_name.haplotype == GBWTGraph::NO_PHASE ? 0 : path_name.haplotype);
  writer.write(haplotype); writer.newfield();
  writer.write(path_name.contig_name); writer.newfield();
  writer.write(path_name.offset); writer.newfield();
  writer.write(path_name.offset + length); writer.newfield();

  for(gbwt::node_type node : path)
  {
    writer.put((gbwt::Node::is_reverse(node) ? '<' : '>'));
    writer.write(gbwt::Node::id(node));
  }

  if(weight != nullptr)
  {
    writer.newfield();
    writer.write(std::string("WT:i:"));
    writer.write(*weight);
  }
  if(cigar != nullptr)
  {
    writer.newfield();
    writer.write(std::string("CG:Z:"));
    writer.write(*cigar);
  }

  writer.newline();
}

//------------------------------------------------------------------------------

/*
  A HandleGraph implementation that stores graph topology without sequences.
*/
class EmptyGraph : public HandleGraph
{
public:
  struct Node
  {
    std::vector<handle_t> predecessors, successors;
  };

  std::unordered_map<nid_t, Node> nodes;
  nid_t min_id, max_id;

  EmptyGraph() : min_id(std::numeric_limits<nid_t>::max()), max_id(0) {}
  virtual ~EmptyGraph() {}

  // Create a new node.
  void create_node(nid_t node_id);

  // Create a new edge and return `true` if the insertion was successful.
  // Throws `std::runtime_error` if the nodes do not exist.
  void create_edge(const handle_t& from, const handle_t& to);

  // Remove all duplicate edges.
  void remove_duplicate_edges();

public:

  // Method to check if a node exists by ID.
  virtual bool has_node(nid_t node_id) const;

  // Look up the handle for the node with the given ID in the given orientation.
  virtual handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;

  // Get the ID from a handle.
  virtual nid_t get_id(const handle_t& handle) const;

  // Get the orientation of a handle.
  virtual bool get_is_reverse(const handle_t& handle) const;

  // Invert the orientation of a handle (potentially without getting its ID).
  virtual handle_t flip(const handle_t& handle) const;

  // Get the length of a node.
  virtual size_t get_length(const handle_t& handle) const;

  // Get the sequence of a node, presented in the handle's local forward
  // orientation.
  virtual std::string get_sequence(const handle_t& handle) const;

  // Returns one base of a handle's sequence, in the orientation of the
  // handle.
  virtual char get_base(const handle_t& handle, size_t index) const;

  // Returns a substring of a handle's sequence, in the orientation of the
  // handle. If the indicated substring would extend beyond the end of the
  // handle's sequence, the return value is truncated to the sequence's end.
  virtual std::string get_subsequence(const handle_t& handle, size_t index, size_t size) const;

  // Return the number of nodes in the graph.
  virtual size_t get_node_count() const;

  // Return the smallest ID in the graph, or some smaller number if the
  // smallest ID is unavailable. Return value is unspecified if the graph is empty.
  virtual nid_t min_node_id() const;

  // Return the largest ID in the graph, or some larger number if the
  // largest ID is unavailable. Return value is unspecified if the graph is empty.
  virtual nid_t max_node_id() const;

protected:

  // Loop over all the handles to next/previous (right/left) nodes. Passes
  // them to a callback which returns false to stop iterating and true to
  // continue. Returns true if we finished and false if we stopped early.
  virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;

  // Loop over all the nodes in the graph in their local forward
  // orientations, in their internal stored order. Stop if the iteratee
  // returns false. Can be told to run in parallel, in which case stopping
  // after a false return value is on a best-effort basis and iteration
  // order is not defined. Returns true if we finished and false if we
  // stopped early.
  virtual bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;

public:

  // Get the number of edges on the right (go_left = false) or left (go_left
  // = true) side of the given handle.
  virtual size_t get_degree(const handle_t& handle, bool go_left) const;

public:

  // Returns an iterator to the specified node.
  std::unordered_map<nid_t, Node>::const_iterator get_node(const handle_t& handle) const
  {
    return this->nodes.find(gbwt::Node::id(handle_to_node(handle)));
  }

  // Returns an iterator to the specified node.
  std::unordered_map<nid_t, Node>::iterator get_node_mut(const handle_t& handle)
  {
    return this->nodes.find(gbwt::Node::id(handle_to_node(handle)));
  }

  // Convert gbwt::node_type to handle_t.
  static handle_t node_to_handle(gbwt::node_type node) { return handlegraph::as_handle(node); }

  // Convert handle_t to gbwt::node_type.
  static gbwt::node_type handle_to_node(const handle_t& handle) { return handlegraph::as_integer(handle); }
};

//------------------------------------------------------------------------------

/*
  A cache that stores GBWT records larger than `bytes` bytes as `DecompressedRecord`
  and supports faster path extraction from the index.
*/
struct LargeRecordCache
{
  LargeRecordCache(const gbwt::GBWT& index, size_t bytes);

  size_t size() const { return this->cache.size(); }
  gbwt::size_type sequences() const { return this->index.sequences(); }

  gbwt::vector_type extract(gbwt::size_type sequence) const;

  const gbwt::GBWT& index;
  std::unordered_map<gbwt::node_type, gbwt::DecompressedRecord> cache;
};

//------------------------------------------------------------------------------

// Sample (sequence offset, GBWT position) at the start of a node approximately
// every `sample_interval` bp, with the first sample at offset 0.
// If `length` is not nullptr, it will be set to the length of the path.
// Sequence offsets are relative to the path, not the full haplotype.
std::vector<std::pair<size_t, gbwt::edge_type>> sample_path_positions(const GBZ& gbz, path_handle_t path, size_t sample_interval, size_t* length = nullptr);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_INTERNAL_H
