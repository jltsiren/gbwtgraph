#ifndef GBWTGRAPH_INTERNAL_H
#define GBWTGRAPH_INTERNAL_H

#include <gbwtgraph/utils.h>

#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/*
  internal.h: Internal support structures.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  A buffered TSV file writer.
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
  A naive HandleGraph implementation that uses sequences from a memory-mapped
  GFA file.
*/
class GFAGraph : public HandleGraph
{
public:
  struct Node
  {
    view_type sequence;
    std::vector<handle_t> predecessors, successors;
  };

  std::unordered_map<nid_t, Node> nodes;
  nid_t min_id, max_id;

  GFAGraph() : min_id(std::numeric_limits<nid_t>::max()), max_id(0) {}
  virtual ~GFAGraph() {}

  // Insert a new node with the given sequence.
  void insert_node(nid_t node_id, view_type sequence);

  // Insert a new edge and return `true` if the insertion was successful.
  // Returns `false` if the endpoints do not exist.
  bool insert_edge(const handle_t& from, const handle_t& to);

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
  
  // Get node sequence as a pointer and length.
  view_type get_sequence_view(const handle_t& handle) const;
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_INTERNAL_H
