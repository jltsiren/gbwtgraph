#ifndef GBWTGRAPH_NAIVE_GRAPH_H
#define GBWTGRAPH_NAIVE_GRAPH_H

#include <gbwtgraph/utils.h>

#include <functional>

/*
  naive_graph.h: A naive HandleGraph implementation used as a helper in
  GBWTGraph construction (e.g. from GFA).
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  A naive HandleGraph implementation based on a hash map of node records and
  concatenated sequences.

  Nodes can be added with `create_node(id, sequence)` or translated from GFA
  segments with `translate_segment(name, sequence, max_bp)`. These two
  approaches must not be mixed. In either case, GFA segment names can be
  translated into half-open ranges of node ids with `translate(name)`.

  The tags object is intended for storing GraphName information for the parent graph.
  If a translation is used, this will be for the translation target of the GBWTGraph
  being built. Otherwise it is either for the supergraph or the same graph,
  depending on whether the paths in the GBWT use all nodes and edges in the graph.
*/
class NaiveGraph : public HandleGraph
{
  /*
    Construction.
  */

public:
  NaiveGraph();
  NaiveGraph(const NaiveGraph& another) = default;
  NaiveGraph(NaiveGraph&& another) = default;
  NaiveGraph& operator=(const NaiveGraph& another) = default;
  NaiveGraph& operator=(NaiveGraph&& another) = default;
  virtual ~NaiveGraph() {}

  // Creates a new node with the given sequence.
  // Cannot be used with `translate_segment()`.
  // Throws `std::runtime_error` if the node already exists or if the sequence is empty.
  void create_node(nid_t node_id, std::string_view sequence);

  // Takes a GFA segment with the given name and sequence.
  // Translates it into one or more nodes of at most `max_length` bp each, assigning them the next unused node ids.
  // Returns the half-open node id range for the translated segment.
  // Cannot be used with `create_node()`.
  // Throws `std::runtime_error` if the segment has already been translated or if the sequence is empty.
  std::pair<nid_t, nid_t> translate_segment(const std::string& name, std::string_view sequence, size_t max_length);

  // Creates a new edge between the two oriented nodes.
  // Throws `std::runtime_error` if the nodes do not exist.
  void create_edge(gbwt::node_type from, gbwt::node_type to);

  // Removes all duplicate edges.
  void remove_duplicate_edges();

//------------------------------------------------------------------------------

  /*
    HandleGraph interface. Methods return default values when the node does not exist.
  */

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

//------------------------------------------------------------------------------

  /*
    Data members.
  */

private:

  struct Node
  {
    std::vector<gbwt::node_type> predecessors, successors;
    size_t sequence_offset; // Offset into NaiveGraph::sequences
    size_t sequence_length;
  };

  std::unordered_map<nid_t, Node> nodes;
  std::vector<char> sequences;
  nid_t min_id, max_id;

  // For each segment name, a half-open range of node ids.
  std::unordered_map<std::string, std::pair<nid_t, nid_t>> segment_translation;

  // Tags for storing GraphName information.
  gbwt::Tags tags;

//------------------------------------------------------------------------------

  /*
    Custom interface.
  */

public:

  // Extension for a file storing the node-to-segment translation.
  const static std::string TRANSLATION_EXTENSION; // ".trans"

  // Returns `true` if the graph contains a node with the given id.
  // Potentially faster than `has_node()`, as this is not virtual.
  bool has_node_impl(nid_t node_id) const
  {
    return (this->nodes.find(node_id) != this->nodes.end());
  }

  // Returns an iterator to the specified node by node id.
  std::unordered_map<nid_t, Node>::const_iterator get_node(nid_t node_id) const
  {
    return this->nodes.find(node_id);
  }

  // Returns an iterator to the specified node by handle.
  std::unordered_map<nid_t, Node>::const_iterator get_node(const handle_t& handle) const
  {
    return this->nodes.find(gbwt::Node::id(handle_to_node(handle)));
  }

  // Returns an iterator to the specified node by node id.
  std::unordered_map<nid_t, Node>::iterator get_node_mut(nid_t node_id)
  {
    return this->nodes.find(node_id);
  }

  // Returns an iterator to the specified node by handle.
  std::unordered_map<nid_t, Node>::iterator get_node_mut(const handle_t& handle)
  {
    return this->nodes.find(gbwt::Node::id(handle_to_node(handle)));
  }

  // Returns a view of the sequence of the specified node in the forward orientation.
  // Returns an empty view if the node does not exist.
  std::string_view get_sequence_view(nid_t node_id) const;

  // Convert gbwt::node_type to handle_t.
  static handle_t node_to_handle(gbwt::node_type node) { return handlegraph::as_handle(node); }

  // Convert handle_t to gbwt::node_type.
  static gbwt::node_type handle_to_node(const handle_t& handle) { return handlegraph::as_integer(handle); }

  // Returns `true` if the graph uses node-to-segment translation.
  // That means that nodes were created using `translate_segment()` rather than `add_node()`.
  bool uses_translation() const { return !(this->segment_translation.empty()); }

  // Returns the number of segments in the translation.
  size_t get_segment_count() const { return this->segment_translation.size(); }

  // Returns `true` if the graph contains a translation for a segment with the given name.
  bool has_segment(const std::string& name) const
  {
    return (this->segment_translation.find(name) != this->segment_translation.end());
  }

  // Returns `true` if the graph contains a node/segment with the given name,
  // depending on whether translation is used.
  bool has_node_or_segment(const std::string& name) const;

  // An empty node id range indicating that the translation failed.
  constexpr static std::pair<nid_t, nid_t> no_translation()
  {
    return std::pair<nid_t, nid_t>(0, 0);
  }

  // Returns an iterator to the translations in an arbitrary order.
  std::unordered_map<std::string, std::pair<nid_t, nid_t>>::const_iterator translation_begin() const
  {
    return this->segment_translation.cbegin();
  }

  // Returns an iterator to the end of the translations.
  std::unordered_map<std::string, std::pair<nid_t, nid_t>>::const_iterator translation_end() const
  {
    return this->segment_translation.cend();
  }

  // Translates the segment name into a half-open range of node identifiers.
  // Returns (id, id + 1) if the graph does not use translation and the segment name
  // is a valid integer node id.
  // Returns `no_translation()` if the translation fails.
  std::pair<nid_t, nid_t> translate(const std::string& segment_name) const;

  // Returns `StringArray` of segment names and `sd_vector<>` mapping node ids to names.
  // If `is_present` returns false, the corresponding segment name will be empty.
  // Uses multiple OpenMP threads.
  std::pair<gbwt::StringArray, sdsl::sd_vector<>> invert_translation(const std::function<bool(std::pair<nid_t, nid_t>)>& is_present) const;

  // Stores the given GraphName information in the tags.
  void set_graph_name(const GraphName& graph_name) { graph_name.set_tags(this->tags); }

  // Retrieves the GraphName information from the tags.
  GraphName graph_name() const { return GraphName(this->tags); }
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph


#endif // GBWTGRAPH_NAIVE_GRAPH_H