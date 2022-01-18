#ifndef GBWTGRAPH_CACHED_GBWTGRAPH_H
#define GBWTGRAPH_CACHED_GBWTGRAPH_H

#include <vector>

#include <gbwtgraph/gbwtgraph.h>

/*
  cached_gbwtgraph.h: A cached variant of GBWTGraph.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  A variant of GBWTGraph intended for algorithms that repeatedly access the edges
  of a small subgraph. Provides an easy way of using the cached GBWTGraph interface
  in HandleGraph algorithms.
  NOTE: The cache is not thread-safe. Use a separate CachedGBWTGraph for each
  thread.
  NOTE: For performance reasons, this implementations replicates much of GBWTGraph
  functionality instead of calling it through virtual functions.
*/

class CachedGBWTGraph : public HandleGraph
{
public:
  CachedGBWTGraph();
  CachedGBWTGraph(const CachedGBWTGraph& source);
  CachedGBWTGraph(CachedGBWTGraph&& source);
  virtual ~CachedGBWTGraph();

  explicit CachedGBWTGraph(const GBWTGraph& graph);

  void swap(CachedGBWTGraph& another);
  CachedGBWTGraph& operator=(const CachedGBWTGraph& source);
  CachedGBWTGraph& operator=(CachedGBWTGraph&& source);

  const GBWTGraph* graph;
  gbwt::CachedGBWT cache;

//------------------------------------------------------------------------------

  /*
    Standard HandleGraph interface.
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

  /*
    More efficient reimplementations.
  */

  // Get the number of edges on the right (go_left = false) or left (go_left
  // = true) side of the given handle.
  virtual size_t get_degree(const handle_t& handle, bool go_left) const;

  // Returns true if there is an edge that allows traversal from the left
  // handle to the right handle.
  virtual bool has_edge(const handle_t& left, const handle_t& right) const;

//------------------------------------------------------------------------------

  /*
    GBWTGraph specific interface.
  */

public:

  // Convert gbwt::node_type to handle_t.
  static handle_t node_to_handle(gbwt::node_type node) { return GBWTGraph::node_to_handle(node); }

  // Convert handle_t to gbwt::node_type.
  static gbwt::node_type handle_to_node(const handle_t& handle) { return GBWTGraph::handle_to_node(handle); }

  // Get node sequence as a pointer and length.
  view_type get_sequence_view(const handle_t& handle) const { return this->graph->get_sequence_view(handle); }

  // Determine if the node sequence starts with the given character.
  bool starts_with(const handle_t& handle, char c) const { return this->graph->starts_with(handle, c); }

  // Determine if the node sequence ends with the given character.
  bool ends_with(const handle_t& handle, char c) const { return this->graph->ends_with(handle, c); }

  // Convert handle_t to gbwt::SearchState.
  // Note that the state may be empty if the handle does not correspond to a real node.
  gbwt::SearchState get_state(const handle_t& handle) const { return this->cache.find(handle_to_node(handle)); }

  // Convert handle_t to gbwt::BidirectionalState.
  // Note that the state may be empty if the handle does not correspond to a real node.
  gbwt::BidirectionalState get_bd_state(const handle_t& handle) const { return this->cache.bdFind(handle_to_node(handle)); }

  // Get the search state corresponding to the vector of handles.
  gbwt::SearchState find(const std::vector<handle_t>& path) const { return this->graph->find(this->cache, path); }

  // Get the bidirectional search state corresponding to the vector of handles.
  gbwt::BidirectionalState bd_find(const std::vector<handle_t>& path) const { return this->graph->bd_find(this->cache, path); }

  // Visit all successor states of this state and call iteratee for the state.
  // Stop and return false if the iteratee returns false.
  // Note that this does not visit empty successor states.
  bool follow_paths(gbwt::SearchState state, const std::function<bool(const gbwt::SearchState&)>& iteratee) const
  {
    return this->graph->follow_paths(this->cache, state, iteratee);
  }

  // Visit all predecessor/successor states of this state and call iteratee for the state.
  // Stop and return false if the iteratee returns false.
  // Note that this does not visit empty predecessor/successor states.
  // Each state corresponds to a path. Going backward extends the path left, while going
  // extends it right. When going backward, the state is for the reverse orientation.
  bool follow_paths(gbwt::BidirectionalState state, bool backward,
                    const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const
  {
    return this->graph->follow_paths(this->cache, state, backward, iteratee); 
  }

//------------------------------------------------------------------------------

private:
  void copy(const CachedGBWTGraph& source);
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_CACHED_GBWTGRAPH_H
