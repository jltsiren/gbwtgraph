#ifndef GBWTGRAPH_GBWTGRAPH_H
#define GBWTGRAPH_GBWTGRAPH_H

#include <vector>

#include <gbwt/cached_gbwt.h>

#include <gbwtgraph/utils.h>

/*
  gbwtgraph.h: GBWT-based Handle Graph.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  A HandleGraph implementation that uses GBWT for graph topology and extracts
  sequences from another HandleGraph. Faster sequence access but slower
  graph navigation than in XG. Also supports a version of follow_edges() that
  takes only paths supported by the indexed haplotypes.

  Graph file format versions:

    1  The initial version.
*/

class GBWTGraph : public HandleGraph, public SerializableHandleGraph
{
public:
  GBWTGraph(); // Call deserialize() and set_gbwt() before using the graph.
  GBWTGraph(const GBWTGraph& source);
  GBWTGraph(GBWTGraph&& source);
  ~GBWTGraph();

  GBWTGraph(const gbwt::GBWT& gbwt_index, const HandleGraph& sequence_source);
  GBWTGraph(const gbwt::GBWT& gbwt_index, const SequenceSource& sequence_source);

  void swap(GBWTGraph& another);
  GBWTGraph& operator=(const GBWTGraph& source);
  GBWTGraph& operator=(GBWTGraph&& source);

  struct Header
  {
    std::uint32_t tag, version;
    std::uint64_t nodes;
    std::uint64_t flags;

    constexpr static std::uint32_t TAG = 0x6B3764AF;
    constexpr static std::uint32_t VERSION = Version::GRAPH_VERSION;

    Header();
    void sanitize();
    bool check() const;

    bool operator==(const Header& another) const;
    bool operator!=(const Header& another) const { return !(this->operator==(another)); }
  };

  const gbwt::GBWT*   index;

  Header              header;
  StringArray         sequences;
  sdsl::bit_vector    real_nodes;

  constexpr static size_t CHUNK_SIZE = 1024; // For parallel for_each_handle().

  const static std::string EXTENSION; // ".gg"

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

  /// Get the number of edges on the right (go_left = false) or left (go_left
  /// = true) side of the given handle.
  virtual size_t get_degree(const handle_t& handle, bool go_left) const;

  /// Returns true if there is an edge that allows traversal from the left
  /// handle to the right handle.
  virtual bool has_edge(const handle_t& left, const handle_t& right) const;

//------------------------------------------------------------------------------

  /*
    SerializableHandleGraph interface.
  */

public:

  // Set the GBWT index used for graph topology.
  // Call deserialize() before using the graph.
  // MUST be called before using the graph if the graph is deserialize()-ed.
  void set_gbwt(const gbwt::GBWT& gbwt_index);
  
  /// Return a magic number to identify serialized GBWTGraphs.
  virtual uint32_t get_magic_number() const;
  
protected:
    
  /// Underlying implementation for "serialize" method.
  /// Serialize the sequences to the ostream.
  virtual void serialize_members(std::ostream& out) const;

  /// Underlying implementation to "deserialize" method.
  /// Load the sequences from the istream.
  /// User must call set_gbwt() before using the graph.
  virtual void deserialize_members(std::istream& in);

//------------------------------------------------------------------------------

  /*
    GBWTGraph specific interface.
  */

public:

  // Convert gbwt::node_type to handle_t.
  static handle_t node_to_handle(gbwt::node_type node) { return handlegraph::as_handle(node); }

  // Convert handle_t to gbwt::node_type.
  static gbwt::node_type handle_to_node(const handle_t& handle) { return handlegraph::as_integer(handle); }

  // Get node sequence as a pointer and length.
  view_type get_sequence_view(const handle_t& handle) const;

  // Determine if the node sequence starts with the given character.
  bool starts_with(const handle_t& handle, char c) const;

  // Determine if the node sequence ends with the given character.
  bool ends_with(const handle_t& handle, char c) const;

  // Convert handle_t to gbwt::SearchState.
  // Note that the state may be empty if the handle does not correspond to a real node.
  gbwt::SearchState get_state(const handle_t& handle) const { return this->index->find(handle_to_node(handle)); }

  // Convert handle_t to gbwt::BidirectionalState.
  // Note that the state may be empty if the handle does not correspond to a real node.
  gbwt::BidirectionalState get_bd_state(const handle_t& handle) const { return this->index->bdFind(handle_to_node(handle)); }

  // Get the search state corresponding to the vector of handles.
  gbwt::SearchState find(const std::vector<handle_t>& path) const;

  // Get the bidirectional search state corresponding to the vector of handles.
  gbwt::BidirectionalState bd_find(const std::vector<handle_t>& path) const;

  // Visit all successor states of this state and call iteratee for the state.
  // Stop and return false if the iteratee returns false.
  // Note that this does not visit empty successor states.
  bool follow_paths(gbwt::SearchState state, const std::function<bool(const gbwt::SearchState&)>& iteratee) const
  {
    return this->follow_paths(this->get_single_cache(), state, iteratee);
  }

  // Visit all predecessor/successor states of this state and call iteratee for the state.
  // Stop and return false if the iteratee returns false.
  // Note that this does not visit empty predecessor/successor states.
  // Each state corresponds to a path. Going backward extends the path left, while going
  // extends it right. When going backward, the state is for the reverse orientation.
  bool follow_paths(gbwt::BidirectionalState state, bool backward,
                    const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const
  {
    return this->follow_paths(this->get_single_cache(), state, backward, iteratee); 
  }

//------------------------------------------------------------------------------

  /*
    Cached GBWTGraph specific interface. Each thread must use a separate cache.
  */

  // Return a cache for the GBWT index. Note: The cache is not thread-safe.
  gbwt::CachedGBWT get_cache() const { return gbwt::CachedGBWT(*(this->index), false); }

  // Return a single-node cache for the GBWT index. Mostly for internal use.
  // Note: The cache is not thread-safe.
  gbwt::CachedGBWT get_single_cache() const { return gbwt::CachedGBWT(*(this->index), true); }

  // Convert handle_t to gbwt::SearchState.
  gbwt::SearchState get_state(const gbwt::CachedGBWT& cache, const handle_t& handle) const
  {
    return cache.find(handle_to_node(handle));
  }

  // Convert handle_t to gbwt::BidirectionalState.
  gbwt::BidirectionalState get_bd_state(const gbwt::CachedGBWT& cache, const handle_t& handle) const
  {
    return cache.bdFind(handle_to_node(handle));
  }

  // Get the search state corresponding to the vector of handles.
  gbwt::SearchState find(const gbwt::CachedGBWT& cache, const std::vector<handle_t>& path) const;

  // Get the bidirectional search state corresponding to the vector of handles.
  gbwt::BidirectionalState bd_find(const gbwt::CachedGBWT& cache, const std::vector<handle_t>& path) const;

  // Visit all successor states of this state and call iteratee for the state.
  // Stop and return false if the iteratee returns false.
  // Note that the state may be empty if no path continues to that node.
  bool follow_paths(const gbwt::CachedGBWT& cache, gbwt::SearchState state,
                    const std::function<bool(const gbwt::SearchState&)>& iteratee) const;

  // Visit all predecessor/successor states of this state and call iteratee for the state.
  // Stop and return false if the iteratee returns false.
  // Note that the state may be empty if no path continues to that node.
  // Each state corresponds to a path. Going backward extends the path left, while going
  // extends it right. When going backward, the state is for the reverse orientation.
  bool follow_paths(const gbwt::CachedGBWT& cache, gbwt::BidirectionalState state, bool backward,
                    const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const;

  // Loop over all the handles to next/previous (right/left) nodes. Passes
  // them to a callback which returns false to stop iterating and true to
  // continue. Returns true if we finished and false if we stopped early.
  bool cached_follow_edges(const gbwt::CachedGBWT& cache, const handle_t& handle, bool go_left,
                           const std::function<bool(const handle_t&)>& iteratee) const;

//------------------------------------------------------------------------------

private:
  friend class CachedGBWTGraph;

  // Construction helpers.
  void determine_real_nodes();

  void copy(const GBWTGraph& source);

  size_t node_offset(gbwt::node_type node) const { return node - this->index->firstNode(); }
  size_t node_offset(const handle_t& handle) const { return this->node_offset(handle_to_node(handle)); }
};

//------------------------------------------------------------------------------

/*
  Traverse all haplotype-consistent windows in the graph and call lambda() for each window.
  Uses multiple threads, so the lambda should be thread-safe.
  A window starts with the sequence of a node and is followed by window_size - 1 bases
  from subsequent nodes. If no extensions are possible, a shorter substring of
  length >= window_size also qualifies as a window.
*/
void for_each_haplotype_window(const GBWTGraph& graph, size_t window_size,
                               const std::function<void(const std::vector<handle_t>&, const std::string&)>& lambda,
                               bool parallel);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_GBWTGRAPH_H
