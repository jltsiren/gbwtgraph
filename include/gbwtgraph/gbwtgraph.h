#ifndef GBWTGRAPH_GBWTGRAPH_H
#define GBWTGRAPH_GBWTGRAPH_H

#include <vector>

#include <gbwt/cached_gbwt.h>
#include <handlegraph/util.hpp>
#include <sdsl/int_vector.hpp>

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

  /*
    The sequence source must implement the following subset of the HandleGraph interface
    for all nodes in forward orientation:
    - get_handle()
    - get_length()
    - get_sequence()
  */
  template<class Source>
  GBWTGraph(const gbwt::GBWT& gbwt_index, const Source& sequence_source);

  void swap(GBWTGraph& another);
  GBWTGraph& operator=(const GBWTGraph& source);
  GBWTGraph& operator=(GBWTGraph&& source);

  struct Header
  {
    std::uint32_t tag, version;
    std::uint64_t nodes;
    std::uint64_t flags;

    constexpr static std::uint32_t TAG = 0x6B3764AF;
    constexpr static std::uint32_t VERSION = 1;
    constexpr static std::uint32_t MIN_VERSION = 1;

    Header();
    void sanitize();
    bool check() const;

    bool operator==(const Header& another) const;
    bool operator!=(const Header& another) const { return !(this->operator==(another)); }
  };

  const gbwt::GBWT*   index;

  Header              header;
  std::vector<char>   sequences;
  sdsl::int_vector<0> offsets;
  sdsl::bit_vector    real_nodes;

  constexpr static size_t CHUNK_SIZE = 1024; // For parallel for_each_handle().
  constexpr static size_t BLOCK_SIZE = 64 * gbwt::MEGABYTE; // For serialization.

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

//------------------------------------------------------------------------------

  /*
    SerializableHandleGraph interface.
  */

public:

  // Serialize the sequences to the ostream.
  virtual void serialize(std::ostream& out) const;

  // Load the sequences from the istream.
  // Call set_gbwt() before using the graph.
  virtual void deserialize(std::istream& in);

  // Set the GBWT index used for graph topology.
  // Call deserialize() before using the graph.
  void set_gbwt(const gbwt::GBWT& gbwt_index);

//------------------------------------------------------------------------------

  /*
    GBWTGraph specific interface.
  */

public:

  // In-place view of the sequence; (start, length).
  typedef std::pair<const char*, size_t> view_type;

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
  gbwt::SearchState get_state(const handle_t& handle) const { return this->index->find(handle_to_node(handle)); }

  // Convert handle_t to gbwt::BidirectionalState.
  gbwt::BidirectionalState get_bd_state(const handle_t& handle) const { return this->index->bdFind(handle_to_node(handle)); }

  // Get the search state corresponding to the vector of handles.
  gbwt::SearchState find(const std::vector<handle_t>& path) const;

  // Get the bidirectional search state corresponding to the vector of handles.
  gbwt::BidirectionalState bd_find(const std::vector<handle_t>& path) const;

  // Visit all successor states of this state and call iteratee for the state.
  // Stop and return false if the iteratee returns false.
  // Note that the state may be empty if no path continues to that node.
  bool follow_paths(gbwt::SearchState state, const std::function<bool(const gbwt::SearchState&)>& iteratee) const;

  // Visit all predecessor/successor states of this state and call iteratee for the state.
  // Stop and return false if the iteratee returns false.
  // Note that the state may be empty if no path continues to that node.
  // Each state corresponds to a path. Going backward extends the path left, while going
  // extends it right. When going backward, the state is for the reverse orientation.
  bool follow_paths(gbwt::BidirectionalState state, bool backward,
                    const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const;

//------------------------------------------------------------------------------

  /*
    Cached GBWTGraph specific interface. Each thread must use a separate cache.
  */

  // Return a cache for the GBWT index. Note: The cache is not thread-safe.
  gbwt::CachedGBWT get_cache() const { return gbwt::CachedGBWT(*(this->index), false); }

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

//------------------------------------------------------------------------------

private:
  // Construction helpers.
  void determine_real_nodes();
  std::vector<handle_t> cache_handles(const std::function<handle_t(gbwt::node_type)>& get_source_handle) const;
  void allocate_arrays(const std::function<size_t(size_t)>& get_source_length);
  void cache_sequences(const std::function<std::string(size_t)>& get_source_sequence);

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

// Implementation of the main GBWTGraph constructor.
template<class Source>
GBWTGraph::GBWTGraph(const gbwt::GBWT& gbwt_index, const Source& sequence_source) :
  index(nullptr), header()
{
  // Set GBWT and do sanity checks.
  this->set_gbwt(gbwt_index);

  // Add the sentinel to the offset vector of an empty graph just in case.
  if(this->index->empty())
  {
    this->offsets = sdsl::int_vector<0>(1, 0);
    return;
  }

  // Determine the real nodes and cache the handles. We do these in parallel, as the first
  // is always a slow operation and the second is equally slow in XG.
  // Node n is real, if real_nodes[node_offset(n) / 2] is true.
  std::vector<handle_t> handle_cache;
  #pragma omp parallel
  {
    #pragma omp single
    {
      #pragma omp task
      {
        this->determine_real_nodes();
      }
      #pragma omp task
      {
        handle_cache = this->cache_handles([&sequence_source](gbwt::node_type node) -> handle_t
        {
          return sequence_source.get_handle(gbwt::Node::id(node), false);
        });
      }
    }
  }

  // Allocate space for the sequence and offset arrays.
  this->allocate_arrays([&sequence_source, &handle_cache](size_t offset) -> size_t
  {
    return sequence_source.get_length(handle_cache[offset]);
  });

  // Store the concatenated sequences and their offset ranges for both orientations of all nodes.
  // Given GBWT node n, the sequence is sequences[node_offset(n)] to sequences[node_offset(n + 1) - 1].
  this->cache_sequences([&sequence_source, &handle_cache](size_t offset) -> std::string
  {
    return sequence_source.get_sequence(handle_cache[offset]);
  });
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_GBWTGRAPH_H
