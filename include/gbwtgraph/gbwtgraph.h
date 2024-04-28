#ifndef GBWTGRAPH_GBWTGRAPH_H
#define GBWTGRAPH_GBWTGRAPH_H

#include <limits>
#include <vector>

#include <gbwt/cached_gbwt.h>
#include <gbwt/metadata.h>

#include <gbwtgraph/utils.h>

/*
  gbwtgraph.h: GBWT-based Handle Graph.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  A HandleGraph implementation that uses GBWT for graph topology and extracts
  sequences from another HandleGraph.

  Main features:

    * Faster sequence access but slower graph navigation than in XG.
    * `follow_paths()`: A version of `follow_edges()` that only uses edges
      consistent with the indexed haplotypes.
    * Exposes haplotype, reference, and generic paths.
    * Translation from (intervals of) node ids to GFA segment names.

  The `PathHandleGraph` interface requires GBWT metadata with sample/contig/path
  names. Information about the named paths is cached using multiple OpenMP
  threads when the graph is created or loaded.

  Graph file format versions:

    3  The compressed version uses simple-sds serialization. Non-compressed
       (SDSL) file format is compatible with versions 1 and 2.

    2  Translation between GFA segment names and (intervals of) node ids.
       Optional compressed serialization format. Compatible with version 1.

    1  The initial version.
*/

class GBWTGraph : public PathHandleGraph, public SerializableHandleGraph, public NamedNodeBackTranslation
{
public:
  GBWTGraph(); // Call (deserialize() and set_gbwt()) or simple_sds_load() before using the graph.
  GBWTGraph(const GBWTGraph& source);
  GBWTGraph(GBWTGraph&& source);
  virtual ~GBWTGraph();

  // Build the graph from another `HandleGraph` and an optional named segment space over it.
  GBWTGraph(const gbwt::GBWT& gbwt_index, const HandleGraph& sequence_source, const NamedNodeBackTranslation* segment_space = nullptr);

  // Build the graph (and possibly the translation) from a `SequenceSource` object.
  // If the translation is present, some parts of the construction are multithreaded.
  GBWTGraph(const gbwt::GBWT& gbwt_index, const SequenceSource& sequence_source);

  // Returns a GBWTGraph for the subgraph defined by the given GBWT index.
  // This is faster than using the graph as a HandleGraph in the constructor.
  GBWTGraph subgraph(const gbwt::GBWT& gbwt_index) const;

  // Makes some sanity checks on the internal consistency of the structure.
  // Requires that the GBWT index has been set.
  // This is called internally during deserialization an is also used in tests.
  // Throws sdsl::simple_sds::InvalidData if the checks fail.
  void sanity_checks();

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

    constexpr static std::uint64_t FLAG_MASK = 0x0003;
    constexpr static std::uint64_t FLAG_TRANSLATION = 0x0001;
    constexpr static std::uint64_t FLAG_SIMPLE_SDS = 0x0002;

    // Old compatible versions.
    constexpr static std::uint32_t TRANS_VERSION = 2;
    constexpr static std::uint64_t TRANS_FLAG_MASK = 0x0001;

    constexpr static std::uint32_t OLD_VERSION = 1;
    constexpr static std::uint64_t OLD_FLAG_MASK = 0x0000;

    Header();

    // Throws `sdsl::simple_sds::InvalidData` if the header is invalid.
    void check() const;

    void set_version() { this->version = VERSION; }

    void set(std::uint64_t flag) { this->flags |= flag; }
    void unset(std::uint64_t flag) { this->flags &= ~flag; }
    bool get(std::uint64_t flag) const { return (this->flags & flag); }

    bool operator==(const Header& another) const;
    bool operator!=(const Header& another) const { return !(this->operator==(another)); }
  };

  const gbwt::GBWT* index;

  Header            header;
  gbwt::StringArray sequences;
  sdsl::bit_vector  real_nodes;

  // Segment to node translation. Node `v` maps to segment `node_to_segment.predecessor(v)->first`.
  gbwt::StringArray segments;
  sdsl::sd_vector<> node_to_segment;

  // Cached named path information.
  std::vector<NamedPath>                      named_paths;
  std::unordered_map<std::string, size_t>     name_to_path; // To offset in `named_paths`.
  std::unordered_map<gbwt::size_type, size_t> id_to_path; // To offset in `named_paths`.
  std::unordered_set<std::string>             reference_samples; // Parsed from tags in the GBWT.
  // Path handles are either indexes into named_paths, or, if larger than
  // named_paths, are an offset of the size of named_paths plus a path number in
  // our metadata object. This syntactically allows for aliasing: cached paths
  // in named_paths also ultimately refer to path numbers in the metadata. So,
  // when iterating, we need to remember to skip path numbers in the metadata
  // that are also cached named paths.

  constexpr static size_t CHUNK_SIZE = 1024; // For parallel for_each_handle().

  // TODO: This should be 0, as in GFA W-lines.
  constexpr static size_t NO_PHASE = std::numeric_limits<gbwt::PathName::path_name_type>::max();

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

  // Get the number of edges on the right (go_left = false) or left (go_left
  // = true) side of the given handle.
  virtual size_t get_degree(const handle_t& handle, bool go_left) const;

  // Returns true if there is an edge that allows traversal from the left
  // handle to the right handle.
  virtual bool has_edge(const handle_t& left, const handle_t& right) const;

//------------------------------------------------------------------------------

  /*
    PathHandleGraph interface.
  */

public:

  /// Returns the number of paths stored in the graph
  virtual size_t get_path_count() const;

  /// Determine if a path name exists and is legal to get a path handle for.
  virtual bool has_path(const std::string& path_name) const;

  /// Look up the path handle for the given path name.
  /// The path with that name must exist.
  virtual path_handle_t get_path_handle(const std::string& path_name) const;

  /// Look up the name of a path from a handle to it
  virtual std::string get_path_name(const path_handle_t& path_handle) const;

  /// Look up whether a path is circular
  virtual bool get_is_circular(const path_handle_t& path_handle) const;

  /// Returns the number of node steps in the path
  virtual size_t get_step_count(const path_handle_t& path_handle) const;

  /// Returns the number of node steps on a handle
  virtual size_t get_step_count(const handle_t& handle) const;

  /// Get a node handle (node ID and orientation) from a handle to an step on a path
  virtual handle_t get_handle_of_step(const step_handle_t& step_handle) const;

  /// Returns a handle to the path that an step is on
  virtual path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;

  /// Get a handle to the first step, which will be an arbitrary step in a circular path
  /// that we consider "first" based on our construction of the path. If the path is empty,
  /// then the implementation must return the same value as path_end().
  virtual step_handle_t path_begin(const path_handle_t& path_handle) const;

  /// Get a handle to a fictitious position past the end of a path. This position is
  /// returned by get_next_step for the final step in a path in a non-circular path.
  /// Note: get_next_step will *NEVER* return this value for a circular path.
  virtual step_handle_t path_end(const path_handle_t& path_handle) const;

  /// Get a handle to the last step, which will be an arbitrary step in a circular path that
  /// we consider "last" based on our construction of the path. If the path is empty
  /// then the implementation must return the same value as path_front_end().
  virtual step_handle_t path_back(const path_handle_t& path_handle) const;

  /// Get a handle to a fictitious position before the beginning of a path. This position is
  /// return by get_previous_step for the first step in a path in a non-circular path.
  /// Note: get_previous_step will *NEVER* return this value for a circular path.
  virtual step_handle_t path_front_end(const path_handle_t& path_handle) const;

  /// Returns true if the step is not the last step in a non-circular path.
  virtual bool has_next_step(const step_handle_t& step_handle) const;

  /// Returns true if the step is not the first step in a non-circular path.
  virtual bool has_previous_step(const step_handle_t& step_handle) const;

  /// Returns a handle to the next step on the path. If the given step is the final step
  /// of a non-circular path, this method has undefined behavior. In a circular path,
  /// the "last" step will loop around to the "first" step.
  virtual step_handle_t get_next_step(const step_handle_t& step_handle) const;

  /// Returns a handle to the previous step on the path. If the given step is the first
  /// step of a non-circular path, this method has undefined behavior. In a circular path,
  /// it will loop around from the "first" step (i.e. the one returned by path_begin) to
  /// the "last" step.
  virtual step_handle_t get_previous_step(const step_handle_t& step_handle) const;

  using PathHandleGraph::for_each_path_handle;
  using PathHandleGraph::for_each_step_on_handle;

protected:

  /// Execute a function on each path in the graph. If it returns false, stop
  /// iteration. Returns true if we finished and false if we stopped early.
  virtual bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;

  /// Execute a function on each step of a handle in any path. If it
  /// returns false, stop iteration. Returns true if we finished and false if
  /// we stopped early.
  virtual bool for_each_step_on_handle_impl(const handle_t& handle,
      const std::function<bool(const step_handle_t&)>& iteratee) const;

//------------------------------------------------------------------------------

  /*
    PathMetadata interface, actually exposing threads.
  */

public:

    /// What is the given path meant to be representing?
    virtual PathSense get_sense(const path_handle_t& handle) const;

    /// Get the name of the sample or assembly asociated with the
    /// path-or-thread, or NO_SAMPLE_NAME if it does not belong to one.
    virtual std::string get_sample_name(const path_handle_t& handle) const;

    /// Get the name of the contig or gene asociated with the path-or-thread,
    /// or NO_LOCUS_NAME if it does not belong to one.
    virtual std::string get_locus_name(const path_handle_t& handle) const;

    /// Get the haplotype number (0 or 1, for diploid) of the path-or-thread,
    /// or NO_HAPLOTYPE if it does not belong to one.
    virtual size_t get_haplotype(const path_handle_t& handle) const;

    /// Get the phase block number (contiguously phased region of a sample,
    /// contig, and haplotype) of the path-or-thread, or NO_PHASE_BLOCK if it
    /// does not belong to one.
    virtual size_t get_phase_block(const path_handle_t& handle) const;

    /// Get the bounds of the path-or-thread that are actually represented
    /// here. Should be NO_SUBRANGE if the entirety is represented here, and
    /// 0-based inclusive start and exclusive end positions of the stored
    /// region on the full path-or-thread if a subregion is stored.
    ///
    /// If no end position is stored, NO_END_POSITION may be returned for the
    /// end position.
    virtual subrange_t get_subrange(const path_handle_t& handle) const;

protected:

    /// Loop through all the paths matching the given query. Query elements
    /// which are null match everything. Returns false and stops if the
    /// iteratee returns false.
    virtual bool for_each_path_matching_impl(const std::unordered_set<PathSense>* senses,
                                             const std::unordered_set<std::string>* samples,
                                             const std::unordered_set<std::string>* loci,
                                             const std::function<bool(const path_handle_t&)>& iteratee) const;

    /// Loop through all steps on the given handle for paths with the given
    /// sense. Returns false and stops if the iteratee returns false.
    virtual bool for_each_step_of_sense_impl(const handle_t& visited, const PathSense& sense, const std::function<bool(const step_handle_t&)>& iteratee) const;

private:

    /// Internal iteration method to find all the GBWT edges and their path
    /// numbers on a node. Only looks at forward sequence for each path, but
    /// looks at both orientations of the node.
    bool for_each_edge_and_path_on_handle(const handle_t& handle, const std::function<bool(const gbwt::edge_type&, const gbwt::size_type&)>& iteratee) const;

    /// Get all the sample numbers that might be relevant for the given
    /// user-visible sample name.
    std::vector<gbwt::size_type> sample_numbers_for_sample_name(const std::unordered_set<PathSense>* senses, const std::string& sample_name) const;

    /// Iterate over all paths of the given senses with the given sample (which
    /// can be NO_SAMPLE_NAME) and locus
    bool for_each_path_matching_sample_and_locus(const std::unordered_set<PathSense>* senses,
                                                 const std::string& sample_name,
                                                 const std::string& locus_name,
                                                 const std::function<bool(const path_handle_t&)>& iteratee) const;

    /// Iterate over all paths of the given senses with the given sample (which
    /// can be NO_SAMPLE_NAME)
    bool for_each_path_matching_sample(const std::unordered_set<PathSense>* senses,
                                       const std::string& sample_name,
                                       const std::function<bool(const path_handle_t&)>& iteratee) const;

    /// Iterate over all paths of the given senses with the given locus
    bool for_each_path_matching_locus(const std::unordered_set<PathSense>* senses,
                                      const std::string& locus_name,
                                      const std::function<bool(const path_handle_t&)>& iteratee) const;

//------------------------------------------------------------------------------

  /*
    SerializableHandleGraph interface. Serialization / deserialization throws
    `std::runtime_error` on failure.
  */

public:

  // Set the GBWT index used for graph topology and cache named path information.
  // Call deserialize() before using the graph.
  // MUST be called before using the graph if the graph is deserialize()-ed.
  void set_gbwt(const gbwt::GBWT& gbwt_index);

  // Set the address of the GBWT index. This assumes that the new index is identical
  // to the one that was already in use. Intended for implementing containers that
  // store both GBWT and GBWTGraph and may get moved around.
  void set_gbwt_address(const gbwt::GBWT& gbwt_index);

  /// Return a magic number to identify serialized GBWTGraphs.
  virtual uint32_t get_magic_number() const;

protected:

  // Underlying implementation for "serialize" method.
  // Serialize the sequences to the ostream in SDSL format.
  virtual void serialize_members(std::ostream& out) const;

  // Underlying implementation to "deserialize" method.
  // Load the sequences from the istream.
  // User must call set_gbwt() before using the graph.
  virtual void deserialize_members(std::istream& in);

//------------------------------------------------------------------------------

  /*
    WIP SegmentHandleGraph interface.

    NOTE: The implementation stores the translations also for segments that were not used
    on any path. Because the corresponding nodes are missing from the graph, the methods
    will treat them like missing segments.
  */

public:

  /// Returns `true` if the graph contains a translation from node ids to segment names.
  virtual bool has_segment_names() const;

  /// Returns (GFA segment name, semiopen node id range) containing the handle.
  /// If there is no such translation, returns ("id", (id, id + 1)).
  virtual std::pair<std::string, std::pair<nid_t, nid_t>> get_segment(const handle_t& handle) const;

  /// Returns (GFA segment name, starting offset in the same orientation) for the handle.
  /// If there is no translation, returns ("id", 0).
  virtual std::pair<std::string, size_t> get_segment_name_and_offset(const handle_t& handle) const;

  /// Returns the name of the original GFA segment corresponding to the handle.
  /// If there is no translation, returns the node id as a string.
  virtual std::string get_segment_name(const handle_t& handle) const;

  /// Returns the starting offset in the original GFA segment corresponding to the handle
  /// in the same orientation as the handle.
  /// If there is no translation, returns 0.
  virtual size_t get_segment_offset(const handle_t& handle) const;

  /// Calls `iteratee` with each segment name as a string, and the semiopen
  /// interval of node ids corresponding to it as a std::pair of nid_t
  /// values. Stops early if the call returns `false`.
  /// Returns false if iteration was stopped, and true otherwise.
  template<typename Iteratee>
  bool for_each_segment(const Iteratee& iteratee, bool parallel = false) const {
    return for_each_segment_impl(handlegraph::BoolReturningWrapper<Iteratee>::wrap(iteratee), parallel);
  }

  /// Calls `iteratee` with each inter-segment edge (as an edge_t) and the
  /// corresponding segment names in the canonical orientation (as two
  /// strings). Stops early if the call returns `false`.
  /// Returns false if iteration was stopped, and true otherwise.
  template<typename Iteratee>
  bool for_each_link(const Iteratee& iteratee, bool parallel = false) const {
    return for_each_link_impl(handlegraph::BoolReturningWrapper<Iteratee>::wrap(iteratee), parallel);
  }

protected:

  // Calls `iteratee` with each segment name and the semiopen interval of node ids
  // corresponding to it. Stops early if the call returns `false`.
  // In GBWTGraph, the segments are visited in sorted order by node ids.
  virtual bool for_each_segment_impl(const std::function<bool(const std::string&, const std::pair<nid_t, nid_t>&)>& iteratee, bool parallel) const;

  // Calls `iteratee` with each inter-segment edge and the corresponding segment names
  // in the canonical orientation. Stops early if the call returns `false`.
  virtual bool for_each_link_impl(const std::function<bool(const edge_t&, const std::string&, const std::string&)>& iteratee, bool parallel) const;

//------------------------------------------------------------------------------

  /*
    NamedNodeBackTranslation interface.
  */

public:

  /// Translate a node range back to segment space.
  /// Return value is not defined if no segments exist or the node does not
  /// exist.
  virtual std::vector<oriented_node_range_t> translate_back(const oriented_node_range_t& range) const;

  /// Get a segment name. Return value is not defined if no segments exist or
  /// if the segment is out of range.
  virtual std::string get_back_graph_node_name(const nid_t& back_node_id) const;

//------------------------------------------------------------------------------

  /*
    GBWTGraph specific interface. Serialization / deserialization throws
    `std::runtime_error` on failure.
  */

public:

  // Serialize the the graph into the output stream in the simple-sds format.
  void simple_sds_serialize(std::ostream& out) const;

  // Deserialize or decompress the graph from the input stream and set the given
  // GBWT index. Note that the GBWT index is essential for loading the structure.
  void simple_sds_load(std::istream& in, const gbwt::GBWT& gbwt_index);

  // Returns the size of the serialized structure in elements.
  size_t simple_sds_size() const;

  // Convert gbwt::node_type to handle_t.
  static handle_t node_to_handle(gbwt::node_type node) { return handlegraph::as_handle(node); }

  // Convert handle_t to gbwt::node_type.
  static gbwt::node_type handle_to_node(const handle_t& handle) { return handlegraph::as_integer(handle); }

  // Convert GBWT path id to path_handle_t.
  path_handle_t path_to_handle(gbwt::size_type path) const;

  // Convert path_handle_t to GBWT path id.
  gbwt::size_type handle_to_path(const path_handle_t& path_handle) const;

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
  void cache_named_paths();

  void copy(const GBWTGraph& source);

  // Copies a translation over this graph's node IDs from a
  // NamedNodeBackTranslation into the internal representation used in a
  // GBWTGraph.
  // Returns `StringArray` of segment names and `sd_vector<>` mapping node ids to names.
  // Runs of nonexistent nodes become segments with empty names.
  // Throws if the translation cannot be represented (i.e. segments aren't
  // forward strands of contiguous ascending node ID ranges).
  std::pair<gbwt::StringArray, sdsl::sd_vector<>>
  copy_translation(const NamedNodeBackTranslation& translation) const;

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

/*
  As above, but prioritize more windows with less overall redundancy. If the length of the
  initial node is at least window_size, it becomes a separate window. Extension windows
  then take only the last window_size - 1 bases from it.
*/
void for_each_nonredundant_window(
  const GBWTGraph& graph, size_t window_size,
  const std::function<void(const std::vector<handle_t>&, const std::string&)>& lambda,
  bool parallel);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_GBWTGRAPH_H
