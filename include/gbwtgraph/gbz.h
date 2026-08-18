#ifndef GBWTGRAPH_GBZ_H
#define GBWTGRAPH_GBZ_H

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/utils.h>

#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>
#endif

/*
  gbz.h: GBZ file format.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  GBZ file format wrapper, as specified in SERIALIZATION.md. The wrapper owns the
  GBWT index and the GBWTGraph.

  Constructors, serialization, and loading throw `std::runtime_error` on failure.

  File format versions:

    2  Contains GBWTGraph version 4 with zstd-compressed sequences.
       Compatible with version 1.

    1  The initial version.
*/

template <typename SAAllocator = std::allocator<char>> class CoreGBZ;

// The allocator only affects `graph.sequences`; every other member and
// operation is the same regardless of it. Callers that don't care about
// shared memory should use this plain-allocator instantiation, not
// `CoreGBZ` directly. Declared here (rather than after the class) because
// `CoreGBZ`'s own merge and supergraph constructors take this type by value
// and reference.
typedef CoreGBZ<std::allocator<char>> GBZ;

template <typename SAAllocator>
class CoreGBZ
{
public:
  // This is a valid graph, unlike the default GBWTGraph.
  CoreGBZ();
#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
  // As above, with the graph's sequences bound to a segment (see
  // CoreGBWTGraph's matching constructor).
  CoreGBZ(bi::managed_shared_memory* shared_memory)
    requires gbwt::SharedMemory<SAAllocator>;
#endif
  CoreGBZ(const CoreGBZ& source);
  CoreGBZ(CoreGBZ&& source);
  ~CoreGBZ();

  /*
    Build GBZ from a GBWT index and a graph, either as the structures
    `gfa_to_gbwt()` returns (which are reset to `nullptr`) or as a plain index
    and graph (which are copied; mostly for testing). Both call
    compute_pggname() internally.

    As with CoreGBWTGraph, the versions without a segment build the graph's
    sequences on the heap and so are absent for an allocator that puts
    characters in shared memory.
  */

  CoreGBZ(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<NaiveGraph>& graph)
    requires (!gbwt::SharedMemory<SAAllocator>);

  CoreGBZ(const gbwt::GBWT& index, const NaiveGraph& graph)
    requires (!gbwt::SharedMemory<SAAllocator>);

#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
  CoreGBZ(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<NaiveGraph>& graph, bi::managed_shared_memory* shared_memory)
    requires gbwt::SharedMemory<SAAllocator>;

  CoreGBZ(const gbwt::GBWT& index, const NaiveGraph& graph, bi::managed_shared_memory* shared_memory)
    requires gbwt::SharedMemory<SAAllocator>;
#endif

private:
  // The body shared by the two `gfa_to_gbwt()` constructors above. Any
  // trailing arguments say where the graph's sequences are to be stored and go
  // straight to CoreGBWTGraph.
  template<typename... StorageArgs>
  void take_index_and_graph(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<NaiveGraph>& graph, StorageArgs&&... storage);

public:

  // Builds GBZ from a set of non-overlapping subgraphs.
  // Moves the subgraphs out of the provided vector.
  // Fails with `std::exit()` if the graphs have overlapping node ids.
  // Sets the union of reference samples in the GBWT indexes as reference
  // samples. Does not create a node-to-segment translation.
  // Calls compute_pggname() internally but does not set any relationships.
  // Subgraphs are typically small and short-lived, so merging is only
  // defined for the plain, heap-allocated GBZ.
  explicit CoreGBZ(std::vector<GBZ>&& subgraphs)
    requires (!gbwt::SharedMemory<SAAllocator>);

  // Builds a GBZ from a GBWT index and a GBZ supergraph.
  // Reference samples in the supergraph that are present in the GBWT
  // index will be set as reference samples. Calls compute_pggname()
  // internally. The provided GBWT index will be moved into the GBZ.
  // Subgraphs are always plain, heap-allocated, regardless of the
  // supergraph's allocator; see CoreGBWTGraph::subgraph().
  CoreGBZ(gbwt::GBWT&& index, const GBZ& supergraph)
    requires (!gbwt::SharedMemory<SAAllocator>);

  // Build GBZ from a GBWT index and a `HandleGraph`, with an optional
  // translation. Because the parent graph does not store GraphName
  // information, compute_pggname() must be called separately after
  // construction. The provided GBWT index will be moved into the GBZ.
  CoreGBZ(gbwt::GBWT&& index, const HandleGraph& graph, const NamedNodeBackTranslation* segment_space)
    requires (!gbwt::SharedMemory<SAAllocator>);

  void swap(CoreGBZ& another);
  CoreGBZ& operator=(const CoreGBZ& source);
  CoreGBZ& operator=(CoreGBZ&& source);

//------------------------------------------------------------------------------

  /*
    Stable graph names (pggname) and known relationships between graphs.
    See `GraphName` documentation in utils.h.
  */

  enum class ParentGraphType {
    // Determine the relationship heuristically.
    HEURISTIC,
    // The parent graph is a supergraph of this graph, unless the names are the same.
    SUPERGRAPH,
    // The parent graph is a translation target of this graph, unless the names are the same.
    TRANSLATION_TARGET
  };

  /*
    Computes the pggname for this graph and stores it in the tags.
    Returns true on success, false on failure.

    If a parent graph is given and it has a set name, adds the corresponding
    relationship and imports all known relationships from the other graph.
    If an explicit relationship type is not given, the following heuristic
    will be used:

    1. If the GBWTGraph has a node-to-segment translation, the relationship
       is a translation to the parent graph.
    2. Otherwise, if the parent graph's pggname is different from the computed
       name for this graph, the relationship is a subgraph relationship.
    3. Otherwise, no relationship is added.

    When the GBZ is built from a SequenceSource or another GBZ, this function
    is called automatically by the constructor. When the parent graph is a
    generic HandleGraph, GraphName information cannot be imported and this
    function must be called separately after construction.
  */
  bool compute_pggname(const GraphName* parent, ParentGraphType relationship = ParentGraphType::HEURISTIC);

  // Returns the graph name object for this graph based on the information
  // stored in the tags.
  GraphName graph_name() const { return GraphName(this->tags); }

  // Returns the pggname for this graph, or an empty string if not set.
  std::string pggname() const { return this->tags.get(GraphName::GBZ_NAME_TAG); }

  // Returns the pggname of the translation target for the node-to-segment
  // translation, or an empty string if not set.
  std::string translation_target() const { return this->tags.get(GraphName::GBZ_TRANSLATION_TARGET_TAG); }

//------------------------------------------------------------------------------

  /*
    Reference samples and named paths. Reference samples cannot be changed
    conveniently within GBWTGraph, because the pointer to the GBWT index is const.
  */

  // Sets the given sample names as reference samples, but only if they are
  // present in the GBWT metadata. Returns the number of reference samples.
  // This is somewhat expensive, as the GBWTGraph must recache named paths.
  size_t set_reference_samples(const sample_name_set& samples);

  // Returns the set of reference samples.
  // Some of these samples may not exist in the GBWT metadata.
  const sample_name_set& get_reference_samples() const
  {
    return this->graph.reference_samples;
  }

  // Returns the number of paths in the graph.
  size_t paths() const { return this->index.metadata.paths(); }

  // Returns the number of reference and generic paths in the graph.
  size_t named_paths() const { return this->graph.named_paths.size(); }

//------------------------------------------------------------------------------

  struct Header
  {
    std::uint32_t tag, version;
    std::uint64_t flags;

    constexpr static std::uint32_t TAG = 0x205A4247; // "GBZ "
    constexpr static std::uint32_t VERSION = Version::GBZ_VERSION;

    constexpr static std::uint64_t FLAG_MASK = 0x0000;

    // Compatible versions.
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

//------------------------------------------------------------------------------

  Header     header;
  gbwt::Tags tags;
  gbwt::GBWT index;
  // There is no `shared_memory` member here to query: check
  // `graph.sequences.shared_memory` instead.
  CoreGBWTGraph<SAAllocator>  graph;

  const static std::string EXTENSION; // ".gbz"

  // Serialize the the GBZ into the output stream in the Simple-SDS format.
  void simple_sds_serialize(std::ostream& out) const;

  // Serialize the graph into the output stream in the Simple-SDS format, version 1.
  // This is for compatibility with tools that do not support version 2 with Zstandard compression.
  void simple_sds_serialize_v1(std::ostream& out) const;

  // Serialize the given GBWT and GBWTGraph objects in the GBZ format.
  // NOTE: GBZ tags will be empty, except for the source tag.
  static void simple_sds_serialize(const gbwt::GBWT& index, const CoreGBWTGraph<SAAllocator>& graph, std::ostream& out);

  // Deserialize or decompress the GBZ from the input stream.
  void simple_sds_load(std::istream& in);

  // Returns the size of the serialized structure in elements.
  size_t simple_sds_size() const;

  // Returns the GBZ tags from the given file, without loading it entirely.
  static gbwt::Tags simple_sds_load_tags(const std::string& filename);

  // Serialize the GBWT (Simple-SDS format) and the GBWTGraph to separate files.
  // Default graph format is libhandlegraph / SDSL.
  // NOTE: GBZ tags are not serialized.
  void serialize_to_files(const std::string& gbwt_name, const std::string& graph_name, bool simple_sds_graph = false) const;

  // Loads the GBWT (Simple-SDS format) and the GBWTGraph from separate files.
  // Graph format is libhandlegraph / SDSL; the Simple-SDS format cannot be read.
  // NOTE: GBZ tags are not loaded.
  void load_from_files(const std::string& gbwt_name, const std::string& graph_name);

private:
  void copy(const CoreGBZ& source);
  void reset_tags();
  void add_source();
  void set_gbwt();
  void set_gbwt_address();
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph


#endif // GBWTGRAPH_GBZ_H
