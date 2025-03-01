#ifndef GBWTGRAPH_GBZ_H
#define GBWTGRAPH_GBZ_H

#include <gbwtgraph/gbwtgraph.h>

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

    1  The initial version.
*/

class GBZ
{
public:
  GBZ();
  GBZ(const GBZ& source);
  GBZ(GBZ&& source);
  ~GBZ();

  // Build GBZ from the structures returned by `gfa_to_gbwt()`.
  // Resets the pointers to `nullptr`.
  GBZ(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<SequenceSource>& source);

  // Build GBZ from a GBWT index and a `HandleGraph`.
  // Resets the GBWT pointer to `nullptr`.
  GBZ(std::unique_ptr<gbwt::GBWT>& index, const HandleGraph& source);

  // Build GBZ from a GBWT index and a sequence source.
  // Note that the GBZ will store a copy of the GBWT index.
  GBZ(const gbwt::GBWT& index, const SequenceSource& source);

  // Build GBZ from a GBWT index and a `HandleGraph`.
  // Note that the GBZ will store a copy of the GBWT index.
  GBZ(const gbwt::GBWT& index, const HandleGraph& source);

  void swap(GBZ& another);
  GBZ& operator=(const GBZ& source);
  GBZ& operator=(GBZ&& source);

//------------------------------------------------------------------------------

  /*
    Reference samples and named paths. Reference samples cannot be changed
    conveniently within GBWTGraph, because the pointer to the GBWT index is const.
  */

  // Sets the given sample names as reference samples, but only if they are
  // present in the GBWT metadata. Returns the number of reference samples.
  // This is somewhat expensive, as the GBWTGraph must recache named paths.
  size_t set_reference_samples(const std::unordered_set<std::string>& samples);

  // Returns the set of reference samples.
  // Some of these samples may not exist in the GBWT metadata.
  const std::unordered_set<std::string>& get_reference_samples() const
  {
    return this->graph.reference_samples;
  }

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
  GBWTGraph  graph;

  const static std::string EXTENSION; // ".gbz"

  // Serialize the the GBZ into the output stream in the simple-sds format.
  void simple_sds_serialize(std::ostream& out) const;

  // Serialize the given GBWT and GBWTGraph objects in the GBZ format.
  static void simple_sds_serialize(const gbwt::GBWT& index, const GBWTGraph& graph, std::ostream& out);

  // Deserialize or decompress the GBZ from the input stream.
  void simple_sds_load(std::istream& in);

  // Returns the size of the serialized structure in elements.
  size_t simple_sds_size() const;

  // Serialize the GBWT (simple-sds format) and the GBWTGraph to separate files.
  // Default graph format is libhandlegraph / SDSL.
  void serialize_to_files(const std::string& gbwt_name, const std::string& graph_name, bool simple_sds_graph = false) const;

  // Loads the GBWT (simple-sds format) and the GBWTGraph from separate files.
  // Graph format is libhandlegraph / SDSL; the simple-sds format cannot be read.
  void load_from_files(const std::string& gbwt_name, const std::string& graph_name);

private:
  void copy(const GBZ& source);
  void reset_tags();
  void add_source();
  void set_gbwt();
  void set_gbwt_address();
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_GBZ_H
