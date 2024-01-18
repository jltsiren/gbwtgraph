#ifndef GBWTGRAPH_INTERNAL_H
#define GBWTGRAPH_INTERNAL_H

#include <gbwt/gbwt.h>
#include <gbwt/metadata.h>
#include <gbwtgraph/utils.h>

#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
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
  A buffered TSV writer for single-threaded situtations.
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
  FIXME: This should be made public.

  A structure for building GBWT metadata.

  Constructor and the methods for handling paths/walks throw `std::runtime_error`
  on failure.
*/
struct MetadataBuilder
{

  constexpr static size_t NO_FIELD = std::numeric_limits<size_t>::max();
  struct PathMetadataBuilder
  {
    std::regex parser;

    // Mapping from regex submatches to GBWT path name components.
    size_t sample_field, contig_field, haplotype_field, fragment_field;

    PathSense sense;

    PathMetadataBuilder(const std::string& path_name_regex, const std::string& path_name_prefix, PathSense path_sense);
  };

  std::vector<PathMetadataBuilder> path_name_formats;

  // GBWT metadata.
  std::map<std::string, size_t> sample_names, contig_names;
  std::set<std::pair<size_t, size_t>> haplotypes; // (sample id, phase id)
  std::vector<std::vector<gbwt::PathName>> path_names;
  std::map<gbwt::PathName, size_t> counts;

  bool ref_path_sample_warning;

  // Construct a MetadataBuilder with no path name formats.
  MetadataBuilder();
  
  // Construct a MetadataBuilder from pre-existing metadata and no path name formats.
  MetadataBuilder(const gbwt::Metadata& source);

  // Construct a MetadataBuilder with one path name format.
  MetadataBuilder(const std::string& path_name_regex, const std::string& path_name_prefix, PathSense path_sense = PathSense::GENERIC);

  // Register a format for parsing path names. Formats are tried in order until one matches.
  void add_path_name_format(const std::string& path_name_regex, const std::string& path_name_prefix, PathSense path_sense);
  
  // Add a path defined by libhandlegraph metadata to the given job.
  // Doesn't create metadata for samples or contigs if the no-name sentinel is
  // passed for a sense that usually has them.
  void add_path(PathSense sense, const std::string& sample_name, const std::string& locus_name, size_t haplotype, size_t phase_block, const handlegraph::subrange_t& subrange, size_t job = 0);
  
  // Parse a path name using a regex to determine metadata, and assign it to the given job.
  void add_path(const std::string& name, size_t job = 0);

  // Add a path based on walk metadata and assign it to the given job.
  void add_walk(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start, size_t job = 0);

  // Add a haplotype path and assign it to the given job.
  void add_haplotype(const std::string& sample, const std::string& contig, size_t haplotype, size_t fragment, size_t job = 0);

  // Add a named path as a generic named path and assign it to the given job.
  void add_generic_path(const std::string& name, size_t job = 0);
  
  bool empty() const { return this->path_names.empty(); }

  // Build GBWT metadata from the current contents.
  // Paths come our ordered by job.
  gbwt::Metadata get_metadata() const;

  void clear()
  {
    this->sample_names = std::map<std::string, size_t>();
    this->contig_names = std::map<std::string, size_t>();
    this->haplotypes = std::set<std::pair<size_t, size_t>>();
    this->path_names = std::vector<std::vector<gbwt::PathName>>();
    this->counts = std::map<gbwt::PathName, size_t>();
  }

  static std::vector<std::string> map_to_vector(const std::map<std::string, size_t>& source)
  {
    std::vector<std::string> result(source.size());
    for(auto& name : source) { result[name.second] = name.first; }
    return result;
  }

private:
  void add_path_name(const gbwt::PathName& path_name, size_t job)
  {
    if(job >= this->path_names.size())
    {
      this->path_names.resize(job + 1);
    }
    this->path_names[job].push_back(path_name);
  }
};

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

} // namespace gbwtgraph

#endif // GBWTGRAPH_INTERNAL_H
