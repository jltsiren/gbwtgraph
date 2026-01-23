#ifndef GBWTGRAPH_INTERNAL_H
#define GBWTGRAPH_INTERNAL_H

#include <gbwtgraph/gbz.h>

#include <iostream>
#include <string>
#include <tuple>
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

  void write(std::string_view view);
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

  void write(std::string_view view) { this->buffer.insert(this->buffer.end(), view.begin(), view.end()); }
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
write_gfa_segment(Writer& writer, nid_t node_id, std::string_view sequence)
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

struct GFAGrammarIterator;

/*
  A collection of grammar rules from GFA Q-lines.
  This class provides iterators over the expansions of any rule, both in
  forward and reverse complement orientations.

  TODO: Should we store the expansion length for each rule?
*/
class GFAGrammar
{
public:
  // The expansion of a rule is a list of (rule / segment name, is_reverse) pairs.
  using expansion_type = std::vector<std::pair<std::string, bool>>;

  // A rule is a mapping from a rule name to its expansion.
  using rule_type = std::unordered_map<std::string, expansion_type>::const_iterator;

  size_t size() const { return this->rules.size(); }
  bool empty() const { return this->rules.empty(); }

  // Inserts a new rule. Returns true if the insertion was successful
  // and false if the rule already existed.
  // No attempt is made to avoid cycles or empty expansions.
  bool insert(std::string&& rule_name, expansion_type&& expansion)
  {
    return this->rules.emplace(rule_name, expansion).second;
  }

  // Returns an iterator over the expansion of the given rule in the given direction.
  GFAGrammarIterator iter(const std::string& rule_name, bool is_reverse) const;

  // Returns an iterator to the expansion of the given rule, or `no_rule()` if the rule does not exist.
  rule_type expand(const std::string& rule_name) const { return this->rules.find(rule_name); }

  // Returns an iterator that indicates that the rule does not exist.
  rule_type no_rule() const { return this->rules.end(); }

  // Returns the first segment name in the expansion of the given symbol in forward orientation.
  // Returns the symbol itself if it is not a rule.
  // This can be used for assigning a rule to a graph component.
  std::string first_segment(const std::string& symbol) const;

  // Validates the grammar:
  // * there are no cycles;
  // * names do not clash with segment names (in the given source);
  // * all names on the right side of productions exist as segments or rules; and
  // * all expansions are nontrivial (length >= 2).
  // Throws `std::runtime_error` if the grammar is invalid.
  void validate(const SequenceSource& source) const;

private:
  std::unordered_map<std::string, expansion_type> rules;
};

struct GFAGrammarIterator
{
  explicit GFAGrammarIterator(const GFAGrammar& grammar, const std::string& rule_name, bool is_reverse);

  // Returns the next segment in the expansion and its orientation,
  // or an empty string if the expansion is complete.
  std::pair<std::string_view, bool> next();

  // Returns `true` if the iterator is empty.
  // An empty iterator either corresponds to a nonexistent rule or has already
  // returned an empty string.
  bool empty() const { return this->stack.empty(); }

  const GFAGrammar& grammar;

  // Stack of (grammar iterator, is_reverse, number of symbols processed).
  std::vector<std::tuple<GFAGrammar::rule_type, bool, size_t>> stack;
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
std::vector<std::pair<size_t, gbwt::edge_type>>
sample_path_positions(const GBZ& gbz, path_handle_t path, size_t sample_interval, size_t* length = nullptr);

//------------------------------------------------------------------------------

// A hasher for gbwt::PathName.
struct PathNameHasher
{
  static size_t hash(gbwt::PathName::path_name_type v)
  {
    return wang_hash_64(v);
  }

  static size_t combine(size_t h, size_t v)
  {
    h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }

  size_t operator() (const gbwt::PathName& key) const
  {
    size_t h = hash(key.sample);
    h = combine(h, hash(key.contig));
    h = combine(h, hash(key.phase));
    h = combine(h, hash(key.count));
    return h;
  }
};

// A hasher for GBWT search states.
struct SearchStateHasher
{
  static size_t hash(gbwt::size_type v)
  {
    return wang_hash_64(v);
  }

  static size_t combine(size_t h, size_t v)
  {
    h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }

  size_t operator() (const gbwt::SearchState& state) const
  {
    size_t h = hash(state.node);
    h = combine(h, hash(state.range.first));
    h = combine(h, hash(state.range.second));
    return h;
  }
};

//------------------------------------------------------------------------------

/*
  Mapping from path names to integers in [0, 64). Depending on the GBWT metadata,
  the mapping may be based on (sample, contig, haplotype), (sample, haplotype),
  or (sample). If there are too many samples, all paths map to 0.

  Path names that map to the same integer are considered parts of the same
  haplotype. We can then use a 64-bit integer for storing subsets of haplotypes.
  The mapping depends on lexicographic order of path names and is therefore
  independent of the order of the paths in the GBWT index.
*/
struct PathIdMap
{
  constexpr static size_t MAX_HAPLOTYPES = 64;

  enum class KeyType { NONE, SAMPLE, SAMPLE_HAPLOTYPE, SAMPLE_CONTIG_HAPLOTYPE };

  explicit PathIdMap(const gbwt::Metadata& metadata);

  // Returns the number of distinct haplotypes in the map.
  // Mapped values are in [0, size()).
  size_t size() const { return this->path_to_id.size(); }

  // Returns the haplotype identifier for the given path name.
  // Returns 0 if the identifier cannot be determined.
  size_t id(const gbwt::PathName& path) const
  {
    gbwt::PathName key = this->mask_name(path);
    auto iter = this->path_to_id.find(key);
    if(iter == this->path_to_id.end()) { return 0; }
    return iter->second;
  }

  // Returns the type of keys used in the mapping.
  KeyType key_type() const { return this->key; }

  // Returns a string representation of the given key type.
  static std::string key_type_str(KeyType key);

  std::unordered_map<gbwt::PathName, size_t, PathNameHasher> path_to_id;
  gbwt::PathName mask;
  KeyType key;

private:
  gbwt::PathName mask_name(const gbwt::PathName& name) const
  {
    return gbwt::PathName
    (
      name.sample & this->mask.sample,
      name.contig & this->mask.contig,
      name.phase & this->mask.phase,
      name.count & this->mask.count
    );
  }

  bool build_map(const std::vector<gbwt::PathName>& sorted_paths);
};

//------------------------------------------------------------------------------

/*
  Extracts the subpath corresponding to the given kmer.
  Assumes that the arguments are valid.

  Offset node_offset of node path[path_offset] is like a position returned by
  KmerIndex. If is_reverse is false, the kmer is on the forward strand, and the
  path extends forward from the position. If it is true, the kmer is on the
  reverse strand, and the path extends backward from the position. In the
  latter case, the returned path is a subpath of the reverse complement of the
  given path.
*/
std::vector<gbwt::node_type>
extract_kmer_path(const GBWTGraph& graph, const std::vector<handle_t>& path, size_t path_offset, size_t node_offset, size_t k, bool is_reverse);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_INTERNAL_H
