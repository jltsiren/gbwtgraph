#include <gbwtgraph/gbwtgraph.h>

#include <algorithm>
#include <cassert>
#include <queue>
#include <stack>
#include <string>

#include <arpa/inet.h>

#include <omp.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t GBWTGraph::CHUNK_SIZE;

constexpr std::uint32_t GBWTGraph::Header::TAG;
constexpr std::uint32_t GBWTGraph::Header::VERSION;

constexpr std::uint64_t GBWTGraph::Header::FLAG_MASK;
constexpr std::uint64_t GBWTGraph::Header::FLAG_TRANSLATION;
constexpr std::uint64_t GBWTGraph::Header::FLAG_COMPRESSED;

constexpr std::uint32_t GBWTGraph::Header::OLD_VERSION;
constexpr std::uint64_t GBWTGraph::Header::OLD_FLAG_MASK;

//------------------------------------------------------------------------------

// Other class variables.

const std::string GBWTGraph::EXTENSION = ".gg";
const std::string GBWTGraph::COMPRESSED_EXTENSION = ".gbz";

//------------------------------------------------------------------------------

GBWTGraph::Header::Header() :
  tag(TAG), version(VERSION),
  nodes(0),
  flags(0)
{
}

bool
GBWTGraph::Header::check() const
{
  if(this->tag != TAG) { return false; }
  switch(this->version)
  {
  case VERSION:
    return ((this->flags & FLAG_MASK) == this->flags);
  case OLD_VERSION:
    return ((this->flags & OLD_FLAG_MASK) == this->flags);
  default:
    return false;
  }
}

bool
GBWTGraph::Header::operator==(const Header& another) const
{
  return (this->tag == another.tag && this->version == another.version &&
          this->nodes == another.nodes &&
          this->flags == another.flags);
}

//------------------------------------------------------------------------------

GBWTGraph::GBWTGraph() :
  index(nullptr), header()
{
}

GBWTGraph::GBWTGraph(const GBWTGraph& source)
{
  this->copy(source);
}

GBWTGraph::GBWTGraph(GBWTGraph&& source)
{
  *this = std::move(source);
}

GBWTGraph::~GBWTGraph()
{
}

void
GBWTGraph::swap(GBWTGraph& another)
{
  if(&another == this) { return; }

  std::swap(this->index, another.index);
  std::swap(this->header, another.header);
  this->sequences.swap(another.sequences);
  this->real_nodes.swap(another.real_nodes);
  this->segments.swap(another.segments);
  this->node_to_segment.swap(another.node_to_segment);
}

GBWTGraph&
GBWTGraph::operator=(const GBWTGraph& source)
{
  if(&source != this) { this->copy(source); }
  return *this;
}

GBWTGraph&
GBWTGraph::operator=(GBWTGraph&& source)
{
  if(&source != this)
  {
    this->index = std::move(source.index);
    this->header = std::move(source.header);
    this->sequences = std::move(source.sequences);
    this->real_nodes = std::move(source.real_nodes);
    this->segments = std::move(source.segments);
    this->node_to_segment = std::move(source.node_to_segment);
  }
  return *this;
}

void
GBWTGraph::copy(const GBWTGraph& source)
{
  this->index = source.index;
  this->header = source.header;
  this->sequences = source.sequences;
  this->real_nodes = source.real_nodes;
  this->segments = source.segments;
  this->node_to_segment = source.node_to_segment;
}

//------------------------------------------------------------------------------

GBWTGraph::GBWTGraph(const gbwt::GBWT& gbwt_index, const HandleGraph& sequence_source) :
  index(nullptr)
{
  // Set GBWT and do sanity checks.
  this->set_gbwt(gbwt_index);
  if(this->index->empty()) { return; }

  // Build real_nodes to support has_node().
  this->determine_real_nodes();

  // Store the sequences.
  this->sequences = StringArray(this->index->sigma() - this->index->firstNode(),
  [&](size_t offset) -> size_t
  {
    gbwt::node_type node = offset + this->index->firstNode();
    nid_t id = gbwt::Node::id(node);
    if(!(this->has_node(id))) { return 0; }
    handle_t handle = sequence_source.get_handle(id, gbwt::Node::is_reverse(node));
    return sequence_source.get_length(handle);
  },
  [&](size_t offset) -> std::string
  {
    gbwt::node_type node = offset + this->index->firstNode();
    nid_t id = gbwt::Node::id(node);
    if(!(this->has_node(id))) { return std::string(); }
    handle_t handle = sequence_source.get_handle(id, gbwt::Node::is_reverse(node));
    return sequence_source.get_sequence(handle);
  });
}

GBWTGraph::GBWTGraph(const gbwt::GBWT& gbwt_index, const SequenceSource& sequence_source) :
  index(nullptr)
{
  // Set GBWT and do sanity checks.
  this->set_gbwt(gbwt_index);
  if(this->index->empty()) { return; }

  // Build real_nodes to support has_node().
  this->determine_real_nodes();

  // Store the sequences.
  this->sequences = StringArray(this->index->sigma() - this->index->firstNode(),
  [&](size_t offset) -> size_t
  {
    gbwt::node_type node = offset + this->index->firstNode();
    nid_t id = gbwt::Node::id(node);
    if(!(this->has_node(id))) { return 0; }
    return sequence_source.get_length(id);
  },
  [&](size_t offset) -> std::string
  {
    gbwt::node_type node = offset + this->index->firstNode();
    nid_t id = gbwt::Node::id(node);
    if(!(this->has_node(id))) { return std::string(); }
    std::string result = sequence_source.get_sequence(id);
    if(gbwt::Node::is_reverse(node)) { reverse_complement_in_place(result); }
    return result;
  });

  // Store the node to segment translation.
  if(sequence_source.uses_translation())
  {
    this->header.set(Header::FLAG_TRANSLATION);
    std::tie(this->segments, this->node_to_segment) = sequence_source.invert_translation();
  }
}

//------------------------------------------------------------------------------

void
GBWTGraph::determine_real_nodes()
{
  // Sometimes (e.g. in `decompress()`) we call this function after the header already
  // contains the correct number of nodes.
  this->header.nodes = 0;
  if(this->index->empty())
  {
    this->real_nodes = sdsl::bit_vector();
    return;
  }

  size_t potential_nodes = this->index->sigma() - this->index->firstNode();
  this->real_nodes = sdsl::bit_vector(potential_nodes / 2, 0);
  for(gbwt::node_type node = this->index->firstNode(); node < this->index->sigma(); node += 2)
  {
    if(!(this->index->empty(node)))
    {
      this->real_nodes[this->node_offset(node) / 2] = 1;
      this->header.nodes++;
    }
  }
}

//------------------------------------------------------------------------------

bool
GBWTGraph::has_node(nid_t node_id) const
{
  size_t offset = this->node_offset(gbwt::Node::encode(node_id, false)) / 2;
  return (offset < this->real_nodes.size() && this->real_nodes[offset]);
}

handle_t
GBWTGraph::get_handle(const nid_t& node_id, bool is_reverse) const
{
  return node_to_handle(gbwt::Node::encode(node_id, is_reverse));
}

nid_t
GBWTGraph::get_id(const handle_t& handle) const
{
  return gbwt::Node::id(handle_to_node(handle));
}

bool
GBWTGraph::get_is_reverse(const handle_t& handle) const
{
  return gbwt::Node::is_reverse(handle_to_node(handle));
}

handle_t
GBWTGraph::flip(const handle_t& handle) const
{
  return node_to_handle(gbwt::Node::reverse(handle_to_node(handle)));
}

size_t
GBWTGraph::get_length(const handle_t& handle) const
{
  size_t offset = this->node_offset(handle);
  return this->sequences.length(offset);
}

std::string
GBWTGraph::get_sequence(const handle_t& handle) const
{
  size_t offset = this->node_offset(handle);
  return this->sequences.str(offset);
}

char
GBWTGraph::get_base(const handle_t& handle, size_t index) const
{
  size_t offset = this->node_offset(handle);
  view_type view = this->sequences.view(offset);
  return *(view.first + index);
}

std::string
GBWTGraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const
{
  size_t offset = this->node_offset(handle);
  view_type view = this->sequences.view(offset);
  index = std::min(index, view.second);
  size = std::min(size, view.second - index);
  return std::string(view.first + index, view.first + index + size);
}

size_t
GBWTGraph::get_node_count() const
{
  return this->header.nodes;
}

nid_t
GBWTGraph::min_node_id() const
{
  return gbwt::Node::id(this->index->firstNode());
}

nid_t
GBWTGraph::max_node_id() const
{
  nid_t next_id = gbwt::Node::id(this->index->sigma());
  return next_id - 1;
}

bool
GBWTGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const
{
  return this->cached_follow_edges(this->get_single_cache(), handle, go_left, iteratee);
}

bool
GBWTGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const
{
  if(parallel)
  {
    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for(gbwt::node_type node = this->index->firstNode(); node < this->index->sigma(); node += 2)
    {
      if(!(this->real_nodes[this->node_offset(node) / 2])) { continue; }
      if(!iteratee(node_to_handle(node)))
      {
        // We should stop early but it's not worth the effort.
      }
    }
  }
  else
  {
    for(gbwt::node_type node = this->index->firstNode(); node < this->index->sigma(); node += 2)
    {
      if(!(this->real_nodes[this->node_offset(node) / 2])) { continue; }
      if(!iteratee(node_to_handle(node))) { return false; }
    }
  }

  return true;
}

size_t
GBWTGraph::get_degree(const handle_t& handle, bool go_left) const
{
  // Cache the node.
  gbwt::node_type curr = handle_to_node(handle);
  if(go_left) { curr = gbwt::Node::reverse(curr); }
  gbwt::CachedGBWT cache = this->get_single_cache();
  gbwt::size_type cache_index = cache.findRecord(curr);

  // The outdegree reported by GBWT might account for the endmarker, which is
  // always the first successor.
  size_t result = cache.outdegree(cache_index);
  if(result > 0 && cache.successor(cache_index, 0) == gbwt::ENDMARKER) { result--; }
  return result;
}

bool
GBWTGraph::has_edge(const handle_t& left, const handle_t& right) const
{
  // Cache the node.
  gbwt::node_type curr = handle_to_node(left);
  gbwt::CachedGBWT cache = this->get_single_cache();
  gbwt::size_type cache_index = cache.findRecord(curr);

  for(gbwt::rank_type outrank = 0; outrank < cache.outdegree(cache_index); outrank++)
  {
    gbwt::node_type next = cache.successor(cache_index, outrank);
    if(node_to_handle(next) == right) { return true; }
  }

  return false;
}

//------------------------------------------------------------------------------

bool
GBWTGraph::has_segment_names() const
{
  return this->header.get(Header::FLAG_TRANSLATION);
}

std::pair<std::string, size_t>
GBWTGraph::get_segment_name_and_offset(const handle_t& handle) const
{
  // If there is no translation, the predecessor is always at the end.
  nid_t id = this->get_id(handle);
  auto iter = this->node_to_segment.predecessor(id);
  if(iter == this->node_to_segment.one_end())
  {
    return std::pair<std::string, size_t>(std::to_string(id), 0);
  }

  // Determine the total length of nodes in this segment that precede `id`
  // in the given orientation.
  size_t start = 0, limit = 0;
  if(this->get_is_reverse(handle))
  {
    auto successor = iter; ++successor;
    start = this->node_offset(gbwt::Node::encode(id + 1, false));
    limit = this->node_offset(gbwt::Node::encode(successor->second, false));
  }
  else
  {
    start = this->node_offset(gbwt::Node::encode(iter->second, false));
    limit = this->node_offset(gbwt::Node::encode(id, false));
  }
  size_t offset = this->sequences.length(start, limit) / 2;

  return std::pair<std::string, size_t>(this->segments.str(iter->first), offset);
}

std::string
GBWTGraph::get_segment_name(const handle_t& handle) const
{
  // If there is no translation, the predecessor is always at the end.
  nid_t id = this->get_id(handle);
  auto iter = this->node_to_segment.predecessor(id);
  if(iter == this->node_to_segment.one_end()) { return std::to_string(id); }
  return this->segments.str(iter->first);
}

size_t
GBWTGraph::get_segment_offset(const handle_t& handle) const
{
  // If there is no translation, the predecessor is always at the end.
  nid_t id = this->get_id(handle);
  auto iter = this->node_to_segment.predecessor(id);
  if(iter == this->node_to_segment.one_end()) { return 0; }

  // Determine the total length of nodes in this segment that precede `id`
  // in the given orientation.
  size_t start = 0, limit = 0;
  if(this->get_is_reverse(handle))
  {
    auto successor = iter; ++successor;
    start = this->node_offset(gbwt::Node::encode(id + 1, false));
    limit = this->node_offset(gbwt::Node::encode(successor->second, false));
  }
  else
  {
    start = this->node_offset(gbwt::Node::encode(iter->second, false));
    limit = this->node_offset(gbwt::Node::encode(id, false));
  }
  size_t offset = this->sequences.length(start, limit) / 2;

  return offset;
}

//------------------------------------------------------------------------------

uint32_t
GBWTGraph::get_magic_number() const {
    // Specify what it should look like on the wire
    const char* bytes = "GBG ";
    // Convert to a host byte order number
    return ntohl(*((const uint32_t*) bytes));
}

void
GBWTGraph::serialize_members(std::ostream& out) const
{
  out.write(reinterpret_cast<const char*>(&(this->header)), sizeof(Header));
  this->sequences.serialize(out);
  this->real_nodes.serialize(out);
  if(this->header.get(Header::FLAG_TRANSLATION))
  {
    this->segments.serialize(out);
    this->node_to_segment.serialize(out);
  }
}

void
GBWTGraph::deserialize_members(std::istream& in)
{
  // Load the header and update to the current version.
  in.read(reinterpret_cast<char*>(&(this->header)), sizeof(Header));
  if(!(this->header.check()))
  {
    std::cerr << "GBWTGraph::deserialize_members(): Invalid or old graph file" << std::endl;
    std::cerr << "GBWTGraph::deserialize_members(): Graph version is " << this->header.version << "; expected " << Header::VERSION << std::endl;
    return;
  }
  if(this->header.get(Header::FLAG_COMPRESSED))
  {
    std::cerr << "GBWTGraph::deserialize_members(): Cannot deserialize compressed file format" << std::endl;
    std::cerr << "GBWTGraph::deserialize_members(): GBWT index must be specified for decompression" << std::endl;
    return;
  }
  this->header.version = Header::VERSION;

  // Load the graph.
  this->sequences.deserialize(in);
  this->real_nodes.load(in);

  // Load the translation.
  if(this->header.get(Header::FLAG_TRANSLATION))
  {
    this->segments.deserialize(in);
    this->node_to_segment.load(in);
  }
}

void
GBWTGraph::set_gbwt(const gbwt::GBWT& gbwt_index)
{
  this->index = &gbwt_index;

  // Sanity checks for the GBWT index.
  assert(this->index->bidirectional());
}

//------------------------------------------------------------------------------

void
GBWTGraph::compress(std::ostream& out) const
{
  // Serialize the header.
  Header copy = this->header;
  copy.set(Header::FLAG_COMPRESSED);
  out.write(reinterpret_cast<const char*>(&copy), sizeof(Header));

  // Compress the sequences. `real_nodes` can be rebuilt from the GBWT.
  {
    StringArray forward_only(this->sequences.size() / 2,
    [&](size_t offset) -> size_t
    {
      return this->sequences.length(2 * offset);
    },
    [&](size_t offset) -> view_type
    {
      return this->sequences.view(2 * offset);
    });
    forward_only.compress(out);
  }

  // Compress the translation.
  if(this->header.get(Header::FLAG_TRANSLATION))
  {
    this->segments.compress(out);
    this->node_to_segment.serialize(out);
  }
}

bool
GBWTGraph::decompress(std::istream& in, const gbwt::GBWT& gbwt_index)
{
  // Set the GBWT so we can rebuild `real_nodes` later.
  this->set_gbwt(gbwt_index);

  // Load the header and update to the current version.
  in.read(reinterpret_cast<char*>(&(this->header)), sizeof(Header));
  if(!(this->header.check()))
  {
    std::cerr << "GBWTGraph::decompress(): Invalid or old graph file" << std::endl;
    std::cerr << "GBWTGraph::decompress(): Graph version is " << this->header.version << "; expected " << Header::VERSION << std::endl;
    return false;
  }
  this->header.version = Header::VERSION;
  bool compressed = this->header.get(Header::FLAG_COMPRESSED);
  this->header.unset(Header::FLAG_COMPRESSED);

  // Load the graph.
  if(compressed)
  {
    {
      StringArray forward_only;
      forward_only.decompress(in);
      this->sequences = StringArray(2 * forward_only.size(),
      [&](size_t offset) -> size_t
      {
        return forward_only.length(offset / 2);
      },
      [&](size_t offset) -> std::string
      {
        std::string result = forward_only.str(offset / 2);
        if(offset & 1) { reverse_complement_in_place(result); }
        return result;
      });
    }
    this->determine_real_nodes();
  }
  else
  {
    this->sequences.deserialize(in);
    this->real_nodes.load(in);
  }

  // Load the translation.
  if(this->header.get(Header::FLAG_TRANSLATION))
  {
    if(compressed)
    {
      this->segments.decompress(in);
      this->node_to_segment.load(in);
    }
    else
    {
      this->segments.deserialize(in);
      this->node_to_segment.load(in);
    }
  }

  return true;
}

//------------------------------------------------------------------------------

view_type
GBWTGraph::get_sequence_view(const handle_t& handle) const
{
  size_t offset = this->node_offset(handle);
  return this->sequences.view(offset);
}

bool
GBWTGraph::starts_with(const handle_t& handle, char c) const
{
  size_t offset = this->node_offset(handle);
  view_type view = this->sequences.view(offset);
  return (view.second > 0 && *view.first == c);
}

bool
GBWTGraph::ends_with(const handle_t& handle, char c) const
{
  size_t offset = this->node_offset(handle);
  view_type view = this->sequences.view(offset);
  return (view.second > 0 && *(view.first + (view.second - 1)) == c);
}

gbwt::SearchState
GBWTGraph::find(const std::vector<handle_t>& path) const
{
  if(path.empty()) { return gbwt::SearchState(); }
  gbwt::SearchState result = this->get_state(path[0]);
  for(size_t i = 1; i < path.size() && !result.empty(); i++)
  {
    result = this->index->extend(result, handle_to_node(path[i]));
  }
  return result;
}

gbwt::BidirectionalState
GBWTGraph::bd_find(const std::vector<handle_t>& path) const
{
  if(path.empty()) { return gbwt::BidirectionalState(); }
  gbwt::BidirectionalState result = this->get_bd_state(path[0]);
  for(size_t i = 1; i < path.size() && !result.empty(); i++)
  {
    result = this->index->bdExtendForward(result, handle_to_node(path[i]));
  }
  return result;
}

//------------------------------------------------------------------------------

gbwt::SearchState
GBWTGraph::find(const gbwt::CachedGBWT& cache, const std::vector<handle_t>& path) const
{
  if(path.empty()) { return gbwt::SearchState(); }
  gbwt::SearchState result = this->get_state(cache, path[0]);
  for(size_t i = 1; i < path.size() && !result.empty(); i++)
  {
    result = cache.extend(result, handle_to_node(path[i]));
  }
  return result;
}

gbwt::BidirectionalState
GBWTGraph::bd_find(const gbwt::CachedGBWT& cache, const std::vector<handle_t>& path) const
{
  if(path.empty()) { return gbwt::BidirectionalState(); }
  gbwt::BidirectionalState result = this->get_bd_state(cache, path[0]);
  for(size_t i = 1; i < path.size() && !result.empty(); i++)
  {
    result = cache.bdExtendForward(result, handle_to_node(path[i]));
  }
  return result;
}

bool
GBWTGraph::follow_paths(const gbwt::CachedGBWT& cache, gbwt::SearchState state,
                        const std::function<bool(const gbwt::SearchState&)>& iteratee) const
{
  gbwt::size_type cache_index = cache.findRecord(state.node);
  for(gbwt::rank_type outrank = 0; outrank < cache.outdegree(cache_index); outrank++)
  {
    if(cache.successor(cache_index, outrank) == gbwt::ENDMARKER) { continue; }
    gbwt::SearchState next_state = cache.cachedExtend(state, cache_index, outrank);
    if(next_state.empty()) { continue; }
    if(!iteratee(next_state)) { return false; }
  }

  return true;
}

bool
GBWTGraph::follow_paths(const gbwt::CachedGBWT& cache, gbwt::BidirectionalState state, bool backward,
                        const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const
{
  gbwt::size_type cache_index = cache.findRecord(backward ? state.backward.node : state.forward.node);
  for(gbwt::rank_type outrank = 0; outrank < cache.outdegree(cache_index); outrank++)
  {
    if(cache.successor(cache_index, outrank) == gbwt::ENDMARKER) { continue; }
    gbwt::BidirectionalState next_state = (backward ? cache.cachedExtendBackward(state, cache_index, outrank) :  cache.cachedExtendForward(state, cache_index, outrank));
    if(next_state.empty()) { continue; }
    if(!iteratee(next_state)) { return false; }
  }

  return true;
}

bool
GBWTGraph::cached_follow_edges(const gbwt::CachedGBWT& cache, const handle_t& handle, bool go_left,
                               const std::function<bool(const handle_t&)>& iteratee) const
{
  // Incoming edges correspond to the outgoing edges of the reverse node.
  gbwt::node_type curr = handle_to_node(handle);
  if(go_left) { curr = gbwt::Node::reverse(curr); }

  // Cache the node.
  gbwt::size_type cache_index = cache.findRecord(curr);

  for(gbwt::rank_type outrank = 0; outrank < cache.outdegree(cache_index); outrank++)
  {
    gbwt::node_type next = cache.successor(cache_index, outrank);
    if(next == gbwt::ENDMARKER) { continue; }

    // If we started from the reverse node, we must reverse the successor nodes to get
    // the predecessor nodes of the original node.
    if(go_left) { next = gbwt::Node::reverse(next); }
    if(!iteratee(node_to_handle(next))) { return false; }
  }

  return true;
}

//------------------------------------------------------------------------------

/*
  Haplotype-consistent traversal of the graph and the corresponding sequence.
  The traversal always starts at the beginning of a node, but it may end in the
  middle of a node.
*/

struct GBWTTraversal
{
  std::vector<handle_t> traversal;
  size_t length;
  gbwt::SearchState state; // GBWT search state at the end of the traversal.

  std::string get_sequence(const HandleGraph& graph) const
  {
    std::string result;
    result.reserve(this->length);
    for(handle_t handle : this->traversal)
    {
      result.append(graph.get_sequence(handle), 0, this->length - result.length());
    }
    return result;
  }

  std::string get_sequence(const GBWTGraph& graph) const
  {
    std::string result;
    result.reserve(this->length);
    for(handle_t handle : this->traversal)
    {
      auto view = graph.get_sequence_view(handle);
      result.append(view.first, std::min(view.second, this->length - result.length()));
    }
    return result;
  }
};

void
for_each_haplotype_window(const GBWTGraph& graph, size_t window_size,
                          const std::function<void(const std::vector<handle_t>&, const std::string&)>& lambda,
                          bool parallel)
{
  // Traverse all starting nodes in parallel.
  graph.for_each_handle([&](const handle_t& h) -> bool
  {
    // Get a GBWT cache.
    gbwt::CachedGBWT cache = graph.get_cache();

    // Initialize the stack with both orientations.
    std::stack<GBWTTraversal> windows;
    size_t node_length = graph.get_length(h);
    for(bool is_reverse : { false, true })
    {
      handle_t handle = (is_reverse ? graph.flip(h) : h);
      gbwt::SearchState state = graph.get_state(cache, handle);
      if(state.empty()) { continue; }
      GBWTTraversal window { { handle }, node_length, state };
      windows.push(window);
    }

    // Extend the windows.
    size_t target_length = node_length + window_size - 1;
    while(!windows.empty())
    {
      GBWTTraversal window = windows.top(); windows.pop();
      // Report the full window.
      if(window.length >= target_length)
      {
        lambda(window.traversal, window.get_sequence(graph));
        continue;
      }

      // Try to extend the window to all successor nodes.
      bool extend_success = false;
      graph.follow_paths(cache, window.state, [&](const gbwt::SearchState& next_state) -> bool
      {
        handle_t next_handle = GBWTGraph::node_to_handle(next_state.node);
        GBWTTraversal next_window = window;
        next_window.traversal.push_back(next_handle);
        next_window.length += std::min(graph.get_length(next_handle), target_length - window.length);
        next_window.state = next_state;
        windows.push(next_window);
        extend_success = true;
        return true;
      });

      // Report sufficiently long kmers that cannot be extended.
      if(!extend_success && window.length >= window_size)
      {
        lambda(window.traversal, window.get_sequence(graph));
      }
    }

    return true;
  }, parallel);
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
