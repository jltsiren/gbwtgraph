#include <gbwtgraph/gbwtgraph.h>

#include <algorithm>
#include <cassert>
#include <queue>
#include <stack>

#include <omp.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t GBWTGraph::CHUNK_SIZE;
constexpr size_t GBWTGraph::BLOCK_SIZE;

constexpr std::uint32_t GBWTGraph::Header::TAG;
constexpr std::uint32_t GBWTGraph::Header::VERSION;

//------------------------------------------------------------------------------

// Other class variables.

const std::string GBWTGraph::EXTENSION = ".gg";

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
  return (this->tag == TAG && this->version == VERSION && this->flags == 0);
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
  this->offsets.swap(another.offsets);
  this->real_nodes.swap(another.real_nodes);
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
    this->offsets = std::move(source.offsets);
    this->real_nodes = std::move(source.real_nodes);
  }
  return *this;
}

void
GBWTGraph::determine_real_nodes()
{
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

std::vector<handle_t>
GBWTGraph::cache_handles(const std::function<handle_t(gbwt::node_type)>& get_source_handle) const
{
  size_t potential_nodes = this->index->sigma() - this->index->firstNode();
  std::vector<handle_t> handle_cache(potential_nodes / 2);
  for(gbwt::node_type node = this->index->firstNode(); node < this->index->sigma(); node += 2)
  {
    handle_cache[this->node_offset(node) / 2] = get_source_handle(node);
  }
  return handle_cache;
}

void
GBWTGraph::allocate_arrays(const std::function<size_t(size_t)>& get_source_length)
{
  size_t total_length = 0;
  for(gbwt::node_type node = this->index->firstNode(); node < this->index->sigma(); node += 2)
  {
    size_t offset = this->node_offset(node) / 2;
    if(this->real_nodes[offset])
    {
      total_length += get_source_length(offset);
    }
  }

  this->sequences.reserve(2 * total_length);
  size_t potential_nodes = this->index->sigma() - this->index->firstNode();
  this->offsets = sdsl::int_vector<0>(potential_nodes + 1, 0, gbwt::bit_length(this->sequences.capacity()));
}

void
GBWTGraph::cache_sequences(const std::function<std::string(size_t)>& get_source_sequence)
{
  for(gbwt::node_type node = this->index->firstNode(); node < this->index->sigma(); node += 2)
  {
    std::string seq;
    size_t offset = this->node_offset(node);
    if(this->real_nodes[offset / 2])
    {
      seq = get_source_sequence(offset / 2);
    }
    this->sequences.insert(this->sequences.end(), seq.begin(), seq.end());
    this->offsets[offset + 1] = this->sequences.size();
    reverse_complement_in_place(seq);
    this->sequences.insert(this->sequences.end(), seq.begin(), seq.end());
    this->offsets[offset + 2] = this->sequences.size();
  }
}

void
GBWTGraph::copy(const GBWTGraph& source)
{
  this->index = source.index;
  this->header = source.header;
  this->sequences = source.sequences;
  this->offsets = source.offsets;
  this->real_nodes = source.real_nodes;
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
  return this->offsets[offset + 1] - this->offsets[offset];
}

std::string
GBWTGraph::get_sequence(const handle_t& handle) const
{
  size_t offset = this->node_offset(handle);
  return std::string(this->sequences.begin() + this->offsets[offset], this->sequences.begin() + this->offsets[offset + 1]);
}

char
GBWTGraph::get_base(const handle_t& handle, size_t index) const
{
  size_t offset = this->node_offset(handle);
  return this->sequences[this->offsets[offset] + index];
}

std::string
GBWTGraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const
{
  size_t offset = this->node_offset(handle);
  size_t start = std::min(static_cast<size_t>(this->offsets[offset] + index), static_cast<size_t>(this->offsets[offset + 1]));
  size = std::min(size, static_cast<size_t>(this->offsets[offset + 1] - start));
  return std::string(this->sequences.begin() + start, this->sequences.begin() + start + size);
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
  // Incoming edges correspond to the outgoing edges of the reverse node.
  gbwt::node_type curr = handle_to_node(handle);
  if(go_left) { curr = gbwt::Node::reverse(curr); }

  // Cache the node.
  gbwt::CachedGBWT cache(*(this->index), true);
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
  gbwt::CachedGBWT cache(*(this->index), true);
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
  gbwt::CachedGBWT cache(*(this->index), true);
  gbwt::size_type cache_index = cache.findRecord(curr);

  for(gbwt::rank_type outrank = 0; outrank < cache.outdegree(cache_index); outrank++)
  {
    gbwt::node_type next = cache.successor(cache_index, outrank);
    if(node_to_handle(next) == right) { return true; }
  }

  return false;
}

//------------------------------------------------------------------------------

void
GBWTGraph::serialize(std::ostream& out) const
{
  // Serialize the header.
  out.write(reinterpret_cast<const char*>(&(this->header)), sizeof(Header));

  // Serialize the sequences manually.
  size_t sequence_size = this->sequences.size();
  sdsl::write_member(sequence_size, out);
  for(size_t offset = 0; offset < this->sequences.size(); offset += BLOCK_SIZE)
  {
    size_t bytes = std::min(BLOCK_SIZE, this->sequences.size() - offset);
    out.write(this->sequences.data() + offset, bytes);
  }

  // Serialize the rest.
  this->offsets.serialize(out);
  this->real_nodes.serialize(out);
}

void
GBWTGraph::deserialize(std::istream& in)
{
  // Load the header.
  in.read(reinterpret_cast<char*>(&(this->header)), sizeof(Header));
  if(!(this->header.check()))
  {
    std::cerr << "GBWTGraph::deserialize(): Invalid or old graph file" << std::endl;
    std::cerr << "GBWTGraph::deserialize(): Graph version is " << this->header.version << "; expected " << Header::VERSION << std::endl;
    return;
  }

  // Load the sequences manually.
  size_t sequence_size = 0;
  sdsl::read_member(sequence_size, in);
  this->sequences = std::vector<char>(sequence_size, 0);
  for(size_t offset = 0; offset < this->sequences.size(); offset += BLOCK_SIZE)
  {
    size_t bytes = std::min(BLOCK_SIZE, this->sequences.size() - offset);
    in.read(this->sequences.data() + offset, bytes);
  }

  // Load the rest.
  this->offsets.load(in);
  this->real_nodes.load(in);
}

void
GBWTGraph::set_gbwt(const gbwt::GBWT& gbwt_index)
{
  this->index = &gbwt_index;

  // Sanity checks for the GBWT index.
  assert(this->index->bidirectional());
}

//------------------------------------------------------------------------------

GBWTGraph::view_type
GBWTGraph::get_sequence_view(const handle_t& handle) const
{
  size_t offset = this->node_offset(handle);
  return std::make_pair(this->sequences.data() + this->offsets[offset], this->offsets[offset + 1] - this->offsets[offset]);
}

bool
GBWTGraph::starts_with(const handle_t& handle, char c) const
{
  size_t offset = this->node_offset(handle);
  if(this->offsets[offset + 1] <= this->offsets[offset]) { return false; }
  return (this->sequences[this->offsets[offset]] == c);
}

bool
GBWTGraph::ends_with(const handle_t& handle, char c) const
{
  size_t offset = this->node_offset(handle);
  if(this->offsets[offset + 1] <= this->offsets[offset]) { return false; }
  return (this->sequences[this->offsets[offset + 1] - 1] == c);
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

bool
GBWTGraph::follow_paths(gbwt::SearchState state, const std::function<bool(const gbwt::SearchState&)>& iteratee) const
{
  gbwt::CachedGBWT cache(*(this->index), true);
  return this->follow_paths(cache, state, iteratee);
}

bool
GBWTGraph::follow_paths(gbwt::BidirectionalState state, bool backward,
                        const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const
{
  gbwt::CachedGBWT cache(*(this->index), true);
  return this->follow_paths(cache, state, backward, iteratee);
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
