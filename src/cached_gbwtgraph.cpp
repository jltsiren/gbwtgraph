#include <gbwtgraph/cached_gbwtgraph.h>

#include <algorithm>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

CachedGBWTGraph::CachedGBWTGraph() :
  graph(nullptr), cache()
{
}

CachedGBWTGraph::CachedGBWTGraph(const CachedGBWTGraph& source)
{
  this->copy(source);
}

CachedGBWTGraph::CachedGBWTGraph(CachedGBWTGraph&& source)
{
  *this = std::move(source);
}

CachedGBWTGraph::~CachedGBWTGraph()
{
}

void
CachedGBWTGraph::swap(CachedGBWTGraph& another)
{
  if(&another == this) { return; }

  std::swap(this->graph, another.graph);
  this->cache.swap(another.cache);
}

CachedGBWTGraph&
CachedGBWTGraph::operator=(const CachedGBWTGraph& source)
{
  if(&source != this) { this->copy(source); }
  return *this;
}

CachedGBWTGraph&
CachedGBWTGraph::operator=(CachedGBWTGraph&& source)
{
  if(&source != this)
  {
    this->graph = std::move(source.graph);
    this->cache = std::move(source.cache);
  }
  return *this;
}

void
CachedGBWTGraph::copy(const CachedGBWTGraph& source)
{
  this->graph = source.graph;
  this->cache = source.cache;
}

//------------------------------------------------------------------------------

CachedGBWTGraph::CachedGBWTGraph(const GBWTGraph& graph) :
  graph(&graph), cache(this->graph->get_cache())
{
}

//------------------------------------------------------------------------------

bool
CachedGBWTGraph::has_node(nid_t node_id) const
{
  size_t offset = this->graph->node_offset(gbwt::Node::encode(node_id, false)) / 2;
  return (offset < this->graph->real_nodes.size() && this->graph->real_nodes[offset]);
}

handle_t
CachedGBWTGraph::get_handle(const nid_t& node_id, bool is_reverse) const
{
  return node_to_handle(gbwt::Node::encode(node_id, is_reverse));
}

nid_t
CachedGBWTGraph::get_id(const handle_t& handle) const
{
  return gbwt::Node::id(handle_to_node(handle));
}

bool
CachedGBWTGraph::get_is_reverse(const handle_t& handle) const
{
  return gbwt::Node::is_reverse(handle_to_node(handle));
}

handle_t
CachedGBWTGraph::flip(const handle_t& handle) const
{
  return node_to_handle(gbwt::Node::reverse(handle_to_node(handle)));
}

size_t
CachedGBWTGraph::get_length(const handle_t& handle) const
{
  size_t offset = this->graph->node_offset(handle);
  return this->graph->sequences.length(offset);
}

std::string
CachedGBWTGraph::get_sequence(const handle_t& handle) const
{
  size_t offset = this->graph->node_offset(handle);
  return this->graph->sequences.str(offset);
}

char
CachedGBWTGraph::get_base(const handle_t& handle, size_t index) const
{
  size_t offset = this->graph->node_offset(handle);
  std::string_view view = this->graph->sequences.view(offset);
  return view[index];
}

std::string
CachedGBWTGraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const
{
  size_t offset = this->graph->node_offset(handle);
  std::string_view view = this->graph->sequences.view(offset);
  index = std::min(index, view.size());
  size = std::min(size, view.size() - index);
  return std::string(view.data() + index, view.data() + index + size);
}

size_t
CachedGBWTGraph::get_node_count() const
{
  return this->graph->header.nodes;
}

nid_t
CachedGBWTGraph::min_node_id() const
{
  return gbwt::Node::id(this->cache.firstNode());
}

nid_t
CachedGBWTGraph::max_node_id() const
{
  nid_t next_id = gbwt::Node::id(this->cache.sigma());
  return next_id - 1;
}

bool
CachedGBWTGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const
{
  return this->graph->cached_follow_edges(this->cache, handle, go_left, iteratee);
}

bool
CachedGBWTGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const
{
  // This is expensive enough that we can afford using virtual functions.
  return this->graph->for_each_handle_impl(iteratee, parallel);
}

size_t
CachedGBWTGraph::get_degree(const handle_t& handle, bool go_left) const
{
  // Cache the node.
  gbwt::node_type curr = handle_to_node(handle);
  if(go_left) { curr = gbwt::Node::reverse(curr); }
  gbwt::size_type cache_index = this->cache.findRecord(curr);

  // The outdegree reported by GBWT might account for the endmarker, which is
  // always the first successor.
  size_t result = this->cache.outdegree(cache_index);
  if(result > 0 && this->cache.successor(cache_index, 0) == gbwt::ENDMARKER) { result--; }
  return result;
}

bool
CachedGBWTGraph::has_edge(const handle_t& left, const handle_t& right) const
{
  // Cache the node.
  gbwt::node_type curr = handle_to_node(left);
  gbwt::size_type cache_index = this->cache.findRecord(curr);

  for(gbwt::rank_type outrank = 0; outrank < this->cache.outdegree(cache_index); outrank++)
  {
    gbwt::node_type next = this->cache.successor(cache_index, outrank);
    if(node_to_handle(next) == right) { return true; }
  }

  return false;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
