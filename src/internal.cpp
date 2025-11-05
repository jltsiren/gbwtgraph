#include <gbwtgraph/internal.h>
#include <gbwtgraph/gbwtgraph.h>

#include <algorithm>
#include <cctype>
#include <functional>
#include <limits>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t TSVWriter::BUFFER_SIZE;

constexpr size_t ManualTSVWriter::BUFFER_SIZE;
constexpr size_t ManualTSVWriter::BUFFER_FULL;

constexpr size_t PathIdMap::MAX_HAPLOTYPES;

//------------------------------------------------------------------------------

TSVWriter::TSVWriter(std::ostream& out) :
  out(out)
{
  this->buffer.reserve(BUFFER_SIZE);
}

TSVWriter::~TSVWriter()
{
  this->flush();
}

void
TSVWriter::write(view_type view)
{
  size_t offset = 0;
  while(offset < view.second)
  {
    size_t length = std::min(view.second - offset, BUFFER_SIZE - this->buffer.size());
    this->buffer.insert(this->buffer.end(), view.first + offset, view.first + offset + length);
    if(this->buffer.size() >= BUFFER_SIZE) { this->flush(); }
    offset += length;
  }
}

void
TSVWriter::flush()
{
  if(!(this->buffer.empty()))
  {
    this->out.write(this->buffer.data(), this->buffer.size());
    this->buffer.clear();
  }
}

//------------------------------------------------------------------------------

ManualTSVWriter::ManualTSVWriter(std::ostream& out) :
  out(out)
{
  this->buffer.reserve(BUFFER_SIZE);
}

void
ManualTSVWriter::flush()
{
  if(!(this->buffer.empty()))
  {
    this->out.write(this->buffer.data(), this->buffer.size());
    this->buffer.clear();
  }
}

//------------------------------------------------------------------------------

void
EmptyGraph::create_node(nid_t node_id)
{
  this->nodes[node_id] = { {}, {} };
  this->min_id = std::min(this->min_id, node_id);
  this->max_id = std::max(this->max_id, node_id);
}

void
EmptyGraph::create_edge(const handle_t& from, const handle_t& to)
{
  auto from_iter = this->get_node_mut(from);
  auto to_iter = this->get_node_mut(to);
  if(from_iter == this->nodes.end() || to_iter == this->nodes.end())
  {
    nid_t from_id = gbwt::Node::id(handle_to_node(from));
    nid_t to_id = gbwt::Node::id(handle_to_node(to));
    throw std::runtime_error("EmptyGraph: Cannot create an edge between nodes " + std::to_string(from_id) + " and " + std::to_string(to_id));
  }

  // from -> to
  if(this->get_is_reverse(from))
  {
    from_iter->second.predecessors.push_back(this->flip(to));
  }
  else
  {
    from_iter->second.successors.push_back(to);
  }

  // to -> from
  if(this->get_is_reverse(to))
  {
    to_iter->second.successors.push_back(this->flip(from));
  }
  else
  {
    to_iter->second.predecessors.push_back(from);
  }
}

void
EmptyGraph::remove_duplicate_edges()
{
  this->for_each_handle([&](const handle_t& handle) {
    auto iter = this->get_node_mut(handle);
    gbwt::removeDuplicates(iter->second.predecessors, false);
    gbwt::removeDuplicates(iter->second.successors, false);
  });
}

bool
EmptyGraph::has_node(nid_t node_id) const
{
  return (this->nodes.find(node_id) != this->nodes.end());
}

handle_t
EmptyGraph::get_handle(const nid_t& node_id, bool is_reverse) const
{
  return node_to_handle(gbwt::Node::encode(node_id, is_reverse));
}

nid_t
EmptyGraph::get_id(const handle_t& handle) const
{
  return gbwt::Node::id(handle_to_node(handle));
}

bool
EmptyGraph::get_is_reverse(const handle_t& handle) const
{
  return gbwt::Node::is_reverse(handle_to_node(handle));
}

handle_t
EmptyGraph::flip(const handle_t& handle) const
{
  return node_to_handle(gbwt::Node::reverse(handle_to_node(handle)));
}

size_t
EmptyGraph::get_length(const handle_t&) const
{
  return 0;
}

std::string
EmptyGraph::get_sequence(const handle_t&) const
{
  return std::string();
}

char
EmptyGraph::get_base(const handle_t&, size_t) const
{
  return 'N';
}

std::string
EmptyGraph::get_subsequence(const handle_t&, size_t, size_t) const
{
  return std::string();
}

size_t
EmptyGraph::get_node_count() const
{
  return this->nodes.size();
}

nid_t
EmptyGraph::min_node_id() const
{
  return this->min_id;
}

nid_t
EmptyGraph::max_node_id() const
{
  return this->max_id;
}

bool
EmptyGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const
{
  auto iter = this->get_node(handle);
  bool flip = this->get_is_reverse(handle);
  const std::vector<handle_t>& edges = (go_left ^ flip ? iter->second.predecessors : iter->second.successors);
  for(const handle_t& next : edges)
  {
    handle_t actual = (flip ? this->flip(next) : next);
    if(!iteratee(actual)) { return false; }
  }
  return true;
}

bool
EmptyGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool) const
{
  for(auto iter = this->nodes.begin(); iter != this->nodes.end(); ++iter)
  {
    if(!iteratee(this->get_handle(iter->first, false))) { return false; }
  }
  return true;
}

size_t
EmptyGraph::get_degree(const handle_t& handle, bool go_left) const
{
  auto iter = this->get_node(handle);
  bool flip = this->get_is_reverse(handle);
  const std::vector<handle_t>& edges = (go_left ^ flip ? iter->second.predecessors : iter->second.successors);
  return edges.size();
}

//------------------------------------------------------------------------------

LargeRecordCache::LargeRecordCache(const gbwt::GBWT& index, size_t bytes) :
  index(index)
{
  for(gbwt::node_type node = this->index.firstNode(); node < this->index.sigma(); node++)
  {
    std::pair<gbwt::size_type, gbwt::size_type> range = this->index.bwt.getRange(this->index.toComp(node));
    if(range.second - range.first > bytes && !(this->index.empty(node)))
    {
      this->cache[node] = gbwt::DecompressedRecord(this->index.record(node));
    }
  }
}

gbwt::vector_type
LargeRecordCache::extract(gbwt::size_type sequence) const
{
  gbwt::vector_type result;
  if(sequence > this->sequences()) { return result; }

  gbwt::edge_type pos = this->index.start(sequence);
  while(pos.first != gbwt::ENDMARKER)
  {
    result.push_back(pos.first);
    auto iter = this->cache.find(pos.first);
    if(iter != this->cache.end()) { pos = iter->second.LF(pos.second); }
    else { pos = this->index.LF(pos); }
  }

  return result;
}

//------------------------------------------------------------------------------

std::vector<std::pair<size_t, gbwt::edge_type>>
sample_path_positions(const GBZ& gbz, path_handle_t path, size_t sample_interval, size_t* length)
{
  std::vector<std::pair<size_t, gbwt::edge_type>> result;
  gbwt::size_type seq_id = gbwt::Path::encode(gbz.graph.handle_to_path(path), false);

  size_t offset = 0, next_sample = 0;
  for(gbwt::edge_type pos = gbz.index.start(seq_id); pos.first != gbwt::ENDMARKER; pos = gbz.index.LF(pos))
  {
    if(offset >= next_sample)
    {
      result.push_back({ offset, pos });
      next_sample = offset + sample_interval;
    }
    offset += gbz.graph.get_length(GBWTGraph::node_to_handle(pos.first));
  }
  if(length != nullptr) { *length = offset; }

  return result;
}

//------------------------------------------------------------------------------

PathIdMap::PathIdMap(const gbwt::Metadata& metadata) :
  mask
  (
    std::numeric_limits<gbwt::PathName::path_name_type>::max(),
    std::numeric_limits<gbwt::PathName::path_name_type>::max(),
    std::numeric_limits<gbwt::PathName::path_name_type>::max(),
    0
  ), key(KeyType::SAMPLE_CONTIG_HAPLOTYPE)
{
  // Sort the paths by (sample, phase, contig, count) instead of the natural order.
  // This is because we fall back from (sample, haplotype, contig) to (sample, haplotype).
  std::vector<gbwt::PathName> sorted_paths;
  sorted_paths.reserve(metadata.paths());
  for(gbwt::size_type i = 0; i < metadata.paths(); i++)
  {
    sorted_paths.push_back(metadata.path(i));
  }
  std::sort(sorted_paths.begin(), sorted_paths.end(), [](const gbwt::PathName& left, const gbwt::PathName& right)
  {
    if(left.sample != right.sample) { return (left.sample < right.sample); }
    if(left.phase != right.phase) { return (left.phase < right.phase); }
    if(left.contig != right.contig) { return (left.contig < right.contig); }
    return (left.count < right.count);
  });

  // First attempt: (sample, contig, haplotype).
  if(this->build_map(sorted_paths)) { return; }

  // Second attempt: (sample, haplotype).
  this->mask.contig = 0;
  this->key = KeyType::SAMPLE_HAPLOTYPE;
  if(this->build_map(sorted_paths)) { return; }

  // Third attempt: (sample). Falls back to an empty map on failure.
  this->mask.phase = 0;
  this->key = KeyType::SAMPLE;
  this->build_map(sorted_paths);
}

bool
PathIdMap::build_map(const std::vector<gbwt::PathName>& sorted_paths)
{
  size_t next = 0;
  for(const gbwt::PathName& path : sorted_paths)
  {
    gbwt::PathName key = this->mask_name(path);
    if(this->path_to_id.find(key) == this->path_to_id.end())
    {
      if(next >= MAX_HAPLOTYPES)
      {
        this->path_to_id.clear();
        this->key = KeyType::NONE;
        return false;
      }
      this->path_to_id[key] = next;
      next++;
    }
  }
  return true;
}

std::string
PathIdMap::key_type_str(KeyType key)
{
  switch(key)
  {
    case KeyType::NONE:
      return "()";
    case KeyType::SAMPLE:
      return "(sample)";
    case KeyType::SAMPLE_HAPLOTYPE:
      return "(sample, haplotype)";
    case KeyType::SAMPLE_CONTIG_HAPLOTYPE:
      return "(sample, contig, haplotype)";
    default:
      return "(unknown)";
  }
}

//------------------------------------------------------------------------------

std::vector<gbwt::node_type>
extract_kmer_path(const GBWTGraph& graph, const std::vector<handle_t>& path, size_t path_offset, size_t node_offset, size_t k, bool is_reverse)
{
  if(is_reverse) { node_offset = graph.get_length(path[path_offset]) - node_offset - 1; }

  std::vector<gbwt::node_type> result;
  size_t path_length = 0;
  while(path_length < k && path_offset < path.size())
  {
    handle_t handle = path[path_offset];
    path_length += graph.get_length(handle) - node_offset;
    node_offset = 0;
    gbwt::node_type node = GBWTGraph::handle_to_node(handle);
    if(is_reverse)
    {
      result.push_back(gbwt::Node::reverse(node));
      if(path_offset == 0) { break; }
      path_offset--;
    }
    else
    {
      result.push_back(node);
      path_offset++;
    }
  }

  return result;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
