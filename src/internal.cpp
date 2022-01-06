#include <gbwtgraph/internal.h>

#include <algorithm>
#include <functional>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

BufferedHashSet::BufferedHashSet()
{
}

BufferedHashSet::~BufferedHashSet()
{
  // Wait for the insertion thread to finish.
  if(this->worker.joinable()) { this->worker.join(); }
}

void
BufferedHashSet::finish()
{
  // Flush the buffer if necessary.
  this->flush();

  // Wait for the insertion thread to finish.
  if(this->worker.joinable()) { this->worker.join(); }
}

void
insert_keys(std::unordered_set<std::string>& data, std::vector<std::string>& buffer)
{
  for(std::string& str : buffer) { data.insert(std::move(str)); }
}

void
BufferedHashSet::flush()
{
  // Wait for the insertion thread to finish.
  if(this->worker.joinable()) { this->worker.join(); }

  // Swap the input buffer and the internal buffer.
  this->internal_buffer.swap(this->input_buffer);
  this->input_buffer.clear();

  // Launch a new construction thread if necessary.
  if(this->internal_buffer.size() > 0)
  {
    this->worker = std::thread(insert_keys, std::ref(this->data), std::ref(this->internal_buffer));
  }
}

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

bool
GFAGraph::insert_edge(const handle_t& from, const handle_t& to)
{
  auto from_iter = this->get_node_mut(from);
  auto to_iter = this->get_node_mut(to);
  if(!(from_iter == this->nodes.end() || to_iter == this->nodes.end())) { return false; }

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

  return true;
}

void
GFAGraph::remove_duplicate_edges()
{
  this->for_each_handle([&](const handle_t& handle) {
    auto iter = this->get_node_mut(handle);
    gbwt::removeDuplicates(iter->second.predecessors, false);
    gbwt::removeDuplicates(iter->second.successors, false);
  });
}

bool
GFAGraph::has_node(nid_t node_id) const
{
  return (this->nodes.find(node_id) != this->nodes.end());
}

handle_t
GFAGraph::get_handle(const nid_t& node_id, bool is_reverse) const
{
  return node_to_handle(gbwt::Node::encode(node_id, is_reverse));
}

nid_t
GFAGraph::get_id(const handle_t& handle) const
{
  return gbwt::Node::id(handle_to_node(handle));
}

bool
GFAGraph::get_is_reverse(const handle_t& handle) const
{
  return gbwt::Node::is_reverse(handle_to_node(handle));
}

handle_t
GFAGraph::flip(const handle_t& handle) const
{
  return node_to_handle(gbwt::Node::reverse(handle_to_node(handle)));
}

size_t
GFAGraph::get_length(const handle_t& handle) const
{
  auto iter = this->get_node(handle);
  return iter->second.sequence.second;
}

std::string
GFAGraph::get_sequence(const handle_t& handle) const
{
  auto iter = this->get_node(handle);
  view_type view = iter->second.sequence;
  return std::string(view.first, view.second);
}

char
GFAGraph::get_base(const handle_t& handle, size_t index) const
{
  auto iter = this->get_node(handle);
  view_type view = iter->second.sequence;
  return *(view.first + index);
}

std::string
GFAGraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const
{
  auto iter = this->get_node(handle);
  view_type view = iter->second.sequence;
  index = std::min(index, view.second);
  size = std::min(size, view.second - index);
  return std::string(view.first + index, view.first + index + size);
}

size_t
GFAGraph::get_node_count() const
{
  return this->nodes.size();
}

nid_t
GFAGraph::min_node_id() const
{
  return this->min_id;
}

nid_t
GFAGraph::max_node_id() const
{
  return this->max_id;
}

bool
GFAGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const
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
GFAGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool) const
{
  for(auto iter = this->nodes.begin(); iter != this->nodes.end(); ++iter)
  {
    if(!iteratee(this->get_handle(iter->first, false))) { return false; }
  }
  return true;
}

size_t
GFAGraph::get_degree(const handle_t& handle, bool go_left) const
{
  auto iter = this->get_node(handle);
  bool flip = this->get_is_reverse(handle);
  const std::vector<handle_t>& edges = (go_left ^ flip ? iter->second.predecessors : iter->second.successors);
  return edges.size();
}

view_type
GFAGraph::get_sequence_view(const handle_t& handle) const
{
  auto iter = this->get_node(handle);
  return iter->second.sequence;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
