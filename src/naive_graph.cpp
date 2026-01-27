#include <gbwtgraph/naive_graph.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Other class variables.

const std::string NaiveGraph::TRANSLATION_EXTENSION = ".trans";

//------------------------------------------------------------------------------

// Construction.

NaiveGraph::NaiveGraph() : min_id(std::numeric_limits<nid_t>::max()), max_id(0)
{
}

void
NaiveGraph::create_node(nid_t node_id, std::string_view sequence)
{
  assert(!this->uses_translation());
  if(sequence.empty())
  {
    throw std::runtime_error("NaiveGraph::create_node: empty sequence for node " + std::to_string(node_id));
  }

  auto result = this->nodes.try_emplace(node_id, Node());
  if(!result.second)
  {
    throw std::runtime_error("NaiveGraph::create_node: duplicate node " + std::to_string(node_id));
  }
  result.first->second.sequence_offset = this->sequences.size();
  result.first->second.sequence_length = sequence.size();
  this->sequences.insert(this->sequences.end(), sequence.begin(), sequence.end());

  if(node_id < this->min_id) { this->min_id = node_id; }
  if(node_id > this->max_id) { this->max_id = node_id; }
}

std::pair<nid_t, nid_t>
NaiveGraph::translate_segment(const std::string& name, std::string_view sequence, size_t max_length)
{
  assert(this->nodes.empty() || this->uses_translation());
  if(sequence.empty())
  {
    throw std::runtime_error("NaiveGraph::translate_segment: empty sequence for segment " + name);
  }

  auto result = this->segment_translation.try_emplace(name, std::pair<nid_t, nid_t>(this->max_id + 1, this->max_id + 1));
  if(!result.second)
  {
    throw std::runtime_error("NaiveGraph::translate_segment: duplicate segment " + name);
  }

  for(size_t offset = 0; offset < sequence.size(); offset += max_length)
  {
    size_t length = std::min(max_length, sequence.size() - offset);
    Node& node = this->nodes[result.first->second.second];
    result.first->second.second++;
    node.sequence_offset = this->sequences.size() + offset;
    node.sequence_length = length;
  }
  this->sequences.insert(this->sequences.end(), sequence.begin(), sequence.end());

  // We always add at least one node beyond the previous max_id.
  this->min_id = std::min(this->min_id, result.first->second.first);
  this->max_id = result.first->second.second - 1;

  return result.first->second;
}

void
NaiveGraph::create_edge(gbwt::node_type from, gbwt::node_type to)
{
  nid_t from_id = gbwt::Node::id(from);
  nid_t to_id = gbwt::Node::id(to);
  auto from_iter = this->get_node_mut(from_id);
  auto to_iter = this->get_node_mut(to_id);
  if(from_iter == this->nodes.end() || to_iter == this->nodes.end())
  {
    throw std::runtime_error("NaiveGraph::create_edge: Cannot create an edge between nodes " + std::to_string(from_id) + " and " + std::to_string(to_id));
  }

  if(gbwt::Node::is_reverse(from))
  {
    from_iter->second.predecessors.push_back(gbwt::Node::reverse(to));
  }
  else
  {
    from_iter->second.successors.push_back(to);
  }

  if(gbwt::Node::is_reverse(to))
  {
    to_iter->second.successors.push_back(gbwt::Node::reverse(from));
  }
  else
  {
    to_iter->second.predecessors.push_back(from);
  }
}

void
NaiveGraph::remove_duplicate_edges()
{
  for(auto& node_pair : this->nodes)
  {
    gbwt::removeDuplicates(node_pair.second.predecessors, false);
    gbwt::removeDuplicates(node_pair.second.successors, false);
  }
}

//------------------------------------------------------------------------------

// HandleGraph interface.

bool
NaiveGraph::has_node(nid_t node_id) const
{
  return this->has_node_impl(node_id);
}

handle_t
NaiveGraph::get_handle(const nid_t& node_id, bool is_reverse) const
{
  return node_to_handle(gbwt::Node::encode(node_id, is_reverse));
}

nid_t
NaiveGraph::get_id(const handle_t& handle) const
{
  return gbwt::Node::id(handle_to_node(handle));
}

bool
NaiveGraph::get_is_reverse(const handle_t& handle) const
{
  return gbwt::Node::is_reverse(handle_to_node(handle));
}

handle_t
NaiveGraph::flip(const handle_t& handle) const
{
  return node_to_handle(gbwt::Node::reverse(handle_to_node(handle)));
}

size_t
NaiveGraph::get_length(const handle_t& handle) const
{
  auto iter = this->get_node(handle);
  if(iter == this->nodes.end()) { return 0; }
  return iter->second.sequence_length;
}

std::string
NaiveGraph::get_sequence(const handle_t& handle) const
{
  auto iter = this->get_node(handle);
  if(iter == this->nodes.end()) { return std::string(); }

  size_t offset = iter->second.sequence_offset;
  size_t length = iter->second.sequence_length;
  std::string sequence(this->sequences.begin() + offset, this->sequences.begin() + offset + length);
  if(this->get_is_reverse(handle))
  {
    reverse_complement_in_place(sequence);
  }
  return sequence;
}

char
NaiveGraph::get_base(const handle_t& handle, size_t index) const
{
  auto iter = this->get_node(handle);
  if(iter == this->nodes.end()) { return 'N'; }

  size_t offset = iter->second.sequence_offset;
  size_t length = iter->second.sequence_length;
  if(index >= length) { return 'N'; }

  char base = this->sequences[offset + index];
  if(this->get_is_reverse(handle))
  {
    base = complement(base);
  }
  return base;
}

std::string
NaiveGraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const
{
  auto iter = this->get_node(handle);
  if(iter == this->nodes.end()) { return std::string(); }

  size_t offset = iter->second.sequence_offset;
  size_t length = iter->second.sequence_length;
  if(index >= length) { return std::string(); }

  size_t actual_size = std::min(size, length - index);
  std::string subsequence(this->sequences.begin() + offset + index, this->sequences.begin() + offset + index + actual_size);
  if(this->get_is_reverse(handle))
  {
    reverse_complement_in_place(subsequence);
  }
  return subsequence;
}

size_t
NaiveGraph::get_node_count() const
{
  return this->nodes.size();
}

nid_t
NaiveGraph::min_node_id() const
{
  return this->min_id;
}

nid_t
NaiveGraph::max_node_id() const
{
  return this->max_id;
}

bool
NaiveGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const
{
  auto iter = this->get_node(handle);
  if(iter == this->nodes.end()) { return true; }

  bool flip = this->get_is_reverse(handle);
  const std::vector<gbwt::node_type>& edges = (go_left ^ flip ? iter->second.predecessors : iter->second.successors);
  for(gbwt::node_type node : edges)
  {
    handle_t next = node_to_handle(flip ? gbwt::Node::reverse(node) : node);
    if(!iteratee(next)) { return false; }
  }
  return true;
}

bool
NaiveGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool) const
{
  for(auto iter = this->nodes.begin(); iter != this->nodes.end(); ++iter)
  {
    handle_t handle = node_to_handle(gbwt::Node::encode(iter->first, false));
    if(!iteratee(handle)) { return false; }
  }
  return true;
}

size_t
NaiveGraph::get_degree(const handle_t& handle, bool go_left) const
{
  auto iter = this->get_node(handle);
  if(iter == this->nodes.end()) { return 0; }

  bool flip = this->get_is_reverse(handle);
  const std::vector<gbwt::node_type>& edges = (go_left ^ flip ? iter->second.predecessors : iter->second.successors);
  return edges.size();
}

//------------------------------------------------------------------------------

// Custom interface.

std::string_view
NaiveGraph::get_sequence_view(nid_t node_id) const
{
  auto iter = this->get_node(node_id);
  if(iter == this->nodes.end()) { return std::string_view(); }

  size_t offset = iter->second.sequence_offset;
  size_t length = iter->second.sequence_length;
  return std::string_view(this->sequences.data() + offset, length);
}

std::pair<nid_t, nid_t>
NaiveGraph::translate(const std::string& segment_name) const
{
  if(this->uses_translation())
  {
    auto iter = this->segment_translation.find(segment_name);
    if(iter == this->segment_translation.end()) { return no_translation(); }
    return iter->second;
  }
  else
  {
    auto parse = parse_unsigned<nid_t>(segment_name);
    if(parse.second && this->has_node_impl(parse.first))
    {
      return std::pair<nid_t, nid_t>(parse.first, parse.first + 1);
    }
    return no_translation();
  }
}

// invert_translation
std::pair<gbwt::StringArray, sdsl::sd_vector<>>
NaiveGraph::invert_translation(const std::function<bool(std::pair<nid_t, nid_t>)>& is_present) const
{
  std::pair<gbwt::StringArray, sdsl::sd_vector<>> result;

  // Invert the translation.
  // This stores half-open ranges of node identifiers corresponding to segments, and views to their segment names.
  std::vector<std::pair<std::pair<nid_t, nid_t>, std::string_view>> inverse;
  inverse.reserve(this->segment_translation.size());
  for(auto iter = this->segment_translation.begin(); iter != this->segment_translation.end(); ++iter)
  {
    inverse.emplace_back(iter->second, std::string_view(iter->first));
  }
  gbwt::parallelQuickSort(inverse.begin(), inverse.end());

  // Store the segment names.
  std::string empty;
  result.first = gbwt::StringArray(inverse.size(),
  [&](size_t offset) -> std::string_view
  {
    // This produces a view to each string to store.
    if(is_present(inverse[offset].first)) { return inverse[offset].second; }
    else { return std::string_view(empty); }
  });

  // Store the mapping.
  sdsl::sd_vector_builder builder(this->max_id + 1, inverse.size());
  for(auto& translation : inverse)
  {
    builder.set_unsafe(translation.first.first);
  }
  result.second = sdsl::sd_vector<>(builder);

  return result;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
