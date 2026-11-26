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
NaiveGraph::normalize_sequence(Node& node)
{
  std::vector<char>& sequence = (node.string_index == 0 ? this->concatenated_sequences : this->superstrings[node.string_index - 1]);
  this->normalize_sequence(sequence, node.sequence_offset, node.sequence_length);
}

void
NaiveGraph::normalize_sequence(std::vector<char>& sequence, size_t from_offset, size_t length)
{
  size_t end = from_offset + length;
  for(size_t i = from_offset; i < end; i++)
  {
    sequence[i] = std::toupper(static_cast<unsigned char>(sequence[i]));
  }
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
  result.first->second.string_index = 0; // We append the sequence to `concatenated_sequences`.
  result.first->second.sequence_offset = this->concatenated_sequences.size();
  result.first->second.sequence_length = sequence.size();
  this->concatenated_sequences.insert(this->concatenated_sequences.end(), sequence.begin(), sequence.end());
  this->normalize_sequence(result.first->second);

  if(node_id < this->min_id) { this->min_id = node_id; }
  if(node_id > this->max_id) { this->max_id = node_id; }
}

void
NaiveGraph::add_superstring(std::vector<char>&& superstring)
{
  this->normalize_sequence(superstring, 0, superstring.size());
  this->superstrings.push_back(std::move(superstring));
}

void
NaiveGraph::create_dedup_node(nid_t node_id, size_t superstring_index, size_t offset, size_t length)
{
  assert(!this->uses_translation());
  if(length == 0)
  {
    throw std::runtime_error("NaiveGraph::create_dedup_node: empty sequence for node " + std::to_string(node_id));
  }
  if(superstring_index >= this->superstrings.size())
  {
    std::string msg = "NaiveGraph::create_dedup_node: superstring " + std::to_string(superstring_index)
      + " does not exist for node " + std::to_string(node_id);
    throw std::runtime_error(msg);
  }
  if
  (
    offset >= this->superstrings[superstring_index].size() ||
    length > this->superstrings[superstring_index].size() - offset
  )
  {
    std::string msg = "NaiveGraph::create_dedup_node: the sequence for node " + std::to_string(node_id)
      + " is out of bounds in superstring " + std::to_string(superstring_index);
    throw std::runtime_error(msg);
  }

  auto result = this->nodes.try_emplace(node_id, Node());
  if(!result.second)
  {
    throw std::runtime_error("NaiveGraph::create_dedup_node: duplicate node " + std::to_string(node_id));
  }
  result.first->second.string_index = superstring_index + 1;
  result.first->second.sequence_offset = offset;
  result.first->second.sequence_length = length;
  // The superstring was already normalized in `add_superstring()`.

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
    node.string_index = 0; // We append the sequence to `concatenated_sequences`.
    node.sequence_offset = this->concatenated_sequences.size() + offset;
    node.sequence_length = length;
  }
  size_t from_offset = this->concatenated_sequences.size();
  this->concatenated_sequences.insert(this->concatenated_sequences.end(), sequence.begin(), sequence.end());
  this->normalize_sequence(this->concatenated_sequences, from_offset, sequence.size());

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

const std::vector<char>&
NaiveGraph::get_source(const Node& node) const
{
  return (node.string_index == 0 ? this->concatenated_sequences : this->superstrings[node.string_index - 1]);
}

std::string
NaiveGraph::get_sequence(const handle_t& handle) const
{
  auto iter = this->get_node(handle);
  if(iter == this->nodes.end()) { return std::string(); }

  const std::vector<char>& source = this->get_source(iter->second);
  size_t offset = iter->second.sequence_offset;
  size_t length = iter->second.sequence_length;
  std::string sequence = get_substring(source, offset, length);
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

  const std::vector<char>& source = this->get_source(iter->second);
  if(this->get_is_reverse(handle))
  {
    return complement(source[offset + length - 1 - index]);
  }
  return source[offset + index];
}

std::string
NaiveGraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const
{
  auto iter = this->get_node(handle);
  if(iter == this->nodes.end()) { return std::string(); }

  const std::vector<char>& source = this->get_source(iter->second);
  size_t offset = iter->second.sequence_offset;
  size_t length = iter->second.sequence_length;
  if(index >= length) { return std::string(); }

  size_t actual_size = std::min(size, length - index);
  if(this->get_is_reverse(handle))
  {
    // Index `index` in the reverse orientation corresponds to the end of the
    // substring in the forward orientation.
    std::string subsequence = get_substring(source, offset + length - index - actual_size, actual_size);
    reverse_complement_in_place(subsequence);
    return subsequence;
  }
  return get_substring(source, offset + index, actual_size);
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

// Custom interface: Sequences.

std::string_view
NaiveGraph::get_sequence_view(nid_t node_id) const
{
  auto iter = this->get_node(node_id);
  if(iter == this->nodes.end()) { return std::string_view(); }

  const std::vector<char>& source = this->get_source(iter->second);
  size_t offset = iter->second.sequence_offset;
  size_t length = iter->second.sequence_length;
  return get_subview(source, offset, length);
}

NaiveGraph::StringSource
NaiveGraph::get_sequence_source(nid_t node_id) const
{
  auto iter = this->get_node(node_id);
  if(iter == this->nodes.end())
  {
    return { std::string_view(), this->get_num_superstrings() };
  }

  const std::vector<char>& source = this->get_source(iter->second);
  std::string_view view = get_subview(source, iter->second.sequence_offset, iter->second.sequence_length);
  size_t index = (iter->second.string_index == 0 ? this->get_num_superstrings() : iter->second.string_index - 1);
  return { view, index };
}

//------------------------------------------------------------------------------

// Custom interface: Translation.

bool
NaiveGraph::has_node_or_segment(const std::string& name) const
{
  if(this->uses_translation())
  {
    return (this->segment_translation.find(name) != this->segment_translation.end());
  }
  else
  {
    auto parse = parse_unsigned<nid_t>(name);
    return (parse.second && this->has_node_impl(parse.first));
  }
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
