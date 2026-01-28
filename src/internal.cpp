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
TSVWriter::write(std::string_view view)
{
  size_t offset = 0;
  while(offset < view.size())
  {
    size_t length = std::min(view.size() - offset, BUFFER_SIZE - this->buffer.size());
    this->buffer.insert(this->buffer.end(), view.begin() + offset, view.begin() + offset + length);
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

std::string
GFAGrammar::first_segment(const std::string& symbol) const
{
  std::string result = symbol;
  for(rule_type rule = this->expand(result); rule != this->no_rule(); rule = this->expand(result))
  {
    if(rule->second.empty()) { break; }
    result = rule->second.front().first;
  }
  return result;
}

void
GFAGrammar::validate(const NaiveGraph& graph) const
{
  // Check the rules individually.
  for(const auto& rule : this->rules)
  {
    if(rule.first.empty())
    {
      throw std::runtime_error("GFAGrammar::validate(): Empty rule name");
    }
    if(graph.has_node_or_segment(rule.first))
    {
      throw std::runtime_error("GFAGrammar::validate(): Rule name " + rule.first + " clashes with a node/segment name");
    }
    if(rule.second.size() < 2)
    {
      throw std::runtime_error("GFAGrammar::validate(): Expansion of rule " + rule.first + " is trivial");
    }
    for(const auto& symbol : rule.second)
    {
      rule_type expansion = this->expand(symbol.first);
      if(graph.has_node_or_segment(symbol.first))
      {
        if(expansion != this->no_rule())
        {
          throw std::runtime_error("GFAGrammar::validate(): Symbol " + symbol.first + " in the expansion of rule " + rule.first + " is both a node/segment and a rule");
        }
      }
      else
      {
        if(expansion == this->no_rule())
        {
          throw std::runtime_error("GFAGrammar::validate(): Symbol " + symbol.first + " in the expansion of rule " + rule.first + " is undefined");
        }
      }
    }
  }

  // TODO: Use a stack to avoid stack overflows on unbalanced grammars.
  // Check for cycles using depth-first search.
  enum class State { ACTIVE, VISITED };
  std::unordered_map<std::string, State> states;
  std::function<void(const std::string&)> dfs = [&](const std::string& rule_name)
  {
    states[rule_name] = State::ACTIVE;
    const expansion_type& expansion = this->rules.at(rule_name);
    for(const auto& symbol : expansion)
    {
      if(this->rules.find(symbol.first) == this->rules.end()) { continue; }
      auto result = states.emplace(symbol.first, State::ACTIVE);
      if(result.second)
      {
        // Unvisited rule.
        dfs(symbol.first);
      }
      else if(result.first->second == State::ACTIVE)
      {
        throw std::runtime_error("GFAGrammar::validate(): Cycle detected at rule " + symbol.first);
      }
    }
    states[rule_name] = State::VISITED;
  };

  for(const auto& rule : this->rules)
  {
    if(states.find(rule.first) == states.end())
    {
      dfs(rule.first);
    }
  }
}

GFAGrammarIterator
GFAGrammar::iter(const std::string& rule_name, bool is_reverse) const
{
  return GFAGrammarIterator(*this, rule_name, is_reverse);
}

GFAGrammarIterator::GFAGrammarIterator(const GFAGrammar& grammar, const std::string& rule_name, bool is_reverse) :
  grammar(grammar)
{
  GFAGrammar::rule_type rule = grammar.expand(rule_name);
  if(rule != grammar.no_rule())
  {
    this->stack.push_back({ rule, is_reverse, 0 });
  }
}

std::pair<std::string_view, bool>
GFAGrammarIterator::next()
{
  while(!(this->stack.empty()))
  {
    auto& top = this->stack.back();
    auto& expansion = std::get<0>(top)->second;
    if(std::get<2>(top) >= expansion.size())
    {
      this->stack.pop_back();
      continue;
    }

    // Determine the next symbol and its orientation.
    bool is_reverse = std::get<1>(top);
    size_t index = (is_reverse ? (expansion.size() - 1 - std::get<2>(top)) : std::get<2>(top));
    std::get<2>(top)++;
    const std::string& symbol = expansion[index].first;
    is_reverse ^= expansion[index].second;
    GFAGrammar::rule_type rule = this->grammar.expand(symbol);
    if(rule != this->grammar.no_rule()) { this->stack.push_back({ rule, is_reverse, 0 }); }
    else { return std::make_pair(std::string_view(symbol), is_reverse); }
  }
  return std::make_pair(std::string_view(), false);
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
