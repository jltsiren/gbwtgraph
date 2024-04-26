#include <gbwtgraph/subgraph.h>
#include <gbwtgraph/internal.h>

#include <sdsl/bit_vectors.hpp>

#include <queue>
#include <unordered_map>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Numerical class constants.
constexpr size_t PathIndex::DEFAULT_SAMPLE_INTERVAL;

//------------------------------------------------------------------------------

PathIndex::PathIndex(const GBZ& gbz, size_t sample_interval) :
  sequence_positions(gbz.graph.named_paths.size()), gbwt_positions(gbz.graph.named_paths.size())
{
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t i = 0; i < gbz.graph.named_paths.size(); i++)
  {
    size_t length = 0;
    auto positions = sample_path_positions(gbz, handlegraph::as_path_handle(i), sample_interval, &length);
    sdsl::sd_vector_builder builder(length, positions.size());
    for(const auto& pos : positions)
    {
      builder.set_unsafe(pos.first);
      this->gbwt_positions[i].push_back(pos.second);
    }
    this->sequence_positions[i] = sdsl::sd_vector<>(builder);
  }
}

std::pair<size_t, gbwt::edge_type>
PathIndex::sampled_position(path_handle_t handle, size_t offset) const
{
  size_t path_id = handlegraph::as_integer(handle);
  if(path_id >= this->paths()) { return std::make_pair(0, gbwt::invalid_edge()); }
  auto iter = this->sequence_positions[path_id].predecessor(offset);
  if(iter == this->sequence_positions[path_id].one_end()) { return std::make_pair(0, gbwt::invalid_edge()); }
  return std::make_pair(iter->second, this->gbwt_positions[path_id][iter->first]);
}

//------------------------------------------------------------------------------

SubgraphQuery
SubgraphQuery::path_offset(path_handle_t path, size_t offset, size_t context, HaplotypeOutput output)
{
  SubgraphQuery result
  {
    QueryType::path_offset_query, output,
    path, offset,
    0,
    context
  };
  return result;
}

SubgraphQuery
SubgraphQuery::path_interval(path_handle_t path, size_t from, size_t to, size_t context, HaplotypeOutput output)
{
  if(from >= to)
  {
    std::string msg = "SubgraphQuery::path_interval(): Empty interval [" + std::to_string(from) + ", " + std::to_string(to) + ")";
    throw std::runtime_error(msg);
  }
  size_t offset = from + (to - from) / 2;

  SubgraphQuery result
  {
    QueryType::path_offset_query, output,
    path, offset,
    0,
    context
  };
  return result;
}

SubgraphQuery
SubgraphQuery::node(nid_t node, size_t context, HaplotypeOutput output)
{
  SubgraphQuery result
  {
    QueryType::node_query, output,
    handlegraph::as_path_handle(std::numeric_limits<size_t>::max()), 0,
    node,
    context
  };
  return result;
}

std::string
output_type(SubgraphQuery::HaplotypeOutput output)
{
  switch(output)
  {
  case SubgraphQuery::HaplotypeOutput::all_haplotypes: return "all";
  case SubgraphQuery::HaplotypeOutput::distinct_haplotypes: return "distinct";
  case SubgraphQuery::HaplotypeOutput::reference_only: return "reference only";
  default: return "invalid output";
  }
}

std::string
SubgraphQuery::to_string(const GBZ& gbz) const
{
  if(this->type == QueryType::path_offset_query)
  {
    std::string path_name = gbz.graph.get_path_name(this->path);
    std::string output = output_type(this->haplotype_output);
    return "(path " + path_name + ", offset " + std::to_string(this->offset) + ", context " + std::to_string(this->context) + ", " + output + ")";
  }
  else if(this->type == QueryType::node_query)
  {
    std::string output = output_type(this->haplotype_output);
    return "(node " + std::to_string(this->node_id) + ", context " + std::to_string(this->context) + ", " + output + ")";
  }
  else
  {
    return "(invalid query)";
  }
}

//------------------------------------------------------------------------------

/*
  Helper functions for the constructor.
*/

std::pair<pos_t, gbwt::edge_type>
find_position(const GBZ& gbz, const PathIndex& path_index, const SubgraphQuery& query)
{
  std::pair<size_t, gbwt::edge_type> position = path_index.sampled_position(query.path, query.offset);
  if(position.second == gbwt::invalid_edge())
  {
    std::string msg = "Subgraph::Subgraph(): Could not find an indexed path position for query " + std::to_string(query.offset);
    throw std::runtime_error(msg);
  }

  while(true)
  {
    size_t node_len = gbz.graph.get_length(GBWTGraph::node_to_handle(position.second.first));
    if(position.first + node_len > query.offset)
    {
      pos_t graph_pos = make_pos_t(gbwt::Node::id(position.second.first), gbwt::Node::is_reverse(position.second.first), query.offset - position.first);
      return std::make_pair(graph_pos, position.second);
    }
    position.first += node_len;
    position.second = gbz.index.LF(position.second);
    if(position.second == gbwt::invalid_edge())
    {
      std::string msg = "Subgraph::Subgraph(): Path " + gbz.graph.get_path_name(query.path) + " does not contain offset " + std::to_string(query.offset);
      throw std::runtime_error(msg);
    }
  }
}

std::set<nid_t>
find_subgraph(const GBWTGraph& graph, pos_t position, size_t context)
{
  std::set<nid_t> result;
  result.insert(id(position));

  struct Comparator
  {
    bool operator()(const std::pair<size_t, pos_t>& a, const std::pair<size_t, pos_t>& b) const
    {
      return a.first > b.first;
    }
  };

  // (distance, position)
  std::priority_queue<std::pair<size_t, pos_t>, std::vector<std::pair<size_t, pos_t>>, Comparator> queue;
  queue.push(std::make_pair(0, position));

  while(!queue.empty())
  {
    std::pair<size_t, pos_t> current = queue.top(); queue.pop();
    if(result.find(id(current.second)) != result.end()) { continue; }
    for(bool is_reverse : { false, true })
    {
      handle_t handle = graph.get_handle(id(current.second), is_reverse);
      size_t distance = current.first;
      if(is_reverse == is_rev(current.second)) { distance += graph.get_length(handle) - offset(current.second); }
      else { distance += offset(current.second) + 1; }
      if(distance <= context)
      {
        graph.follow_edges(handle, false, [&](const handle_t& next)
        {
          queue.push(std::make_pair(distance, make_pos_t(graph.get_id(next), graph.get_is_reverse(next), 0)));
        });
      }
      result.insert(id(current.second));
    }
  }

  return result;
}

void
append_edit(std::vector<std::pair<char, size_t>>& edits, char op, size_t length)
{
  if(!edits.empty() && edits.back().first == op)
  {
    edits.back().second += length;
  }
  else
  {
    edits.emplace_back(op, length);
  }
}

size_t
path_interval_length(const GBWTGraph& graph, const gbwt::vector_type& path, std::pair<size_t, size_t> interval)
{
  size_t result = 0;
  for(size_t i = interval.first; i < interval.second; i++)
  {
    result += graph.get_length(GBWTGraph::node_to_handle(path[i]));
  }
  return result;
}

// Create edits for two path intervals, which are assumed to be diverging.
void
append_edits(
  std::vector<std::pair<char, size_t>>& edits, const GBWTGraph& graph,
  const gbwt::vector_type& reference, std::pair<size_t, size_t> ref_interval,
  const gbwt::vector_type& path, std::pair<size_t, size_t> path_interval
)
{
  size_t ref_length = path_interval_length(graph, reference, ref_interval);
  size_t path_length = path_interval_length(graph, path, path_interval);
  if(path_length == ref_length && path_length > 0 && path_length < 5)
  {
    append_edit(edits, 'M', path_length); // Mismatch.
  }
  else
  {
    if(path_length > 0) { append_edit(edits, 'I', path_length); } // Insertion.
    if(ref_length > 0) { append_edit(edits, 'D', ref_length); } // Deletion.
  }
}

// Returns a CIGAR string.
std::string
align_paths(const GBWTGraph& graph, const gbwt::vector_type& reference, const gbwt::vector_type& path)
{
  // TODO: This uses the basic quadratic longest common subsequence algorithm
  // over the node sequences. If the subgraph is large, we could use something
  // more efficient.
  std::vector<std::vector<std::uint32_t>> dp(path.size() + 1, std::vector<std::uint32_t>(reference.size() + 1, 0));
  for(size_t path_offset = 0; path_offset < path.size(); path_offset++)
  {
    for(size_t ref_offset = 0; ref_offset < reference.size(); ref_offset++)
    {
      if(path[path_offset] == reference[ref_offset])
      {
        dp[path_offset + 1][ref_offset + 1] = dp[path_offset][ref_offset] + 1;
      }
      else
      {
        dp[path_offset + 1][ref_offset + 1] = std::max(dp[path_offset][ref_offset + 1], dp[path_offset + 1][ref_offset]);
      }
    }
  }

  // Trace back the LCS.
  std::vector<std::pair<size_t, size_t>> lcs;
  size_t path_offset = path.size(), ref_offset = reference.size();
  while(path_offset > 0 && ref_offset > 0)
  {
    if(path[path_offset - 1] == reference[ref_offset - 1])
    {
      lcs.push_back(std::make_pair(path_offset - 1, ref_offset - 1));
      path_offset--; ref_offset--;
    }
    else if(dp[path_offset - 1][ref_offset] > dp[path_offset][ref_offset - 1])
    {
      path_offset--;
    }
    else
    {
      ref_offset--;
    }
  }
  std::reverse(lcs.begin(), lcs.end());

  // Convert the LCS a sequence of edits.
  std::vector<std::pair<char, size_t>> edits;
  path_offset = 0; ref_offset = 0;
  for(const std::pair<size_t, size_t>& match : lcs)
  {
    append_edits(edits, graph, reference, std::make_pair(ref_offset, match.second), path, std::make_pair(path_offset, match.first));
    append_edit(edits, 'M', graph.get_length(GBWTGraph::node_to_handle(path[match.first])));
    path_offset = match.first + 1;
    ref_offset = match.second + 1;
  }
  append_edits(edits, graph, reference, std::make_pair(ref_offset, reference.size()), path, std::make_pair(path_offset, path.size()));

  // Convert the edits to a CIGAR string.
  std::string result;
  for(const std::pair<char, size_t>& edit : edits)
  {
    result += std::to_string(edit.second) + edit.first;
  }
  return result;
}

//------------------------------------------------------------------------------

/*
  Private member functions used in the constructor.
*/

void
Subgraph::extract_paths(const GBZ& gbz, size_t query_offset, const std::pair<pos_t, gbwt::edge_type>& ref_pos)
{
  const GBWTGraph& graph = gbz.graph;
  const gbwt::GBWT& index = gbz.index;

  // For each GBWT node in the subgraph, mark which GBWT offsets have a
  // predecessor in the subgraph.
  std::unordered_map<gbwt::node_type, sdsl::bit_vector> has_predecessor;
  has_predecessor.reserve(this->nodes.size() * 2);
  for(nid_t node : this->nodes)
  {
    for(bool is_reverse : { false, true })
    {
      gbwt::node_type gbwt_node = gbwt::Node::encode(node, is_reverse);
      gbwt::CompressedRecord record = index.record(gbwt_node);
      has_predecessor[gbwt_node] = sdsl::bit_vector(record.size(), 0);
    }
  }

  // Mark the positions that have predecessors in the subgraph.
  for(nid_t node : this->nodes)
  {
    for(bool is_reverse : { false, true })
    {
      gbwt::node_type gbwt_node = gbwt::Node::encode(node, is_reverse);
      gbwt::CompressedRecord record = index.record(gbwt_node);
      gbwt::DecompressedRecord decompressed(record);
      for(size_t i = 0; i < decompressed.size(); i++)
      {
        gbwt::edge_type successor = decompressed.LF(i);
        auto iter = has_predecessor.find(successor.first);
        if(iter != has_predecessor.end())
        {
          iter->second[successor.second] = true;
        }
      }
    }
  }

  // FIXME: If the GBWT index is invalid, we could have infinite loops.
  // Extract all paths and determine reference path information if necessary.
  for(const auto& start_node : has_predecessor)
  {
    for(size_t gbwt_offset = 0; gbwt_offset < start_node.second.size(); gbwt_offset++)
    {
      if(start_node.second[gbwt_offset]) { continue; } // There is a predecessor.
      gbwt::edge_type curr(start_node.first, gbwt_offset);
      bool is_ref = false;
      gbwt::vector_type path;
      size_t path_length = 0;
      while(curr.first != gbwt::ENDMARKER)
      {
        if(curr == ref_pos.second)
        {
          this->reference_path = this->paths.size();
          this->reference_start = query_offset - path_length - offset(ref_pos.first);
          is_ref = true;
        }
        path.push_back(curr.first);
        path_length += graph.get_length(GBWTGraph::node_to_handle(curr.first));
        curr = index.LF(curr);
      }
      if(is_ref)
      {
        if(!path_is_canonical(path))
        {
          // TODO: Do we want this?
          std::cerr << "Subgraph::Subgraph(): Warning: Reference path is not in canonical orientation" << std::endl;
        }
        this->paths.push_back(path);
        this->path_lengths.push_back(path_length);
      }
      else if(path_is_canonical(path))
      {
        this->paths.push_back(path);
        this->path_lengths.push_back(path_length);
      }
    }
  }

  if(ref_pos.second != gbwt::invalid_edge() && this->reference_path == std::numeric_limits<size_t>::max())
  {
    std::string msg = "Subgraph::Subgraph(): Reference path not found in the subgraph";
    throw std::runtime_error(msg);
  }
}

void
Subgraph::update_paths(const SubgraphQuery& query)
{
  if(query.haplotype_output == SubgraphQuery::HaplotypeOutput::all_haplotypes) { return; }

  if(query.haplotype_output == SubgraphQuery::HaplotypeOutput::reference_only)
  {
    std::vector<gbwt::vector_type> new_paths { this->paths[this->reference_path] };
    this->paths.swap(new_paths);
    this->reference_path = 0;
    return;
  }

  // TODO: This could be more efficient by using normal GBWT operations.
  if(query.haplotype_output == SubgraphQuery::HaplotypeOutput::distinct_haplotypes)
  {
    gbwt::vector_type ref_path = this->paths[this->reference_path];
    std::sort(this->paths.begin(), this->paths.end());
    std::vector<gbwt::vector_type> new_paths;
    std::vector<size_t> new_lengths;
    for(size_t i = 0; i < this->paths.size(); i++)
    {
      const gbwt::vector_type& path = this->paths[i];
      if(new_paths.empty() || path != new_paths.back())
      {
        if(path == ref_path) { this->reference_path = new_paths.size(); }
        new_paths.push_back(path);
        new_lengths.push_back(this->path_lengths[i]);
        this->path_weights.push_back(1);
      }
      else { this->path_weights.back()++; }
    }
    this->paths.swap(new_paths);
    this->path_lengths.swap(new_lengths);
  }
}

//------------------------------------------------------------------------------

Subgraph::Subgraph(const GBZ& gbz, const PathIndex* path_index, const SubgraphQuery& query) :
  reference_path(std::numeric_limits<size_t>::max()),
  reference_handle(handlegraph::as_path_handle(std::numeric_limits<size_t>::max())),
  reference_start(0)
{
  std::pair<pos_t, gbwt::edge_type> position;
  size_t context = query.context;

  switch(query.type)
  {
  case SubgraphQuery::QueryType::path_offset_query:
    {
      if(path_index == nullptr)
      {
        throw std::runtime_error("Subgraph::Subgraph(): Path index required for path queries");
      }
      position = find_position(gbz, *path_index, query);
      this->reference_handle = query.path;
      break;
    }
  case SubgraphQuery::QueryType::node_query:
    {
      // TODO: This has the same issue as in GBZ-base. If node length is even,
      // the right context is one bp too long. What about odd lengths?
      size_t node_len = gbz.graph.get_length(gbz.graph.get_handle(query.node_id));
      position.first = make_pos_t(query.node_id, false, node_len / 2);
      position.second = gbwt::invalid_edge();
      context += node_len / 2;
      break;
    }
  default:
    {
      throw std::runtime_error("Subgraph::Subgraph(): Invalid query type");
    }
  }

  this->nodes = find_subgraph(gbz.graph, position.first, context);
  this->extract_paths(gbz, query.offset, position);
  this->update_paths(query);

  if(this->reference_path < this->paths.size())
  {
    for(size_t i = 0; i < this->paths.size(); i++)
    {
      if(i == this->reference_path) { this->path_cigars.push_back(""); }
      else { this->path_cigars.push_back(align_paths(gbz.graph, this->paths[this->reference_path], this->paths[i])); }
    }
  }
}

//------------------------------------------------------------------------------

gbwt::FullPathName
Subgraph::reference_path_name(const GBZ& gbz) const
{
  gbwt::size_type path_id = gbz.graph.handle_to_path(this->reference_handle);
  gbwt::FullPathName path_name = gbz.index.metadata.full_path(path_id);
  path_name.offset = this->reference_start;
  return path_name;
}

const size_t*
Subgraph::weight(size_t path_id) const
{
  if(path_id >= this->path_weights.size()) { return nullptr; }
  return &(this->path_weights[path_id]);
}

const std::string*
Subgraph::cigar(size_t path_id) const
{
  if(path_id >= this->path_cigars.size()) { return nullptr; }
  if(path_id == this->reference_path) { return nullptr; }
  return &(this->path_cigars[path_id]);
}

void
Subgraph::to_gfa(const GBZ& gbz, std::ostream& out) const
{
  TSVWriter writer(out);

  // GFA header.
  if(this->reference_path < this->paths.size())
  {
    std::string sample_name = gbz.graph.get_sample_name(this->reference_handle);
    write_gfa_header(writer, &sample_name);
  }
  else { write_gfa_header(writer, nullptr); }

  // Nodes / segments.
  for(nid_t node : this->nodes)
  {
    handle_t handle = gbz.graph.get_handle(node);
    write_gfa_segment(writer, node, gbz.graph.get_sequence_view(handle));
  }

  // Edges / links.
  for(nid_t node : this->nodes)
  {
    for(bool is_reverse : { false, true })
    {
      gbwt::node_type from = gbwt::Node::encode(node, is_reverse);
      handle_t handle = GBWTGraph::node_to_handle(from);
      gbz.graph.follow_edges(handle, false, [&](const handle_t& next)
      {
        gbwt::node_type to = GBWTGraph::handle_to_node(next);
        nid_t successor = gbwt::Node::id(to);
        if(this->nodes.find(successor) != this->nodes.end() && edge_is_canonical(from, to))
        {
          write_gfa_link(writer, from, to);
        }
      });
    }
  }

  // W-lines: reference path
  std::string contig_name = "unknown";
  if(this->reference_path < this->paths.size())
  {
    const gbwt::vector_type& path = this->paths[this->reference_path];
    gbwt::FullPathName path_name = this->reference_path_name(gbz);
    contig_name = path_name.contig_name;
    size_t length = this->path_lengths[this->reference_path];
    const size_t* weight = this->weight(this->reference_path);
    const std::string* cigar = this->cigar(this->reference_path);
    write_gfa_walk(writer, path, path_name, length, weight, cigar);
  }

  // W-lines: anonymous haplotypes
  size_t haplotype_id = 1;
  for(size_t i = 0; i < this->paths.size(); i++)
  {
    if(i == this->reference_path) { continue; }
    const gbwt::vector_type& path = this->paths[i];
    gbwt::FullPathName path_name { "unknown", contig_name, haplotype_id, 0 };
    size_t length = this->path_lengths[i];
    const size_t* weight = this->weight(i);
    const std::string* cigar = this->cigar(i);
    write_gfa_walk(writer, path, path_name, length, weight, cigar);
    haplotype_id++;
  }

  writer.flush();
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
