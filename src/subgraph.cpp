#include <gbwtgraph/subgraph.h>
#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/internal.h>

#include <sdsl/bit_vectors.hpp>

#include <queue>
#include <map>

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
SubgraphQuery::path_offset(const gbwt::FullPathName& path_name, size_t offset, size_t context, HaplotypeOutput output)
{
  gbwt::FullPathName path_query = path_name;
  path_query.offset = offset;

  SubgraphQuery result
  {
    QueryType::path_offset_query, output,
    path_query, 0,
    context
  };
  return result;
}

SubgraphQuery
SubgraphQuery::path_interval(const gbwt::FullPathName& path_name, size_t from, size_t to, size_t context, HaplotypeOutput output)
{
  if(from >= to)
  {
    std::string msg = "SubgraphQuery::path_interval(): Empty interval [" + std::to_string(from) + ", " + std::to_string(to) + ")";
    throw std::runtime_error(msg);
  }
  gbwt::FullPathName path_query = path_name;
  path_query.offset = from + (to - from) / 2;

  SubgraphQuery result
  {
    QueryType::path_offset_query, output,
    path_query, 0,
    context
  };
  return result;
}

SubgraphQuery
SubgraphQuery::node(nid_t node, size_t context, HaplotypeOutput output)
{
  gbwt::FullPathName empty { "", "", 0, 0 };
  SubgraphQuery result
  {
    QueryType::node_query, output,
    empty, node,
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
SubgraphQuery::to_string() const
{
  if(this->type == QueryType::path_offset_query)
  {
    std::string path_name =
      "sample " + this->path_query.sample_name +
      ", contig " + this->path_query.contig_name +
      ", haplotype " + std::to_string(this->path_query.haplotype);
    std::string output = output_type(this->haplotype_output);
    return "(" + path_name + ", offset " + std::to_string(this->path_query.offset) + ", context " + std::to_string(this->context) + ", " + output + ")";
  }
  else if(this->type == QueryType::node_query)
  {
    std::string output = output_type(this->haplotype_output);
    return "(node " + std::to_string(this->node_id) + ", context " + std::to_string(this->context) + ", " + output + ")";
  }
  else
  {
    // NOTE: `path_interval_query` is currently implemented as a `path_offset_query` around the midpoint.
    return "(invalid query)";
  }
}

//------------------------------------------------------------------------------

/*
  Helper functions for the constructor.
*/

std::pair<pos_t, gbwt::edge_type>
find_position(const GBZ& gbz, const PathIndex& path_index, path_handle_t path, size_t offset)
{
  std::pair<size_t, gbwt::edge_type> position = path_index.sampled_position(path, offset);
  if(position.second == gbwt::invalid_edge())
  {
    std::string msg = "Subgraph::Subgraph(): Could not find an indexed path position for path " +
      gbz.graph.get_path_name(path) + ", offset " + std::to_string(offset);
    throw std::runtime_error(msg);
  }

  while(true)
  {
    size_t node_len = gbz.graph.get_length(GBWTGraph::node_to_handle(position.second.first));
    if(position.first + node_len > offset)
    {
      pos_t graph_pos = make_pos_t(gbwt::Node::id(position.second.first), gbwt::Node::is_reverse(position.second.first), offset - position.first);
      return std::make_pair(graph_pos, position.second);
    }
    position.first += node_len;
    position.second = gbz.index.LF(position.second);
    if(position.second.first == gbwt::ENDMARKER)
    {
      std::string msg = "Subgraph::Subgraph(): Path " + gbz.graph.get_path_name(path) + " does not contain offset " + std::to_string(offset);
      throw std::runtime_error(msg);
    }
  }
}

std::set<nid_t>
find_subgraph(const GBWTGraph& graph, pos_t position, size_t context)
{
  std::set<nid_t> result;

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
    result.insert(id(current.second));
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
    }
  }

  return result;
}

//------------------------------------------------------------------------------

/*
  Helper functions for determining CIGAR strings relative to the reference.
*/

void
append_edit(std::vector<std::pair<char, size_t>>& edits, char op, size_t length)
{
  if(length == 0) { return; }
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
subpath_length(const GBWTGraph& graph, subpath_type subpath)
{
  size_t result = 0;
  for(size_t i = 0; i < subpath.second; i++)
  {
    result += graph.get_length(GBWTGraph::node_to_handle(subpath.first[i]));
  }
  return result;
}

// Returns the number of matching bases at the start of the two subpaths.
size_t
prefix_matches(const GBWTGraph& graph, subpath_type a, subpath_type b)
{
  size_t result = 0;
  size_t a_offset = 0, b_offset = 0;
  size_t a_base = 0, b_base = 0;
  while(a_offset < a.second && b_offset < b.second)
  {
    std::string_view a_seq = graph.get_sequence_view(GBWTGraph::node_to_handle(a.first[a_offset]));
    std::string_view b_seq = graph.get_sequence_view(GBWTGraph::node_to_handle(b.first[b_offset]));
    while(a_base < a_seq.size() && b_base < b_seq.size())
    {
      if(a_seq[a_base] != b_seq[b_base]) { return result; }
      result++;
      a_base++; b_base++;
    }
    if(a_base == a_seq.size()) { a_offset++; a_base = 0; }
    if(b_base == b_seq.size()) { b_offset++; b_base = 0; }
  }
  return result;
}

// Returns the number of matching bases at the end of the two subpaths.
size_t
suffix_matches(const GBWTGraph& graph, subpath_type a, subpath_type b)
{
  size_t result = 0;
  size_t a_offset = a.second, b_offset = b.second;
  size_t a_base = 0, b_base = 0;
  while(a_offset > 0 && b_offset > 0)
  {
    std::string_view a_seq = graph.get_sequence_view(GBWTGraph::node_to_handle(a.first[a_offset - 1]));
    std::string_view b_seq = graph.get_sequence_view(GBWTGraph::node_to_handle(b.first[b_offset - 1]));
    while(a_base < a_seq.size() && b_base < b_seq.size())
    {
      if(a_seq[a_seq.size() - 1 - a_base] != b_seq[b_seq.size() - 1 - b_base]) { return result; }
      result++;
      a_base++; b_base++;
    }
    if(a_base == a_seq.size()) { a_offset--; a_base = 0; }
    if(b_base == b_seq.size()) { b_offset--; b_base = 0; }
  }
  return result;
}

size_t
mismatch_penalty(size_t length)
{
  return 4 * length;
}

size_t
gap_penalty(size_t length)
{
  if(length == 0) { return 0; }
  return 6 + (length - 1);
}

void
align_diverging(
  const GBWTGraph& graph,
  subpath_type path, subpath_type reference,
  std::vector<std::pair<char, size_t>>& edits
)
{
  size_t path_length = subpath_length(graph, path);
  size_t ref_length = subpath_length(graph, reference);
  size_t prefix = prefix_matches(graph, reference, path);
  size_t suffix = suffix_matches(graph, reference, path);
  if(prefix + suffix > path_length) { suffix = path_length - prefix; }
  if(prefix + suffix > ref_length) { suffix = ref_length - prefix; }

  append_edit(edits, 'M', prefix); // Match.
  size_t path_middle = path_length - prefix - suffix;
  size_t ref_middle = ref_length - prefix - suffix;
  if(path_middle == 0)
  {
    append_edit(edits, 'D', ref_middle); // Deletion.
  }
  else if(ref_middle == 0)
  {
    append_edit(edits, 'I', path_middle); // Insertion.
  }
  else
  {
    size_t mismatch = std::min(ref_middle, path_middle);
    size_t mismatch_indel =
      mismatch_penalty(mismatch) +
      gap_penalty(ref_middle - mismatch) +
      gap_penalty(path_middle - mismatch);
    size_t insertion_deletion = gap_penalty(ref_middle) + gap_penalty(path_middle);
    if(mismatch_indel <= insertion_deletion)
    {
      append_edit(edits, 'M', mismatch); // Mismatch.
      append_edit(edits, 'I', path_middle - mismatch); // Insertion.
      append_edit(edits, 'D', ref_middle - mismatch); // Deletion.
    }
    else
    {
      append_edit(edits, 'I', path_middle); // Insertion.
      append_edit(edits, 'D', ref_middle); // Deletion.
    }
  }
  append_edit(edits, 'M', suffix); // Match.
}

std::string
to_cigar(const std::vector<std::pair<char, size_t>>& edits)
{
  std::string result;
  for(std::pair<char, size_t> edit : edits)
  {
    result += std::to_string(edit.second);
    result.push_back(edit.first);
  }
  return result;
}

// Returns a CIGAR string.
std::string
align_path(const GBWTGraph& graph, const gbwt::vector_type& reference, const gbwt::vector_type& path)
{
  // Find the LCS weighted by sequence length.
  // The values are (path offset, reference offset).
  std::vector<std::pair<size_t, size_t>> lcs = path_lcs(graph, path, reference);

  // Convert the LCS to a sequence of edits.
  std::vector<std::pair<char, size_t>> edits;
  size_t path_offset = 0, ref_offset = 0;
  for(const std::pair<size_t, size_t>& match : lcs)
  {
    align_diverging(
      graph,
      get_subpath(path, path_offset, match.first),
      get_subpath(reference, ref_offset, match.second),
      edits
    );
    append_edit(edits, 'M', graph.get_length(GBWTGraph::node_to_handle(path[match.first])));
    path_offset = match.first + 1;
    ref_offset = match.second + 1;
  }
  align_diverging(
    graph,
    get_subpath(path, path_offset, path.size()),
    get_subpath(reference, ref_offset, reference.size()),
    edits
  );

  return to_cigar(edits);
}

//------------------------------------------------------------------------------

/*
  Private member functions used in the constructor.
*/

void
Subgraph::extract_paths(const GBZ& gbz, const SubgraphQuery& query, size_t query_offset, const std::pair<pos_t, gbwt::edge_type>& ref_pos)
{
  const GBWTGraph& graph = gbz.graph;
  const gbwt::GBWT& index = gbz.index;

  // For each GBWT node in the subgraph, mark which GBWT offsets have a
  // predecessor in the subgraph.
  std::map<gbwt::node_type, sdsl::bit_vector> has_predecessor;
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
      while(curr.first != gbwt::ENDMARKER && has_predecessor.find(curr.first) != has_predecessor.end())
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
          std::cerr << "Subgraph::Subgraph(): Warning: Reference path for query " << query.to_string() << " is not in canonical orientation" << std::endl;
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
    gbwt::vector_type ref_path = (this->reference_path < this->paths.size() ? this->paths[this->reference_path] : gbwt::vector_type(1, gbwt::ENDMARKER));
    std::vector<size_t> order; order.reserve(this->paths.size());
    for(size_t i = 0; i < this->paths.size(); i++) { order.push_back(i); }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) -> bool
    {
      return (this->paths[a] < this->paths[b]);
    });
    std::vector<gbwt::vector_type> new_paths;
    std::vector<size_t> new_lengths;
    for(size_t i = 0; i < this->paths.size(); i++)
    {
      const gbwt::vector_type& path = this->paths[order[i]];
      if(new_paths.empty() || path != new_paths.back())
      {
        if(path == ref_path) { this->reference_path = new_paths.size(); }
        new_paths.push_back(path);
        new_lengths.push_back(this->path_lengths[order[i]]);
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
  reference_start(0)
{
  std::pair<pos_t, gbwt::edge_type> position;
  size_t context = query.context;
  size_t query_offset = 0;

  switch(query.type)
  {
  case SubgraphQuery::QueryType::path_offset_query:
    {
      if(path_index == nullptr)
      {
        throw std::runtime_error("Subgraph::Subgraph(): Path index required for path queries");
      }
      const gbwt::Metadata& metadata = gbz.index.metadata;
      gbwt::size_type path_id = metadata.findFragment(query.path_query);
      if(path_id >= metadata.paths())
      {
        std::string msg =
          "Subgraph::Subgraph(): Could not find a path fragment for sample " + query.path_query.sample_name +
          ", contig " + query.path_query.contig_name +
          ", haplotype " + std::to_string(query.path_query.haplotype) +
          ", offset " + std::to_string(query.path_query.offset);
        throw std::runtime_error(msg);
      }
      path_handle_t reference_handle = gbz.graph.path_to_handle(path_id);
      this->reference_name = metadata.fullPath(path_id);
      query_offset = query.path_query.offset - this->reference_name.offset;
      position = find_position(gbz, *path_index, reference_handle, query_offset);
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
      // NOTE: `path_interval_query` is currently implemented as a `path_offset_query` around the midpoint.
      throw std::runtime_error("Subgraph::Subgraph(): Invalid query type");
    }
  }

  this->nodes = find_subgraph(gbz.graph, position.first, context);
  this->extract_paths(gbz, query, query_offset, position);
  this->update_paths(query);

  if(this->reference_path < this->paths.size())
  {
    for(size_t i = 0; i < this->paths.size(); i++)
    {
      if(i == this->reference_path) { this->path_cigars.push_back(""); }
      else { this->path_cigars.push_back(align_path(gbz.graph, this->paths[this->reference_path], this->paths[i])); }
    }
  }
}

//------------------------------------------------------------------------------

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
    write_gfa_header(writer, &(this->reference_name.sample_name));
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
    gbwt::FullPathName path_name = this->reference_name;
    // Offset is for the path relative to the entire haplotype.
    // We transform it into the starting offset of the path fragment within the subgraph relative to the entire haplotype.
    path_name.offset += this->reference_start;
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
