#include <gbwtgraph/gbwtgraph.h>

#include <algorithm>
#include <queue>
#include <stack>
#include <string>
#include <utility>

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
constexpr std::uint64_t GBWTGraph::Header::FLAG_SIMPLE_SDS;

constexpr std::uint32_t GBWTGraph::Header::TRANS_VERSION;
constexpr std::uint64_t GBWTGraph::Header::TRANS_FLAG_MASK;

constexpr std::uint32_t GBWTGraph::Header::OLD_VERSION;
constexpr std::uint64_t GBWTGraph::Header::OLD_FLAG_MASK;

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

void
GBWTGraph::Header::check() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("GBWTGraph: Invalid tag");
  }

  if(this->version > VERSION || this->version < OLD_VERSION)
  {
    std::string msg = "GBWTGraph: Expected v" + std::to_string(OLD_VERSION) + " to v" + std::to_string(VERSION) + ", got v" + std::to_string(this->version);
    throw sdsl::simple_sds::InvalidData(msg);
  }

  std::uint64_t mask = 0;
  switch(this->version)
  {
  case VERSION:
    mask = FLAG_MASK; break;
  case TRANS_VERSION:
    mask = TRANS_FLAG_MASK; break;
  case OLD_VERSION:
    mask = OLD_FLAG_MASK; break;
  }
  if((this->flags & mask) != this->flags)
  {
    throw sdsl::simple_sds::InvalidData("GBWTGraph: Invalid flags");
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
  this->named_paths.swap(another.named_paths);
  this->name_to_path.swap(another.name_to_path);
  this->id_to_path.swap(another.id_to_path);
  this->reference_samples.swap(another.reference_samples);
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
    this->named_paths = std::move(source.named_paths);
    this->name_to_path = std::move(source.name_to_path);
    this->id_to_path = std::move(source.id_to_path);
    this->reference_samples = std::move(source.reference_samples);
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
  this->named_paths = source.named_paths;
  this->name_to_path = source.name_to_path;
  this->id_to_path = source.id_to_path;
  this->reference_samples = source.reference_samples;
}

void
GBWTGraph::sanity_checks()
{
  size_t nodes = sdsl::util::cnt_one_bits(this->real_nodes);
  if(nodes != this->header.nodes)
  {
    throw sdsl::simple_sds::InvalidData("GBWTGraph: Invalid number of set bits in real_nodes");
  }

  size_t potential_nodes = this->sequences.size();
  if(this->index != nullptr && !(this->index->empty())) { potential_nodes =  this->index->sigma() - this->index->firstNode(); }
  if(this->sequences.size() != potential_nodes)
  {
    throw sdsl::simple_sds::InvalidData("GBWTGraph: Node range / sequence count mismatch");
  }
  if(this->real_nodes.size() != potential_nodes / 2)
  {
    throw sdsl::simple_sds::InvalidData("GBWTGraph: Node range / real_nodes size mismatch");
  }

  if(this->node_to_segment.ones() != this->segments.size())
  {
    throw sdsl::simple_sds::InvalidData("GBWTGraph: Segment count / node_to_segment mapping mismatch");
  }
  if(this->segments.size() > 0 && this->index != nullptr && this->node_to_segment.size() != this->index->sigma() / 2)
  {
    throw sdsl::simple_sds::InvalidData("GBWTGraph: GBWT alphabet / node_to_segment size mismatch");
  }
}

std::pair<gbwt::StringArray, sdsl::sd_vector<>>
GBWTGraph::copy_translation(const NamedNodeBackTranslation& translation) const
{
  // This will hold an array of segment name strings, and a bit vector of
  // flags, where the flag is set at node IDs that
  // are the first in their segments.
  // The representation assumes that all nodes in a segment have successive IDs.
  std::pair<gbwt::StringArray, sdsl::sd_vector<>> result;

  // NamedNodeBackTranslation copies out segment names, but to make a StringArray we need a place to store them.
  // We also need to know the first node ID in every segment, for the bit vector, but we need the total set bit count to start the builder.
  // So we store pairs of segment names and first node IDs, and end up re-numbering our segments by rank here.
  std::vector<std::pair<std::string, nid_t>> segment_names_and_starts;

  // We need to track consistency with our contiguous node ID range segment model.
  nid_t prev_segment_number = std::numeric_limits<nid_t>::max();
  bool prev_node_existed = false;
  size_t next_offset_along_segment = std::numeric_limits<size_t>::max();

  if(this->index->firstNode() > 1)
  {
    // Have a segment for the first run of nonexistent nodes
    segment_names_and_starts.emplace_back("", 1);
  }

  for(gbwt::node_type node = this->index->firstNode(); node < this->index->sigma(); node += 2)
  {
    // Get each node in the graph, in ID order.
    nid_t node_id = gbwt::Node::id(node);
    if(!(this->has_node(node_id)))
    {
      // This node doesn't exist.
      if(prev_node_existed)
      {
        // We're starting a nonexistent segment.
        // We need to remember its empty name and start ID.
        segment_names_and_starts.emplace_back("", node_id);
        // Remember we were in a nonexistent segment.
        prev_node_existed = false;
      }
    }
    else
    {
      // Translate it back to a segment range.
      size_t node_length = this->sequences.length(this->node_offset(node));
      oriented_node_range_t range{node_id, false, 0, node_length};
      auto translated_back = translation.translate_back(range);
      if(translated_back.size() != 1)
      {
        // This node didn't come from a segment, or spans multiple segments
        throw InvalidGBWT("GBWTGraph: Node " + std::to_string(node_id) + " did not come from exactly one segment");
      }
      nid_t& segment_number = std::get<0>(translated_back[0]);
      bool& reverse_on_segment = std::get<1>(translated_back[0]);
      size_t& offset_along_segment = std::get<2>(translated_back[0]);
      if(reverse_on_segment)
      {
        // We can only deal with nodes on the forward strands of their segments.
        throw InvalidGBWT("GBWTGraph: Node " + std::to_string(node_id) + " came from the reverse strand of its segment");
      }
      if(!prev_node_existed || prev_segment_number != segment_number)
      {
        // This is a new segment!
        prev_segment_number = segment_number;
        prev_node_existed = true;
        next_offset_along_segment = 0;
        // We need to remember its name and start ID.
        segment_names_and_starts.emplace_back(translation.get_back_graph_node_name(segment_number), node_id);
      }
      if(offset_along_segment != next_offset_along_segment)
      {
        // Actually we're not at the right place in the segment, so we can't store this translation.
        throw InvalidGBWT("GBWTGraph: Node " + std::to_string(node_id) + " not at expected position in segment");
      }
      next_offset_along_segment += node_length;
    }
  }

  // Store the segment names.
  std::string empty;
  result.first = gbwt::StringArray(segment_names_and_starts.size(),
  [&](size_t offset) -> size_t
  {
    // This produces the length of each string to store
    return segment_names_and_starts[offset].first.size();
  },
  [&](size_t offset) -> view_type
  {
    // This produces a view to each string to store.
    return str_to_view(segment_names_and_starts[offset].first);
  });

  // Store the mapping.
  sdsl::sd_vector_builder builder(gbwt::Node::id(this->index->sigma()), segment_names_and_starts.size());
  for(auto& record : segment_names_and_starts)
  {
    builder.set_unsafe(record.second);
  }
  result.second = sdsl::sd_vector<>(builder);

  return result;
}

//------------------------------------------------------------------------------

GBWTGraph::GBWTGraph(const gbwt::GBWT& gbwt_index,
                     const HandleGraph& sequence_source,
                     const NamedNodeBackTranslation* segment_space) :
  index(nullptr)
{
  // Set GBWT, cache named paths, and do sanity checks.
  this->set_gbwt(gbwt_index);
  if(this->index->empty()) { return; }

  // Build real_nodes to support has_node().
  this->determine_real_nodes();

  // Store the sequences.
  this->sequences = gbwt::StringArray(this->index->sigma() - this->index->firstNode(),
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

  // Store the node to segment translation
  if(segment_space)
  {
    this->header.set(Header::FLAG_TRANSLATION);
    std::tie(this->segments, this->node_to_segment) =
      this->copy_translation(*segment_space);
  }
}

GBWTGraph::GBWTGraph(const gbwt::GBWT& gbwt_index, const SequenceSource& sequence_source) :
  index(nullptr)
{
  // Set GBWT, cache named paths, and do sanity checks.
  this->set_gbwt(gbwt_index);
  if(this->index->empty()) { return; }

  // Build real_nodes to support has_node().
  this->determine_real_nodes();

  // Store the sequences.
  this->sequences = gbwt::StringArray(this->index->sigma() - this->index->firstNode(),
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

  // Store the node to segment translation but leave the names of unused segments empty.
  if(sequence_source.uses_translation())
  {
    this->header.set(Header::FLAG_TRANSLATION);
    std::tie(this->segments, this->node_to_segment) =
      sequence_source.invert_translation([&](std::pair<nid_t, nid_t> interval) -> bool
    {
      return this->has_node(interval.first);
    });
  }
}

//------------------------------------------------------------------------------

void
GBWTGraph::determine_real_nodes()
{
  // Sometimes (e.g. in `deserialize()`) we call this function after the header already
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

void
GBWTGraph::cache_named_paths()
{
  this->named_paths.clear();
  this->name_to_path.clear();
  this->id_to_path.clear();

  // There cannot be named paths without sufficient metadata.
  if(this->index == nullptr || !(this->index->hasMetadata()) ||
    !(this->index->metadata.hasSampleNames()) ||
    !(this->index->metadata.hasContigNames()) ||
    !(this->index->metadata.hasPathNames()))
  {
    return;
  }

  // Determine the named paths.
  for(gbwt::size_type sample = 0; sample < this->index->metadata.sample_names.size(); sample++)
  {
    // TODO: It should be efficient to find all the names in the dictionary
    // with the given prefix, but the dictionary doesn't expose that method.
    // When GBWT Dictionaries get that functionality, use it here.
    // For now just scan all sample names.
    std::string sample_name = this->index->metadata.sample(sample);
    if(sample_name.size() >= REFERENCE_PATH_SAMPLE_NAME.size() &&
       std::equal(REFERENCE_PATH_SAMPLE_NAME.begin(), REFERENCE_PATH_SAMPLE_NAME.end(), sample_name.begin()))
    {
      // This sample name starts with the prefix, so index it.

      // Determine what sample name we should present to the user. Either part
      // after the prefix, or the no-name sentinel
      std::string exposed_sample_name;
      // This also determines the sense of the path
      PathSense sense;
      if(sample_name.size() > REFERENCE_PATH_SAMPLE_NAME.size())
      {
        exposed_sample_name = sample_name.substr(REFERENCE_PATH_SAMPLE_NAME.size());
        sense = PathSense::REFERENCE;
      }
      else
      {
        exposed_sample_name = NO_SAMPLE_NAME;
        sense = PathSense::GENERIC;
      }

      std::vector<gbwt::size_type> sample_ids = this->index->metadata.pathsForSample(sample);
      for(size_t i = 0; i < sample_ids.size(); i++)
      {
        // For each path in this reference sample
        const gbwt::PathName& path = this->index->metadata.path(sample_ids[i]);

        // Compose a name for it using the sense, sample and conting name we want to present.
        // Make sure to only send a haplotype number if we can handle one.
        // And make sure to send a subrange start if we are hiding one in the count.
        std::string composed_path_name = PathMetadata::create_path_name(
          sense,
          exposed_sample_name,
          this->index->metadata.contig(path.contig),
          sense == PathSense::REFERENCE ? path.phase : NO_HAPLOTYPE,
          NO_PHASE_BLOCK,
          path.count == 0 ? NO_SUBRANGE : subrange_t(path.count, NO_END_POSITION)
        );

        // Store a mapping from name to index this will appear at in named_paths
        this->name_to_path[composed_path_name] = named_paths.size();
        // And from internal path ID to index this will appear at in named_paths;
        this->id_to_path[sample_ids[i]] = named_paths.size();

        // Create the NamedPath
        this->named_paths.emplace_back();
        // And store the ID
        this->named_paths.back().id = sample_ids[i];
        // And whether we're reference sense
        this->named_paths.back().sense = sense;
      }
    }
  }
  if(this->name_to_path.size() != this->named_paths.size())
  {
     throw InvalidGBWT("GBWTGraph: Named path names are not unique");
  }

  // Cache named path information we get from traversing the paths.
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t i = 0; i < this->named_paths.size(); i++)
  {
    NamedPath& path = this->named_paths[i];
    gbwt::edge_type curr = this->index->start(gbwt::Path::encode(path.id, false));
    path.from = (curr.first == gbwt::ENDMARKER ? gbwt::invalid_edge() : curr);
    path.to = gbwt::invalid_edge();
    path.length = 0;
    while(curr.first != gbwt::ENDMARKER)
    {
      path.to = curr;
      path.length++;
      curr = this->index->LF(curr);
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

size_t
GBWTGraph::get_path_count() const
{
  return this->named_paths.size();
}

bool
GBWTGraph::has_path(const std::string& path_name) const
{
  // Polling for path names that look like haplotypes can be kind of hard, so
  // we get the handle and look for the sentinel.
  return handlegraph::as_integer(this->get_path_handle(path_name)) != std::numeric_limits<size_t>::max();
}

path_handle_t
GBWTGraph::get_path_handle(const std::string& path_name) const
{
  // We define this to return a sentinel path handle when no path is found
  path_handle_t to_return = handlegraph::as_path_handle(std::numeric_limits<size_t>::max());

  auto found = this->name_to_path.find(path_name);
  if(found != this->name_to_path.end())
  {
    // This is the name of a stored named path.
    // Pack up its cache offset as a handle.
    to_return = handlegraph::as_path_handle(found->second);
  }
  else
  {
    // Parse out the path name.
    PathSense sense;
    std::string sample_name;
    std::string contig_name;
    size_t haplotype;
    size_t phase_block;
    subrange_t subrange;
    PathMetadata::parse_path_name(path_name,
                                  sense,
                                  sample_name,
                                  contig_name,
                                  haplotype,
                                  phase_block,
                                  subrange);

    if(sample_name == REFERENCE_PATH_SAMPLE_NAME ||
       this->reference_samples.count(sample_name))
    {
      // We aren't allowed to expose named paths through this mechanism.
      return to_return;
    }

    if(subrange != NO_SUBRANGE)
    {
      // We don't store haplotype subranges.
      return to_return;
    }

    if(haplotype == NO_HAPLOTYPE)
    {
      // We need a haplotype.
      return to_return;
    }

    if(phase_block == NO_PHASE_BLOCK)
    {
      // We need a phase block.
      return to_return;
    }

    auto sample_number = this->index->metadata.sample(sample_name);
    if(sample_number == this->index->metadata.sample_names.size())
    {
      // This sample doesn't exist.
      return to_return;
    }

    auto contig_number = this->index->metadata.contig(contig_name);
    if(contig_number == this->index->metadata.contig_names.size())
    {
      // This contig doesn't exist.
      return to_return;
    }

    // Now we have to look up this path in the metadata and get the path number.
    for(auto& path_id : this->index->metadata.findPaths(sample_number, contig_number))
    {
      // Paths are only indexed by sample and contig, so we have to scan for the
      // right haplotype and phase block.
      auto& structured_name = this->index->metadata.path(path_id);
      if(structured_name.phase == haplotype && structured_name.count == phase_block)
      {
        // This is the right path. Turn it into a haplotype path handle.
        to_return = handlegraph::as_path_handle(this->named_paths.size() + path_id);
        break;
      }
    }
  }
  // Now return what we found, or the sentinel if we found nothing.
  return to_return;
}

std::string
GBWTGraph::get_path_name(const path_handle_t& path_handle) const
{

  auto sense = this->get_sense(path_handle);
  gbwt::size_type path_id = this->get_metadata_index(path_handle);
  // Get the name fields from the metadata.
  auto& structured_name = this->index->metadata.path(path_id);
  return gbwtgraph::compose_path_name(this->index->metadata, structured_name, sense);
}

bool
GBWTGraph::get_is_circular(const path_handle_t&) const
{
  // TODO: We don't track circular paths
  return false;
}

size_t
GBWTGraph::get_step_count(const path_handle_t& path_handle) const
{
  switch(this->get_sense(path_handle))
  {
  case PathSense::GENERIC: // Fall-through
  case PathSense::REFERENCE:
    // This information is cached
    return this->named_paths[handlegraph::as_integer(path_handle)].length;;
    break;
  case PathSense::HAPLOTYPE:
    // This information is not cached.
    {
      size_t count = 0;
      for(step_handle_t here = this->path_begin(path_handle);
          here != this->path_end(path_handle);
          here = this->get_next_step(here))
      {
        // Trace the path and count all the steps.
        count++;
      }
      return count;
    }
    break;
  default:
    throw std::runtime_error("Unimplemented sense!");
  }

}

size_t
GBWTGraph::get_step_count(const handle_t& handle) const
{
  // Use the brute force approach where we total up the number of step handles
  // we iterate over.
  size_t count = 0;
  for_each_step_on_handle_impl(handle, [&](const step_handle_t&)
  {
    count++;
    return true;
  });
  return count;
}

handle_t
GBWTGraph::get_handle_of_step(const step_handle_t& step_handle) const {
  // The GBWT node number is the first field of the step handle, so grab that
  // and turn it into a graph handle.
  return node_to_handle(handlegraph::as_integers(step_handle)[0]);
}

path_handle_t
GBWTGraph::get_path_handle_of_step(const step_handle_t& step_handle) const {
  // To find the thread number we will need to locate the selected visit. So
  // turn it into an edge.
  gbwt::edge_type here;
  here.first = handlegraph::as_integers(step_handle)[0];
  here.second = handlegraph::as_integers(step_handle)[1];
  gbwt::size_type path_id = gbwt::Path::id(index->locate(here));

  // Convert path id number to path handle.
  return this->from_metadata_index(path_id);
}

step_handle_t
GBWTGraph::path_begin(const path_handle_t& path_handle) const {
  // Step handles correspond to GBWT edges.
  gbwt::edge_type from;

  switch(this->get_sense(path_handle))
  {
  case PathSense::GENERIC: // Fall-through
  case PathSense::REFERENCE:
    // This information is cached
    from = this->named_paths[handlegraph::as_integer(path_handle)].from;
    break;
  case PathSense::HAPLOTYPE:
    {
      size_t path_number = this->get_metadata_index(path_handle);
      // We can use the gbwt's start() to get the start of a GBWT sequence.
      // And we're always interested in the forward orientation of the path
      from = this->index->start(gbwt::Path::encode(path_number, false));

      if(from.first == gbwt::ENDMARKER)
      {
        // Path must be empty. Use the past-end sentinel.
        from = gbwt::invalid_edge();
      }
    }
    break;
  default:
    throw std::runtime_error("Unimplemented sense!");
  }

  step_handle_t step;
  // Edges and steps have GBWT node numbers first
  as_integers(step)[0] = from.first;
  // And then offsets into the node's visits
  as_integers(step)[1] = from.second;

  return step;
}

step_handle_t
GBWTGraph::path_end(const path_handle_t&) const {
  // path_end can just be invalid_edge() since it's a past-end.
  gbwt::edge_type past_last_edge = gbwt::invalid_edge();

  step_handle_t step;
  // Edges and steps have GBWT node numbers first
  as_integers(step)[0] = past_last_edge.first;
  // And then offsets into the node's visits
  as_integers(step)[1] = past_last_edge.second;

  return step;
}

step_handle_t
GBWTGraph::path_back(const path_handle_t& path_handle) const {
  gbwt::edge_type to;

  switch(this->get_sense(path_handle))
  {
  case PathSense::GENERIC: // Fall-through
  case PathSense::REFERENCE:
    // This information is cached
    to = this->named_paths[handlegraph::as_integer(path_handle)].to;
    break;
  case PathSense::HAPLOTYPE:
    {
      // This information isn't cached.
      size_t path_number = this->get_metadata_index(path_handle);
      // We need to get the final step along the path.
      //
      // inverseLF() deosnt' work *from* the end marker, because the end marker
      // edge values aren't actually necessarily unique.
      //
      // So we need to find the first step on the reverse version of the path, and
      // scan its node for the last step on the forward version of the path.
      gbwt::edge_type first_reverse_edge = this->index->start(gbwt::Path::encode(path_number, true));
      if(first_reverse_edge.first == gbwt::ENDMARKER)
      {
        // Path is empty.
        return path_front_end(path_handle);
      }
      handlegraph::handle_t last_forward_handle = flip(node_to_handle(first_reverse_edge.first));

      // Look up the GBWT node we're using, in the path-forward orientaton
      gbwt::SearchState node_state = get_state(last_forward_handle);

      // Fill in the node part of the to edge, which doesn't change
      to.first = node_state.node;

      // Get an edge for each single haplotype selected by the search state.
      // Remember that ranges are inclusive at both ends.
      for(to.second = node_state.range.first; to.second <= node_state.range.second; ++to.second)
      {
        gbwt::edge_type next_edge = this->index->LF(to);
        if(next_edge.first != gbwt::ENDMARKER)
        {
          // We're only interested in final visits, which we can ID before
          // locate()
          continue;
        }

        // This candidate is the last visit along some path, so it might be the
        // last visit along the path we're looking for.
        auto sequence_number = this->index->locate(to);
        if(gbwt::Path::is_reverse(sequence_number))
        {
          // We're only interested in forward versions of paths
          continue;
        }

        auto path_number_here = gbwt::Path::id(sequence_number);
        if(path_number_here != path_number)
        {
          // We're only interested in the one path.
          continue;
        }
        // This is the visit! Stop scanning!
        break;
      }
    }
    break;
  default:
    throw std::runtime_error("Unimplemented sense!");
  }

  step_handle_t step;
  // Edges and steps have GBWT node numbers first
  as_integers(step)[0] = to.first;
  // And then offsets into the node's visits
  as_integers(step)[1] = to.second;

  return step;
}

step_handle_t
GBWTGraph::path_front_end(const path_handle_t&) const {
  // path_front_end can just be invalid_edge() since it's a past-end.
  gbwt::edge_type before_first_edge = gbwt::invalid_edge();

  step_handle_t step;
  // Edges and steps have GBWT node numbers first
  as_integers(step)[0] = before_first_edge.first;
  // And then offsets into the node's visits
  as_integers(step)[1] = before_first_edge.second;

  return step;
}

bool
GBWTGraph::has_next_step(const step_handle_t& step_handle) const {
  // Just look ahead and see if we get a past-end
  step_handle_t would_be_next = this->get_next_step(step_handle);
  gbwt::edge_type past_last_edge = gbwt::invalid_edge();
  return as_integers(would_be_next)[0] != past_last_edge.first || as_integers(would_be_next)[1] != past_last_edge.second;
}

bool
GBWTGraph::has_previous_step(const step_handle_t& step_handle) const {
  // Just look back and see if we get a sentinel
  step_handle_t would_be_prev = this->get_previous_step(step_handle);
  gbwt::edge_type before_first_edge = gbwt::invalid_edge();
  return as_integers(would_be_prev)[0] != before_first_edge.first || as_integers(would_be_prev)[1] != before_first_edge.second;
}

step_handle_t
GBWTGraph::get_next_step(const step_handle_t& step_handle) const {
  // Convert into a GBWT edge
  gbwt::edge_type here;
  here.first = as_integers(step_handle)[0];
  here.second = as_integers(step_handle)[1];

  // Follow it, and get either a real edge, or an edge where the node is
  // gbwt::ENDMARKER.
  here = this->index->LF(here);

  if(here.first == gbwt::ENDMARKER)
  {
    // We've hit the end marker. Use an invalid_edge() sentinel instead.
    here = gbwt::invalid_edge();
  }

  // Convert back
  step_handle_t next;
  as_integers(next)[0] = here.first;
  as_integers(next)[1] = here.second;

  return next;
}

step_handle_t
GBWTGraph::get_previous_step(const step_handle_t& step_handle) const {
  // Convert into a GBWT edge
  gbwt::edge_type here;
  here.first = as_integers(step_handle)[0];
  here.second = as_integers(step_handle)[1];

  // Follow it, backward and get either a real edge, or an edge where the
  // node is gbwt::ENDMARKER.
  here = this->index->inverseLF(here);

  if(here.first == gbwt::ENDMARKER)
  {
    // We've hit the end marker. Use an invalid_edge() sentinel instead.
    here = gbwt::invalid_edge();
  }

  // Convert back
  step_handle_t prev;
  as_integers(prev)[0] = here.first;
  as_integers(prev)[1] = here.second;

  return prev;
}

bool
GBWTGraph::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const
{
  for(size_t i = 0; i < this->named_paths.size(); i++)
  {
    // Show the iteratee each path
    bool should_continue = iteratee(handlegraph::as_path_handle(i));
    if(!should_continue)
    {
      // We were told to stop
      return false;
    }
  }

  // We got to the end
  return true;
}

bool
GBWTGraph::for_each_step_on_handle_impl(const handle_t& handle,
  const std::function<bool(const step_handle_t&)>& iteratee) const
{
  // Nothing to do without named paths.
  if(this->get_path_count() == 0) { return true; }

  return this->for_each_edge_and_path_on_handle(handle, [&](const gbwt::edge_type& candidate_edge, const gbwt::size_type& path_number)
  {
    if(this->id_to_path.count(path_number))
    {
      // This is a path on an indexed sample. We elide the haplotypes.
      // Prepare a step.
      step_handle_t step;
      handlegraph::as_integers(step)[0] = candidate_edge.first;
      handlegraph::as_integers(step)[1] = candidate_edge.second;
      // And show it to the iteratee
      return iteratee(step);
    }
    // Otherwise, continue.
    return true;
  });
}

//------------------------------------------------------------------------------

PathSense
GBWTGraph::get_sense(const path_handle_t& handle) const
{
  if(handlegraph::as_integer(handle) < this->named_paths.size())
  {
    // This is a cached named path.

    // To avoid pulling out the whole sample name every time we want to check
    // path sense, we cache the sense in the NamedPath.
    auto& named_path = this->named_paths[handlegraph::as_integer(handle)];
    return named_path.sense;
  }
  // Otherwise it's a haolotype
  return PathSense::HAPLOTYPE;
}

std::string
GBWTGraph::get_sample_name(const path_handle_t& handle) const
{
  PathSense sense = this->get_sense(handle);
  auto& structured_name = this->index->metadata.path(this->get_metadata_index(handle));
  return gbwtgraph::get_path_sample_name(this->index->metadata, structured_name, sense);
}

std::string
GBWTGraph::get_locus_name(const path_handle_t& handle) const
{
  PathSense sense = this->get_sense(handle);
  auto& structured_name = this->index->metadata.path(this->get_metadata_index(handle));
  return gbwtgraph::get_path_locus_name(this->index->metadata, structured_name, sense);
}

size_t
GBWTGraph::get_haplotype(const path_handle_t& handle) const
{
  PathSense sense = this->get_sense(handle);
  auto& structured_name = this->index->metadata.path(this->get_metadata_index(handle));
  return gbwtgraph::get_path_haplotype(this->index->metadata, structured_name, sense);
}

size_t
GBWTGraph::get_phase_block(const path_handle_t& handle) const
{
  PathSense sense = this->get_sense(handle);
  auto& structured_name = this->index->metadata.path(this->get_metadata_index(handle));
  return gbwtgraph::get_path_phase_block(this->index->metadata, structured_name, sense);
}

subrange_t
GBWTGraph::get_subrange(const path_handle_t& handle) const
{
  PathSense sense = this->get_sense(handle);
  auto& structured_name = this->index->metadata.path(this->get_metadata_index(handle));
  return gbwtgraph::get_path_subrange(this->index->metadata, structured_name, sense);
}

std::vector<gbwt::size_type>
GBWTGraph::sample_numbers_for_sample_name(const std::unordered_set<PathSense>* senses, const std::string& sample_name) const
{
  std::vector<gbwt::size_type> sample_numbers;
  if(sample_name == NO_SAMPLE_NAME)
  {
    // Don't try and look up the sentinel
    if(!senses || senses->count(PathSense::GENERIC))
    {
      // But we might have the non-sample sample
      gbwt::size_type sample_number = this->index->metadata.sample(REFERENCE_PATH_SAMPLE_NAME);
      if(sample_number < this->index->metadata.sample_names.size())
      {
        sample_numbers.push_back(sample_number);
      }
    }
    return sample_numbers;
  }
  // Otherwise we aren't working with NO_SAMPLE_NAME.
  if(!senses || senses->count(PathSense::HAPLOTYPE))
  {
    // Include just the same as the user-visible sample name
    gbwt::size_type sample_number = this->index->metadata.sample(sample_name);
    if(sample_number < this->index->metadata.sample_names.size())
    {
      sample_numbers.push_back(sample_number);
    }
  }
  if(!senses || senses->count(PathSense::REFERENCE))
  {
    // Include the sample name with the reference prefix
    gbwt::size_type sample_number = this->index->metadata.sample(REFERENCE_PATH_SAMPLE_NAME + sample_name);
    if(sample_number < this->index->metadata.sample_names.size())
    {
      sample_numbers.push_back(sample_number);
    }
  }
  return sample_numbers;
}

bool
GBWTGraph::for_each_path_matching_impl(const std::unordered_set<PathSense>* senses,
                                       const std::unordered_set<std::string>* samples,
                                       const std::unordered_set<std::string>* loci,
                                       const std::function<bool(const path_handle_t&)>& iteratee) const
{


  // First try to divert queries we have efficient implementations for.
  // TODO: Can we look up all the sample numbers and all the locus numbers once
  // first? That would be faster.

  if((senses && senses->empty()) || (samples && samples->empty()) || (loci && loci->empty()))
  {
    // Nothing to do!
    return true;
  }

  if(samples && samples->size() == 1 && loci && loci->size() == 1)
  {
    // We can look up one sample and locus, and we don't even need to filter ourselves.
    return this->for_each_path_matching_sample_and_locus(senses, *samples->begin(), *loci->begin(), iteratee);
  }

  if (senses && senses->size() == 1 && senses->count(PathSense::GENERIC)) {
    // We only want to look at generic paths, so we know the sample can only be one thing.
    if(!samples || samples->count(NO_SAMPLE_NAME))
    {
      // We can use the only allowed sample.
      if(loci)
      {
        // And we have loci to restrict to.
        for(auto& locus_name : *loci)
        {
          if(!this->for_each_path_matching_sample_and_locus(senses, NO_SAMPLE_NAME, locus_name, iteratee))
          {
            return false;
          }
        }
        return true;
      }
      else
      {
        // Any locus is acceptable in this sample.
        return this->for_each_path_matching_sample(senses, NO_SAMPLE_NAME, iteratee);
      }
    }
    else
    {
      // We banned the only allowed sample.
      return true;
    }
  }

  if(samples)
  {
    // We can look at just these particular samples.
    for(auto& sample_name : *samples)
    {
      if(loci && loci->size() == 1)
      {
        // We can look up the one locus for every sample
        if(!this->for_each_path_matching_sample_and_locus(senses, sample_name, *loci->begin(), iteratee))
        {
          return false;
        }
      }
      else
      {
        // We need to scan and filter on locus for each sample
        bool keep_going = this->for_each_path_matching_sample(senses, sample_name, [&](const path_handle_t& path_handle)
        {
          if(loci && !loci->count(this->get_locus_name(path_handle)))
          {
            return true;
          }
          return iteratee(path_handle);
        });
        if(!keep_going)
        {
          return false;
        }
      }
    }
    return true;
  }

  if(loci)
  {
    // We have no sample names but we have locus names
    for(auto& locus_name : *loci)
    {
      // We don't actually need to filter, because samples is unset, so all samples are allowed.
      if(!this->for_each_path_matching_locus(senses, locus_name, iteratee))
      {
        return false;
      }
    }
    return true;
  }

  // If we get here we have no sample or locus names.

  if(senses && !senses->count(PathSense::HAPLOTYPE))
  {
    // We just have to go through all the named paths
    for(size_t i = 0; i < this->named_paths.size(); i++)
    {
      auto& named_path = this->named_paths[i];
      if(senses->count(named_path.sense))
      {
        // This path is a sense we want
        if(!iteratee(handlegraph::as_path_handle(i)))
        {
          return false;
        }
      }
    }
    return true;
  }

  // We want haplotype paths at least, and we have no samples or loci to select
  // on, so we might as well go through everything.
  for(size_t i = 0; i < this->index->metadata.paths(); i++)
  {
    path_handle_t path_handle = this->from_metadata_index(i);
    if(senses && !senses->count(this->get_sense(path_handle)))
    {
      // THis sense is unwanted.
      continue;
    }
    if(!iteratee(path_handle))
    {
      return false;
    }
  }
  return true;
}

bool
GBWTGraph::for_each_path_matching_sample_and_locus(const std::unordered_set<PathSense>* senses,
                                                   const std::string& sample_name,
                                                   const std::string& locus_name,
                                                   const std::function<bool(const path_handle_t&)>& iteratee) const
{
  auto contig_number = this->index->metadata.contig(locus_name);
  if(contig_number == this->index->metadata.contig_names.size())
  {
    // Looking for a nonexistent locus
    return true;
  }
  for(auto& sample_number : this->sample_numbers_for_sample_name(senses, sample_name))
  {
    // For each sample number that might belong to this sample for one of the
    // senses we want
    for(auto& path_id : this->index->metadata.findPaths(sample_number, contig_number))
    {
      // For each path in that sample-sense and locus, try it.
      if(!iteratee(this->from_metadata_index(path_id)))
      {
        return false;
      }
    }
  }
  return true;
}

bool
GBWTGraph::for_each_path_matching_sample(const std::unordered_set<PathSense>* senses,
                                         const std::string& sample_name,
                                         const std::function<bool(const path_handle_t&)>& iteratee) const
{
  for(auto& sample_number : this->sample_numbers_for_sample_name(senses, sample_name))
  {
    // For each sample number that might belong to this sample for one of the
    // senses we want
    for(auto& path_id : this->index->metadata.pathsForSample(sample_number))
    {
      // For each path in that sample-sense, try it
      if(!iteratee(this->from_metadata_index(path_id)))
      {
        return false;
      }
    }
  }
  return true;
}

bool
GBWTGraph::for_each_path_matching_locus(const std::unordered_set<PathSense>* senses,
                                        const std::string& locus_name,
                                        const std::function<bool(const path_handle_t&)>& iteratee) const
{

  auto contig_number = this->index->metadata.contig(locus_name);
  if(contig_number == this->index->metadata.contig_names.size())
  {
    // Looking for a nonexistent locus
    return true;
  }
  for(auto& path_id : this->index->metadata.pathsForContig(contig_number))
  {
    // For each path in that locus, get the handle
    path_handle_t path_handle = this->from_metadata_index(path_id);
    if(!senses || senses->count(this->get_sense(path_handle)))
    {
      // This is a sense we want.
      // We have to check because locus number doesn't determine sense like
      // sample number does.
      if(!iteratee(path_handle))
      {
        return false;
      }
    }
  }
  return true;
}

bool
GBWTGraph::for_each_step_of_sense_impl(const handle_t& visited, const PathSense& sense, const std::function<bool(const step_handle_t&)>& iteratee) const
{
    return this->for_each_edge_and_path_on_handle(visited, [&](const gbwt::edge_type& candidate_edge, const gbwt::size_type& path_number)
    {
      auto found = this->id_to_path.find(path_number);
      if(found == this->id_to_path.end() && sense != PathSense::HAPLOTYPE)
      {
        // We found a haplotype we don't want
        return true;
      }
      else if(found != this->id_to_path.end() && this->named_paths[found->second].sense != sense)
      {
        // We are looking for reference paths but this is a generic path, or visa versa.
        return true;
      }

      // Make the step handle.
      step_handle_t step;
      handlegraph::as_integers(step)[0] = candidate_edge.first;
      handlegraph::as_integers(step)[1] = candidate_edge.second;
      // And show it to the iteratee
      return iteratee(step);
  });
}

bool
GBWTGraph::for_each_edge_and_path_on_handle(const handle_t& handle, const std::function<bool(const gbwt::edge_type&, const gbwt::size_type&)>& iteratee) const
{

  for(const handle_t oriented_handle : {handle, flip(handle)})
  {
    // We need to look at both orientations, because we only want steps on the
    // forward versions of their paths. And there's no good way to go from an
    // edge on a path's reverse sequence to an edge on the path's forward one.

    // Look up the GBWT node
    gbwt::SearchState node_state = get_state(oriented_handle);

    gbwt::edge_type candidate_edge;
    candidate_edge.first = node_state.node;
    for(candidate_edge.second = node_state.range.first;
        candidate_edge.second <= node_state.range.second;
        ++candidate_edge.second)
    {
      // Get the edge for each haplotype in the start-and-end-inclusive range

      // Get the sequence number the edge is on.
      // TODO: We could have a version of locate() that decompresses the entire DA for the node.
      auto sequence_number = this->index->locate(candidate_edge);

      if(gbwt::Path::is_reverse(sequence_number))
      {
        // We're looking at the reverse version of the path, which doesn't
        // have a corresponding libhandlegraph step, because step handles don't
        // help you at all with path orientation. Skip this one and come up
        // with the forward version of the thread when we look at the other
        // orientation of the handle.
        continue;
      }

      auto path_number = gbwt::Path::id(sequence_number);

      bool should_continue = iteratee(candidate_edge, path_number);
      if(!should_continue)
      {
        // We are supposed to stop now.
        return false;
      }
    }
  }

  // We made it to the end.
  return true;

}

size_t
GBWTGraph::get_metadata_index(const path_handle_t& handle) const
{
  size_t scratch = handlegraph::as_integer(handle);
  if(scratch < this->named_paths.size())
  {
    // Look up the metadata object path number for this cache entry
    return this->named_paths[scratch].id;
  }
  else
  {
    // There's no cache entry; remove the offset and get the metadata object path number.
    return scratch - this->named_paths.size();
  }
}

path_handle_t
GBWTGraph::from_metadata_index(const size_t& metadata_index) const
{
  // This might be a named path or a haplotype path.
  auto found = this->id_to_path.find(metadata_index);
  if (found == this->id_to_path.end()) {
    // This isn't referenced by a stored NamedPath. Must be a haplotype path.
    return handlegraph::as_path_handle(this->named_paths.size() + metadata_index);
  }
  // Otherwise just use the number of the stored NamedPath
  return handlegraph::as_path_handle(found->second);
}

//------------------------------------------------------------------------------

bool
GBWTGraph::has_segment_names() const
{
  return this->header.get(Header::FLAG_TRANSLATION);
}

std::pair<std::string, std::pair<nid_t, nid_t>>
GBWTGraph::get_segment(const handle_t& handle) const
{
  // If there is no translation, the predecessor is always at the end.
  nid_t id = this->get_id(handle);
  auto iter = this->node_to_segment.predecessor(id);
  if(!(this->has_node(id)) || iter == this->node_to_segment.one_end())
  {
    return std::pair<std::string, std::pair<nid_t, nid_t>>(std::to_string(id), std::make_pair(id, id + 1));
  }

  nid_t start = iter->second;
  std::string name = this->segments.str(iter->first);
  ++iter;
  nid_t limit = iter->second;

  return std::make_pair(name, std::make_pair(start, limit));
}

std::pair<std::string, size_t>
GBWTGraph::get_segment_name_and_offset(const handle_t& handle) const
{
  // If there is no translation, the predecessor is always at the end.
  nid_t id = this->get_id(handle);
  auto iter = this->node_to_segment.predecessor(id);
  if(!(this->has_node(id)) || iter == this->node_to_segment.one_end())
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
  if(!(this->has_node(id)) || iter == this->node_to_segment.one_end()) { return 0; }

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

bool
GBWTGraph::for_each_segment_impl(const std::function<bool(const std::string&, const std::pair<nid_t, nid_t>&)>& iteratee, bool parallel) const
{
  if(!(this->has_segment_names())) { return true; }

  bool keep_going = true;
  if(parallel)
  {
    // The parallel version iterates over all handles, converts them to segments, and
    // calls the iteratee if the handle is the first node of the segment.
    keep_going = this->for_each_handle([&](const handle_t& handle) -> bool
    {
      std::pair<std::string, std::pair<nid_t, nid_t>> segment = this->get_segment(handle);
      if(this->get_id(handle) == segment.second.first)
      {
        return iteratee(segment.first, segment.second);
      }
      return true;
    }, parallel);
  }
  else
  {
    auto iter = this->node_to_segment.one_begin();
    while(iter != this->node_to_segment.one_end())
    {
      nid_t start = iter->second;
      std::string name = this->segments.str(iter->first);
      ++iter;
      nid_t limit = iter->second;
      // The translation may include segments that were not used on any path.
      // The corresponding nodes are missing from the graph.
      if(this->has_node(start))
      {
        if(!iteratee(name, std::make_pair(start, limit))) { keep_going = false; break; }
      }
    }
  }

  return keep_going;
}

bool
GBWTGraph::for_each_link_impl(const std::function<bool(const edge_t&, const std::string&, const std::string&)>& iteratee, bool parallel) const
{
  if(!(this->has_segment_names())) { return true; }

  return this->for_each_segment([&](const std::string& from_segment, std::pair<nid_t, nid_t> nodes) -> bool
  {
    bool keep_going = true;
    // Right edges from forward orientation are canonical if the destination node
    // has a greater id or if the edge is a self-loop.
    handle_t last = this->get_handle(nodes.second - 1, false);
    keep_going = this->follow_edges(last, false, [&](const handle_t& next) -> bool
    {
      nid_t next_id = this->get_id(next);
      if(next_id >= nodes.second - 1)
      {
        std::string to_segment = this->get_segment_name(next);
        if(!iteratee(edge_t(last, next), from_segment, to_segment)) { return false; }
      }
      return true;
    });
    if(!keep_going) { return false; }

    // Right edges from reverse orientation are canonical if the destination node
    // has a greater id or if the edge is a self-loop to forward orientation of
    // this node.
    handle_t first = this->get_handle(nodes.first, true);
    keep_going = this->follow_edges(first, false, [&](const handle_t& next) -> bool
    {
      nid_t next_id = this->get_id(next);
      if(next_id > nodes.first || (next_id == nodes.first && !(this->get_is_reverse(next))))
      {
        std::string to_segment = this->get_segment_name(next);
        if(!iteratee(edge_t(first, next), from_segment, to_segment)) { return false; }
      }
      return true;
    });
    return keep_going;
  }, parallel);
}

//------------------------------------------------------------------------------

/// Translate a node range back to segment space
std::vector<oriented_node_range_t>
GBWTGraph::translate_back(const oriented_node_range_t& range) const {
  // If there is no translation, the predecessor is always at the end.
  nid_t id = std::get<0>(range);
  auto iter = this->node_to_segment.predecessor(id);
  if(!(this->has_node(id)) || iter == this->node_to_segment.one_end())
  {
    // No segments or nonexistent node.
    return {};
  }

  // Determine the total length of nodes in this segment that precede `id`
  // in the given orientation.
  size_t start = 0, limit = 0;
  if(std::get<1>(range))
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

  return {oriented_node_range_t(iter->first, std::get<1>(range), offset + std::get<2>(range), std::get<3>(range))};

}

/// Get a segment name
std::string
GBWTGraph::get_back_graph_node_name(const nid_t& back_node_id) const {
    return this->segments.str(back_node_id);
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
  // Read the header.
  Header h = sdsl::simple_sds::load_value<Header>(in);
  h.check();
  bool simple_sds = h.get(Header::FLAG_SIMPLE_SDS);
  h.unset(Header::FLAG_SIMPLE_SDS); // We only set this flag in the serialized header.
  h.set_version(); // Update to the current version.
  this->header = h;

  // Load the graph.
  if(simple_sds)
  {
    if(this->index == nullptr)
    {
      throw InvalidGBWT("GBWTGraph: A GBWT index is required for loading simple-sds format");
    }
    {
      gbwt::StringArray forward_only;
      forward_only.simple_sds_load(in);
      this->sequences = gbwt::StringArray(2 * forward_only.size(),
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
    this->sequences.load(in);
    this->real_nodes.load(in);
  }

  // Load the translation.
  if(simple_sds)
  {
    // The translation may be empty but not absent.
    this->segments.simple_sds_load(in);
    this->node_to_segment.simple_sds_load(in);
    if(this->header.get(Header::FLAG_TRANSLATION) != (this->segments.size() > 0))
    {
      throw sdsl::simple_sds::InvalidData("GBWTGraph: Invalid translation flag in the header");
    }
  }
  else if(this->header.get(Header::FLAG_TRANSLATION))
  {
    this->segments.load(in);
    this->node_to_segment.load(in);
  }

  this->sanity_checks();
}

void
GBWTGraph::set_gbwt(const gbwt::GBWT& gbwt_index)
{
  this->index = &gbwt_index;

  if(!(this->index->bidirectional()))
  {
    throw InvalidGBWT("GBWTGraph: The GBWT index must be bidirectional");
  }

  this->reference_samples = parse_reference_samples_tag(this->index.tags.get(REFERENCE_SAMPLE_LIST_GBWT_TAG));
  this->cache_named_paths();
}

//------------------------------------------------------------------------------

void
GBWTGraph::simple_sds_serialize(std::ostream& out) const
{
  // Serialize the header.
  Header copy = this->header;
  copy.set(Header::FLAG_SIMPLE_SDS); // We only set this flag in the serialized header.
  sdsl::simple_sds::serialize_value(copy, out);

  // Compress the sequences. `real_nodes` can be rebuilt from the GBWT.
  {
    gbwt::StringArray forward_only(this->sequences.size() / 2,
    [&](size_t offset) -> size_t
    {
      return this->sequences.length(2 * offset);
    },
    [&](size_t offset) -> view_type
    {
      return this->sequences.view(2 * offset);
    });
    forward_only.simple_sds_serialize(out);
  }

  // Compress the translation.
  this->segments.simple_sds_serialize(out);
  this->node_to_segment.simple_sds_serialize(out);
}

void
GBWTGraph::simple_sds_load(std::istream& in, const gbwt::GBWT& gbwt_index)
{
  // Set the GBWT so we can rebuild `real_nodes` later.
  this->set_gbwt(gbwt_index);

  // The same deserialize() function can handle the SDSL and simple-sds formats.
  this->deserialize_members(in);
}

size_t
GBWTGraph::simple_sds_size() const
{
  size_t result = sdsl::simple_sds::value_size(this->header);

  // Compress the sequences.
  {
    gbwt::StringArray forward_only(this->sequences.size() / 2,
    [&](size_t offset) -> size_t
    {
      return this->sequences.length(2 * offset);
    },
    [&](size_t offset) -> view_type
    {
      return this->sequences.view(2 * offset);
    });
    result += forward_only.simple_sds_size();
  }

  // Compress the translation.
  result += this->segments.simple_sds_size();
  result += this->node_to_segment.simple_sds_size();

  return result;
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
