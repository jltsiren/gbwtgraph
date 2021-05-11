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
const std::string GBWTGraph::COMPRESSED_EXTENSION = ".gbz";

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
  if(this->tag != TAG) { return false; }
  switch(this->version)
  {
  case VERSION:
    return ((this->flags & FLAG_MASK) == this->flags);
  case TRANS_VERSION:
    return ((this->flags & TRANS_FLAG_MASK) == this->flags);
  case OLD_VERSION:
    return ((this->flags & OLD_FLAG_MASK) == this->flags);
  default:
    return false;
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
  this->tags.swap(another.tags);
  this->sequences.swap(another.sequences);
  this->real_nodes.swap(another.real_nodes);
  this->segments.swap(another.segments);
  this->node_to_segment.swap(another.node_to_segment);
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
    this->tags = std::move(source.tags);
    this->sequences = std::move(source.sequences);
    this->real_nodes = std::move(source.real_nodes);
    this->segments = std::move(source.segments);
    this->node_to_segment = std::move(source.node_to_segment);
  }
  return *this;
}

void
GBWTGraph::copy(const GBWTGraph& source)
{
  this->index = source.index;
  this->header = source.header;
  this->tags = source.tags;
  this->sequences = source.sequences;
  this->real_nodes = source.real_nodes;
  this->segments = source.segments;
  this->node_to_segment = source.node_to_segment;
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

void
GBWTGraph::reset_tags()
{
  this->tags.clear();
  this->add_source();
}

void
GBWTGraph::add_source()
{
  this->tags[Version::SOURCE_KEY] = Version::SOURCE_VALUE;
}

//------------------------------------------------------------------------------

GBWTGraph::GBWTGraph(const gbwt::GBWT& gbwt_index, const HandleGraph& sequence_source) :
  index(nullptr)
{
  // Set GBWT and do sanity checks.
  this->set_gbwt(gbwt_index);
  if(this->index->empty()) { return; }

  // Set the source tag.
  this->add_source();

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
}

GBWTGraph::GBWTGraph(const gbwt::GBWT& gbwt_index, const SequenceSource& sequence_source) :
  index(nullptr)
{
  // Set GBWT and do sanity checks.
  this->set_gbwt(gbwt_index);
  if(this->index->empty()) { return; }

  // Set the source tag.
  this->add_source();

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

  // Store the node to segment translation.
  if(sequence_source.uses_translation())
  {
    this->header.set(Header::FLAG_TRANSLATION);
    std::tie(this->segments, this->node_to_segment) = sequence_source.invert_translation();
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

void
GBWTGraph::for_each_segment(const std::function<bool(const std::string&, std::pair<nid_t, nid_t>)>& iteratee) const
{
  if(!(this->has_segment_names())) { return; }

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
      if(!iteratee(name, std::make_pair(start, limit))) { return; }
    }
  }
}

void
GBWTGraph::for_each_link(const std::function<bool(const edge_t&, const std::string&, const std::string&)>& iteratee) const
{
  if(!(this->has_segment_names())) { return; }

  this->for_each_segment([&](const std::string& from_segment, std::pair<nid_t, nid_t> nodes) -> bool
  {
    bool keep_going = true;
    // Right edges from forward orienation are canonical if the destination node
    // has a greater id or if the edge is a self-loop.
    handle_t last = this->get_handle(nodes.second - 1, false);
    this->follow_edges(last, false, [&](const handle_t& next) -> bool
    {
      nid_t next_id = this->get_id(next);
      if(next_id >= nodes.second - 1)
      {
        std::string to_segment = this->get_segment_name(next);
        keep_going = iteratee(edge_t(last, next), from_segment, to_segment);
      }
      return keep_going;
    });
    if(!keep_going) { return false; }

    // Right edges from reverse orientation are canonical if the destination node
    // has a greater id or if the edge is a self-loop to forward orientation of
    // this node.
    handle_t first = this->get_handle(nodes.first, true);
    this->follow_edges(first, false, [&](const handle_t& next) -> bool
    {
      nid_t next_id = this->get_id(next);
      if(next_id > nodes.first || (next_id == nodes.first && !(this->get_is_reverse(next))))
      {
        std::string to_segment = this->get_segment_name(next);
        keep_going = iteratee(edge_t(first, next), from_segment, to_segment);
      }
      return keep_going;
    });
    return keep_going;
  });
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

  {
    gbwt::StringArray linearized(this->tags);
    linearized.serialize(out);
  }

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
  if(!(h.check()))
  {
    throw sdsl::simple_sds::InvalidData("GBWTGraph: Invalid header");
  }
  bool simple_sds = h.get(Header::FLAG_SIMPLE_SDS);
  bool has_tags = h.version >= 3; // FIXME Replace with symbolic constant.
  h.unset(Header::FLAG_SIMPLE_SDS); // We only set this flag in the serialized header.
  h.set_version(); // Update to the current version.
  this->header = h;

  // Read the tags and set the source to Version::SOURCE_VALUE.
  if(has_tags)
  {
    gbwt::StringArray linearized;
    if(simple_sds) { linearized.simple_sds_load(in); }
    else { linearized.load(in); }
    if(linearized.size() % 2 != 0)
    {
      throw sdsl::simple_sds::InvalidData("GBWTGraphs: Tag without a value");
    }
    this->tags.clear();
    for(size_t i = 0; i < linearized.size(); i += 2)
    {
      std::string key = linearized.str(i);
      for(auto iter = key.begin(); iter != key.end(); ++iter) { *iter = std::tolower(*iter); }
      this->tags[key] = linearized.str(i + 1);
    }
    if(this->tags.size() != linearized.size() / 2)
    {
      throw sdsl::simple_sds::InvalidData("GBWTGraph: Duplicate tags");
    }
    // If we need the original source tag, we should check it here.
    this->add_source();
  }
  else { this->reset_tags(); }

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
}

//------------------------------------------------------------------------------

void
GBWTGraph::simple_sds_serialize(std::ostream& out) const
{
  // Serialize the header.
  Header copy = this->header;
  copy.set(Header::FLAG_SIMPLE_SDS); // We only set this flag in the serialized header.
  sdsl::simple_sds::serialize_value(copy, out);

  {
    gbwt::StringArray linearized(this->tags);
    linearized.simple_sds_serialize(out);
  }

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

  {
    gbwt::StringArray linearized(this->tags);
    result += linearized.simple_sds_size();
  }

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
