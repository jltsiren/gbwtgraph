#include <gbwtgraph/error_handling.h>
#include <gbwt/gbwt.h>
#include <gbwt/utils.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/utils.h>
#include <sdsl/simple_sds.hpp>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Numerical class constants.
template <typename CharAllocatorType>
constexpr std::uint32_t GBZ<CharAllocatorType>::Header::TAG;

template <typename CharAllocatorType>
constexpr std::uint32_t GBZ<CharAllocatorType>::Header::VERSION;

template <typename CharAllocatorType>
constexpr std::uint64_t GBZ<CharAllocatorType>::Header::FLAG_MASK;

//------------------------------------------------------------------------------

// Other class variables.

template <typename CharAllocatorType>
const std::string GBZ<CharAllocatorType>::EXTENSION = ".gbz";

//------------------------------------------------------------------------------

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::Header::Header() :
  tag(TAG), version(VERSION),
  flags(0)
{
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::Header::check() const
{
  if(this->tag != TAG)
  {
    GBWTGRAPH_THROW(sdsl::simple_sds::InvalidData("GBZ: Invalid tag"));
  }

  if(this->version > VERSION || this->version < OLD_VERSION)
  {
    std::string msg = "GBZ: Expected version " + std::to_string(OLD_VERSION) + " to " + std::to_string(VERSION) + ", got version " + std::to_string(this->version);
    GBWTGRAPH_THROW(sdsl::simple_sds::InvalidData(msg));
  }

  std::uint64_t mask = 0;
  switch(this->version)
  {
  case VERSION:
    mask = FLAG_MASK; break;
  case OLD_VERSION:
    mask = OLD_FLAG_MASK; break;
  }
  if((this->flags & mask) != this->flags)
  {
    GBWTGRAPH_THROW(sdsl::simple_sds::InvalidData("GBZ: Invalid flags"));
  }
}

template <typename CharAllocatorType>
bool
GBZ<CharAllocatorType>::Header::operator==(const Header& another) const
{
  return (this->tag == another.tag && this->version == another.version &&
          this->flags == another.flags);
}

//------------------------------------------------------------------------------

template <typename CharAllocatorType>
size_t
GBZ<CharAllocatorType>::set_reference_samples(const sample_name_set& samples)
{
  sample_name_set present_samples = present_sample_names(samples, this->index);

  std::string tag_value = compose_reference_samples_tag(present_samples);
  this->index.tags.set(REFERENCE_SAMPLE_LIST_GBWT_TAG, tag_value);

  // By resetting the pointer to the GBWT index, we tell the GBWTGraph object
  // to pull the tag and recache named paths.
  this->set_gbwt();

  return present_samples.size();
}

//------------------------------------------------------------------------------

#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(bi::managed_shared_memory* shared_memory) :
  graph(shared_memory)
{
  this->add_source();
  this->set_gbwt();
  this->shared_memory = shared_memory;
  this->compute_pggname(nullptr);
}

#else

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ()
{
  this->add_source();
  this->set_gbwt();
  this->compute_pggname(nullptr);
}

#endif

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(const GBZ& source)
#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
  // `graph` would otherwise default-construct via GBWTGraph's default
  // constructor with a null shared memory segment, which is only valid for
  // the plain-allocator instantiation; passing the source's segment here
  // keeps this safe (a no-op) for CharAllocatorType == std::allocator<char>
  // and correct for the shared-memory instantiation.
  : graph(source.shared_memory)
#endif
{
  this->copy(source);
}

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(GBZ&& source)
#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
  : graph(source.shared_memory)
#endif
{
  *this = std::move(source);
}

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::~GBZ()
{}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::swap(GBZ& another)
{
  if(&another == this) { return; }

  std::swap(this->header, another.header);
  this->tags.swap(another.tags);
  this->index.swap(another.index);
  this->graph.swap(another.graph);

  // GBWTGraph did not know that we also swapped the GBWTs.
  this->set_gbwt_address();
  another.set_gbwt_address();
}

template <typename CharAllocatorType>
GBZ<CharAllocatorType>&
GBZ<CharAllocatorType>::operator=(const GBZ& source)
{
  if(&source != this) { this->copy(source); }
  return *this;
}

template <typename CharAllocatorType>
GBZ<CharAllocatorType>&
GBZ<CharAllocatorType>::operator=(GBZ&& source)
{
  if(&source != this)
  {
    this->header = std::move(source.header);
    this->tags = std::move(source.tags);
    this->index = std::move(source.index);
    this->graph = std::move(source.graph);

    // GBWTGraph did not know that we also moved the GBWT.
    this->set_gbwt_address();
  }
  return *this;
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::copy(const GBZ& source)
{
  this->header = source.header;
  this->tags = source.tags;
  this->index = source.index;
  this->graph = source.graph;

  // Use the local copy of the GBWT.
  this->set_gbwt_address();
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::reset_tags()
{
  this->tags.clear();
  this->add_source();
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::add_source()
{
  this->tags.set(Version::SOURCE_KEY, Version::SOURCE_VALUE);
}

//------------------------------------------------------------------------------

#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<NaiveGraph>& graph, bi::managed_shared_memory* shared_memory)
{
  if(index == nullptr || graph == nullptr)
  {
    GBWTGRAPH_THROW(std::runtime_error("GBZ: Index and graph must be non-null"));
  }

  this->add_source();
  this->index = std::move(*index); index.reset();
  GraphName parent = graph->graph_name();
  this->graph = GBWTGraph<CharAllocatorType>(this->index, *graph, shared_memory);
  this->shared_memory = shared_memory;
  graph.reset();
  this->compute_pggname(&parent);
}

#else

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<NaiveGraph>& graph)
{
  if(index == nullptr || graph == nullptr)
  {
    GBWTGRAPH_THROW(std::runtime_error("GBZ: Index and graph must be non-null"));
  }

  this->add_source();
  this->index = std::move(*index); index.reset();
  GraphName parent = graph->graph_name();
  this->graph = GBWTGraph<CharAllocatorType>(this->index, *graph);
  graph.reset();
  this->compute_pggname(&parent);
}

#endif

#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(const gbwt::GBWT& index, const NaiveGraph& graph, bi::managed_shared_memory* shared_memory) :
  index(index),
  graph(this->index, graph, shared_memory)
{
  this->add_source();
  GraphName parent = graph.graph_name();
  this->shared_memory = shared_memory;
  this->compute_pggname(&parent);
}

#else

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(const gbwt::GBWT& index, const NaiveGraph& graph) :
  index(index),
  graph(this->index, graph)
{
  this->add_source();
  GraphName parent = graph.graph_name();
  this->compute_pggname(&parent);
}

#endif

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(std::vector<GBZ>&& subgraphs)
{
  if(subgraphs.empty())
  {
    this->add_source();
    this->set_gbwt();
    this->compute_pggname(nullptr);
    return;
  }

  // Move the GBWTs into a vector for the GBWT merge constructor.
  std::vector<gbwt::GBWT> indexes;
  indexes.reserve(subgraphs.size());
  for(GBZ& subgraph : subgraphs)
  {
    indexes.push_back(std::move(subgraph.index));
    subgraph.graph.set_gbwt_address(indexes.back());
  }

  // Determine the reference samples as a union over the subgraphs.
  sample_name_set reference_samples;
  for(const gbwt::GBWT& index : indexes)
  {
    auto current_ref_samples = parse_reference_samples_tag(index);
    reference_samples.insert(current_ref_samples.begin(), current_ref_samples.end());
  }

  // Merge the GBWTs and set reference samples.
  // This may fail if the subgraphs are not disjoint.
  this->index = gbwt::GBWT(indexes);
  if(!reference_samples.empty())
  {
    std::string reference_samples_tag = compose_reference_samples_tag(reference_samples);
    this->index.tags.set(REFERENCE_SAMPLE_LIST_GBWT_TAG, reference_samples_tag);
  }

  // Create the union of the graphs.
  NaiveGraph sequence_source;
  for(const GBZ& subgraph : subgraphs)
  {
    subgraph.graph.for_each_handle([&](const handle_t& handle)
    {
      std::string_view sequence = subgraph.graph.get_sequence_view(handle);
      sequence_source.create_node(subgraph.graph.get_id(handle), sequence);
    });
  }
  indexes.clear(); subgraphs.clear(); // Free memory before building the GBWTGraph.
  this->graph = GBWTGraph<CharAllocatorType>(this->index, sequence_source);

  // Determine GBZ tags.
  this->add_source();
  this->compute_pggname(nullptr);
}

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(gbwt::GBWT&& index, const GBZ& supergraph) :
  index(std::move(index))
#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
  // Not a real attach: just enough to avoid default-constructing `graph`
  // with a null segment (which the shared-memory instantiation cannot do;
  // see GBWTGraph's default constructor). The `if constexpr` below throws
  // for this instantiation before `graph` is used for anything real.
  , graph(supergraph.shared_memory)
#endif
{
  // `GBWTGraph::subgraph()` always returns a plain, heap-allocated graph (see
  // its declaration in gbwtgraph.h), so it can only be assigned directly into
  // `this->graph` when this GBZ is also using the default allocator. Building
  // a subgraph of a shared-memory GBZ is not a case we support.
  if constexpr (std::is_same<CharAllocatorType, std::allocator<char>>::value)
  {
    this->graph = supergraph.graph.subgraph(this->index);
  }
  else
  {
    GBWTGRAPH_THROW(std::runtime_error("GBZ: Building a subgraph of a shared-memory GBZ is not supported"));
  }
  this->add_source();
  GraphName parent = supergraph.graph_name();
  this->compute_pggname(&parent, ParentGraphType::SUPERGRAPH);
}

template <typename CharAllocatorType>
GBZ<CharAllocatorType>::GBZ(gbwt::GBWT&& index, const HandleGraph& graph, const NamedNodeBackTranslation* segment_space) :
  index(index)
{
  // This constructor takes an arbitrary, already-in-memory `HandleGraph`,
  // which is never itself shared-memory-backed, so there is nothing to gain
  // by supporting the shared-memory instantiation here; see the GBWT&&,
  // GBZ supergraph constructor above for the same reasoning.
  if constexpr (std::is_same<CharAllocatorType, std::allocator<char>>::value)
  {
    this->graph = GBWTGraph<CharAllocatorType>(this->index, graph, segment_space);
  }
  else
  {
    GBWTGRAPH_THROW(std::runtime_error("GBZ: Building from an arbitrary HandleGraph into a shared-memory GBZ is not supported"));
  }
  // Unlike the other constructors, this one does not call compute_pggname():
  // an arbitrary HandleGraph carries no GraphName information to import, so
  // (per this constructor's doc comment in gbz.h) that is left to the caller.
  this->add_source();
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::set_gbwt()
{
  this->graph.set_gbwt(this->index);
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::set_gbwt_address()
{
  this->graph.set_gbwt_address(this->index);
}

//------------------------------------------------------------------------------

template <typename CharAllocatorType>
bool
GBZ<CharAllocatorType>::compute_pggname(const GraphName* parent, ParentGraphType relationship)
{
  // `gbwt_to_canonical_gfa()` only accepts a plain, default-allocator graph
  // (see gfa.h); computing a pggname hashes the entire canonical GFA
  // representation, which is not worth generalizing to shared-memory graphs
  // just for this. A shared-memory GBZ still has fully working sequence and
  // path data; it simply does not get an automatically computed pggname tag.
  if constexpr (std::is_same<CharAllocatorType, std::allocator<char>>::value)
  {
    // Compute the name.
    DigestStream digest_stream(EVP_sha256());
    gbwt_to_canonical_gfa(this->graph, digest_stream);
    std::string digest = digest_stream.finish();
    if(digest.empty()) { return false; }

    // Set the name and copy existing relationships.
    GraphName name(digest);
    name.add_relationships(this->graph_name());

    // Determine the relationship to the parent graph, if given,
    // and copy relationships from it.
    if(parent != nullptr && parent->has_name())
    {
      if(relationship == ParentGraphType::HEURISTIC)
      {
        relationship = (this->graph.has_segment_names() ? ParentGraphType::TRANSLATION_TARGET : ParentGraphType::SUPERGRAPH);
      }
      if(relationship == ParentGraphType::TRANSLATION_TARGET)
      {
        if(!name.same(*parent))
        {
          name.add_translation(name.name(), parent->name());
          this->tags.set(GraphName::GBZ_TRANSLATION_TARGET_TAG, parent->name());
        }
      }
      else
      {
        // This does nothing if the names are the same.
        name.add_subgraph(name.name(), parent->name());
      }
      name.add_relationships(*parent);
    }

    // Store the information back into the tags.
    name.set_tags(this->tags);

    return true;
  }
  else
  {
    return false;
  }
}

//------------------------------------------------------------------------------

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::simple_sds_serialize(std::ostream& out) const
{
  sdsl::simple_sds::serialize_value(this->header, out);
  this->tags.simple_sds_serialize(out);
  this->index.simple_sds_serialize(out);
  this->graph.simple_sds_serialize(out);
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::simple_sds_serialize_v1(std::ostream& out) const
{
  // Only change the version number in the serialized header.
  Header header = this->header;
  header.version = Header::OLD_VERSION;
  sdsl::simple_sds::serialize_value(header, out);

  this->tags.simple_sds_serialize(out);
  this->index.simple_sds_serialize(out);
  this->graph.simple_sds_serialize_v3(out);
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::simple_sds_serialize(const gbwt::GBWT& index, const GBWTGraph<CharAllocatorType>& graph, std::ostream& out)
{
  GBZ<std::allocator<char>> empty;
  sdsl::simple_sds::serialize_value(empty.header, out);
  empty.tags.simple_sds_serialize(out);
  index.simple_sds_serialize(out);
  graph.simple_sds_serialize(out);
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::simple_sds_load(std::istream& in)
{
  this->header = sdsl::simple_sds::load_value<Header>(in);
  this->header.check();
  this->header.set_version();

  // Load the tags and update the source to this library.
  // We could also check if the source was already this library, but we have no
  // uses for that information at the moment.
  this->tags.simple_sds_load(in);
  this->add_source();

  this->index.simple_sds_load(in);
  this->graph.simple_sds_load(in, this->index);
}

template <typename CharAllocatorType>
size_t
GBZ<CharAllocatorType>::simple_sds_size() const
{
  size_t result = sdsl::simple_sds::value_size(this->header);
  result += this->tags.simple_sds_size();
  result += this->index.simple_sds_size();
  result += this->graph.simple_sds_size();
  return result;
}

template <typename CharAllocatorType>
gbwt::Tags
GBZ<CharAllocatorType>::simple_sds_load_tags(const std::string& filename)
{
  std::ifstream in(filename, std::ios_base::binary);
  if(!in)
  {
    GBWTGRAPH_THROW(sdsl::simple_sds::CannotOpenFile(filename, false));
  }
  in.exceptions(std::ios::eofbit | std::ios::badbit | std::ios::failbit);

  Header header = sdsl::simple_sds::load_value<Header>(in);
  header.check();
  gbwt::Tags tags;
  tags.simple_sds_load(in);

  in.close();
  return tags;
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::serialize_to_files(const std::string& gbwt_name, const std::string& graph_name, bool simple_sds_graph) const
{
  sdsl::simple_sds::serialize_to(this->index, gbwt_name);
  if(simple_sds_graph) { sdsl::simple_sds::serialize_to(this->graph, graph_name); }
  else { this->graph.serialize(graph_name); }
}

template <typename CharAllocatorType>
void
GBZ<CharAllocatorType>::load_from_files(const std::string& gbwt_name, const std::string& graph_name)
{
  this->header.set_version();
  this->tags.clear();
  this->add_source();
  sdsl::simple_sds::load_from(this->index, gbwt_name);
  this->set_gbwt();
  this->graph.deserialize(graph_name);
}

template class GBZ<std::allocator<char>>;
#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
template class GBZ<SharedMemCharAllocatorType>;
#endif

//------------------------------------------------------------------------------

} // namespace gbwtgraph
