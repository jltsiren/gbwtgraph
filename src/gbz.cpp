#include <stdexcept>
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
template <typename SAAllocator>
constexpr std::uint32_t CoreGBZ<SAAllocator>::Header::TAG;

template <typename SAAllocator>
constexpr std::uint32_t CoreGBZ<SAAllocator>::Header::VERSION;

template <typename SAAllocator>
constexpr std::uint64_t CoreGBZ<SAAllocator>::Header::FLAG_MASK;

//------------------------------------------------------------------------------

// Other class variables.

template <typename SAAllocator>
const std::string CoreGBZ<SAAllocator>::EXTENSION = ".gbz";

//------------------------------------------------------------------------------

template <typename SAAllocator>
CoreGBZ<SAAllocator>::Header::Header() :
  tag(TAG), version(VERSION),
  flags(0)
{
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::Header::check() const
{
  if(this->tag != TAG)
  {
    throw (sdsl::simple_sds::InvalidData("GBZ: Invalid tag"));
  }

  if(this->version > VERSION || this->version < OLD_VERSION)
  {
    std::string msg = "GBZ: Expected version " + std::to_string(OLD_VERSION) + " to " + std::to_string(VERSION) + ", got version " + std::to_string(this->version);
    throw (sdsl::simple_sds::InvalidData(msg));
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
    throw (sdsl::simple_sds::InvalidData("GBZ: Invalid flags"));
  }
}

template <typename SAAllocator>
bool
CoreGBZ<SAAllocator>::Header::operator==(const Header& another) const
{
  return (this->tag == another.tag && this->version == another.version &&
          this->flags == another.flags);
}

//------------------------------------------------------------------------------

template <typename SAAllocator>
size_t
CoreGBZ<SAAllocator>::set_reference_samples(const sample_name_set& samples)
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

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(bi::managed_shared_memory* shared_memory) :
  graph(shared_memory)
{
  this->add_source();
  this->set_gbwt();
  this->compute_pggname(nullptr);
}

#else

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ()
{
  this->add_source();
  this->set_gbwt();
  this->compute_pggname(nullptr);
}

#endif

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(const CoreGBZ& source)
#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
  // `graph` would otherwise default-construct via GBWTGraph's default
  // constructor with a null shared memory segment, which is only valid for
  // the plain-allocator instantiation; passing the source's segment here
  // keeps this safe (a no-op) for SAAllocator == std::allocator<char>
  // and correct for the shared-memory instantiation.
  : graph(source.graph.sequences.shared_memory)
#endif
{
  this->copy(source);
}

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(CoreGBZ&& source)
#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
  : graph(source.graph.sequences.shared_memory)
#endif
{
  *this = std::move(source);
}

template <typename SAAllocator>
CoreGBZ<SAAllocator>::~CoreGBZ()
{}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::swap(CoreGBZ& another)
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

template <typename SAAllocator>
CoreGBZ<SAAllocator>&
CoreGBZ<SAAllocator>::operator=(const CoreGBZ& source)
{
  if(&source != this) { this->copy(source); }
  return *this;
}

template <typename SAAllocator>
CoreGBZ<SAAllocator>&
CoreGBZ<SAAllocator>::operator=(CoreGBZ&& source)
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

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::copy(const CoreGBZ& source)
{
  this->header = source.header;
  this->tags = source.tags;
  this->index = source.index;
  this->graph = source.graph;

  // Use the local copy of the GBWT.
  this->set_gbwt_address();
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::reset_tags()
{
  this->tags.clear();
  this->add_source();
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::add_source()
{
  this->tags.set(Version::SOURCE_KEY, Version::SOURCE_VALUE);
}

//------------------------------------------------------------------------------

#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<NaiveGraph>& graph, bi::managed_shared_memory* shared_memory)
{
  if(index == nullptr || graph == nullptr)
  {
    throw (std::runtime_error("GBZ: Index and graph must be non-null"));
  }

  this->add_source();
  this->index = std::move(*index); index.reset();
  GraphName parent = graph->graph_name();
  this->graph = CoreGBWTGraph<SAAllocator>(this->index, *graph, shared_memory);
  graph.reset();
  this->compute_pggname(&parent);
}

#else

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<NaiveGraph>& graph)
{
  if(index == nullptr || graph == nullptr)
  {
    throw (std::runtime_error("GBZ: Index and graph must be non-null"));
  }

  this->add_source();
  this->index = std::move(*index); index.reset();
  GraphName parent = graph->graph_name();
  this->graph = CoreGBWTGraph<SAAllocator>(this->index, *graph);
  graph.reset();
  this->compute_pggname(&parent);
}

#endif

#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(const gbwt::GBWT& index, const NaiveGraph& graph, bi::managed_shared_memory* shared_memory) :
  index(index),
  graph(this->index, graph, shared_memory)
{
  this->add_source();
  GraphName parent = graph.graph_name();
  this->compute_pggname(&parent);
}

#else

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(const gbwt::GBWT& index, const NaiveGraph& graph) :
  index(index),
  graph(this->index, graph)
{
  this->add_source();
  GraphName parent = graph.graph_name();
  this->compute_pggname(&parent);
}

#endif

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(std::vector<GBZ>&& subgraphs)
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
  this->graph = CoreGBWTGraph<SAAllocator>(this->index, sequence_source);

  // Determine GBZ tags.
  this->add_source();
  this->compute_pggname(nullptr);
}

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(gbwt::GBWT&& index, const GBZ& supergraph) :
  index(std::move(index))
{
  // GBWTGraph::subgraph() always returns a plain, heap-allocated graph (see
  // its declaration in gbwtgraph.h), so this constructor only exists for the
  // plain-allocator CoreGBZ: assigning that result into `this->graph` does
  // not type-check for any other SAAllocator, and this constructor is never
  // explicitly instantiated for one (see the bottom of this file).
  this->graph = supergraph.graph.subgraph(this->index);
  this->add_source();
  GraphName parent = supergraph.graph_name();
  this->compute_pggname(&parent, ParentGraphType::SUPERGRAPH);
}

template <typename SAAllocator>
CoreGBZ<SAAllocator>::CoreGBZ(gbwt::GBWT&& index, const HandleGraph& graph, const NamedNodeBackTranslation* segment_space) :
  index(index)
{
  // This constructor's input, an arbitrary HandleGraph, is never itself
  // shared-memory-backed, so it always builds a plain graph; like the
  // supergraph constructor above, it is only instantiated for the plain
  // allocator.
  this->graph = CoreGBWTGraph<SAAllocator>(this->index, graph, segment_space);
  // Unlike the other constructors, this one does not call compute_pggname():
  // an arbitrary HandleGraph carries no GraphName information to import, so
  // (per this constructor's doc comment in gbz.h) that is left to the caller.
  this->add_source();
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::set_gbwt()
{
  this->graph.set_gbwt(this->index);
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::set_gbwt_address()
{
  this->graph.set_gbwt_address(this->index);
}

//------------------------------------------------------------------------------

template <typename SAAllocator>
bool
CoreGBZ<SAAllocator>::compute_pggname(const GraphName* parent, ParentGraphType relationship)
{
  // gbwt_to_canonical_gfa() reads sequences and topology through the
  // graph's public interface, so this works the same regardless of
  // SAAllocator; called by every constructor, including the shared-memory
  // ones, so it must not be locked to the plain allocator.

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

//------------------------------------------------------------------------------

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::simple_sds_serialize(std::ostream& out) const
{
  sdsl::simple_sds::serialize_value(this->header, out);
  this->tags.simple_sds_serialize(out);
  this->index.simple_sds_serialize(out);
  this->graph.simple_sds_serialize(out);
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::simple_sds_serialize_v1(std::ostream& out) const
{
  // Only change the version number in the serialized header.
  Header header = this->header;
  header.version = Header::OLD_VERSION;
  sdsl::simple_sds::serialize_value(header, out);

  this->tags.simple_sds_serialize(out);
  this->index.simple_sds_serialize(out);
  this->graph.simple_sds_serialize_v3(out);
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::simple_sds_serialize(const gbwt::GBWT& index, const CoreGBWTGraph<SAAllocator>& graph, std::ostream& out)
{
  CoreGBZ<std::allocator<char>> empty;
  sdsl::simple_sds::serialize_value(empty.header, out);
  empty.tags.simple_sds_serialize(out);
  index.simple_sds_serialize(out);
  graph.simple_sds_serialize(out);
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::simple_sds_load(std::istream& in)
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

template <typename SAAllocator>
size_t
CoreGBZ<SAAllocator>::simple_sds_size() const
{
  size_t result = sdsl::simple_sds::value_size(this->header);
  result += this->tags.simple_sds_size();
  result += this->index.simple_sds_size();
  result += this->graph.simple_sds_size();
  return result;
}

template <typename SAAllocator>
gbwt::Tags
CoreGBZ<SAAllocator>::simple_sds_load_tags(const std::string& filename)
{
  std::ifstream in(filename, std::ios_base::binary);
  if(!in)
  {
    throw (sdsl::simple_sds::CannotOpenFile(filename, false));
  }
  in.exceptions(std::ios::eofbit | std::ios::badbit | std::ios::failbit);

  Header header = sdsl::simple_sds::load_value<Header>(in);
  header.check();
  gbwt::Tags tags;
  tags.simple_sds_load(in);

  in.close();
  return tags;
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::serialize_to_files(const std::string& gbwt_name, const std::string& graph_name, bool simple_sds_graph) const
{
  sdsl::simple_sds::serialize_to(this->index, gbwt_name);
  if(simple_sds_graph) { sdsl::simple_sds::serialize_to(this->graph, graph_name); }
  else { this->graph.serialize(graph_name); }
}

template <typename SAAllocator>
void
CoreGBZ<SAAllocator>::load_from_files(const std::string& gbwt_name, const std::string& graph_name)
{
  this->header.set_version();
  this->tags.clear();
  this->add_source();
  sdsl::simple_sds::load_from(this->index, gbwt_name);
  this->set_gbwt();
  this->graph.deserialize(graph_name);
}

template class CoreGBZ<std::allocator<char>>;

// CoreGBZ<SharedMemCharAllocatorType> is not explicitly instantiated as a
// whole class: the merge, supergraph, and HandleGraph-import constructors
// only type-check for the plain allocator (they build via operations that
// return a plain CoreGBWTGraph). Only the members actually needed for a
// shared-memory CoreGBZ are instantiated below, so the rest are simply
// never compiled for that allocator -- a compile-time (link-time) lock,
// not a runtime one.
#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY
template CoreGBZ<SharedMemCharAllocatorType>::CoreGBZ(bi::managed_shared_memory*);
template CoreGBZ<SharedMemCharAllocatorType>::CoreGBZ(std::unique_ptr<gbwt::GBWT>&, std::unique_ptr<NaiveGraph>&, bi::managed_shared_memory*);
template CoreGBZ<SharedMemCharAllocatorType>::CoreGBZ(const gbwt::GBWT&, const NaiveGraph&, bi::managed_shared_memory*);
template CoreGBZ<SharedMemCharAllocatorType>::CoreGBZ(const CoreGBZ&);
template CoreGBZ<SharedMemCharAllocatorType>::CoreGBZ(CoreGBZ&&);
template CoreGBZ<SharedMemCharAllocatorType>::~CoreGBZ();
template void CoreGBZ<SharedMemCharAllocatorType>::swap(CoreGBZ&);
template CoreGBZ<SharedMemCharAllocatorType>& CoreGBZ<SharedMemCharAllocatorType>::operator=(const CoreGBZ&);
template CoreGBZ<SharedMemCharAllocatorType>& CoreGBZ<SharedMemCharAllocatorType>::operator=(CoreGBZ&&);
template bool CoreGBZ<SharedMemCharAllocatorType>::compute_pggname(const GraphName*, ParentGraphType);
#endif

//------------------------------------------------------------------------------

} // namespace gbwtgraph
