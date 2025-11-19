#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr std::uint32_t GBZ::Header::TAG;
constexpr std::uint32_t GBZ::Header::VERSION;

constexpr std::uint64_t GBZ::Header::FLAG_MASK;

//------------------------------------------------------------------------------

// Other class variables.

const std::string GBZ::EXTENSION = ".gbz";

//------------------------------------------------------------------------------

GBZ::Header::Header() :
  tag(TAG), version(VERSION),
  flags(0)
{
}

void
GBZ::Header::check() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("GBZ: Invalid tag");
  }

  if(this->version != VERSION)
  {
    std::string msg = "GBZ: Expected v" + std::to_string(VERSION) + ", got v" + std::to_string(this->version);
    throw sdsl::simple_sds::InvalidData(msg);
  }

  std::uint64_t mask = 0;
  switch(this->version)
  {
  case VERSION:
    mask = FLAG_MASK; break;
  }
  if((this->flags & mask) != this->flags)
  {
    throw sdsl::simple_sds::InvalidData("GBZ: Invalid flags");
  }
}

bool
GBZ::Header::operator==(const Header& another) const
{
  return (this->tag == another.tag && this->version == another.version &&
          this->flags == another.flags);
}

//------------------------------------------------------------------------------

size_t
GBZ::set_reference_samples(const std::unordered_set<std::string>& samples)
{
  const gbwt::Metadata& metadata = this->index.metadata;
  std::unordered_set<std::string> present_samples;
  for(const std::string& sample : samples)
  {
    if(metadata.sample(sample) < metadata.samples())
    {
      present_samples.insert(sample);
    }
  }

  std::string tag_value = compose_reference_samples_tag(present_samples);
  this->index.tags.set(REFERENCE_SAMPLE_LIST_GBWT_TAG, tag_value);

  // By resetting the pointer to the GBWT index, we tell the GBWTGraph object
  // to pull the tag and recache named paths.
  this->set_gbwt();

  return present_samples.size();
}

//------------------------------------------------------------------------------

GBZ::GBZ()
{
  this->add_source();
  this->set_gbwt();
}

GBZ::GBZ(const GBZ& source)
{
  this->copy(source);
}

GBZ::GBZ(GBZ&& source)
{
  *this = std::move(source);
}

GBZ::~GBZ()
{
}

void
GBZ::swap(GBZ& another)
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

GBZ&
GBZ::operator=(const GBZ& source)
{
  if(&source != this) { this->copy(source); }
  return *this;
}

GBZ&
GBZ::operator=(GBZ&& source)
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

void
GBZ::copy(const GBZ& source)
{
  this->header = source.header;
  this->tags = source.tags;
  this->index = source.index;
  this->graph = source.graph;

  // Use the local copy of the GBWT.
  this->set_gbwt_address();
}

void
GBZ::reset_tags()
{
  this->tags.clear();
  this->add_source();
}

void
GBZ::add_source()
{
  this->tags.set(Version::SOURCE_KEY, Version::SOURCE_VALUE);
}

//------------------------------------------------------------------------------

GBZ::GBZ(std::unique_ptr<gbwt::GBWT>& index, std::unique_ptr<SequenceSource>& source)
{
  if(index == nullptr || source == nullptr)
  {
    throw std::runtime_error("GBZ: Index and sequence source must be non-null");
  }

  this->add_source();
  this->index = std::move(*index); index.reset();
  GraphName parent = source->graph_name();
  this->graph = GBWTGraph(this->index, *source); source.reset();
  this->compute_pggname(&parent);
}

GBZ::GBZ(const gbwt::GBWT& index, const SequenceSource& source) :
  index(index), graph(this->index, source)
{
  this->add_source();
  GraphName parent = source.graph_name();
  this->compute_pggname(&parent);
}

GBZ::GBZ(gbwt::GBWT&& index, const GBZ& supergraph) :
  index(std::move(index)), graph(supergraph.graph.subgraph(this->index))
{
  this->add_source();
  GraphName parent = supergraph.graph_name();
  this->compute_pggname(&parent, ParentGraphType::SUPERGRAPH);

}

GBZ::GBZ(gbwt::GBWT&& index, const HandleGraph& source) :
  index(index), graph(this->index, source)
{
  this->add_source();
}

void
GBZ::set_gbwt()
{
  this->graph.set_gbwt(this->index);
}

void
GBZ::set_gbwt_address()
{
  this->graph.set_gbwt_address(this->index);
}

//------------------------------------------------------------------------------

bool
GBZ::compute_pggname(const GraphName* parent, ParentGraphType relationship)
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
        name.add_translation(view_type(name.name()), view_type(parent->name()));
        this->tags.set(GraphName::GBZ_TRANSLATION_TARGET_TAG, parent->name());
      }
    }
    else
    {
      // This does nothing if the names are the same.
      name.add_subgraph(view_type(name.name()), view_type(parent->name()));
    }
    name.add_relationships(*parent);
  }

  // Store the information back into the tags.
  name.set_tags(this->tags);

  return true;
}

//------------------------------------------------------------------------------

void
GBZ::simple_sds_serialize(std::ostream& out) const
{
  sdsl::simple_sds::serialize_value(this->header, out);
  this->tags.simple_sds_serialize(out);
  this->index.simple_sds_serialize(out);
  this->graph.simple_sds_serialize(out);
}

void
GBZ::simple_sds_serialize(const gbwt::GBWT& index, const GBWTGraph& graph, std::ostream& out)
{
  GBZ empty;
  sdsl::simple_sds::serialize_value(empty.header, out);
  empty.tags.simple_sds_serialize(out);
  index.simple_sds_serialize(out);
  graph.simple_sds_serialize(out);
}

void
GBZ::simple_sds_load(std::istream& in)
{
  this->header = sdsl::simple_sds::load_value<Header>(in);
  this->header.check();

  // Load the tags and update the source to this library.
  // We could also check if the source was already this library, but we have no
  // uses for that information at the moment.
  this->tags.simple_sds_load(in);
  this->add_source();

  this->index.simple_sds_load(in);
  this->graph.simple_sds_load(in, this->index);
}

size_t
GBZ::simple_sds_size() const
{
  size_t result = sdsl::simple_sds::value_size(this->header);
  result += this->tags.simple_sds_size();
  result += this->index.simple_sds_size();
  result += this->graph.simple_sds_size();
  return result;
}

void
GBZ::serialize_to_files(const std::string& gbwt_name, const std::string& graph_name, bool simple_sds_graph) const
{
  sdsl::simple_sds::serialize_to(this->index, gbwt_name);
  if(simple_sds_graph) { sdsl::simple_sds::serialize_to(this->graph, graph_name); }
  else { this->graph.serialize(graph_name); }
}

void
GBZ::load_from_files(const std::string& gbwt_name, const std::string& graph_name)
{
  this->tags.clear();
  this->add_source();
  sdsl::simple_sds::load_from(this->index, gbwt_name);
  this->set_gbwt();
  this->graph.deserialize(graph_name);
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
