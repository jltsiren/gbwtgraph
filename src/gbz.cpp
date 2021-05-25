#include <gbwtgraph/gbz.h>

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

bool
GBZ::Header::check() const
{
  if(this->tag != TAG) { return false; }
  switch(this->version)
  {
  case VERSION:
    return ((this->flags & FLAG_MASK) == this->flags);
  default:
    return false;
  }
}

bool
GBZ::Header::operator==(const Header& another) const
{
  return (this->tag == another.tag && this->version == another.version &&
          this->flags == another.flags);
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
  this->set_gbwt();
  another.set_gbwt();
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
    this->set_gbwt();
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
  this->set_gbwt();
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
    throw std::invalid_argument("GBZ: Index and sequence source must be non-null");
  }

  this->add_source();
  this->index = std::move(*index); index.reset();
  this->graph = GBWTGraph(this->index, *source); source.reset();
}

GBZ::GBZ(std::unique_ptr<gbwt::GBWT>& index, const HandleGraph& source)
{
  if(index == nullptr)
  {
    throw std::invalid_argument("GBZ: Index must be non-null");
  }

  this->add_source();
  this->index = std::move(*index); index.reset();
  this->graph = GBWTGraph(this->index, source);
}

GBZ::GBZ(const gbwt::GBWT& index, const SequenceSource& source) :
  index(index), graph(this->index, source)
{
  this->add_source();
}

void
GBZ::set_gbwt()
{
  this->graph.set_gbwt(this->index);
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
GBZ::simple_sds_load(std::istream& in)
{
  this->header = sdsl::simple_sds::load_value<Header>(in);
  if(!(this->header.check()))
  {
    throw sdsl::simple_sds::InvalidData("GBZ: Invalid header");
  }

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
GBZ::sdsl_serialize(const std::string& gbwt_name, const std::string& graph_name) const
{
  if(!sdsl::store_to_file(this->index, gbwt_name))
  {
    std::cerr << "GBZ: Cannot write GBWT to " << gbwt_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->graph.serialize(graph_name);
}

void
GBZ::sdsl_load(const std::string& gbwt_name, const std::string& graph_name)
{
  if(!sdsl::load_from_file(this->index, gbwt_name))
  {
    std::cerr << "GBZ: Cannot load GBWT from " << gbwt_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->set_gbwt();
  this->graph.deserialize(graph_name);
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
