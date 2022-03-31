#include <gbwtgraph/internal.h>

#include <algorithm>
#include <cctype>
#include <functional>

namespace gbwtgraph
{

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
TSVWriter::write(view_type view)
{
  size_t offset = 0;
  while(offset < view.second)
  {
    size_t length = std::min(view.second - offset, BUFFER_SIZE - this->buffer.size());
    this->buffer.insert(this->buffer.end(), view.first + offset, view.first + offset + length);
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

MetadataBuilder::MetadataBuilder(const std::string& path_name_regex, const std::string& path_name_fields) :
  sample_field(NO_FIELD), contig_field(NO_FIELD), haplotype_field(NO_FIELD), fragment_field(NO_FIELD),
  ref_path_sample_warning(false)
{
  // Initialize the regex.
  try { this->parser = std::regex(path_name_regex); }
  catch(std::regex_error& e)
  {
    throw std::runtime_error("MetadataBuilder: Invalid regex: " + path_name_regex);
  }
  if(path_name_fields.size() > this->parser.mark_count() + 1)
  {
    throw std::runtime_error("MetadataBuilder: Field string too long: " + path_name_fields);
  }

  // Initialize the fields.
  for(size_t i = 0; i < path_name_fields.size(); i++)
  {
    switch(std::tolower(path_name_fields[i]))
    {
      case 's':
        if(this->sample_field != NO_FIELD)
        {
          throw std::runtime_error("MetadataBuilder: Duplicate sample field");
        }
        this->sample_field = i;
        break;
      case 'c':
        if(this->contig_field != NO_FIELD)
        {
          throw std::runtime_error("MetadataBuilder: Duplicate contig field");
        }
        this->contig_field = i;
        break;
      case 'h':
        if(this->haplotype_field != NO_FIELD)
        {
          throw std::runtime_error("MetadataBuilder: Duplicate haplotype field");
        }
        this->haplotype_field = i;
        break;
      case 'f':
        if(this->fragment_field != NO_FIELD)
        {
          throw std::runtime_error("MetadataBuilder: Duplicate fragment field");
        }
        this->fragment_field = i;
        break;
    }
  }
}

void
MetadataBuilder::parse(const std::string& name, size_t job)
{
  std::smatch fields;
  if(!std::regex_match(name, fields, this->parser))
  {
    throw std::runtime_error("MetadataBuilder: Cannot parse path name " + name);
  }

  gbwt::PathName path_name =
  {
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0)
  };

  if(this->sample_field != NO_FIELD)
  {
    std::string sample_name = fields[this->sample_field];
    if(!(this->ref_path_sample_warning) && sample_name == REFERENCE_PATH_SAMPLE_NAME)
    {
      std::cerr << "MetadataBuilder::parse(): Warning: Sample name " << REFERENCE_PATH_SAMPLE_NAME << " is reserved for named paths" << std::endl;
      this->ref_path_sample_warning = true;
    }
    auto iter = this->sample_names.find(sample_name);
    if(iter == this->sample_names.end())
    {
      path_name.sample = this->sample_names.size();
      this->sample_names[sample_name] = path_name.sample;
    }
    else { path_name.sample = iter->second; }
  }

  if(this->contig_field != NO_FIELD)
  {
    std::string contig_name = fields[this->contig_field];
    auto iter = this->contig_names.find(contig_name);
    if(iter == this->contig_names.end())
    {
      path_name.contig = this->contig_names.size();
      this->contig_names[contig_name] = path_name.contig;
    }
    else { path_name.contig = iter->second; }
  }

  if(this->haplotype_field != NO_FIELD)
  {
    try { path_name.phase = std::stoul(fields[this->haplotype_field]); }
    catch(const std::invalid_argument&)
    {
      throw std::runtime_error("MetadataBuilder: Invalid haplotype field " + fields[this->haplotype_field].str());
    }
  }
  this->haplotypes.insert(std::pair<size_t, size_t>(path_name.sample, path_name.phase));

  if(this->fragment_field != NO_FIELD)
  {
    try { path_name.count = std::stoul(fields[this->fragment_field]); }
    catch(const std::invalid_argument&)
    {
      throw std::runtime_error("MetadataBuilder: Invalid fragment field " + fields[this->fragment_field].str());
    }
    if(this->counts.find(path_name) != this->counts.end())
    {
      throw std::runtime_error("MetadataBuilder: Duplicate path name " + name);
    }
    this->counts[path_name] = 1;
  }
  else
  {
    auto iter = this->counts.find(path_name);
    if(iter == this->counts.end()) { this->counts[path_name] = 1; }
    else { path_name.count = iter->second; iter->second++; }
  }

  this->add_path_name(path_name, job);
}

void
MetadataBuilder::add_walk(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start, size_t job)
{
  gbwt::PathName path_name =
  {
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0)
  };

  // Sample.
  {
    if(!(this->ref_path_sample_warning) && sample == REFERENCE_PATH_SAMPLE_NAME)
    {
      std::cerr << "MetadataBuilder::add_walk(): Warning: Sample name " << REFERENCE_PATH_SAMPLE_NAME << " is reserved for named paths" << std::endl;
      this->ref_path_sample_warning = true;
    }
    auto iter = this->sample_names.find(sample);
    if(iter == this->sample_names.end())
    {
      path_name.sample = this->sample_names.size();
      this->sample_names[sample] = path_name.sample;
    }
    else { path_name.sample = iter->second; }
  }

  // Contig.
  {
    auto iter = this->contig_names.find(contig);
    if(iter == this->contig_names.end())
    {
      path_name.contig = this->contig_names.size();
      this->contig_names[contig] = path_name.contig;
    }
    else { path_name.contig = iter->second; }
  }

  // Haplotype.
  {
    try { path_name.phase = std::stoul(haplotype); }
    catch(const std::invalid_argument&)
    {
      throw std::runtime_error("MetadataBuilder: Invalid haplotype field " + haplotype);
    }
  }
  this->haplotypes.insert(std::pair<size_t, size_t>(path_name.sample, path_name.phase));

  // Start position as fragment identifier.
  {
    try { path_name.count = std::stoul(start); }
    catch(const std::invalid_argument&)
    {
      throw std::runtime_error("MetadataBuilder: Invalid start position " + start);
    }
    if(this->counts.find(path_name) != this->counts.end())
    {
      throw std::runtime_error("MetadataBuilder: Duplicate walk " + sample + "\t" + haplotype + "\t" + contig + "\t" + start);
    }
    this->counts[path_name] = 1;
  }

  this->add_path_name(path_name, job);
}

void
MetadataBuilder::add_named_path(const std::string& name, size_t job)
{
  gbwt::PathName path_name =
  {
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0)
  };

  // Sample.
  {
    auto iter = this->sample_names.find(REFERENCE_PATH_SAMPLE_NAME);
    if(iter == this->sample_names.end())
    {
      path_name.sample = this->sample_names.size();
      this->sample_names[REFERENCE_PATH_SAMPLE_NAME] = path_name.sample;
    }
    else { path_name.sample = iter->second; }
  }

  // Contig.
  {
    auto iter = this->contig_names.find(name);
    if(iter == this->contig_names.end())
    {
      path_name.contig = this->contig_names.size();
      this->contig_names[name] = path_name.contig;
    }
    else { path_name.contig = iter->second; }
  }

  // Haplotype.
  this->haplotypes.insert(std::pair<size_t, size_t>(path_name.sample, path_name.phase));

  if(this->counts.find(path_name) != this->counts.end())
  {
    throw std::runtime_error("MetadataBuilder: Duplicate named path " + name);
  }
  this->counts[path_name] = 1;

  this->add_path_name(path_name, job);
}

gbwt::Metadata
MetadataBuilder::get_metadata() const
{
  gbwt::Metadata metadata;

  if(this->sample_names.empty()) { metadata.setSamples(1); }
  else { metadata.setSamples(map_to_vector(this->sample_names)); }
  metadata.setHaplotypes(this->haplotypes.size());
  if(this->contig_names.empty()) { metadata.setContigs(1); }
  else { metadata.setContigs(map_to_vector(this->contig_names)); }
  for(auto& names_by_job : this->path_names)
  {
    for(auto& path_name : names_by_job) { metadata.addPath(path_name); }
  }

  return metadata;
}

//------------------------------------------------------------------------------

void
EmptyGraph::create_node(nid_t node_id)
{
  this->nodes[node_id] = { {}, {} };
  this->min_id = std::min(this->min_id, node_id);
  this->max_id = std::max(this->max_id, node_id);
}

void
EmptyGraph::create_edge(const handle_t& from, const handle_t& to)
{
  auto from_iter = this->get_node_mut(from);
  auto to_iter = this->get_node_mut(to);
  if(from_iter == this->nodes.end() || to_iter == this->nodes.end())
  {
    nid_t from_id = gbwt::Node::id(handle_to_node(from));
    nid_t to_id = gbwt::Node::id(handle_to_node(to));
    throw std::runtime_error("EmptyGraph: Cannot create an edge between nodes " + std::to_string(from_id) + " and " + std::to_string(to_id));
  }

  // from -> to
  if(this->get_is_reverse(from))
  {
    from_iter->second.predecessors.push_back(this->flip(to));
  }
  else
  {
    from_iter->second.successors.push_back(to);
  }

  // to -> from
  if(this->get_is_reverse(to))
  {
    to_iter->second.successors.push_back(this->flip(from));
  }
  else
  {
    to_iter->second.predecessors.push_back(from);
  }
}

void
EmptyGraph::remove_duplicate_edges()
{
  this->for_each_handle([&](const handle_t& handle) {
    auto iter = this->get_node_mut(handle);
    gbwt::removeDuplicates(iter->second.predecessors, false);
    gbwt::removeDuplicates(iter->second.successors, false);
  });
}

bool
EmptyGraph::has_node(nid_t node_id) const
{
  return (this->nodes.find(node_id) != this->nodes.end());
}

handle_t
EmptyGraph::get_handle(const nid_t& node_id, bool is_reverse) const
{
  return node_to_handle(gbwt::Node::encode(node_id, is_reverse));
}

nid_t
EmptyGraph::get_id(const handle_t& handle) const
{
  return gbwt::Node::id(handle_to_node(handle));
}

bool
EmptyGraph::get_is_reverse(const handle_t& handle) const
{
  return gbwt::Node::is_reverse(handle_to_node(handle));
}

handle_t
EmptyGraph::flip(const handle_t& handle) const
{
  return node_to_handle(gbwt::Node::reverse(handle_to_node(handle)));
}

size_t
EmptyGraph::get_length(const handle_t&) const
{
  return 0;
}

std::string
EmptyGraph::get_sequence(const handle_t&) const
{
  return std::string();
}

char
EmptyGraph::get_base(const handle_t&, size_t) const
{
  return 'N';
}

std::string
EmptyGraph::get_subsequence(const handle_t&, size_t, size_t) const
{
  return std::string();
}

size_t
EmptyGraph::get_node_count() const
{
  return this->nodes.size();
}

nid_t
EmptyGraph::min_node_id() const
{
  return this->min_id;
}

nid_t
EmptyGraph::max_node_id() const
{
  return this->max_id;
}

bool
EmptyGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const
{
  auto iter = this->get_node(handle);
  bool flip = this->get_is_reverse(handle);
  const std::vector<handle_t>& edges = (go_left ^ flip ? iter->second.predecessors : iter->second.successors);
  for(const handle_t& next : edges)
  {
    handle_t actual = (flip ? this->flip(next) : next);
    if(!iteratee(actual)) { return false; }
  }
  return true;
}

bool
EmptyGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool) const
{
  for(auto iter = this->nodes.begin(); iter != this->nodes.end(); ++iter)
  {
    if(!iteratee(this->get_handle(iter->first, false))) { return false; }
  }
  return true;
}

size_t
EmptyGraph::get_degree(const handle_t& handle, bool go_left) const
{
  auto iter = this->get_node(handle);
  bool flip = this->get_is_reverse(handle);
  const std::vector<handle_t>& edges = (go_left ^ flip ? iter->second.predecessors : iter->second.successors);
  return edges.size();
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

} // namespace gbwtgraph
