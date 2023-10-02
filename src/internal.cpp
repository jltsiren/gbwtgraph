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

MetadataBuilder::PathMetadataBuilder::PathMetadataBuilder(const std::string& path_name_regex, const std::string& path_name_fields, PathSense path_sense) :
  sample_field(NO_FIELD), contig_field(NO_FIELD), haplotype_field(NO_FIELD), fragment_field(NO_FIELD), sense(path_sense)
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

MetadataBuilder::MetadataBuilder() :
  ref_path_sample_warning(false)
{
}

MetadataBuilder::MetadataBuilder(const gbwt::Metadata& metadata) :
  path_names{metadata.path_names},
  ref_path_sample_warning(false)
{
  for(size_t i = 0; i < metadata.sample_names.size(); i++)
  {
    this->sample_names.emplace(metadata.sample(i), i);
  }

  for(size_t i = 0; i < metadata.contig_names.size(); i++)
  {
    this->contig_names.emplace(metadata.contig(i), i);
  }

  for(gbwt::PathName copy : metadata.path_names)
  {
    // Record all the phases
    this->haplotypes.emplace(copy.sample, copy.phase);

    {
      // Count the actual observed path name
      auto& self_count = this->counts[copy];
      self_count = std::max<size_t>(self_count, 1);
    }

    {
      // Make sure the count stored for the no-count version of the path name
      // will not generate this path name again when inferring counts.
      size_t observed_count = copy.count;
      copy.count = 0;
      auto& guess_count = this->counts[copy];
      guess_count = std::max<size_t>(guess_count, observed_count + 1);
    }
  }
}

MetadataBuilder::MetadataBuilder(const std::string& path_name_regex, const std::string& path_name_fields, PathSense path_sense) :
  ref_path_sample_warning(false)
{
  this->add_path_name_format(path_name_regex, path_name_fields, path_sense);
}

void
MetadataBuilder::add_path_name_format(const std::string& path_name_regex, const std::string& path_name_fields, PathSense path_sense)
{
  this->path_name_formats.emplace_back(path_name_regex, path_name_fields, path_sense);
}

void
MetadataBuilder::add_path(PathSense sense, const std::string& sample_name, const std::string& locus_name, size_t haplotype, size_t phase_block, const handlegraph::subrange_t& subrange, size_t job)
{
  gbwt::PathName path_name =
  {
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0),
    static_cast<gbwt::PathName::path_name_type>(0)
  };

  if(sample_name != PathMetadata::NO_SAMPLE_NAME || sense == PathSense::GENERIC)
  {
    // We need sample name metadata.

    // If using generic sense, use the magic sample name.
    auto& sample_name_to_store = (sense == PathSense::GENERIC) ? REFERENCE_PATH_SAMPLE_NAME : sample_name;
    // Apply the sample name.
    auto iter = this->sample_names.find(sample_name_to_store);
    if(iter == this->sample_names.end())
    {
      path_name.sample = this->sample_names.size();
      this->sample_names[sample_name_to_store] = path_name.sample;
    }
    else { path_name.sample = iter->second; }
  }

  if(locus_name != PathMetadata::NO_LOCUS_NAME)
  {
    // Apply the locus as the contig name
    auto iter = this->contig_names.find(locus_name);
    if(iter == this->contig_names.end())
    {
      path_name.contig = this->contig_names.size();
      this->contig_names[locus_name] = path_name.contig;
    }
    else { path_name.contig = iter->second; }
  }

  if(haplotype == PathMetadata::NO_HAPLOTYPE)
  {
    // Record a sentinel phase number
    path_name.phase = std::numeric_limits<gbwt::PathName::path_name_type>::max();
  }
  else
  {
    // Apply the phase number
    path_name.phase = haplotype;
  }

  // Remember that this phase of this sample exists.
  this->haplotypes.insert(std::pair<size_t, size_t>(path_name.sample, path_name.phase));

  if(phase_block != PathMetadata::NO_PHASE_BLOCK)
  {
    // Use the phase block as the count
    path_name.count = phase_block;
  }
  else if(subrange != PathMetadata::NO_SUBRANGE)
  {
    // Use the subrange start as the count
    path_name.count = subrange.first;
  }

  {
    auto iter = this->counts.find(path_name);
    if(iter != this->counts.end())
    {
      // This path is duplicate or ambiguous.

      // We need to describe it for reporting.
      std::stringstream ss;
      if(sample_name != PathMetadata::NO_SAMPLE_NAME) { ss << " sample " << sample_name; }
      if(locus_name != PathMetadata::NO_LOCUS_NAME) { ss << " contig " << locus_name; }
      if(haplotype != PathMetadata::NO_HAPLOTYPE) { ss << " phase " << haplotype; }
      if(phase_block != PathMetadata::NO_PHASE_BLOCK) { ss << " count " << phase_block; }
      if(subrange != PathMetadata::NO_SUBRANGE)
      {
        ss << "range [" << subrange.first;
        if(subrange.second != PathMetadata::NO_END_POSITION) { ss << "-" << subrange.second; }
        ss << "]";
      }

      if(phase_block != PathMetadata::NO_PHASE_BLOCK || subrange != PathMetadata::NO_SUBRANGE)
      {
        // The count was user-specified, so bail out.
        throw std::runtime_error("MetadataBuilder: Duplicate path for " + ss.str());
      }
      else
      {
        // We are going to infer a count; just warn.
        std::cerr << "MetadataBuilder::add_path(): Warning: Path for " << ss.str() << " is not unique" << std::endl;
        std::cerr << "MetadataBuilder::add_path(): Warning: Using fragment/count field to disambiguate" << std::endl;
        std::cerr << "MetadataBuilder::add_path(): Warning: Decompression may not produce valid GFA" << std::endl;
        path_name.count = iter->second;
        iter->second++;
      }
    }
    else
    {
      // This path is unique but the next one that looks like it won't be.
      this->counts.emplace_hint(iter, path_name, 1);
    }
  }

  // Record the structured path name in the metadata, referencing the strings
  // we have saved.
  this->add_path_name(path_name, job);
}

void
MetadataBuilder::add_path(const std::string& name, size_t job)
{
  for(auto& format : this->path_name_formats)
  {
    std::smatch fields;
    if(!std::regex_match(name, fields, format.parser))
    {
      continue;
    }

    std::string sample_name = PathMetadata::NO_SAMPLE_NAME;
    if(format.sample_field != NO_FIELD)
    {
      sample_name = fields[format.sample_field];
      if(!(this->ref_path_sample_warning) && sample_name == REFERENCE_PATH_SAMPLE_NAME)
      {
        std::cerr << "MetadataBuilder::add_path(): Warning: Sample name " << REFERENCE_PATH_SAMPLE_NAME << " is reserved" << std::endl;
        this->ref_path_sample_warning = true;
      }
    }

    std::string locus_name = PathMetadata::NO_LOCUS_NAME;
    if(format.contig_field != NO_FIELD)
    {
      locus_name = fields[format.contig_field];
    }

    size_t haplotype = PathMetadata::NO_HAPLOTYPE;
    if(format.haplotype_field != NO_FIELD)
    {
      try { haplotype = std::stoul(fields[format.haplotype_field]); }
      catch(const std::invalid_argument&)
      {
        throw std::runtime_error("MetadataBuilder: Invalid haplotype field " + fields[format.haplotype_field].str());
      }
    }

    size_t phase_block = PathMetadata::NO_PHASE_BLOCK;
    if(format.fragment_field != NO_FIELD)
    {
      try { phase_block = std::stoul(fields[format.fragment_field]); }
      catch(const std::invalid_argument&)
      {
        throw std::runtime_error("MetadataBuilder: Invalid fragment field " + fields[format.fragment_field].str());
      }
    }

    this->add_path(format.sense, sample_name, locus_name, haplotype, phase_block, PathMetadata::NO_SUBRANGE, job);
    return;
  }
  throw std::runtime_error("MetadataBuilder: Cannot parse path name " + name);
}

void
MetadataBuilder::add_walk(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start, size_t job)
{
  // Check sample name.
  if(!(this->ref_path_sample_warning) && sample == REFERENCE_PATH_SAMPLE_NAME)
  {
    std::cerr << "MetadataBuilder::add_walk(): Warning: Sample name " << REFERENCE_PATH_SAMPLE_NAME << " is reserved for generic paths" << std::endl;
    this->ref_path_sample_warning = true;
  }


  if(sample == "*")
  {
    // Treat this as an elided sample name and a generic path.

    subrange_t subrange = PathMetadata::NO_SUBRANGE;
    try { subrange.first = std::stoul(start); }
    catch(const std::invalid_argument&)
    {
      throw std::runtime_error("MetadataBuilder: Invalid start position " + start);
    }

    this->add_path(PathSense::GENERIC, PathMetadata::NO_SAMPLE_NAME, contig, PathMetadata::NO_HAPLOTYPE, PathMetadata::NO_PHASE_BLOCK, subrange, job);
  }
  else
  {
    // Parse the haplotype
    size_t haplotype_number = PathMetadata::NO_HAPLOTYPE;
    try { haplotype_number = std::stoul(haplotype); }
    catch(const std::invalid_argument&)
    {
      throw std::runtime_error("MetadataBuilder: Invalid haplotype field " + haplotype);
    }

    // Start position as fragment identifier.
    size_t phase_block = PathMetadata::NO_PHASE_BLOCK;
    try { phase_block = std::stoul(start); }
    catch(const std::invalid_argument&)
    {
      throw std::runtime_error("MetadataBuilder: Invalid start position " + start);
    }

    // Add as a haplotype
    this->add_path(PathSense::HAPLOTYPE, sample, contig, haplotype_number, phase_block, PathMetadata::NO_SUBRANGE, job);
  }
}

void
MetadataBuilder::add_generic_path(const std::string& name, size_t job)
{
  this->add_path(PathSense::GENERIC, PathMetadata::NO_SAMPLE_NAME, name, PathMetadata::NO_HAPLOTYPE, PathMetadata::NO_PHASE_BLOCK, PathMetadata::NO_SUBRANGE, job);
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
