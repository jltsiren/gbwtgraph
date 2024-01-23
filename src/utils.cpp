#include <gbwtgraph/utils.h>
#include <gbwtgraph/gbwtgraph.h>

#include <algorithm>
#include <sstream>

#include <gbwt/utils.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t Version::MAJOR_VERSION;
constexpr size_t Version::MINOR_VERSION;
constexpr size_t Version::PATCH_VERSION;
constexpr size_t Version::GBZ_VERSION;
constexpr size_t Version::GRAPH_VERSION;
constexpr size_t Version::MINIMIZER_VERSION;

constexpr size_t MetadataBuilder::NO_FIELD;

//------------------------------------------------------------------------------

// Global variables.

const std::string REFERENCE_PATH_SAMPLE_NAME = "_gbwt_ref";

const std::string REFERENCE_SAMPLE_LIST_GBWT_TAG = "reference_samples";
const std::string REFERENCE_SAMPLE_LIST_GFA_TAG = "RS";

// Since spaces are allowed in GFA tags, we can use them as sample name separators.
const char REFERENCE_SAMPLE_LIST_SEPARATOR = ' ';

//------------------------------------------------------------------------------

// Other class variables.

const std::string Version::SOURCE_KEY = "source";
const std::string Version::SOURCE_VALUE = "jltsiren/gbwtgraph";

const std::string SequenceSource::TRANSLATION_EXTENSION = ".trans";

//------------------------------------------------------------------------------

std::unordered_set<std::string>
parse_reference_samples_tag(const char* cursor, const char* end)
{
  std::unordered_set<std::string> reference_samples;

  while(cursor != end)
  {
    // Until we run out of sting, parse out a sample name.
    auto name_start = cursor;
    while (cursor != end && *cursor != REFERENCE_SAMPLE_LIST_SEPARATOR)
    {
      ++cursor;
    }
    auto name_end = cursor;

    // Put it in the set
    reference_samples.emplace(name_start, name_end);

    if(cursor != end)
    {
      // Advance past the comma
      ++cursor;
    }
  }

  return reference_samples;
}

std::unordered_set<std::string>
parse_reference_samples_tag(const std::string& tag_value)
{
  const char* cursor = tag_value.c_str();
  const char* end = cursor + tag_value.size();
  return parse_reference_samples_tag(cursor, end);
}

std::unordered_set<std::string>
parse_reference_samples_tag(view_type tag_value)
{
  return parse_reference_samples_tag(tag_value.first, tag_value.first + tag_value.second);
}

std::unordered_set<std::string>
parse_reference_samples_tag(const gbwt::GBWT& index)
{
  return parse_reference_samples_tag(index.tags.get(REFERENCE_SAMPLE_LIST_GBWT_TAG));
}

std::string
compose_reference_samples_tag(const std::unordered_set<std::string>& reference_samples)
{
  std::stringstream ss;
  for(auto it = reference_samples.begin(); it != reference_samples.end(); ++it)
  {
    // Put the name of evcery reference sample
    ss << *it;

    auto next = it;
    ++next;
    if(next != reference_samples.end())
    {
      // And if it isn't the last one, put a separator after it.
      ss << REFERENCE_SAMPLE_LIST_SEPARATOR;
    }
  }
  return ss.str();
}

PathSense
get_path_sense(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, const std::unordered_set<std::string>& reference_samples)
{
  return get_sample_sense(metadata, path_name.sample, reference_samples);
}

PathSense
get_path_sense(const gbwt::GBWT& index, gbwt::size_type path_number, const std::unordered_set<std::string>& reference_samples)
{
  if(!index.hasMetadata() || !index.metadata.hasPathNames() || path_number >= index.metadata.paths())
  {
    return PathSense::HAPLOTYPE;
  }
  return get_path_sense(index.metadata, index.metadata.path(path_number), reference_samples);
}

PathSense
get_sample_sense(const gbwt::Metadata& metadata, gbwt::size_type sample, const std::unordered_set<std::string>& reference_samples)
{
  if(!metadata.hasSampleNames() || sample >= metadata.sample_names.size())
  {
    // If there are no sample names, everything is a haplotype.
    return PathSense::HAPLOTYPE;
  }
  return get_sample_sense(metadata.sample(sample), reference_samples);
}

PathSense
get_sample_sense(const std::string& sample_name, const std::unordered_set<std::string>& reference_samples)
{
  if(sample_name == REFERENCE_PATH_SAMPLE_NAME)
  {
    // Paths with the magic sample are generic named paths.
    return PathSense::GENERIC;
  }
  if(reference_samples.count(sample_name))
  {
    // Paths with these samples are reference paths
    return PathSense::REFERENCE;
  }
  // Other paths are haplotypes.
  return PathSense::HAPLOTYPE;
}

std::string
get_path_sample_name(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense)
{
  if(sense == PathSense::GENERIC)
  {
    // The libhandlegraph sample name should be the no-sample sentinel
    return PathMetadata::NO_SAMPLE_NAME;
  }
  if(!metadata.hasSampleNames())
  {
    // If there are no sample names, use the sample number.
    return std::to_string(path_name.sample);
  }
  // Othwrwise return what we have stored.
  return metadata.sample(path_name.sample);
}

std::string
get_path_sample_name(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense)
{
  if(!index.hasMetadata() || !index.metadata.hasPathNames() || path_number >= index.metadata.paths())
  {
    return PathMetadata::NO_SAMPLE_NAME;
  }
  return get_path_sample_name(index.metadata, index.metadata.path(path_number), sense);
}

// maybe_unused is a C++17 attribute, but it doesn't hurt to have it in lower stnadards, and it might work.
// We want to pass all the fields to all of these in case we need to evolve storage format later.

std::string
get_path_locus_name(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, [[maybe_unused]] PathSense sense)
{
  if(!metadata.hasContigNames() || path_name.contig > metadata.contig_names.size())
  {
    // If there are no contig names, use the contig number.
    return std::to_string(path_name.contig);
  }
  // This is the same for all senses.
  return metadata.contig(path_name.contig);
}

std::string
get_path_locus_name(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense)
{
  if(!index.hasMetadata() || !index.metadata.hasPathNames() || path_number >= index.metadata.paths())
  {
    return PathMetadata::NO_LOCUS_NAME;
  }
  return get_path_locus_name(index.metadata, index.metadata.path(path_number), sense);
}

size_t
get_path_haplotype([[maybe_unused]] const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense)
{
  if(sense == PathSense::GENERIC)
  {
    // Generic paths aren't allowed haplotypes
    return PathMetadata::NO_HAPLOTYPE;
  }
  // Otherwise it's just stored, but we need to detect the sentinel
  return path_name.phase == gbwtgraph::GBWTGraph::NO_PHASE ? PathMetadata::NO_HAPLOTYPE : path_name.phase;
}

size_t
get_path_haplotype(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense)
{
  if(!index.hasMetadata() || !index.metadata.hasPathNames() || path_number >= index.metadata.paths())
  {
    return PathMetadata::NO_HAPLOTYPE;
  }
  return get_path_haplotype(index.metadata, index.metadata.path(path_number), sense);
}

size_t
get_path_phase_block([[maybe_unused]] const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense)
{
  if(sense == PathSense::HAPLOTYPE)
  {
    // Only haplotype paths have phase blocks
    return path_name.count;
  }
  return PathMetadata::NO_PHASE_BLOCK;
}

size_t
get_path_phase_block(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense)
{
  if(!index.hasMetadata() || !index.metadata.hasPathNames() || path_number >= index.metadata.paths())
  {
    return PathMetadata::NO_PHASE_BLOCK;
  }
  return get_path_phase_block(index.metadata, index.metadata.path(path_number), sense);
}

subrange_t
get_path_subrange([[maybe_unused]] const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense)
{
  subrange_t subrange = PathMetadata::NO_SUBRANGE;
  if(sense != PathSense::HAPLOTYPE && path_name.count != 0)
  {
    // If we aren't a haplotype and we have a nonzero count, we use count to
    // store the subrange start.
    subrange.first = path_name.count;
  }
  return subrange;
}

subrange_t
get_path_subrange(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense)
{
  if(!index.hasMetadata() || !index.metadata.hasPathNames() || path_number >= index.metadata.paths())
  {
    return PathMetadata::NO_SUBRANGE;
  }
  return get_path_subrange(index.metadata, index.metadata.path(path_number), sense);
}

std::string
compose_path_name(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense)
{
  // Put everything together into a string path name the way PathMetadata suggests.
  return PathMetadata::create_path_name(
    sense,
    get_path_sample_name(metadata, path_name, sense),
    get_path_locus_name(metadata, path_name, sense),
    get_path_haplotype(metadata, path_name, sense),
    get_path_phase_block(metadata, path_name, sense),
    get_path_subrange(metadata, path_name, sense)
  );
}

std::string
compose_path_name(const gbwt::GBWT& index, gbwt::size_type path_number, PathSense sense)
{
  if(!index.hasMetadata() || !index.metadata.hasPathNames() || path_number >= index.metadata.paths())
  {
    return "";
  }
  return compose_path_name(index.metadata, index.metadata.path(path_number), sense);
}

void
set_sample_path_senses(gbwt::Tags& tags, const std::unordered_map<std::string, PathSense>& senses)
{
  // Find all the current reference samples
  std::unordered_set<std::string> reference_sample_names = parse_reference_samples_tag(tags.get(REFERENCE_SAMPLE_LIST_GBWT_TAG));

  for(auto& kv : senses)
  {
    if(kv.first == PathMetadata::NO_SAMPLE_NAME)
    {
      // If the sample isn't set, it must be a generic named path.
      if(kv.second != PathSense::GENERIC)
      {
        throw std::runtime_error("Cannot store a sense other than generic in GBWT tags for the no-sample sample.");
      }
      continue;
    }

    if(kv.second == PathSense::REFERENCE)
    {
      // Add new reference samples
      reference_sample_names.insert(kv.first);
    }
    else if(kv.second == PathSense::HAPLOTYPE)
    {
      // And remove any that change sense back to haplotype
      reference_sample_names.erase(kv.first);
    } else {
      // We can't actually set this sense.
      throw std::runtime_error("Cannot store sense " + std::to_string((int)kv.second) + " in GBWT tags for sample " + kv.first);
    }
  }

  tags.set(REFERENCE_SAMPLE_LIST_GBWT_TAG, compose_reference_samples_tag(reference_sample_names));
}

//------------------------------------------------------------------------------

std::string
Version::str(bool verbose)
{
  std::ostringstream ss;
  if(verbose) { ss << "GBWTGraph version "; }
  else { ss << "v"; }
  ss << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_VERSION;
  if(verbose) { ss << " (file format version " << GRAPH_VERSION << ")"; }
  return ss.str();
}

void
Version::print(std::ostream& out, const std::string& tool_name, bool verbose, size_t new_lines)
{
  out << tool_name;
  if(verbose) { out << std::endl; }
  else { out << " "; }
  out << str(verbose);
  for(size_t i = 0; i < new_lines; i++) { out << std::endl; }
}

//------------------------------------------------------------------------------

const std::vector<char> complement =
{
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', '$', '#', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', '-', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',

  'N', 'T', 'V', 'G', 'H', 'N', 'N', 'C',   'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N',
  'N', 'Q', 'Y', 'W', 'A', 'A', 'B', 'S',   'N', 'R', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 't', 'v', 'g', 'h', 'N', 'N', 'c',   'd', 'N', 'N', 'm', 'N', 'k', 'n', 'N',
  'N', 'q', 'y', 'w', 'a', 'a', 'b', 's',   'N', 'r', 'N', 'N', 'N', 'N', 'N', 'N',

  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',

  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

std::string
reverse_complement(const std::string& seq)
{
  std::string result = seq;
  reverse_complement_in_place(result);
  return result;
}

void
reverse_complement_in_place(std::string& seq)
{
  size_t swap_size = seq.size() / 2;
  for(size_t i = 0, j = seq.size() - 1; i < swap_size; i++, j--)
  {
    char tmp = seq[i];
    seq[i] = complement[seq[j]];
    seq[j] = complement[tmp];
  }

  if(seq.size() % 2 != 0)
  {
    seq[swap_size] = complement[seq[swap_size]];
  }
}

//------------------------------------------------------------------------------

void
SequenceSource::swap(SequenceSource& another)
{
  if(&another == this) { return; }

  this->nodes.swap(another.nodes);
  this->sequences.swap(another.sequences);
  this->segment_translation.swap(another.segment_translation);
  std::swap(this->next_id, another.next_id);
}

void
SequenceSource::add_node(nid_t id, const std::string& sequence)
{
  if(sequence.empty()) { return; }
  if(this->nodes.find(id) != this->nodes.end()) { return; }
  size_t offset = this->sequences.size();
  this->sequences.insert(this->sequences.end(), sequence.begin(), sequence.end());
  this->nodes[id] = std::pair<size_t, size_t>(offset, sequence.length());
}

void
SequenceSource::add_node(nid_t id, view_type sequence)
{
  if(sequence.second == 0) { return; }
  if(this->nodes.find(id) != this->nodes.end()) { return; }
  size_t offset = this->sequences.size();
  this->sequences.insert(this->sequences.end(), sequence.first, sequence.first + sequence.second);
  this->nodes[id] = std::pair<size_t, size_t>(offset, sequence.second);
}

std::pair<nid_t, nid_t>
SequenceSource::translate_segment(const std::string& name, view_type sequence, size_t max_length)
{
  if(sequence.second == 0) { return invalid_translation(); }
  auto iter = this->segment_translation.find(name);
  if(iter != this->segment_translation.end())
  {
    return iter->second;
  }

  std::pair<nid_t, nid_t> translation(this->next_id, this->next_id + (sequence.second + max_length - 1) / max_length);
  for(nid_t id = translation.first; id < translation.second; id++)
  {
    size_t offset = (id - translation.first) * max_length;
    size_t length = std::min(max_length, sequence.second - offset);
    this->add_node(id, view_type(sequence.first + offset, length));
  }

  this->segment_translation[name] = translation;
  this->next_id = translation.second;
  return translation;
}

std::pair<nid_t, nid_t>
SequenceSource::force_translate(const std::string& segment_name) const
{
  if(this->uses_translation())
  {
    return this->get_translation(segment_name);
  }
  else
  {
    try
    {
      nid_t id = std::stoul(segment_name);
      return std::pair<nid_t, nid_t>(id, id + 1);
    }
    catch(const std::logic_error&)
    {
      return invalid_translation();
    }
  }
}

std::pair<gbwt::StringArray, sdsl::sd_vector<>>
SequenceSource::invert_translation(const std::function<bool(std::pair<nid_t, nid_t>)>& is_present) const
{
  std::pair<gbwt::StringArray, sdsl::sd_vector<>> result;

  // Invert the translation.
  // This stores semiopen ranges of node identifiers corresponding to segments, and views to their segment names.
  std::vector<std::pair<std::pair<nid_t, nid_t>, view_type>> inverse;
  inverse.reserve(this->segment_translation.size());
  for(auto iter = this->segment_translation.begin(); iter != this->segment_translation.end(); ++iter)
  {
    inverse.emplace_back(iter->second, str_to_view(iter->first));
  }
  gbwt::parallelQuickSort(inverse.begin(), inverse.end());

  // Store the segment names.
  std::string empty;
  result.first = gbwt::StringArray(inverse.size(),
  [&](size_t offset) -> size_t
  {
    // This produces the length of each string to store
    if(is_present(inverse[offset].first)) { return inverse[offset].second.second; }
    else { return 0; }
  },
  [&](size_t offset) -> view_type
  {
    // This produces a view to each string to store.
    if(is_present(inverse[offset].first)) { return inverse[offset].second; }
    else { return str_to_view(empty); }
  });

  // Store the mapping.
  sdsl::sd_vector_builder builder(this->next_id, inverse.size());
  for(auto& translation : inverse)
  {
    builder.set_unsafe(translation.first.first);
  }
  result.second = sdsl::sd_vector<>(builder);

  return result;
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
  // Sanity checks.
  if(!metadata.hasSampleNames() && metadata.samples() > 0) {
    throw std::runtime_error("MetadataBuilder(): Cannot use metadata without sample names");
  }
  if(!metadata.hasContigNames() && metadata.contigs() > 0) {
    throw std::runtime_error("MetadataBuilder(): Cannot use metadata without contig names");
  }
  if(!metadata.hasPathNames() && metadata.paths() > 0) {
    throw std::runtime_error("MetadataBuilder(): Cannot use metadata without path names");
  }

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
    path_name.phase = gbwtgraph::GBWTGraph::NO_PHASE;
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
MetadataBuilder::add_haplotype(const std::string& sample, const std::string& contig, size_t haplotype, size_t fragment, size_t job)
{
  this->add_path(PathSense::HAPLOTYPE, sample, contig, haplotype, fragment, PathMetadata::NO_SUBRANGE, job);
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

void
MetadataBuilder::clear()
{
  this->sample_names = std::map<std::string, size_t>();
  this->contig_names = std::map<std::string, size_t>();
  this->haplotypes = std::set<std::pair<size_t, size_t>>();
  this->path_names = std::vector<std::vector<gbwt::PathName>>();
  this->counts = std::map<gbwt::PathName, size_t>();
}

void
MetadataBuilder::add_path_name(const gbwt::PathName& path_name, size_t job)
{
  if(job >= this->path_names.size())
  {
    this->path_names.resize(job + 1);
  }
  this->path_names[job].push_back(path_name);
}

std::vector<std::string>
MetadataBuilder::map_to_vector(const std::map<std::string, size_t>& source)
{
  std::vector<std::string> result(source.size());
  for(auto& name : source) { result[name.second] = name.first; }
  return result;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph

