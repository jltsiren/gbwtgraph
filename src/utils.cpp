#include <gbwtgraph/utils.h>

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
get_path_sense(const std::unordered_set<std::string>& reference_samples, const gbwt::Metadata& metadata, const gbwt::PathName& path_name)
{
  return get_path_sense(reference_samples, metadata, path_name.sample);
}

PathSense
get_path_sense(const std::unordered_set<std::string>& reference_samples, const gbwt::Metadata& metadata, gbwt::size_type sample)
{
  return get_path_sense(reference_samples, metadata.sample(sample));
}

PathSense
get_path_sense(const std::unordered_set<std::string>& reference_samples, const std::string& sample_name)
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
  // Othwrwise return what we have stored.
  return metadata.sample(path_name.sample);
}

// maybe_unused is a C++17 attribute, but it doesn't hurt to have it in lower stnadards, and it might work.
// We want to pass all the fields to all of these in case we need to evolve storage format later.

std::string
get_path_locus_name(const gbwt::Metadata& metadata, const gbwt::PathName& path_name, [[maybe_unused]] PathSense sense)
{
  // This is the same for all senses.
  return metadata.contig(path_name.contig);
}

size_t
get_path_haplotype([[maybe_unused]] const gbwt::Metadata& metadata, const gbwt::PathName& path_name, PathSense sense)
{
  if(sense == PathSense::GENERIC)
  {
    // Generic paths aren't allowed haplotypes
    return PathMetadata::NO_HAPLOTYPE;
  }
  // Otherwise it's just stored
  return path_name.phase;
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

} // namespace gbwtgraph

