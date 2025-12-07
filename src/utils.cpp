#include <gbwtgraph/utils.h>
#include <gbwtgraph/gbwtgraph.h>

#include <algorithm>
#include <deque>
#include <sstream>
#include <unordered_map>

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

constexpr char GraphName::GFA_HEADER_PREFIX;
constexpr char GraphName::GAF_HEADER_PREFIX;
constexpr char GraphName::GFA_GAF_FIELD_SEPARATOR;
constexpr char GraphName::GFA_GAF_TAG_SEPARATOR;
constexpr char GraphName::GFA_GAF_TAG_STR_TYPE;
constexpr char GraphName::RELATIONSHIP_SEPARATOR;
constexpr char GraphName::RELATIONSHIP_LIST_SEPARATOR;

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

const std::string GraphName::GBZ_NAME_TAG = "pggname";
const std::string GraphName::GBZ_SUBGRAPH_TAG = "subgraph";
const std::string GraphName::GBZ_TRANSLATION_TAG = "translation";
const std::string GraphName::GBZ_TRANSLATION_TARGET_TAG = "translation_target";
const std::string GraphName::GFA_NAME_TAG = "NM";
const std::string GraphName::GAF_NAME_TAG = "RN";
const std::string GraphName::GFA_GAF_SUBGRAPH_TAG = "SG";
const std::string GraphName::GFA_GAF_TRANSLATION_TAG = "TL";

//------------------------------------------------------------------------------

std::vector<view_type>
split_view(view_type str_view, char separator)
{
  std::vector<view_type> result;
  const char* cursor = str_view.first;
  const char* end = str_view.first + str_view.second;

  while(cursor != end)
  {
    // Find the next occurrence of the separator.
    const char* next = std::find(cursor, end, separator);
    result.emplace_back(cursor, next - cursor);
    cursor = next;
    if(cursor != end)
    {
      // Advance past the separator.
      ++cursor;
    }
  }

  return result;
}

//------------------------------------------------------------------------------

std::unordered_set<std::string>
parse_reference_samples_tag(const char* cursor, const char* end)
{
  return parse_reference_samples_tag(view_type(cursor, end - cursor));
}

std::unordered_set<std::string>
parse_reference_samples_tag(const std::string& tag_value)
{
  return parse_reference_samples_tag(view_type(tag_value));
}

std::unordered_set<std::string>
parse_reference_samples_tag(view_type tag_value)
{
  std::vector<view_type> sample_names = split_view(tag_value, REFERENCE_SAMPLE_LIST_SEPARATOR);
  std::unordered_set<std::string> reference_samples;
  for(const auto& sample_view : sample_names)
  {
    reference_samples.emplace(sample_view.first, sample_view.second);
  }
  return reference_samples;
}

std::unordered_set<std::string>
parse_reference_samples_tag(const gbwt::GBWT& index)
{
  return parse_reference_samples_tag(index.tags.get(REFERENCE_SAMPLE_LIST_GBWT_TAG));
}

std::string
compose_reference_samples_tag(const std::unordered_set<std::string>& reference_samples)
{
  // We sort the sample names to make the output deterministic across
  // standard library implementations.
  std::vector<view_type> sorted_sample_names;
  sorted_sample_names.reserve(reference_samples.size());
  for(const auto& sample_name : reference_samples)
  {
    sorted_sample_names.push_back(view_type(sample_name));
  }
  std::sort(sorted_sample_names.begin(), sorted_sample_names.end());

  std::string result;
  for(size_t i = 0; i < sorted_sample_names.size(); i++)
  {
    const auto& sample_view = sorted_sample_names[i];
    result.append(sample_view.first, sample_view.second);
    if(i + 1 < sorted_sample_names.size())
    {
      result.push_back(REFERENCE_SAMPLE_LIST_SEPARATOR);
    }
  }

  return result;
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

DigestBuf::DigestBuf(const EVP_MD* algorithm) :
  std::streambuf(),
  context(EVP_MD_CTX_new())
{
  if(this->context != nullptr && EVP_DigestInit_ex(this->context, algorithm, nullptr) == -1)
  {
    EVP_MD_CTX_free(this->context);
    this->context = nullptr;
  }
}

DigestBuf::~DigestBuf()
{
  if(this->context != nullptr)
  {
    EVP_MD_CTX_free(this->context);
    this->context = nullptr;
  }
}

int
DigestBuf::overflow(int_type ch)
{
  if(this->context == nullptr) { return traits_type::eof(); }
  if(ch != traits_type::eof())
  {
    unsigned char buf = static_cast<unsigned char>(ch);
    if(EVP_DigestUpdate(this->context, &buf, sizeof(buf)) == -1) { return traits_type::eof(); }
  }
  return ch;
}

std::streamsize
DigestBuf::xsputn(const char* s, std::streamsize n)
{
  if(this->context == nullptr) { return 0; }
  if(EVP_DigestUpdate(this->context, s, static_cast<size_t>(n)) == -1) { return 0; }
  return n;
}

std::string
DigestBuf::finish()
{
  if(this->context == nullptr) { return ""; }

  unsigned char md_value[EVP_MAX_MD_SIZE];
  unsigned int md_len = 0;

  if(EVP_DigestFinal_ex(this->context, md_value, &md_len) == -1) { return ""; }
  EVP_MD_CTX_free(this->context);
  this->context = nullptr;

  std::ostringstream ss;
  ss << std::hex;
  for(unsigned int i = 0; i < md_len; i++)
  {
    ss.width(2);
    ss.fill('0');
    ss << static_cast<unsigned int>(md_value[i]);
  }

  return ss.str();
}

DigestStream::DigestStream(const EVP_MD* algorithm) :
  std::ostream(&buffer),
  buffer(algorithm)
{
  if(!this->buffer.good())
  {
    this->setstate(std::ios::badbit);
  }
}

std::string
DigestStream::finish()
{
  std::string digest = this->buffer.finish();
  this->setstate(std::ios::badbit);
  return digest;
}

//------------------------------------------------------------------------------

std::vector<view_type>
parse_relationship(view_type entry, const std::string& type)
{
  std::vector<view_type> names = split_view(entry, GraphName::RELATIONSHIP_SEPARATOR);
  if(names.size() != 2 || names[0].second == 0 || names[1].second == 0)
  {
    std::string msg = "Cannot parse " + type + " relationship: " + entry.to_string();
    throw std::runtime_error(msg);
  }
  return names;
}

GraphName::GraphName(const gbwt::Tags& tags) :
  pggname(tags.get(GBZ_NAME_TAG))
{
  std::string subgraph_tag = tags.get(GBZ_SUBGRAPH_TAG);
  std::vector<view_type> subgraph_entries = split_view(view_type(subgraph_tag), RELATIONSHIP_LIST_SEPARATOR);
  for(const auto& entry_view : subgraph_entries)
  {
    std::vector<view_type> names = parse_relationship(entry_view, "subgraph");
    this->add_subgraph(names[0], names[1]);
  }

  std::string translation_tag = tags.get(GBZ_TRANSLATION_TAG);
  std::vector<view_type> translation_entries = split_view(view_type(translation_tag), RELATIONSHIP_LIST_SEPARATOR);
  for(const auto& entry_view : translation_entries)
  {
    std::vector<view_type> names = parse_relationship(entry_view, "translation");
    this->add_translation(names[0], names[1]);
  }
}

// TODO: This should be in gfa.cpp/.h.
std::tuple<view_type, char, view_type>
parse_typed_field(view_type field)
{
  if(field.second < 5 || field.first[2] != ':' || field.first[4] != ':')
  {
    std::string msg = "Cannot parse typed field: " + field.to_string();
    throw std::runtime_error(msg);
  }
  view_type tag(field.first, 2);
  char type = field.first[3];
  view_type value(field.first + 5, field.second - 5);
  return std::make_tuple(tag, type, value);
}

GraphName::GraphName(const std::vector<std::string>& header_lines)
{
  for(const std::string& line : header_lines)
  {
    std::vector<view_type> fields = split_view(view_type(line), GFA_GAF_FIELD_SEPARATOR);
    if(fields.empty()) { continue; }

    if(fields[0].second == 1 && fields[0].first[0] == GFA_HEADER_PREFIX)
    {
      for(size_t i = 1; i < fields.size(); i++)
      {
        auto parsed = parse_typed_field(fields[i]);
        if(std::get<0>(parsed) == GFA_NAME_TAG && std::get<1>(parsed) == GFA_GAF_TAG_STR_TYPE)
        {
          this->pggname = std::get<2>(parsed).to_string();
        }
        else if(std::get<0>(parsed) == GFA_GAF_SUBGRAPH_TAG && std::get<1>(parsed) == GFA_GAF_TAG_STR_TYPE)
        {
          std::vector<view_type> names = parse_relationship(std::get<2>(parsed), "subgraph");
          this->add_subgraph(names[0], names[1]);
        }
        else if(std::get<0>(parsed) == GFA_GAF_TRANSLATION_TAG && std::get<1>(parsed) == GFA_GAF_TAG_STR_TYPE)
        {
          std::vector<view_type> names = parse_relationship(std::get<2>(parsed), "translation");
          this->add_translation(names[0], names[1]);
        }
      }
    }
    else if(fields[0].second > 1 && fields[0].first[0] == GAF_HEADER_PREFIX)
    {
      fields[0].first++; fields[0].second--; // Skip '@'.
      if(fields[0] == GAF_NAME_TAG)
      {
        if(fields.size() < 2)
        {
          std::string msg = "Cannot parse GAF name line: " + line;
          throw std::runtime_error(msg);
        }
        this->pggname = fields[1].to_string();
      }
      else if(fields[0] == GFA_GAF_SUBGRAPH_TAG)
      {
        if(fields.size() < 3)
        {
          std::string msg = "Cannot parse GAF subgraph line: " + line;
          throw std::runtime_error(msg);
        }
        this->add_subgraph(fields[1], fields[2]);
      }
      else if(fields[0] == GFA_GAF_TRANSLATION_TAG)
      {
        if(fields.size() < 3)
        {
          std::string msg = "Cannot parse GAF translation line: " + line;
          throw std::runtime_error(msg);
        }
        this->add_translation(fields[1], fields[2]);
      }
    }
    else
    {
      std::string msg = "Cannot parse header line: " + line;
      throw std::runtime_error(msg);
    }
  }
}

void
GraphName::add_subgraph(view_type subgraph, view_type supergraph)
{
  if(subgraph == supergraph || subgraph.second == 0 || supergraph.second == 0) { return; }
  this->subgraph[subgraph.to_string()].insert(supergraph.to_string());
}

void
GraphName::add_translation(view_type from, view_type to)
{
  if(from == to || from.second == 0 || to.second == 0) { return; }
  this->translation[from.to_string()].insert(to.to_string());
}

void
GraphName::add_relationships(const GraphName& another)
{
  for(const auto& kv : another.subgraph)
  {
    for(const auto& supergraph_name : kv.second)
    {
      this->add_subgraph(view_type(kv.first), view_type(supergraph_name));
    }
  }

  for(const auto& kv : another.translation)
  {
    for(const auto& to_name : kv.second)
    {
      this->add_translation(view_type(kv.first), view_type(to_name));
    }
  }
}

void
compose_relationship_list(const std::map<std::string, std::set<std::string>>& relationships, gbwt::Tags& tags, const std::string& tag_name)
{
  std::string relationship_tag;
  for(const auto& kv : relationships)
  {
    for(const auto& related_name : kv.second)
    {
      if(!relationship_tag.empty()) { relationship_tag.push_back(GraphName::RELATIONSHIP_LIST_SEPARATOR); }
      relationship_tag.append(kv.first);
      relationship_tag.push_back(GraphName::RELATIONSHIP_SEPARATOR);
      relationship_tag.append(related_name);
    }
  }
  if(relationship_tag.empty())
  {
    tags.unset(tag_name);
  }
  else
  {
    tags.set(tag_name, relationship_tag);
  }
}

void
GraphName::set_tags(gbwt::Tags& tags) const
{
  if(this->pggname.empty())
  {
    tags.unset(GBZ_NAME_TAG);
  }
  else
  {
    tags.set(GBZ_NAME_TAG, this->pggname);
  }
  compose_relationship_list(this->subgraph, tags, GBZ_SUBGRAPH_TAG);
  compose_relationship_list(this->translation, tags, GBZ_TRANSLATION_TAG);
}

void
append_typed_field(std::string& line, const std::string& tag, char type, const std::string& value)
{
  line.push_back(GraphName::GFA_GAF_FIELD_SEPARATOR);
  line.append(tag);
  line.push_back(GraphName::GFA_GAF_TAG_SEPARATOR);
  line.push_back(type);
  line.push_back(GraphName::GFA_GAF_TAG_SEPARATOR);
  line.append(value);
}

std::vector<std::string>
GraphName::gfa_header_lines() const
{
  std::vector<std::string> result;

  if(!this->pggname.empty())
  {
    std::string line;
    line.push_back(GFA_HEADER_PREFIX);
    append_typed_field(line, GFA_NAME_TAG, GFA_GAF_TAG_STR_TYPE, this->pggname);
    result.push_back(line);
  }

  for(const auto& kv : this->subgraph)
  {
    for(const auto& supergraph_name : kv.second)
    {
      std::string line;
      line.push_back(GFA_HEADER_PREFIX);
      append_typed_field(line, GFA_GAF_SUBGRAPH_TAG, GFA_GAF_TAG_STR_TYPE, kv.first + GraphName::RELATIONSHIP_SEPARATOR + supergraph_name);
      result.push_back(line);
    }
  }

  for(const auto& kv : this->translation)
  {
    for(const auto& to_name : kv.second)
    {
      std::string line;
      line.push_back(GFA_HEADER_PREFIX);
      append_typed_field(line, GFA_GAF_TRANSLATION_TAG, GFA_GAF_TAG_STR_TYPE, kv.first + GraphName::RELATIONSHIP_SEPARATOR + to_name);
      result.push_back(line);
    }
  }

  return result;
}

std::vector<std::string>
GraphName::gaf_header_lines() const
{
  std::vector<std::string> result;

  if(!this->pggname.empty())
  {
    std::string line;
    line.push_back(GAF_HEADER_PREFIX);
    line.append(GAF_NAME_TAG);
    line.push_back(GraphName::GFA_GAF_FIELD_SEPARATOR);
    line.append(this->pggname);
    result.push_back(line);
  }

  for(const auto& kv : this->subgraph)
  {
    for(const auto& supergraph_name : kv.second)
    {
      std::string line;
      line.push_back(GAF_HEADER_PREFIX);
      line.append(GFA_GAF_SUBGRAPH_TAG);
      line.push_back(GraphName::GFA_GAF_FIELD_SEPARATOR);
      line.append(kv.first);
      line.push_back(GraphName::GFA_GAF_FIELD_SEPARATOR);
      line.append(supergraph_name);
      result.push_back(line);
    }
  }

  for(const auto& kv : this->translation)
  {
    for(const auto& to_name : kv.second)
    {
      std::string line;
      line.push_back(GAF_HEADER_PREFIX);
      line.append(GFA_GAF_TRANSLATION_TAG);
      line.push_back(GraphName::GFA_GAF_FIELD_SEPARATOR);
      line.append(kv.first);
      line.push_back(GraphName::GFA_GAF_FIELD_SEPARATOR);
      line.append(to_name);
      result.push_back(line);
    }
  }

  return result;
}

bool
GraphName::subgraph_of(const GraphName& another) const
{
  GraphName combined = *this;
  combined.add_relationships(another);
  return !(combined.find_subgraph_path(*this, another).empty());
}

bool
GraphName::translates_to(const GraphName& another) const
{
  GraphName combined = *this;
  combined.add_relationships(another);
  return !(combined.find_path(*this, another).empty());
}

void
append_description(std::string& result, size_t step, const std::string& description)
{
  result.append("Name ");
  result.append(std::to_string(step));
  result.append(" is for ");
  result.append(description);
  result.push_back('\n');
}

void
append_relationship(std::string& result, size_t step, bool is_translation)
{
  result.append("Graph ");
  result.append(std::to_string(step));
  if(is_translation)
  {
    result.append(" translates to ");
  }
  else
  {
    result.append(" is a subgraph of ");
  }
  result.append("graph ");;
  result.append(std::to_string(step + 1));
  result.push_back('\n');
}

void
append_graph(std::string& result, size_t step, const std::string& name)
{
  result.append(std::to_string(step));
  result.push_back('\t');
  result.append(name);
  result.push_back('\n');
}

std::string
GraphName::describe_relationship(const GraphName& another, const std::string& this_desc, const std::string& another_desc) const
{
  GraphName combined = *this;
  combined.add_relationships(another);

  // Try to find a path from this to another, or the other way around.
  auto path = combined.find_path(*this, another);
  std::pair<std::string, std::string> from(this->pggname, this_desc);
  std::pair<std::string, std::string> to(another.pggname, another_desc);
  if(path.empty())
  {
    path = combined.find_path(another, *this);
    if(!path.empty())
    {
      std::swap(from, to);
    }
  }

  // Graph descriptions and relationships.
  std::string result;
  append_description(result, 1, from.second);
  for(size_t i = 0; i + 1 < path.size(); i++)
  {
    append_relationship(result, i + 1, path[i].second);
  }
  if(path.empty())
  {
    append_description(result, 2, to.second);
  }
  else
  {
    append_description(result, path.size(), to.second);
  }

  // Graph names involved in the relationships.
  result.append("With graph names:\n");
  if(path.empty())
  {
    std::string missing_name = "(no name)";
    append_graph(result, 1, (from.first.empty() ? missing_name : from.first));
    append_graph(result, 2, (to.first.empty() ? missing_name : to.first));
  }
  else
  {
    for(size_t i = 0; i < path.size(); i++)
    {
      append_graph(result, i + 1, path[i].first);
    }
  }

  return result;
}

std::vector<std::string>
GraphName::find_subgraph_path(const GraphName& from, const GraphName& to) const
{
  std::vector<std::string> result;
  if(from.pggname.empty() || to.pggname.empty()) { return result; }

  std::unordered_map<std::string, std::string> predecessor;
  predecessor[from.pggname] = "";
  std::deque<std::string> queue;
  queue.push_back(from.pggname);
  while(!queue.empty())
  {
    std::string curr = queue.front(); queue.pop_front();
    if(curr == to.pggname) { break; }

    auto it = this->subgraph.find(curr);
    if(it == this->subgraph.end()) { continue; }
    for(const auto& neighbor : it->second)
    {
      if(predecessor.find(neighbor) == predecessor.end())
      {
        predecessor[neighbor] = curr;
        queue.push_back(neighbor);
      }
    }
  }

  if(predecessor.find(to.pggname) == predecessor.end()) { return result; }
  for(std::string curr = to.pggname; !curr.empty(); curr = predecessor[curr]) { result.push_back(curr); }
  std::reverse(result.begin(), result.end());

  return result;
}

std::vector<std::pair<std::string, bool>>
GraphName::find_path(const GraphName& from, const GraphName& to) const
{
  std::vector<std::pair<std::string, bool>> result;
  if(from.pggname.empty() || to.pggname.empty()) { return result; }

  std::unordered_map<std::string, std::pair<std::string, bool>> predecessor;
  predecessor[from.pggname] = std::make_pair("", false);
  std::deque<std::string> queue;
  queue.push_back(from.pggname);
  while(!queue.empty())
  {
    std::string curr = queue.front(); queue.pop_front();
    if(curr == to.pggname) { break; }

    // Prioritize subgraph edges.
    auto it = this->subgraph.find(curr);
    if(it != this->subgraph.end())
    {
      for(const auto& neighbor : it->second)
      {
        if(predecessor.find(neighbor) == predecessor.end())
        {
          predecessor[neighbor] = std::make_pair(curr, false);
          queue.push_back(neighbor);
        }
      }
    }

    // Then try translation edges.
    it = this->translation.find(curr);
    if(it != this->translation.end())
    {
      for(const auto& neighbor : it->second)
      {
        if(predecessor.find(neighbor) == predecessor.end())
        {
          predecessor[neighbor] = std::make_pair(curr, true);
          queue.push_back(neighbor);
        }
      }
    }
  }

  auto iter = predecessor.find(to.pggname);
  if(iter == predecessor.end()) { return result; }
  result.push_back(std::make_pair(to.pggname, false));
  while(iter->second.first != "")
  {
    result.push_back(iter->second);
    iter = predecessor.find(iter->second.first);
  }
  std::reverse(result.begin(), result.end());

  return result;
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

bool edge_is_canonical(gbwt::node_type from, gbwt::node_type to)
{
  nid_t from_id = gbwt::Node::id(from), to_id = gbwt::Node::id(to);
  if(!gbwt::Node::is_reverse(from)) { return (to_id >= from_id); }
  else
  {
    return ((to_id > from_id) || ((to_id == from_id) && !gbwt::Node::is_reverse(to)));
  }
}

bool path_is_canonical(const gbwt::vector_type& path)
{
  if(path.empty()) { return true; }
  if(gbwt::Node::is_reverse(path.front()) == gbwt::Node::is_reverse(path.back())) { return !gbwt::Node::is_reverse(path.front()); }
  return edge_is_canonical(path.front(), path.back());
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
    inverse.emplace_back(iter->second, view_type(iter->first));
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
    else { return view_type(empty); }
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
  // For error handling we need to be able to describe the path
  // This returns a string starting with " ".
  auto describe_path = [&]() {
    std::stringstream ss;
    if(sample_name != PathMetadata::NO_SAMPLE_NAME) { ss << " sample " << sample_name; }
    if(locus_name != PathMetadata::NO_LOCUS_NAME) { ss << " contig " << locus_name; }
    if(haplotype != PathMetadata::NO_HAPLOTYPE) { ss << " phase " << haplotype; }
    if(phase_block != PathMetadata::NO_PHASE_BLOCK) { ss << " count " << phase_block; }
    if(subrange != PathMetadata::NO_SUBRANGE)
    {
      ss << " range [" << subrange.first;
      if(subrange.second != PathMetadata::NO_END_POSITION) { ss << "-" << subrange.second; }
      ss << "]";
    }
    return ss.str();
  };

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

  if(subrange != PathMetadata::NO_SUBRANGE)
  {
    if(phase_block != PathMetadata::NO_PHASE_BLOCK && phase_block != 0 && subrange.first != 0)
    {
      // TODO: We can't represent both a nonzero phase block and a nonzero subrange.
      throw std::runtime_error("MetadataBuilder: Can't represent a nonzero phase block and a nonzero subrange in the same path" + describe_path());
    }
    // Use the subrange start as the count
    path_name.count = subrange.first;
  }
  else if(phase_block != PathMetadata::NO_PHASE_BLOCK)
  {
    // Use the phase block as the count
    path_name.count = phase_block;
  }

  {
    auto iter = this->counts.find(path_name);
    if(iter != this->counts.end())
    {
      // This path is duplicate or ambiguous.

      if(phase_block != PathMetadata::NO_PHASE_BLOCK || subrange != PathMetadata::NO_SUBRANGE)
      {
        // The count was user-specified, so bail out.
        throw std::runtime_error("MetadataBuilder: Duplicate path for" + describe_path());
      }
      else
      {
        // We are going to infer a count; just warn.
        std::cerr << "MetadataBuilder::add_path(): Warning: Path for" << describe_path() << " is not unique" << std::endl;
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
    if (start == "*")
    {
      phase_block = 0;
    }
    else
    {
      try { phase_block = std::stoul(start); }
      catch(const std::invalid_argument&)
      {
        throw std::runtime_error("MetadataBuilder: Invalid start position " + start);
      }
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

  size_t path_count = 0;
  for(auto& names_by_job : this->path_names)
  {
    for(auto& path_name : names_by_job)
    {
      metadata.addPath(path_name);
      path_count++;
    }
  }

  if(this->sample_names.empty() && path_count > 0) { metadata.setSamples(1); }
  else { metadata.setSamples(map_to_vector(this->sample_names)); }
  metadata.setHaplotypes(this->haplotypes.size());
  if(this->contig_names.empty() && path_count > 0) { metadata.setContigs(1); }
  else { metadata.setContigs(map_to_vector(this->contig_names)); }

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

