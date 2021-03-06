#include <gbwtgraph/utils.h>

#include <algorithm>
#include <sstream>

#include <gbwt/utils.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Global variables.

const std::string REFERENCE_PATH_SAMPLE_NAME = "_gbwt_ref";

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t StringArray::BLOCK_SIZE;

//------------------------------------------------------------------------------

// Other class variables.

const std::string SequenceSource::TRANSLATION_EXTENSION = ".trans";

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

void
SequenceSource::translate_segment(const std::string& name, view_type sequence, size_t max_length)
{
  if(this->segment_translation.find(name) != this->segment_translation.end() || sequence.second == 0) { return; }

  nid_t start = this->next_id;
  nid_t limit = start + (sequence.second + max_length - 1) / max_length;
  for(nid_t id = start; id < limit; id++)
  {
    size_t offset = (id - start) * max_length;
    size_t length = std::min(max_length, sequence.second - offset);
    this->add_node(id, view_type(sequence.first + offset, length));
  }

  this->segment_translation[name] = std::make_pair(start, limit);
  this->next_id = limit;
}

std::pair<StringArray, sdsl::sd_vector<>>
SequenceSource::invert_translation() const
{
  std::pair<StringArray, sdsl::sd_vector<>> result;

  // Invert the translation.
  std::vector<std::pair<std::pair<nid_t, nid_t>, view_type>> inverse;
  inverse.reserve(this->segment_translation.size());
  for(auto iter = this->segment_translation.begin(); iter != this->segment_translation.end(); ++iter)
  {
    inverse.emplace_back(iter->second, str_to_view(iter->first));
  }
  gbwt::parallelQuickSort(inverse.begin(), inverse.end());

  // Store the segment names.
  result.first = StringArray(inverse.size(),
  [&](size_t offset) -> size_t
  {
    return inverse[offset].second.second;
  },
  [&](size_t offset) -> view_type
  {
    return inverse[offset].second;
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

StringArray::StringArray(const std::vector<std::string>& source)
{
  *this = StringArray(source.size(),
  [&source](size_t i) -> size_t
  {
    return source[i].length();
  },
  [&source](size_t i) -> view_type
  {
    return str_to_view(source[i]);
  });
}

StringArray::StringArray(size_t n, const std::function<size_t(size_t)>& length, const std::function<view_type(size_t)>& sequence)
{
  size_t total_length = 0;
  for(size_t i = 0; i < n; i++) { total_length += length(i); }
  this->sequences.reserve(total_length);
  this->offsets = sdsl::int_vector<0>(n + 1, 0, sdsl::bits::length(total_length));

  size_t total = 0;
  for(size_t i = 0; i < n; i++)
  {
    view_type view = sequence(i);
    this->sequences.insert(this->sequences.end(), view.first, view.first + view.second);
    this->offsets[i] = total;
    total += view.second;
  }
  this->offsets[n] = total;
}

StringArray::StringArray(size_t n, const std::function<size_t(size_t)>& length, const std::function<std::string(size_t)>& sequence)
{
  size_t total_length = 0;
  for(size_t i = 0; i < n; i++) { total_length += length(i); }
  this->sequences.reserve(total_length);
  this->offsets = sdsl::int_vector<0>(n + 1, 0, sdsl::bits::length(total_length));

  size_t total = 0;
  for(size_t i = 0; i < n; i++)
  {
    std::string str = sequence(i);
    this->sequences.insert(this->sequences.end(), str.begin(), str.end());
    this->offsets[i] = total;
    total += str.length();
  }
  this->offsets[n] = total;
}

void
StringArray::swap(StringArray& another)
{
  this->sequences.swap(another.sequences);
  this->offsets.swap(another.offsets);
}

void
StringArray::serialize(std::ostream& out) const
{
  size_t sequence_size = this->sequences.size();
  sdsl::write_member(sequence_size, out);
  for(size_t offset = 0; offset < this->sequences.size(); offset += BLOCK_SIZE)
  {
    size_t bytes = std::min(BLOCK_SIZE, this->sequences.size() - offset);
    out.write(this->sequences.data() + offset, bytes);
  }
  this->offsets.serialize(out);
}

void
StringArray::deserialize(std::istream& in)
{
  size_t sequence_size = 0;
  sdsl::read_member(sequence_size, in);
  this->sequences = std::vector<char>(sequence_size, 0);
  for(size_t offset = 0; offset < this->sequences.size(); offset += BLOCK_SIZE)
  {
    size_t bytes = std::min(BLOCK_SIZE, this->sequences.size() - offset);
    in.read(this->sequences.data() + offset, bytes);
  }
  this->offsets.load(in);
}

void
StringArray::compress(std::ostream& out) const
{
  // Determine and serialize the alphabet.
  std::vector<std::uint8_t> char_to_comp(256, 0);
  for(char c : this->sequences) { char_to_comp[static_cast<std::uint8_t>(c)] = 1; }
  size_t sigma = 0;
  for(auto c : char_to_comp) { if(c != 0) { sigma++; } }
  sigma = std::max(sigma, size_t(1));
  sdsl::int_vector<8> comp_to_char(sigma, 0);
  for(size_t i = 0, found = 0; i < char_to_comp.size(); i++)
  {
    if(char_to_comp[i] != 0)
    {
      comp_to_char[found] = i;
      char_to_comp[i] = found;
      found++;
    }
  }
  comp_to_char.serialize(out);

  // Compress the sequences.
  {
    sdsl::int_vector<> compressed(this->sequences.size(), 0, sdsl::bits::length(sigma - 1));
    for(size_t i = 0; i < this->sequences.size(); i++)
    {
      compressed[i] = char_to_comp[static_cast<uint8_t>(this->sequences[i])];
    }
    compressed.serialize(out);
  }

  // Compress the offsets.
  {
    sdsl::sd_vector<> v(this->offsets.begin(), this->offsets.end());
    v.serialize(out);
  }
}

void
StringArray::decompress(std::istream& in)
{
  // Load the alphabet.
  sdsl::int_vector<8> comp_to_char;
  comp_to_char.load(in);

  // Decompress the sequences.
  {
    sdsl::int_vector<> compressed; compressed.load(in);
    this->sequences = std::vector<char>(); this->sequences.reserve(compressed.size());
    for(auto c : compressed) { this->sequences.push_back(comp_to_char[c]); }
  }

  // Decompress the offsets.
  {
    sdsl::sd_vector<> v; v.load(in);
    this->offsets = sdsl::int_vector<>(v.ones(), 0);
    size_t i = 0;
    for(auto iter = v.one_begin(); iter != v.one_end(); ++iter, i++) { this->offsets[i] = iter->second; }
    sdsl::util::bit_compress(this->offsets);
  }
}

bool StringArray::operator==(const StringArray& another) const
{
  return (this->sequences == another.sequences && this->offsets == another.offsets);
}

bool StringArray::operator!=(const StringArray& another) const
{
  return !(this->operator==(another));
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph

