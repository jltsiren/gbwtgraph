#include <gbwtgraph/gfa.h>
#include <gbwtgraph/internal.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <regex>
#include <string>
#include <unordered_set>
#include <utility>

#include <cctype>

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Global variables.

const std::string GFA_EXTENSION = ".gfa";

// Class constants.
const std::string GFAParsingParameters::DEFAULT_REGEX = ".*";
const std::string GFAParsingParameters::DEFAULT_FIELDS = "S";

//------------------------------------------------------------------------------

// Parse a nonnegative integer, assuming that the string is valid.
std::uint64_t
stoul_unsafe(const std::string& str)
{
  std::uint64_t result = 0;
  for(char c : str) { result = 10 * result + (c - '0'); }
  return result;
}

//------------------------------------------------------------------------------

struct GFAFile
{
  // Memory mapped file.
  int    fd;
  size_t file_size;
  char*  ptr;

  // GFA information.
  bool valid_gfa;
  bool translate_segment_ids;
  size_t max_segment_length, max_path_length;

  // Bit masks for field separators.
  size_t field_end[4];
  size_t subfield_end[4];
  size_t walk_subfield_end[4];

  // Pointers to line starts.
  std::vector<const char*> s_lines;
  std::vector<const char*> p_lines;
  std::vector<const char*> w_lines;

  struct field_type
  {
    const char* begin;
    const char* end;
    size_t      line_num;
    char        type;
    bool        has_next;

    size_t size() const { return this->end - this->begin; }
    bool empty() const { return (this->size() == 0); }

    std::string str() const { return std::string(this->begin, this->end); }
    view_type view() const { return view_type(this->begin, this->end - this->begin); }
    char front() const { return *(this->begin); }
    char back() const { return *(this->end - 1); }

    // For path segment subfields.
    bool valid_path_segment() const { return (this->size() >= 2 && (this->back() == '-' || this->back() == '+')); }
    std::string path_segment() const { return std::string(this->begin, this->end - 1); }
    view_type path_segment_view() const { return view_type(this->begin, (this->end - 1) - this->begin); }
    bool is_reverse_path_segment() const { return (this->back() == '-'); }

    // Usually the next field/subfield starts at `end + 1`, because `end` points
    // to the separator. Walk subfields include the separator as a part of the field,
    // so they start at `end` instead. Before we go to the first subfield, we must
    // increment `end` (which points to the preceding '\t' before the call).
    void start_walk() { this->end++; }

    // For walk segment subfields.
    bool valid_walk_segment() const { return (this->size() >= 2 && (this->front() == '<' || this->front() == '>')); }
    std::string walk_segment() const { return std::string(this->begin + 1, this->end); }
    view_type walk_segment_view() const { return view_type(this->begin + 1, this->end - (this->begin + 1)); }
    bool is_reverse_walk_segment() const { return (this->front() == '<'); }
  };

  // Memory map and validate a GFA file. The constructor checks the following:
  //
  // * The mandatory fields used for GBWTGraph construction exist and are nonempty.
  // * There are no duplicate segment names.
  // * All segments used on paths/walks exist.
  //
  // There is no check for duplicate path/walk names, because that check is cheaper
  // using the parsed metadata.
  GFAFile(const std::string& filename, bool show_progress);

  ~GFAFile();

  bool ok() const { return (this->fd >= 0 && this->file_size > 0 && this->ptr != nullptr && this->valid_gfa); }
  size_t size() const { return this->file_size; }
  size_t segments() const { return this->s_lines.size(); }
  size_t paths() const { return this->p_lines.size(); }
  size_t walks() const { return this->w_lines.size(); }

private:
  // Preprocess a new S-line. Returns an iterator at the start of the next line or
  // nullptr if the parse failed.
  const char* add_s_line(const char* iter, std::unordered_set<std::string>& found_segments, size_t line_num);

  // Preprocess a new P-line. Returns an iterator at the start of the next line or
  // nullptr if the parse failed.
  const char* add_p_line(const char* iter, BufferedHashSet& required_segments, size_t line_num);

  // Preprocess a new W-line. Returns an iterator at the start of the next line or
  // nullptr if the parse failed.
  const char* add_w_line(const char* iter, BufferedHashSet& required_segments, size_t line_num);

  // Returns true if the field is valid. Otherwise marks the GFA file invalid,
  // prints an error message, and returns false.
  bool check_field(const field_type& field, const std::string& field_name, bool should_have_next);

  const char* begin() const { return this->ptr; }
  const char* end() const { return this->ptr + this->size(); }

  // Return an iterator to the beginning of the next line.
  const char* next_line(const char* iter) const
  {
    while(iter != this->end() && *iter != '\n') { ++iter; }
    if(iter != this->end()) { ++iter; }
    return iter;
  }

  // Return the first tab-separated field of the line.
  field_type first_field(const char* line_start, size_t line_num = 0) const
  {
    const char* limit = line_start;
    while(limit != this->end() && !(this->is_field_end(limit))) { ++limit; }
    return { line_start, limit, line_num, *line_start, (limit != this->end() && *limit == '\t') };
  }

  // Return the next tab-separated field, assuming there is one.
  field_type next_field(const field_type& field) const
  {
    const char* limit = field.end + 1;
    while(limit != this->end() && !(this->is_field_end(limit))) { ++limit; }
    return { field.end + 1, limit, field.line_num, field.type, (limit != this->end() && *limit == '\t') };
  }

  // Return the next comma-separated subfield, assuming there is one.
  field_type next_subfield(const field_type& field) const
  {
    const char* limit = field.end + 1;
    while(limit != this->end() && !(this->is_subfield_end(limit))) { ++limit; }
    return { field.end + 1, limit, field.line_num, field.type, (limit != this->end() && *limit == ',') };
  }

  // Return the next walk subfield, assuming there is one.
  // The orientation symbol at the start of the segment is also used as subfield separator.
  field_type next_walk_subfield(const field_type& field) const
  {
    const char* limit = field.end;
    if(limit != this->end() && (*limit == '<' || *limit == '>'))
    {
      do { limit++; }
      while(limit != this->end() && !(this->is_walk_subfield_end(limit)));
    }
    return { field.end, limit, field.line_num, field.type, (limit != this->end() && (*limit == '<' || *limit == '>')) };
  }

  bool is_field_end(const char* iter) const
  {
    unsigned char c = *iter;
    return (this->field_end[c / 64] & (size_t(1) << (c & 0x3F)));
  }

  bool is_subfield_end(const char* iter) const
  {
    unsigned char c = *iter;
    return (this->subfield_end[c / 64] & (size_t(1) << (c & 0x3F)));
  }

  bool is_walk_subfield_end(const char* iter) const
  {
    unsigned char c = *iter;
    return (this->walk_subfield_end[c / 64] & (size_t(1) << (c & 0x3F)));
  }

  void add_field_end(unsigned char c)
  {
    this->field_end[c / 64] |= size_t(1) << (c & 0x3F);
  }

  void add_subfield_end(unsigned char c)
  {
    this->subfield_end[c / 64] |= size_t(1) << (c & 0x3F);
  }

  void add_walk_subfield_end(unsigned char c)
  {
    this->walk_subfield_end[c / 64] |= size_t(1) << (c & 0x3F);
  }

public:
  /*
    Iterate over the S-lines, calling segment() for all segments. Stops early if segment()
    returns false.
  */
  void for_each_segment(const std::function<bool(const std::string& name, view_type sequence)>& segment) const;

  /*
    Iterate over the file, calling path() for each path. Stops early if path() returns false.
  */
  void for_each_path_name(const std::function<bool(const std::string& name)>& path) const;

  /*
    Iterate over the file, calling path() for each path, path_segment() for
    each path segment, and finish_path() after parsing each path. Stops early
    if any call returns false.
  */
  void for_each_path(const std::function<bool(const std::string& name)>& path,
                     const std::function<bool(const std::string& name, bool is_reverse)>& path_segment,
                     const std::function<bool()>& finish_path) const;

  /*
    Iterate over the file, calling walk() for each walk. Stops early if walk() returns false.
  */
  void for_each_walk_name(const std::function<bool(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)>& walk) const;

  /*
    Iterate over the file, calling walk() for each walk, walk_segment() for
    each walk segment, and finish_walk() after parsing each walk. Stops early
    if any call returns false.
  */
  void for_each_walk(const std::function<bool(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)>& walk,
                     const std::function<bool(const std::string& name, bool is_reverse)>& walk_segment,
                     const std::function<bool()>& finish_walk) const;
};

//------------------------------------------------------------------------------

GFAFile::GFAFile(const std::string& filename, bool show_progress) :
  fd(-1), file_size(0), ptr(nullptr),
  valid_gfa(true), translate_segment_ids(false),
  max_segment_length(0), max_path_length(0)
{
  if(show_progress)
  {
    std::cerr << "Opening GFA file " << filename << std::endl;
  }

  // Open the file.
  this->fd = ::open(filename.c_str(), O_RDONLY);
  if(this->fd < 0)
  {
    std::cerr << "GFAFile::GFAFile(): Cannot open file " << filename << std::endl;
    return;
  }

  // Memory map the file.
  struct stat st;
  if(::fstat(this->fd, &st) < 0)
  {
    std::cerr << "GFAFile::GFAFile(): Cannot stat file " << filename << std::endl;
    return;
  }
  this->file_size = st.st_size;

  void* temp_ptr = ::mmap(nullptr, file_size, PROT_READ, MAP_FILE | MAP_SHARED, this->fd, 0);
  if(temp_ptr == MAP_FAILED)
  {
    std::cerr << "GFAFile::GFAFile(): Cannot memory map file " << filename << std::endl;
    return;
  }
  ::madvise(temp_ptr, file_size, MADV_SEQUENTIAL); // We will be making sequential passes over the data.
  this->ptr = static_cast<char*>(temp_ptr);

  // Mark characters indicating field/subfield end. This could depend on the GFA version.
  // TODO: If that happens, we need variables for field/subfield separators.
  this->field_end[0] = 0; this->field_end[1] = 0;
  this->field_end[2] = 0; this->field_end[3] = 0;
  this->add_field_end('\n'); this->add_field_end('\t');
  this->subfield_end[0] = 0; this->subfield_end[1] = 0;
  this->subfield_end[2] = 0; this->subfield_end[3] = 0;
  this->add_subfield_end('\n'); this->add_subfield_end('\t'); this->add_subfield_end(',');
  this->walk_subfield_end[0] = 0; this->walk_subfield_end[1] = 0;
  this->walk_subfield_end[2] = 0; this->walk_subfield_end[3] = 0;
  this->add_walk_subfield_end('\n'); this->add_walk_subfield_end('\t');
  this->add_walk_subfield_end('<'); this->add_walk_subfield_end('>');

  // Preprocess and validate the file.
  double start = gbwt::readTimer();
  if(show_progress)
  {
    std::cerr << "Validating GFA file " << filename << std::endl;
  }
  const char* iter = this->begin();
  size_t line_num = 0;
  std::unordered_set<std::string> found_segments;
  BufferedHashSet required_segments;
  while(iter != this->end())
  {
    switch(*iter)
    {
    case 'S':
      iter = this->add_s_line(iter, found_segments, line_num);
      break;
    case 'P':
      iter = this->add_p_line(iter, required_segments, line_num);
      break;
    case 'W':
      iter = this->add_w_line(iter, required_segments, line_num);
      break;
    default:
      iter = this->next_line(iter);
      break;
    }
    if(iter == nullptr) { return; }
    line_num++;
  }

  // Flush the buffers and check that we have found all the required segments.
  required_segments.finish();
  for(const std::string& segment_name : required_segments.data)
  {
    if(found_segments.find(segment_name) == found_segments.end())
    {
      std::cerr << "GFAFile::GFAFile(): Segment " << segment_name << " used on paths/walks but missing from the file" << std::endl;
      this->valid_gfa = false;
      return;
    }
  }

  if(show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Found " << this->segments() << " segments, " << this->paths() << " paths, and " << this->walks() << " walks in " << seconds << " seconds" << std::endl;
  }
}

const char*
GFAFile::add_s_line(const char* iter, std::unordered_set<std::string>& found_segments, size_t line_num)
{
  this->s_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->first_field(iter, line_num);
  if(!(this->check_field(field, "record type", true))) { return nullptr; }

  // Segment name field.
  field = this->next_field(field);
  if(!(this->check_field(field, "segment name", true))) { return nullptr; }
  std::string name = field.str();
  if(found_segments.find(name) != found_segments.end())
  {
    std::cerr << "GFAFile::add_s_line(): Duplicate segment name " << name << " on line " << line_num << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  found_segments.insert(name);
  if(!(this->translate_segment_ids))
  {
    try
    {
      nid_t id = std::stoul(name);
      if (id == 0) { this->translate_segment_ids = true; }
    }
    catch(const std::invalid_argument&) { this->translate_segment_ids = true; }
  }

  // Sequence field.
  field = this->next_field(field);
  if(!(this->check_field(field, "sequence", false))) { return nullptr; }
  this->max_segment_length = std::max(this->max_segment_length, field.size());

  return this->next_line(field.end);
}

const char*
GFAFile::add_p_line(const char* iter, BufferedHashSet& required_segments, size_t line_num)
{
  this->p_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->first_field(iter, line_num);
  if(!(this->check_field(field, "record type", true))) { return nullptr; }

  // Path name field.
  field = this->next_field(field);
  if(!(this->check_field(field, "path_name", true))) { return nullptr; }

  // Segment names field.
  size_t path_length = 0;
  do
  {
    field = this->next_subfield(field);
    if(!(field.valid_path_segment()))
    {
      std::cerr << "GFAFile::add_p_line(): Invalid path segment " << field.str() << " on line " << line_num << std::endl;
      this->valid_gfa = false;
      return nullptr;
    }
    required_segments.insert(field.path_segment_view());
    path_length++;
  }
  while(field.has_next);
  if(path_length == 0)
  {
    std::cerr << "GFAFile::add_p_line(): The path on line " << line_num << " is empty" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  this->max_path_length = std::max(this->max_path_length, path_length);

  return this->next_line(field.end);
}

const char*
GFAFile::add_w_line(const char* iter, BufferedHashSet& required_segments, size_t line_num)
{
  this->w_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->first_field(iter, line_num);
  if(!(this->check_field(field, "record type", true))) { return nullptr; }

  // Sample name field.
  field = this->next_field(field);
  if(!(this->check_field(field, "sample name", true))) { return nullptr; }

  // Haplotype index field.
  field = this->next_field(field);
  if(!(this->check_field(field, "haplotype index", true))) { return nullptr; }

  // Contig name field.
  field = this->next_field(field);
  if(!(this->check_field(field, "contig name", true))) { return nullptr; }

  // Start position field.
  field = this->next_field(field);
  if(!(this->check_field(field, "start position", true))) { return nullptr; }

  // Skip the end position field.
  field = this->next_field(field);
  if(!(this->check_field(field, "end position", true))) { return nullptr; }

  // Segment names field.
  size_t path_length = 0;
  field.start_walk();
  do
  {
    field = this->next_walk_subfield(field);
    if(!(field.valid_walk_segment()))
    {
      std::cerr << "GFAFile::add_w_line(): Invalid walk segment " << field.str() << " on line " << line_num << std::endl;
      this->valid_gfa = false;
      return nullptr;
    }
    required_segments.insert(field.walk_segment_view());
    path_length++;
  }
  while(field.has_next);
  if(path_length == 0)
  {
    std::cerr << "GFAFile::add_w_line(): The walk on line " << line_num << " is empty" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  this->max_path_length = std::max(this->max_path_length, path_length);

  return this->next_line(field.end);
}

bool
GFAFile::check_field(const field_type& field, const std::string& field_name, bool should_have_next)
{
  if(field.empty())
  {
    std::cerr << "GFAFile::check_field(): " << field.type << "-line " << field.line_num << " has no " << field_name << std::endl;
    this->valid_gfa = false;
    return false;
  }
  if(should_have_next && !(field.has_next))
  {
    std::cerr << "GFAFile::check_field(): " << field.type << "-line " << field.line_num << " ended after " << field_name << std::endl;
    this->valid_gfa = false;
    return false;
  }
  return true;
}

GFAFile::~GFAFile()
{
  if(this->ptr != nullptr)
  {
    ::munmap(static_cast<void*>(this->ptr), this->file_size);
    this->file_size = 0;
    this->ptr = nullptr;
  }
  if(this->fd >= 0)
  {
    ::close(this->fd);
    this->fd = -1;
  }
}

//------------------------------------------------------------------------------

void
GFAFile::for_each_segment(const std::function<bool(const std::string& name, view_type sequence)>& segment) const
{
  if(!(this->ok())) { return; }

  for(const char* iter : this->s_lines)
  {
    // Skip the record type field.
    field_type field = this->first_field(iter);

    // Segment name field.
    field = this->next_field(field);
    std::string name = field.str();

    // Sequence field.
    field = this->next_field(field);
    view_type sequence = field.view();
    if(!segment(name, sequence)) { return; }
  }
}

void
GFAFile::for_each_path_name(const std::function<bool(const std::string& name)>& path) const
{
  if(!(this->ok())) { return; }

  for(const char* iter : this->p_lines)
  {
    // Skip the record type field.
    field_type field = this->first_field(iter);

    // Path name field.
    field = this->next_field(field);
    std::string path_name = field.str();
    if(!path(path_name)) { return; }
  }
}

void
GFAFile::for_each_path(const std::function<bool(const std::string& name)>& path,
                       const std::function<bool(const std::string& name, bool is_reverse)>& path_segment,
                       const std::function<bool()>& finish_path) const
{
  if(!(this->ok())) { return; }

  for(const char* iter : this->p_lines)
  {
    // Skip the record type field.
    field_type field = this->first_field(iter);

    // Path name field.
    field = this->next_field(field);
    if(!path(field.str())) { return; }

    // Segment names field.
    do
    {
      field = this->next_subfield(field);
      std::string segment_name = field.path_segment();
      if(!path_segment(segment_name, field.is_reverse_path_segment())) { return; }
    }
    while(field.has_next);

    if(!finish_path()) { return; }
  }
}

void
GFAFile::for_each_walk_name(const std::function<bool(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)>& walk) const
{
  if(!(this->ok())) { return; }

  for(const char* iter : this->w_lines)
  {
    // Skip the record type field.
    field_type field = this->first_field(iter);

    // Sample field.
    field = this->next_field(field);
    std::string sample = field.str();

    // Haplotype field.
    field = this->next_field(field);
    std::string haplotype = field.str();

    // Contig field.
    field = this->next_field(field);
    std::string contig = field.str();

    // Start field.
    field = this->next_field(field);
    std::string start = field.str();

    if(!walk(sample, haplotype, contig, start)) { return; }
  }
}

void
GFAFile::for_each_walk(const std::function<bool(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)>& walk,
                       const std::function<bool(const std::string& name, bool is_reverse)>& walk_segment,
                       const std::function<bool()>& finish_walk) const
{
  if(!(this->ok())) { return; }

  for(const char* iter : this->w_lines)
  {
    // Skip the record type field.
    field_type field = this->first_field(iter);

    // Sample field.
    field = this->next_field(field);
    std::string sample = field.str();

    // Haplotype field.
    field = this->next_field(field);
    std::string haplotype = field.str();

    // Contig field.
    field = this->next_field(field);
    std::string contig = field.str();

    // Start field.
    field = this->next_field(field);
    std::string start = field.str();

    if(!walk(sample, haplotype, contig, start)) { return; }

    // Skip the end field.
    field = this->next_field(field);

    // Segment names field.
    field.start_walk();
    do
    {
      field = this->next_walk_subfield(field);
      std::string segment_name = field.walk_segment();
      if(!walk_segment(segment_name, field.is_reverse_walk_segment())) { return; }
    }
    while(field.has_next);

    if(!finish_walk()) { return; }
  }
}

//------------------------------------------------------------------------------

struct PathNameParser
{
  std::regex parser;

  // Mapping from regex submatches to GBWT path name components.
  constexpr static size_t NO_FIELD = std::numeric_limits<size_t>::max();
  size_t sample_field, contig_field, haplotype_field, fragment_field;

  // GBWT metadata.
  std::map<std::string, size_t> sample_names, contig_names;
  std::set<std::pair<size_t, size_t>> haplotypes; // (sample id, phase id)
  std::vector<gbwt::PathName> path_names;
  std::map<gbwt::PathName, size_t> counts;

  bool ref_path_sample_warning;
  bool ok;

  explicit PathNameParser(const GFAParsingParameters& parameters);

  // Parse a path name using a regex. Returns true if successful.
  // This should not be used with add_walk() or add_reference_path().
  bool parse(const std::string& name);

  // Add a path based on walk metadata. Returns true if successful.
  // This should not be used with parse().
  bool add_walk(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start);

  // Add a reference path. Returns true if successful.
  // This should not be used with parse().
  bool add_reference_path(const std::string& name);

  bool empty() const { return this->path_names.empty(); }

  void set_metadata(gbwt::Metadata& metadata)
  {
    if(this->sample_names.empty()) { metadata.setSamples(1); }
    else { metadata.setSamples(map_to_vector(this->sample_names)); }
    metadata.setHaplotypes(this->haplotypes.size());
    if(this->contig_names.empty()) { metadata.setContigs(1); }
    else { metadata.setContigs(map_to_vector(this->contig_names)); }
    for(auto& path_name : this->path_names) { metadata.addPath(path_name); }
  }

  void clear()
  {
    this->sample_names = std::map<std::string, size_t>();
    this->contig_names = std::map<std::string, size_t>();
    this->haplotypes = std::set<std::pair<size_t, size_t>>();
    this->path_names = std::vector<gbwt::PathName>();
    this->counts = std::map<gbwt::PathName, size_t>();
  }

  static std::vector<std::string> map_to_vector(const std::map<std::string, size_t>& source)
  {
    std::vector<std::string> result(source.size());
    for(auto& name : source) { result[name.second] = name.first; }
    return result;
  }
};

//------------------------------------------------------------------------------

PathNameParser::PathNameParser(const GFAParsingParameters& parameters) :
  sample_field(NO_FIELD), contig_field(NO_FIELD), haplotype_field(NO_FIELD), fragment_field(NO_FIELD),
  ref_path_sample_warning(false), ok(true)
{
  // Initialize the regex.
  try { this->parser = std::regex(parameters.path_name_regex); }
  catch(std::regex_error& e)
  {
    std::cerr << "PathNameParser::PathNameParser(): Invalid regex " << parameters.path_name_regex << std::endl;
    std::cerr << "PathNameParser::PathNameParser(): Error was: " << e.what() << std::endl;
    this->ok = false; return;
  }
  if(parameters.path_name_fields.size() > this->parser.mark_count() + 1)
  {
    std::cerr << "PathNameParser::PathNameParser(): Field string too long: " << parameters.path_name_fields << std::endl;
    this->ok = false; return;
  }

  // Initialize the fields.
  for(size_t i = 0; i < parameters.path_name_fields.size(); i++)
  {
    switch(std::tolower(parameters.path_name_fields[i]))
    {
      case 's':
        if(this->sample_field != NO_FIELD)
        {
          std::cerr << "PathNameParser::PathNameParser(): Duplicate sample field" << std::endl;
          this->ok = false; return;
        }
        this->sample_field = i;
        break;
      case 'c':
        if(this->contig_field != NO_FIELD)
        {
          std::cerr << "PathNameParser::PathNameParser(): Duplicate contig field" << std::endl;
          this->ok = false; return;
        }
        this->contig_field = i;
        break;
      case 'h':
        if(this->haplotype_field != NO_FIELD)
        {
          std::cerr << "PathNameParser::PathNameParser(): Duplicate haplotype field" << std::endl;
          this->ok = false; return;
        }
        this->haplotype_field = i;
        break;
      case 'f':
        if(this->fragment_field != NO_FIELD)
        {
          std::cerr << "PathNameParser::PathNameParser(): Duplicate fragment field" << std::endl;
          this->ok = false; return;
        }
        this->fragment_field = i;
        break;
    }
  }
}

bool
PathNameParser::parse(const std::string& name)
{
  std::smatch fields;
  if(!std::regex_match(name, fields, this->parser))
  {
    std::cerr << "PathNameParser::parse(): Invalid path name " << name << std::endl;
    return false;
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
      std::cerr << "PathNameParser::parse(): Warning: Sample name " << REFERENCE_PATH_SAMPLE_NAME << " is reserved for reference paths" << std::endl;
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
      std::cerr << "PathNameParser::parse(): Invalid haplotype field " << fields[this->haplotype_field] << std::endl;
      return false;
    }
  }
  this->haplotypes.insert(std::pair<size_t, size_t>(path_name.sample, path_name.phase));

  if(this->fragment_field != NO_FIELD)
  {
    try { path_name.count = std::stoul(fields[this->fragment_field]); }
    catch(const std::invalid_argument&)
    {
      std::cerr << "PathNameParser::parse(): Invalid fragment field " << fields[this->fragment_field] << std::endl;
      return false;
    }
    if(this->counts.find(path_name) != this->counts.end())
    {
      std::cerr << "PathNameParser::parse(): Duplicate path name " << name << std::endl;
      return false;
    }
    this->counts[path_name] = 1;
  }
  else
  {
    auto iter = this->counts.find(path_name);
    if(iter == this->counts.end()) { this->counts[path_name] = 1; }
    else { path_name.count = iter->second; iter->second++; }
  }

  this->path_names.push_back(path_name);
  return true;
}

bool
PathNameParser::add_walk(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)
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
      std::cerr << "PathNameParser::add_walk(): Warning: Sample name " << REFERENCE_PATH_SAMPLE_NAME << " is reserved for reference paths" << std::endl;
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
      std::cerr << "PathNameParser::add_walk(): Invalid haplotype field " << haplotype << std::endl;
      return false;
    }
  }
  this->haplotypes.insert(std::pair<size_t, size_t>(path_name.sample, path_name.phase));

  // Start position as fragment identifier.
  {
    try { path_name.count = std::stoul(start); }
    catch(const std::invalid_argument&)
    {
      std::cerr << "PathNameParser::add_walk(): Invalid start_position " << start << std::endl;
      return false;
    }
    if(this->counts.find(path_name) != this->counts.end())
    {
      std::cerr << "PathNameParser::add_walk(): Duplicate walk " << sample << "\t" << haplotype << "\t" << contig << "\t" << start << ")" << std::endl;
      return false;
    }
    this->counts[path_name] = 1;
  }

  this->path_names.push_back(path_name);
  return true;
}

bool
PathNameParser::add_reference_path(const std::string& name)
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
    std::cerr << "PathNameParser::add_reference_path(): Duplicate path " << name << std::endl;
    return false;
  }
  this->counts[path_name] = 1;

  this->path_names.push_back(path_name);
  return true;
}

//------------------------------------------------------------------------------

bool
check_gfa_file(const GFAFile& gfa_file, const GFAParsingParameters& parameters)
{
  if(!(gfa_file.ok())) { return false; }

  if(gfa_file.segments() == 0)
  {
    std::cerr << "check_gfa_file(): No segments in the GFA file" << std::endl;
    return false;
  }
  if(gfa_file.paths() > 0 && gfa_file.walks() > 0)
  {
    if(parameters.show_progress)
    {
      std::cerr << "Storing reference paths as sample " << REFERENCE_PATH_SAMPLE_NAME << std::endl;
    }
  }
  if(gfa_file.paths() == 0 && gfa_file.walks() == 0)
  {
    std::cerr << "check_gfa_file(): No paths or walks in the GFA file" << std::endl;
    return false;
  }

  return true;
}

gbwt::size_type
determine_batch_size(const GFAFile& gfa_file, const GFAParsingParameters& parameters)
{
  gbwt::size_type batch_size = parameters.batch_size;
  if(parameters.automatic_batch_size)
  {
    gbwt::size_type min_size = gbwt::DynamicGBWT::MIN_SEQUENCES_PER_BATCH * (gfa_file.max_path_length + 1);
    batch_size = std::max(min_size, batch_size);
    batch_size = std::min(static_cast<gbwt::size_type>(gfa_file.size()), batch_size);
  }
  if(parameters.show_progress)
  {
    std::cerr << "GBWT insertion batch size: " << batch_size << " nodes" << std::endl;
  }
  return batch_size;
}

std::unique_ptr<SequenceSource>
parse_segments(const GFAFile& gfa_file, const GFAParsingParameters& parameters)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Parsing segments" << std::endl;
  }

  // Determine if we need translation.
  bool translate = false;
  size_t max_node_length = (parameters.max_node_length == 0 ? std::numeric_limits<size_t>::max() : parameters.max_node_length);
  if(gfa_file.max_segment_length > max_node_length)
  {
    translate = true;
    if(parameters.show_progress)
    {
      std::cerr << "Breaking segments into " << max_node_length << " bp nodes" << std::endl;
    }
  }
  else if(gfa_file.translate_segment_ids)
  {
    translate = true;
    if(parameters.show_progress)
    {
      std::cerr << "Translating segment ids into valid node ids" << std::endl;
    }
  }

  std::unique_ptr<SequenceSource> result(new SequenceSource());
  gfa_file.for_each_segment([&](const std::string& name, view_type sequence) -> bool
  {
    if(translate)
    {
      result->translate_segment(name, sequence, max_node_length);
    }
    else
    {
      result->add_node(stoul_unsafe(name), sequence);
    }
    return true;
  });

  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Parsed " << result->get_node_count() << " nodes in " << seconds << " seconds" << std::endl;
  }
  return result;
}

bool
parse_metadata(const GFAFile& gfa_file, const GFAParsingParameters& parameters, PathNameParser& parser, gbwt::GBWTBuilder& builder)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Parsing metadata" << std::endl;
  }
  builder.index.addMetadata();

  // Parse walks.
  if(gfa_file.walks() > 0)
  {
    // Parse reference paths.
    bool failed = false;
    if(gfa_file.paths() > 0)
    {
      gfa_file.for_each_path_name([&](const std::string& name) -> bool
      {
        if(!(parser.add_reference_path(name))) { failed = true; return false; }
        return true;
      });
      if(failed)
      {
        std::cerr << "parse_metadata(): Could not parse GBWT metadata from reference path names" << std::endl;
        return false;
      }
    }
    // Parse walks.
    gfa_file.for_each_walk_name([&](const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start) -> bool
    {
      if(!(parser.add_walk(sample, haplotype, contig, start))) { failed = true; return false; }
      return true;
    });
    if(failed)
    {
      std::cerr << "parse_metadata(): Could not parse GBWT metadata from walks" << std::endl;
      return false;
    }
  }

  // Parse paths.
  else if(gfa_file.paths() > 0)
  {
    bool failed = false;
    gfa_file.for_each_path_name([&](const std::string& name) -> bool
    {
      if(!(parser.parse(name))) { failed = true; return false; }
      return true;
    });
    if(failed)
    {
      std::cerr << "parse_metadata(): Could not parse GBWT metadata from path names" << std::endl;
      return false;
    }
  }

  parser.set_metadata(builder.index.metadata);
  parser.clear();
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Parsed metadata in " << seconds << " seconds" << std::endl;
    std::cerr << "Metadata: "; gbwt::operator<<(std::cerr, builder.index.metadata) << std::endl;
  }
  return true;
}

void
parse_paths(const GFAFile& gfa_file, const GFAParsingParameters& parameters, const SequenceSource& source, gbwt::GBWTBuilder& builder)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Indexing paths/walks" << std::endl;
  }

  gbwt::vector_type current_path;
  auto add_segment = [&](const std::string& name, bool is_reverse) -> bool
  {
    if(source.uses_translation())
    {
      std::pair<nid_t, nid_t> range = source.get_translation(name);
      if(is_reverse)
      {
        for(nid_t id = range.second; id > range.first; id--)
        {
          current_path.push_back(gbwt::Node::encode(id - 1, is_reverse));
        }
      }
      else
      {
        for(nid_t id = range.first; id < range.second; id++)
        {
          current_path.push_back(gbwt::Node::encode(id, is_reverse));
        }
      }
    }
    else
    {
      current_path.push_back(gbwt::Node::encode(stoul_unsafe(name), is_reverse));
    }
    return true;
  };

  // Parse paths.
  gfa_file.for_each_path([&](const std::string&) -> bool
  {
    return true;
  }, add_segment, [&]() -> bool
  {
    builder.insert(current_path, true); current_path.clear();
    return true;
  });

  // Parse walks
  gfa_file.for_each_walk([&](const std::string&, const std::string&, const std::string&, const std::string&) -> bool
  {
    return true;
  }, add_segment, [&]() -> bool
  {
    builder.insert(current_path, true); current_path.clear();
    return true;
  });

  // Finish construction.
  builder.finish();
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Indexed " << gfa_file.paths() << " paths and " << gfa_file.walks() << " walks in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<SequenceSource>>
gfa_to_gbwt(const std::string& gfa_filename, const GFAParsingParameters& parameters)
{
  // Path name parsing.
  PathNameParser parser(parameters);
  if(!parser.ok)
  {
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }

  // GFA parsing.
  GFAFile gfa_file(gfa_filename, parameters.show_progress);
  if(!check_gfa_file(gfa_file, parameters))
  {
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }

  // Adjust batch size by GFA size and maximum path length.
  gbwt::size_type batch_size = determine_batch_size(gfa_file, parameters);

  // Parse segments.
  std::unique_ptr<SequenceSource> source = parse_segments(gfa_file, parameters);

  // Parse metadata from path names and walks.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::GBWTBuilder builder(parameters.node_width, batch_size, parameters.sample_interval);
  if(!parse_metadata(gfa_file, parameters, parser, builder))
  {
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }

  // Build GBWT from the paths and the walks.
  parse_paths(gfa_file, parameters, *source, builder);

  return std::make_pair(std::unique_ptr<gbwt::GBWT>(new gbwt::GBWT(builder.index)), std::move(source));
}

//------------------------------------------------------------------------------

// Cache segment names and lengths (in nodes). Assume that segment names are short enough
// that small string optimization avoids unnecessary memory allocations.
struct SegmentCache
{
  explicit SegmentCache(const GBWTGraph& graph) :
    graph(graph), segments((graph.index->sigma() - graph.index->firstNode()) / 2)
  {
    if(graph.has_segment_names())
    {
      graph.for_each_segment([&](const std::string& name, std::pair<nid_t, nid_t> nodes) -> bool
      {
        size_t relative = (gbwt::Node::encode(nodes.first, false) - graph.index->firstNode()) / 2;
        size_t length = nodes.second - nodes.first;
        for(size_t i = relative; i < relative + length; i++)
        {
          this->segments[i] = std::pair<size_t, size_t>(this->names.size(), length);
        }
        this->names.emplace_back(name);
        return true;
      });
    }
    else
    {
      graph.for_each_handle([&](const handle_t& handle)
      {
        size_t relative = (GBWTGraph::handle_to_node(handle) - graph.index->firstNode()) / 2;
        this->segments[relative] = std::pair<size_t, size_t>(this->names.size(), 1);
        this->names.emplace_back(std::to_string(graph.get_id(handle)));
      });
    }
  }

  size_t size() const { return this->names.size(); }

  std::pair<view_type, size_t> get(const handle_t& handle) const
  {
    return this->get(GBWTGraph::handle_to_node(handle));
  }

  std::pair<view_type, size_t> get(gbwt::node_type node) const
  {
    size_t relative = (node - this->graph.index->firstNode()) / 2;
    size_t offset = this->segments[relative].first;
    return std::make_pair(str_to_view(this->names[offset]), this->segments[relative].second);
  }

  const GBWTGraph& graph;

  // This vector goes over the same range as `graph.real_nodes`. The first component
  // is offset in `names` and the second is the length of the segment in nodes.
  std::vector<std::pair<size_t, size_t>> segments;
  std::vector<std::string> names;
};

//------------------------------------------------------------------------------

void
write_segments(const GBWTGraph& graph, const SegmentCache& cache, TSVWriter& writer, bool show_progress)
{
  double start = gbwt::readTimer();
  size_t segments = 0;
  if(show_progress)
  {
    std::cerr << "Writing segments" << std::endl;
  }

  view_type prev(nullptr, 0);
  graph.for_each_handle([&](const handle_t& handle)
  {
    auto segment = cache.get(handle);
    if(segment.first != prev)
    {
      if(prev.first != nullptr) { writer.newline(); }
      prev = segment.first;
      writer.put('S'); writer.newfield();
      writer.write(segment.first); writer.newfield();
      segments++;
    }
    writer.write(graph.get_sequence_view(handle));
  });
  if(prev.first != nullptr) { writer.newline(); }

  if(show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Wrote " << segments << " segments in " << seconds << " seconds" << std::endl;
  }
}

void
write_links(const GBWTGraph& graph, const SegmentCache& cache, TSVWriter& writer, bool show_progress)
{
  double start = gbwt::readTimer();
  size_t links = 0;
  if(show_progress)
  {
    std::cerr << "Writing links" << std::endl;
  }

  if(graph.has_segment_names())
  {
    // TODO this could be faster with for_each_edge and the cache.
    graph.for_each_link([&](const edge_t& edge, const std::string& from, const std::string& to) -> bool
    {
      writer.put('L'); writer.newfield();
      writer.write(from); writer.newfield();
      writer.put((graph.get_is_reverse(edge.first) ? '-' : '+')); writer.newfield();
      writer.write(to); writer.newfield();
      writer.put((graph.get_is_reverse(edge.second) ? '-' : '+')); writer.newfield();
      writer.put('*'); writer.newline();
      links++;
      return true;
    });
  }
  else
  {
    graph.for_each_edge([&](const edge_t& edge)
    {
      writer.put('L'); writer.newfield();
      writer.write(cache.get(edge.first).first); writer.newfield();
      writer.put((graph.get_is_reverse(edge.first) ? '-' : '+')); writer.newfield();
      writer.write(cache.get(edge.second).first); writer.newfield();
      writer.put((graph.get_is_reverse(edge.second) ? '-' : '+')); writer.newfield();
      writer.put('*'); writer.newline();
      links++;
    });
  }

  if(show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Wrote " << links << " links in " << seconds << " seconds" << std::endl;
  }
}

void
write_paths(const GBWTGraph& graph, const SegmentCache& cache, TSVWriter& writer, gbwt::size_type ref_sample, bool show_progress)
{
  double start = gbwt::readTimer();
  if(show_progress)
  {
    std::cerr << "Writing paths" << std::endl;
  }

  const gbwt::GBWT& index = *(graph.index);
  std::vector<gbwt::size_type> ref_paths = index.metadata.pathsForSample(ref_sample);
  for(gbwt::size_type path_id : ref_paths)
  {
    gbwt::vector_type path = index.extract(gbwt::Path::encode(path_id, false));
    writer.put('P'); writer.newfield();
    writer.write(index.metadata.contig(index.metadata.path(path_id).contig)); writer.newfield();
    size_t segments = 0, offset = 0;
    while(offset < path.size())
    {
      auto segment = cache.get(path[offset]);
      writer.write(segment.first);
      writer.put((gbwt::Node::is_reverse(path[offset]) ? '-' : '+'));
      segments++; offset += segment.second;
      if(offset < path.size()) { writer.put(','); }
    }
    writer.newfield();
    for(size_t i = 1; i < segments; i++)
    {
      writer.put('*');
      if(i + 1 < segments) { writer.put(','); }
    }
    writer.newline();
  }

  if(show_progress && !(ref_paths.empty()))
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Wrote " << ref_paths.size() << " paths in " << seconds << " seconds" << std::endl;
  }
}

void
write_walks(const GBWTGraph& graph, const SegmentCache& cache, TSVWriter& writer, gbwt::size_type ref_sample, bool show_progress)
{
  double start = gbwt::readTimer();
  size_t walks = 0;
  if(show_progress)
  {
    std::cerr << "Writing walks" << std::endl;
  }

  const gbwt::GBWT& index = *(graph.index);
  for(gbwt::size_type path_id = 0; path_id < index.metadata.paths(); path_id++)
  {
    const gbwt::PathName& path_name = index.metadata.path(path_id);
    if(path_name.sample == ref_sample) { continue; }
    walks++;
    gbwt::vector_type path = index.extract(gbwt::Path::encode(path_id, false));
    size_t length = 0;
    for(auto node : path) { length += graph.get_length(GBWTGraph::node_to_handle(node)); }
    writer.put('W'); writer.newfield();
    writer.write(index.metadata.sample(path_name.sample)); writer.newfield();
    writer.write(path_name.phase); writer.newfield();
    writer.write(index.metadata.contig(path_name.contig)); writer.newfield();
    writer.write(path_name.count); writer.newfield();
    writer.write(path_name.count + length); writer.newfield();
    size_t offset = 0;
    while(offset < path.size())
    {
      auto segment = cache.get(path[offset]);
      writer.put((gbwt::Node::is_reverse(path[offset]) ? '<' : '>'));
      writer.write(segment.first);
      offset += segment.second;
    }
    writer.newline();
  }

  if(show_progress && walks > 0)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Wrote " << walks << " walks in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

void
gbwt_to_gfa(const GBWTGraph& graph, std::ostream& out, bool show_progress)
{
  if(!(graph.index->hasMetadata() && graph.index->metadata.hasPathNames()))
  {
    throw InvalidGBWT("gbwt_to_gfa: The GBWT index must contain metadata with path names");
  }

  // Cache segment names.
  double start = gbwt::readTimer();
  if(show_progress)
  {
    std::cerr << "Caching segments" << std::endl;
  }
  SegmentCache cache(graph);
  if(show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Cached " << cache.size() << " segments in " << seconds << " seconds" << std::endl;
  }

  // GFA header.
  TSVWriter writer(out);
  writer.put('H'); writer.newfield();
  writer.write(std::string("VN:Z:1.0")); writer.newline();

  // Write the graph.
  write_segments(graph, cache, writer, show_progress);
  write_links(graph, cache, writer, show_progress);

  // Write the paths.
  gbwt::size_type ref_sample = graph.index->metadata.sample(REFERENCE_PATH_SAMPLE_NAME);
  write_paths(graph, cache, writer, ref_sample, show_progress);
  write_walks(graph, cache, writer, ref_sample, show_progress);
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
