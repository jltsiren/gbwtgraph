#include <gbwtgraph/gfa.h>

#include <algorithm>
#include <functional>
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
  size_t segments, paths;

  // Bit masks for field separators.
  size_t field_end[4];
  size_t subfield_end[4];

  // Pointers to line starts.
  std::vector<const char*> s_lines;
  std::vector<const char*> p_lines;

  struct field_type
  {
    const char* begin;
    const char* end;
    bool        has_next;

    field_type(const char* begin, const char* end, bool has_next) : begin(begin), end(end), has_next(has_next) {}

    size_t size() const { return this->end - this->begin; }
    bool empty() const { return (this->size() == 0); }

    std::string str() const { return std::string(this->begin, this->end); }
    view_type view() const { return view_type(this->begin, this->end - this->begin); }
    std::string oriented_str() const { return std::string(this->begin, this->end - 1); }
    char front() const { return *(this->begin); }
    char back() const { return *(this->end - 1); }
    bool valid_orientation() const { return (this->back() == '-' || this->back() == '+'); }
    bool is_reverse() const { return (this->back() == '-'); }
  };

  // Memory map and validate a GFA file.
  GFAFile(const std::string& filename, bool show_progress);

  ~GFAFile();

  bool ok() const { return (this->fd >= 0 && this->file_size > 0 && this->ptr != nullptr && this->valid_gfa); }
  size_t size() const { return this->file_size; }

private:
  // Preprocess a new S-line. Returns an iterator at the start of the next line or
  // nullptr if the parse failed.
  const char* add_s_line(const char* iter, std::unordered_set<std::string>& found_segments, size_t line_num);

  // Preprocess a new P-line. Returns an iterator at the start of the next line or
  // nullptr if the parse failed.
  const char* add_p_line(const char* iter, std::unordered_set<std::string>& found_paths, std::unordered_set<std::string>& required_segments, size_t line_num);

  const char* begin() const { return this->ptr; }
  const char* end() const { return this->ptr + this->size(); }

  // Return an iterator to the beginning of the next line.
  const char* next_line(const char* iter) const
  {
    while(iter != this->end() && *iter != '\n') { ++iter; }
    if(iter != this->end()) { ++iter; }
    return iter;
  }

  // Return the tab-separated field starting at the iterator.
  field_type get_field(const char* iter) const
  {
    const char* limit = iter;
    while(limit != this->end() && !(this->is_field_end(limit))) { ++limit; }
    return field_type(iter, limit, (limit != this->end() && *limit == '\t'));
  }

  // Return the comma-separated subfield starting at the iterator.
  field_type get_subfield(const char* iter) const
  {
    const char* limit = iter;
    while(limit != this->end() && !(this->is_subfield_end(limit))) { ++limit; }
    return field_type(iter, limit, (limit != this->end() && *limit == ','));
  }

  bool is_segment(const char* iter) const { return *iter == 'S'; }
  bool is_path(const char* iter) const { return *iter == 'P'; }

  bool is_field_end(const char* iter) const {
    unsigned char c = *iter;
    return (this->field_end[c / 64] & (size_t(1) << (c & 0x3F)));
  }

  bool is_subfield_end(const char* iter) const {
    unsigned char c = *iter;
    return (this->subfield_end[c / 64] & (size_t(1) << (c & 0x3F)));
  }

  void add_field_end(unsigned char c)
  {
    this->field_end[c / 64] |= size_t(1) << (c & 0x3F);
  }

  void add_subfield_end(unsigned char c)
  {
    this->subfield_end[c / 64] |= size_t(1) << (c & 0x3F);
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
                     const std::function<bool(const std::string& name)>& finish_path) const;
};

//------------------------------------------------------------------------------

GFAFile::GFAFile(const std::string& filename, bool show_progress) :
  fd(-1), file_size(0), ptr(nullptr),
  valid_gfa(true), translate_segment_ids(false),
  max_segment_length(0), max_path_length(0),
  segments(0), paths(0)
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
  this->ptr = static_cast<char*>(temp_ptr);

  // Mark characters indicating field/subfield end. This could depend on the GFA version.
  // TODO: If that happens, we need variables for field/subfield separators.
  this->field_end[0] = 0; this->field_end[1] = 0;
  this->field_end[2] = 0; this->field_end[3] = 0;
  this->add_field_end('\n'); this->add_field_end('\t');
  this->subfield_end[0] = 0; this->subfield_end[1] = 0;
  this->subfield_end[2] = 0; this->subfield_end[3] = 0;
  this->add_subfield_end('\n'); this->add_subfield_end('\t'); this->add_subfield_end(',');

  // Preprocess and validate the file.
  double start = gbwt::readTimer();
  if(show_progress)
  {
    std::cerr << "Validating GFA file " << filename << std::endl;
  }
  const char* iter = this->begin();
  size_t line_num = 0;
  std::unordered_set<std::string> found_paths, found_segments, required_segments;
  while(iter != nullptr && iter != this->end())
  {
    if(this->is_segment(iter))
    {
      iter = this->add_s_line(iter, found_segments, line_num);
    }
    else if(this->is_path(iter))
    {
      iter = this->add_p_line(iter, found_paths, required_segments, line_num);
    }
    else
    {
      iter = this->next_line(iter);
    }
    line_num++;
  }

  for(std::string segment_name : required_segments)
  {
    if(found_segments.find(segment_name) == found_segments.end())
    {
      std::cerr << "GFAFile::GFAFile(): Segment " << segment_name << " used on paths but missing from the file" << std::endl;
      this->valid_gfa = false;
      return;
    }
  }

  if(show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Found " << this->segments << " segments and " << this->paths << " paths in " << seconds << " seconds" << std::endl;
  }
}

const char*
GFAFile::add_s_line(const char* iter, std::unordered_set<std::string>& found_segments, size_t line_num)
{
  this->s_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->get_field(iter);
  if(!(field.has_next))
  {
    std::cerr << "GFAFile::GFAFile(): Segment line " << line_num << " ended after record type field" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }

  // Segment name field.
  field = this->get_field(field.end + 1);
  if(field.empty())
  {
    std::cerr << "GFAFile::GFAFile(): Segment line " << line_num << " has no name" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  if(!(field.has_next))
  {
    std::cerr << "GFAFile::GFAFile(): Segment line " << line_num << " ended after name field" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  std::string name = field.str();
  if(found_segments.find(name) != found_segments.end())
  {
    std::cerr << "GFAFile::GFAFile(): Duplicate segment name " << name << " on line " << line_num << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  found_segments.insert(name);
  if(!(this->translate_segment_ids))
  {
    try
    {
      nid_t id = std::stoul(name);
      if(id == 0) { this->translate_segment_ids = true; }
    }
    catch(const std::invalid_argument&) { this->translate_segment_ids = true; }
  }

  // Sequence field.
  field = this->get_field(field.end + 1);
  if(field.empty())
  {
    std::cerr << "GFAFile::GFAFile(): Segment line " << line_num << " has no sequence" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  this->max_segment_length = std::max(this->max_segment_length, field.size());
  this->segments++;

  return this->next_line(field.end);
}

const char*
GFAFile::add_p_line(const char* iter, std::unordered_set<std::string>& found_paths, std::unordered_set<std::string>& required_segments, size_t line_num)
{
  this->p_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->get_field(iter);
  if(!(field.has_next))
  {
    std::cerr << "GFAFile::GFAFile(): Path line " << line_num << " ended after record type field" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }

  // Path name field.
  field = this->get_field(field.end + 1);
  if(field.empty())
  {
    std::cerr << "GFAFile::GFAFile(): Path line " << line_num << " has no name" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  if(!(field.has_next))
  {
    std::cerr << "GFAFile::GFAFile(): Path line " << line_num << " ended after name field" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  std::string name = field.str();
  if(found_paths.find(name) != found_paths.end())
  {
    std::cerr << "GFAFile::GFAFile(): Duplicate path name " << name << " on line " << line_num << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  found_paths.insert(name);

  // Segment names field.
  size_t path_length = 0;
  do
  {
    field = this->get_subfield(field.end + 1);
    if(field.size() < 2 || !(field.valid_orientation()))
    {
      std::cerr << "GFAFile::GFAFile(): Invalid path segment " << field.str() << " on line " << line_num << std::endl;
      this->valid_gfa = false;
      return nullptr;
    }
    required_segments.insert(field.oriented_str());
    path_length++;
  }
  while(field.has_next);
  if(path_length == 0)
  {
    std::cerr << "GFAFile::GFAFile(): Path " << name << " on line " << line_num << " is empty" << std::endl;
    this->valid_gfa = false;
    return nullptr;
  }
  this->max_path_length = std::max(this->max_path_length, path_length);
  this->paths++;

  return this->next_line(field.end);
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
    field_type field = this->get_field(iter);

    // Segment name field.
    field = this->get_field(field.end + 1);
    std::string name = field.str();

    // Sequence field.
    field = this->get_field(field.end + 1);
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
    field_type field = this->get_field(iter);

    // Path name field.
    field = this->get_field(field.end + 1);
    std::string path_name = field.str();
    if(!path(path_name)) { return; }
  }
}

void
GFAFile::for_each_path(const std::function<bool(const std::string& name)>& path,
                       const std::function<bool(const std::string& name, bool is_reverse)>& path_segment,
                       const std::function<bool(const std::string& name)>& finish_path) const
{
  if(!(this->ok())) { return; }

  for(const char* iter : this->p_lines)
  {
    // Skip the record type field.
    field_type field = this->get_field(iter);

    // Path name field.
    field = this->get_field(field.end + 1);
    std::string path_name = field.str();
    if(!path(path_name)) { return; }

    // Segment names field.
    do
    {
      field = this->get_subfield(field.end + 1);
      std::string segment_name = field.oriented_str();
      if(!path_segment(segment_name, field.is_reverse())) { return; }
    }
    while(field.has_next);

    if(!finish_path(path_name)) { return; }
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

  bool ok;

  explicit PathNameParser(const GFAParsingParameters& parameters);

  // Parse a path name using a regex. Returns true if successful.
  // This should not be used with add_path().
  bool parse(const std::string& name);

  // Add a path based on walk metadata. Returns true if successful.
  // This should not be used with parse().
  // FIXME implement
  bool add_path(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start);

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
  ok(true)
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
PathNameParser::add_path(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)
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
      std::cerr << "PathNameParser::add_path(): Invalid haplotype field " << haplotype << std::endl;
      return false;
    }
  }
  this->haplotypes.insert(std::pair<size_t, size_t>(path_name.sample, path_name.phase));

  // Start position as fragment identifier.
  {
    try { path_name.count = std::stoul(start); }
    catch(const std::invalid_argument&)
    {
      std::cerr << "PathNameParser::add_path(): Invalid start_position " << start << std::endl;
      return false;
    }
    if(this->counts.find(path_name) != this->counts.end())
    {
      std::cerr << "PathNameParser::add_path(): Duplicate walk: " << sample << "\t" << haplotype << "\t" << contig << "\t" << start << ")" << std::endl;
      return false;
    }
    this->counts[path_name] = 1;
  }

  this->path_names.push_back(path_name);
  return true;
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
  if(!(gfa_file.ok()))
  {
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }
  if(gfa_file.segments == 0)
  {
    std::cerr << "gfa_to_gbwt(): No segments in the GFA file" << std::endl;
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }
  if(gfa_file.paths == 0)
  {
    std::cerr << "gfa_to_gbwt(): No paths in the GFA file" << std::endl;
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }

  // Adjust batch size by GFA size and maximum path length.
  gbwt::size_type batch_size = parameters.batch_size;
  if(parameters.automatic_batch_size)
  {
    // FIXME Source for the multiplier? Should it be a parameter?
    batch_size = std::max(static_cast<gbwt::size_type>(20 * (gfa_file.max_path_length + 1)), parameters.batch_size);
    batch_size = std::min(static_cast<gbwt::size_type>(gfa_file.size()), parameters.batch_size);
  }
  if(parameters.show_progress)
  {
    std::cerr << "GBWT insertion batch size: " << batch_size << " nodes" << std::endl;
  }

  bool translate = false;
  if(gfa_file.max_segment_length > parameters.max_node_length)
  {
    translate = true;
    if(parameters.show_progress)
    {
      std::cerr << "Breaking segments into " << parameters.max_node_length << " bp nodes" << std::endl;
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

  // Parse segments.
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Parsing segments" << std::endl;
  }
  std::unique_ptr<SequenceSource> source(new SequenceSource());
  gfa_file.for_each_segment([&](const std::string& name, view_type sequence) -> bool
  {
    if(translate)
    {
      source->translate_segment(name, sequence, parameters.max_node_length);
    }
    else
    {
      source->add_node(std::stoul(name), sequence);
    }
    return true;
  });
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Parsed " << source->get_node_count() << " nodes in " << seconds << " seconds" << std::endl;
  }

  // Parse metadata from path names.
  start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Parsing metadata" << std::endl;
  }
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::GBWTBuilder builder(parameters.node_width, batch_size, parameters.sample_interval);
  builder.index.addMetadata();
  bool failed = false;
  gfa_file.for_each_path_name([&](const std::string& name) -> bool
  {
    if(!(parser.parse(name))) { failed = true; return false; }
    return true;
  });
  if(failed)
  {
    std::cerr << "gfa_to_gbwt(): Could not parse GBWT metadata from path names" << std::endl;
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }
  parser.set_metadata(builder.index.metadata);
  parser.clear();
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Parsed metadata in " << seconds << " seconds" << std::endl;
    std::cerr << "Metadata: "; gbwt::operator<<(std::cerr, builder.index.metadata) << std::endl;
  }

  // Build GBWT from the paths.
  start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Indexing paths" << std::endl;
  }
  gbwt::vector_type current_path;
  gfa_file.for_each_path([&](const std::string&) -> bool
  {
    return true;
  }, [&](const std::string& name, bool is_reverse) -> bool
  {
    if(translate)
    {
      std::pair<nid_t, nid_t> range = source->get_translation(name);
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
      current_path.push_back(gbwt::Node::encode(std::stoul(name), is_reverse));
    }
    return true;
  }, [&](const std::string&) -> bool
  {
    builder.insert(current_path, true); current_path.clear();
    return true;
  });
  builder.finish();
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Indexed " << gfa_file.paths << " paths in " << seconds << " seconds" << std::endl;
  }

  return std::make_pair(std::unique_ptr<gbwt::GBWT>(new gbwt::GBWT(builder.index)), std::move(source));
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
