#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/gfa.h>
#include <gbwtgraph/internal.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <string>
#include <utility>

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

constexpr size_t GFAParsingParameters::APPROXIMATE_NUM_JOBS;
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
  bool translate_segment_ids;
  size_t max_segment_length, max_path_length;

  // Bit masks for field separators.
  size_t field_end[4];
  size_t subfield_end[4];
  size_t walk_subfield_end[4];

  // Pointers to line starts.
  std::vector<const char*> s_lines;
  std::vector<const char*> l_lines;
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

    // For segment orientations in links.
    bool valid_orientation() const { return (this->size() == 1 && (this->back() == '-' || this->back() == '+')); }
    bool is_reverse_orientation() const { return (this->back() == '-'); }

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

  // Memory map and validate a GFA file. The constructor checks that all mandatory
  // fields used for GBWTGraph construction exist and are nonempty.
  // There are no checks for duplicates.
  // Throws `std::runtime_error` on failure.
  GFAFile(const std::string& filename, bool show_progress);

  ~GFAFile();

  size_t size() const { return this->file_size; }
  size_t segments() const { return this->s_lines.size(); }
  size_t links() const { return this->l_lines.size(); }
  size_t paths() const { return this->p_lines.size(); }
  size_t walks() const { return this->w_lines.size(); }

private:
  // Preprocess a new S-line. Returns an iterator at the start of the next line or
  // throws `std::runtime_error` if the parse failed.
  const char* add_s_line(const char* iter, size_t line_num);

  // Preprocess a new S-line. Returns an iterator at the start of the next line or
  // throws `std::runtime_error` if the parse failed.
  const char* add_l_line(const char* iter, size_t line_num);

  // Preprocess a new P-line. Returns an iterator at the start of the next line or
  // throws `std::runtime_error` if the parse failed.
  const char* add_p_line(const char* iter, size_t line_num);

  // Preprocess a new W-line. Returns an iterator at the start of the next line or
  // throws `std::runtime_error` if the parse failed.
  const char* add_w_line(const char* iter, size_t line_num);

  // Throws `std::runtime_error` if the field is invalid.
  void check_field(const field_type& field, const std::string& field_name, bool should_have_next);

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
    Iterate over the S-lines, calling segment() for all segments.
  */
  void for_each_segment(const std::function<void(const std::string& name, view_type sequence)>& segment) const;

  /*
    Iterate over the L-lines, calling link() for all segments.
  */
 void for_each_link(const std::function<void(const std::string& from, bool from_is_reverse, const std::string& to, bool to_is_reverse)>& link) const;

  /*
    Iterate over the file, calling path_start() for each path.
  */
  void for_each_path_start(const std::function<void(const char* line_start, const std::string& first_segment)>& path_start) const;

  /*
    Iterate over the file, calling path() for the selected paths.
  */
  void for_these_path_names(const std::vector<const char*>& selected_paths, const std::function<void(const std::string& name)>& path) const;

  /*
    Iterate over the file, calling path() for the selected paths,
    path_segment() for each path segment, and finish_path() after
    parsing each path.
  */
  void for_these_paths(const std::vector<const char*>& selected_paths,
                       const std::function<void(const std::string& name)>& path,
                       const std::function<void(const std::string& name, bool is_reverse)>& path_segment,
                       const std::function<void()>& finish_path) const;

  /*
    Iterate over the file, calling walk_start() for each walk.
  */
  void for_each_walk_start(const std::function<void(const char* line_start, const std::string& first_segment)>& walk_start) const;

  /*
    Iterate over the file, calling walk() for the selected walks.
  */
  void for_these_walk_names(const std::vector<const char*>& selected_walks,
                            const std::function<void(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)>& walk) const;

  /*
    Iterate over the file, calling walk() for the selected walks,
    walk_segment() for each walk segment, and finish_walk() after
    parsing each walk.
  */
  void for_these_walks(const std::vector<const char*>& selected_walks,
                       const std::function<void(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)>& walk,
                       const std::function<void(const std::string& name, bool is_reverse)>& walk_segment,
                       const std::function<void()>& finish_walk) const;
};

//------------------------------------------------------------------------------

GFAFile::GFAFile(const std::string& filename, bool show_progress) :
  fd(-1), file_size(0), ptr(nullptr),
  translate_segment_ids(false),
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
    throw std::runtime_error("GFAFile: Cannot open file " + filename);
  }

  // Memory map the file.
  struct stat st;
  if(::fstat(this->fd, &st) < 0)
  {
    throw std::runtime_error("GFAFile: Cannot stat file " + filename);
  }
  this->file_size = st.st_size;

  void* temp_ptr = ::mmap(nullptr, file_size, PROT_READ, MAP_FILE | MAP_SHARED, this->fd, 0);
  if(temp_ptr == MAP_FAILED)
  {
    throw std::runtime_error("GFAFile: Cannot memory map file " + filename);
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
  while(iter != this->end())
  {
    switch(*iter)
    {
    case 'S':
      iter = this->add_s_line(iter, line_num);
      break;
    case 'L':
      iter = this->add_l_line(iter, line_num);
      break;
    case 'P':
      iter = this->add_p_line(iter, line_num);
      break;
    case 'W':
      iter = this->add_w_line(iter, line_num);
      break;
    default:
      iter = this->next_line(iter);
      break;
    }
    if(iter == nullptr) { return; }
    line_num++;
  }

  if(show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Found " << this->segments() << " segments, " << this->links() << " links, " << this->paths() << " paths, and " << this->walks() << " walks in " << seconds << " seconds" << std::endl;
  }
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

const char*
GFAFile::add_s_line(const char* iter, size_t line_num)
{
  this->s_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->first_field(iter, line_num);
  this->check_field(field, "record type", true);

  // Segment name field.
  field = this->next_field(field);
  this->check_field(field, "segment name", true);
  std::string name = field.str();
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
  this->check_field(field, "sequence", false);
  this->max_segment_length = std::max(this->max_segment_length, field.size());

  return this->next_line(field.end);
}

const char*
GFAFile::add_l_line(const char* iter, size_t line_num)
{
  this->l_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->first_field(iter, line_num);
  this->check_field(field, "record type", true);

  // Source segment field.
  field = this->next_field(field);
  this->check_field(field, "source segment", true);

  // Source orientation field.
  field = this->next_field(field);
  this->check_field(field, "source orientation", true);
  if(!(field.valid_orientation()))
  {
    throw std::runtime_error("GFAFile: Invalid source orientation " + field.str() + " on line " + std::to_string(line_num));
  }

  // Destination segment field.
  field = this->next_field(field);
  this->check_field(field, "destination segment", true);

  // Destination orientation field.
  field = this->next_field(field);
  this->check_field(field, "destination orientation", false);
  if(!(field.valid_orientation()))
  {
    throw std::runtime_error("GFAFile: Invalid destination orientation " + field.str() + " on line " + std::to_string(line_num));
  }

  return this->next_line(field.end);
}

const char*
GFAFile::add_p_line(const char* iter, size_t line_num)
{
  this->p_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->first_field(iter, line_num);
  this->check_field(field, "record type", true);

  // Path name field.
  field = this->next_field(field);
  this->check_field(field, "path name", true);

  // Segment names field.
  size_t path_length = 0;
  do
  {
    field = this->next_subfield(field);
    if(!(field.valid_path_segment()))
    {
      throw std::runtime_error("GFAFile: Invalid path segment " + field.str() + " on line " + std::to_string(line_num));
    }
    path_length++;
  }
  while(field.has_next);
  if(path_length == 0)
  {
    throw std::runtime_error("GFAFile: The path on line " + std::to_string(line_num) + " is empty");
  }
  this->max_path_length = std::max(this->max_path_length, path_length);

  return this->next_line(field.end);
}

const char*
GFAFile::add_w_line(const char* iter, size_t line_num)
{
  this->w_lines.push_back(iter);

  // Skip the record type field.
  field_type field = this->first_field(iter, line_num);
  this->check_field(field, "record type", true);

  // Sample name field.
  field = this->next_field(field);
  this->check_field(field, "sample name", true);

  // Haplotype index field.
  field = this->next_field(field);
  this->check_field(field, "haplotype index", true);

  // Contig name field.
  field = this->next_field(field);
  this->check_field(field, "contig name", true);

  // Start position field.
  field = this->next_field(field);
  this->check_field(field, "start position", true);

  // Skip the end position field.
  field = this->next_field(field);
  this->check_field(field, "end position", true);

  // Segment names field.
  size_t path_length = 0;
  field.start_walk();
  do
  {
    field = this->next_walk_subfield(field);
    if(!(field.valid_walk_segment()))
    {
      throw std::runtime_error("GFAFile: Invalid walk segment " + field.str() + " on line " + std::to_string(line_num));
    }
    path_length++;
  }
  while(field.has_next);
  if(path_length == 0)
  {
    throw std::runtime_error("GFAFile: The walk on line " + std::to_string(line_num) + " is empty");
  }
  this->max_path_length = std::max(this->max_path_length, path_length);

  return this->next_line(field.end);
}

void
GFAFile::check_field(const field_type& field, const std::string& field_name, bool should_have_next)
{
  if(field.empty())
  {
    throw std::runtime_error("GFAFile: " + std::string(field.type, 1) + "-line " + std::to_string(field.line_num) + " has no " + field_name);
  }
  if(should_have_next && !(field.has_next))
  {
    throw std::runtime_error("GFAFile: " + std::string(field.type, 1) + "-line " + std::to_string(field.line_num) + " ended after " + field_name);
  }
}

//------------------------------------------------------------------------------

void
GFAFile::for_each_segment(const std::function<void(const std::string& name, view_type sequence)>& segment) const
{
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
    segment(name, sequence);
  }
}

void
GFAFile::for_each_link(const std::function<void(const std::string& from, bool from_is_reverse, const std::string& to, bool to_is_reverse)>& link) const
{
  for(const char* iter : this->l_lines)
  {
    // Skip the record type field.
    field_type field = this->first_field(iter);

    // Source segment field.
    field = this->next_field(field);
    std::string from = field.str();

    // Source orientation field.
    field = this->next_field(field);
    bool from_is_reverse = field.is_reverse_orientation();

    // Destination segment field.
    field = this->next_field(field);
    std::string to = field.str();

    // Destination orientation field.
    field = this->next_field(field);
    bool to_is_reverse = field.is_reverse_orientation();

    link(from, from_is_reverse, to, to_is_reverse);
  }
}

//------------------------------------------------------------------------------

void
GFAFile::for_each_path_start(const std::function<void(const char* line_start, const std::string& first_segment)>& path_start) const
{
  for(const char* iter : this->p_lines)
  {
    const char* line_start = iter;

    // Skip the record type field and the path name field.
    field_type field = this->first_field(iter);
    field = this->next_field(field);

    // First segment.
    field = this->next_subfield(field);
    std::string first_segment = field.path_segment();
    path_start(line_start, first_segment);
  }
}

void
GFAFile::for_these_path_names(const std::vector<const char*>& selected_paths, const std::function<void(const std::string& name)>& path) const
{
  for(const char* iter : selected_paths)
  {
    // Skip the record type field.
    field_type field = this->first_field(iter);

    // Path name field.
    field = this->next_field(field);
    std::string path_name = field.str();
    path(path_name);
  }
}

void
GFAFile::for_these_paths(const std::vector<const char*>& selected_paths,
                         const std::function<void(const std::string& name)>& path,
                         const std::function<void(const std::string& name, bool is_reverse)>& path_segment,
                         const std::function<void()>& finish_path) const
{
  for(const char* iter : selected_paths)
  {
    // Skip the record type field.
    field_type field = this->first_field(iter);

    // Path name field.
    field = this->next_field(field);
    path(field.str());

    // Segment names field.
    do
    {
      field = this->next_subfield(field);
      std::string segment_name = field.path_segment();
      path_segment(segment_name, field.is_reverse_path_segment());
    }
    while(field.has_next);

    finish_path();
  }
}

//------------------------------------------------------------------------------

void
GFAFile::for_each_walk_start(const std::function<void(const char* line_start, const std::string& first_segment)>& walk_start) const
{
  for(const char* iter : this->w_lines)
  {
    const char* line_start = iter;

    // Skip the record type field and metadata fields.
    field_type field = this->first_field(iter);
    field = this->next_field(field); // Sample.
    field = this->next_field(field); // Haplotype.
    field = this->next_field(field); // Contig.
    field = this->next_field(field); // Start.
    field = this->next_field(field); // End.

    // First segment.
    field.start_walk();
    field = this->next_walk_subfield(field);
    std::string first_segment = field.walk_segment();
    walk_start(line_start, first_segment);
  }
}

void
GFAFile::for_these_walk_names(const std::vector<const char*>& selected_walks,
                              const std::function<void(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)>& walk) const
{
  for(const char* iter : selected_walks)
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

    walk(sample, haplotype, contig, start);
  }
}

void
GFAFile::for_these_walks(const std::vector<const char*>& selected_walks,
                         const std::function<void(const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)>& walk,
                         const std::function<void(const std::string& name, bool is_reverse)>& walk_segment,
                         const std::function<void()>& finish_walk) const
{
  for(const char* iter : selected_walks)
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

    walk(sample, haplotype, contig, start);

    // Skip the end field.
    field = this->next_field(field);

    // Segment names field.
    field.start_walk();
    do
    {
      field = this->next_walk_subfield(field);
      std::string segment_name = field.walk_segment();
      walk_segment(segment_name, field.is_reverse_walk_segment());
    }
    while(field.has_next);

    finish_walk();
  }
}

//------------------------------------------------------------------------------

void
check_gfa_file(const GFAFile& gfa_file, const GFAParsingParameters& parameters)
{
  if(gfa_file.segments() == 0)
  {
    throw std::runtime_error("No segments in the GFA file");
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
    throw std::runtime_error("No paths or walks in the GFA file");
  }
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

struct ConstructionJob
{
  size_t num_nodes;
  size_t id;
  std::vector<const char*> p_lines;
  std::vector<const char*> w_lines;

  // Largest jobs first.
  bool operator<(const ConstructionJob& another) const
  {
    return (this->num_nodes > another.num_nodes);
  }
};

//------------------------------------------------------------------------------

std::pair<std::unique_ptr<SequenceSource>, std::unique_ptr<EmptyGraph>>
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

  std::pair<std::unique_ptr<SequenceSource>, std::unique_ptr<EmptyGraph>> result(new SequenceSource(), new EmptyGraph());
  gfa_file.for_each_segment([&](const std::string& name, view_type sequence)
  {
    if(translate)
    {
      std::pair<nid_t, nid_t> translation = result.first->translate_segment(name, sequence, max_node_length);
      for(nid_t id = translation.first; id < translation.second; id++)
      {
        result.second->create_node(id);
      }
    }
    else
    {
      nid_t id = stoul_unsafe(name);
      result.first->add_node(id, sequence);
      result.second->create_node(id);
    }
  });

  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Parsed " << result.first->get_node_count() << " nodes in " << seconds << " seconds" << std::endl;
  }
  return result;
}

void
parse_links(const GFAFile& gfa_file, const SequenceSource& source, EmptyGraph& graph, const GFAParsingParameters& parameters)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Parsing links" << std::endl;
  }

  size_t edge_count = 0;
  gfa_file.for_each_link([&](const std::string& from, bool from_is_reverse, const std::string& to, bool to_is_reverse)
  {
    std::pair<nid_t, nid_t> from_nodes = source.force_translate(from);
    if(from_nodes == SequenceSource::invalid_translation())
    {
      throw std::runtime_error("Invalid source segment " + from);
    }
    std::pair<nid_t, nid_t> to_nodes = source.force_translate(to);
    if(to_nodes == SequenceSource::invalid_translation())
    {
      throw std::runtime_error("Invalid destination segment " + to);
    }
    nid_t from_node = (from_is_reverse ? from_nodes.first : from_nodes.second - 1);
    nid_t to_node = (to_is_reverse ? to_nodes.second - 1 : to_nodes.first);
    graph.create_edge(graph.get_handle(from_node, from_is_reverse), graph.get_handle(to_node, to_is_reverse));
    edge_count++;
  });

  // Add edges inside segments if necessary.
  if(source.uses_translation())
  {
    for(auto iter = source.segment_translation.begin(); iter != source.segment_translation.end(); ++iter)
    {
      for(nid_t id = iter->second.first; id + 1 < iter->second.second; id++)
      {
        graph.create_edge(graph.get_handle(id, false), graph.get_handle(id + 1, false));
        edge_count++;
      }
    }
  }

  // Get rid of duplicate edges if the GFA graph had them.
  graph.remove_duplicate_edges();

  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Parsed " << edge_count << " edges in " << seconds << " seconds" << std::endl;
  }
}

gbwt::Metadata
parse_metadata(const GFAFile& gfa_file, const std::vector<ConstructionJob>& jobs, MetadataBuilder& metadata, const GFAParsingParameters& parameters)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Parsing metadata" << std::endl;
  }

  for(size_t i = 0; i < jobs.size(); i++)
  {
    if(gfa_file.walks() > 0)
    {
      // Parse reference paths.
      gfa_file.for_these_path_names(jobs[i].p_lines, [&](const std::string& name)
      {
        metadata.add_reference_path(name, jobs[i].id);
      });
      // Parse walks.
      gfa_file.for_these_walk_names(jobs[i].w_lines, [&](const std::string& sample, const std::string& haplotype, const std::string& contig, const std::string& start)
      {
        metadata.add_walk(sample, haplotype, contig, start, jobs[i].id);
      });
    }
    else
    {
      // Parse path names.
      gfa_file.for_these_path_names(jobs[i].p_lines, [&](const std::string& name)
      {
        metadata.parse(name, jobs[i].id);
      });
    }
  }

  gbwt::Metadata result = metadata.get_metadata();
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Metadata: "; gbwt::operator<<(std::cerr, result) << std::endl;
    std::cerr << "Parsed metadata in " << seconds << " seconds" << std::endl;
  }
  return result;
}

std::unique_ptr<gbwt::GBWT>
parse_paths(const GFAFile& gfa_file, const std::vector<ConstructionJob>& jobs, const SequenceSource& source, const GFAParsingParameters& parameters, gbwt::size_type batch_size)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Indexing paths/walks" << std::endl;
  }

  // Prepare for GBWT construction.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  omp_set_num_threads(std::max(parameters.parallel_jobs, size_t(1)));
  std::vector<gbwt::GBWT> partial_indexes(jobs.size());
  std::vector<gbwt::vector_type> current_paths(parameters.parallel_jobs);

  auto add_segment = [&](const std::string& name, bool is_reverse)
  {
    std::pair<nid_t, nid_t> range = source.force_translate(name);
    if(range == SequenceSource::invalid_translation())
    {
      throw std::runtime_error("Invalid segment " + name);
    }
    gbwt::vector_type& current_path = current_paths[omp_get_thread_num()];
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
  };

  // Build the partial indexes in parallel.
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t i = 0; i < jobs.size(); i++)
  {
    double job_start = gbwt::readTimer();
    if(parameters.show_progress)
    {
      #pragma omp critical
      {
        std::cerr << "Starting job " << i << " (" << jobs[i].num_nodes << " nodes, " << jobs[i].p_lines.size() << " paths, " << jobs[i].w_lines.size() << " walks)" << std::endl;
      }
    }
    gbwt::GBWTBuilder builder(parameters.node_width, batch_size, parameters.sample_interval);
    size_t thread_num = omp_get_thread_num();
    try
    {
      gfa_file.for_these_paths(jobs[i].p_lines, [&](const std::string&) {}, add_segment, [&]()
      {
        builder.insert(current_paths[thread_num], true);
        current_paths[thread_num].clear();
      });
      gfa_file.for_these_walks(jobs[i].w_lines, [&](const std::string&, const std::string&, const std::string&, const std::string&) {}, add_segment, [&]()
      {
        builder.insert(current_paths[thread_num], true);
        current_paths[thread_num].clear();
      });
    }
    catch(const std::runtime_error& e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }    
    builder.finish();
    // Order the partial indexes by original ids (minimum node ids) to make the construction deterministic.
    partial_indexes[jobs[i].id] = gbwt::GBWT(builder.index);
    // Deleting a dynamic GBWT is a bit expensive, so we do it manually to include it in the measured time.
    builder.index = gbwt::DynamicGBWT();
    if(parameters.show_progress)
    {
      double seconds = gbwt::readTimer() - job_start;
      #pragma omp critical
      {
        std::cerr << "Finished job " << i << " in " << seconds << " seconds" << std::endl;
      }
    }
  }

  // Merge the indexes.
  if(parameters.show_progress)
  {
    std::cerr << "Merging partial indexes" << std::endl;
  }
  std::unique_ptr<gbwt::GBWT> result(new gbwt::GBWT(partial_indexes));
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Indexed " << gfa_file.paths() << " paths and " << gfa_file.walks() << " walks in " << seconds << " seconds" << std::endl;
  }

  return result;
}

//------------------------------------------------------------------------------

// Graph will be cleared, as we do not need graph topology after this.
std::vector<ConstructionJob>
determine_jobs(const GFAFile& gfa_file, const SequenceSource& source, std::unique_ptr<EmptyGraph>& graph, const GFAParsingParameters& parameters)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Creating jobs" << std::endl;
  }

  // Find weakly connected components.
  std::vector<std::vector<nid_t>> components = weakly_connected_components(*graph);
 
  // Assign graph components to jobs.
  size_t nodes = graph->get_node_count();
  size_t target_size = nodes / std::max(size_t(1), parameters.approximate_num_jobs);
  std::unordered_map<nid_t, size_t> node_to_job;
  node_to_job.reserve(nodes);
  std::vector<ConstructionJob> jobs;
  for(size_t i = 0; i < components.size(); i++)
  {
    if(jobs.empty() || jobs.back().num_nodes + components[i].size() > target_size)
    {
      jobs.push_back({ 0, jobs.size(), {}, {} });
    }
    jobs.back().num_nodes += components[i].size();
    for(nid_t id : components[i]) { node_to_job[id] = jobs.back().id; }
  }

  // Assign P-lines and W-lines to jobs.
  gfa_file.for_each_path_start([&](const char* line_start, const std::string& first_segment)
  {
    nid_t id = source.force_translate(first_segment).first; // 0 on failure.
    auto iter = node_to_job.find(id);
    if(iter != node_to_job.end())
    {
      jobs[iter->second].p_lines.push_back(line_start);
    }
    else
    {
      throw std::runtime_error("Invalid path segment " + first_segment);
    }
  });
  gfa_file.for_each_walk_start([&](const char* line_start, const std::string& first_segment)
  {
    nid_t id = source.force_translate(first_segment).first; // 0 on failure.
    auto iter = node_to_job.find(id);
    if(iter != node_to_job.end())
    {
      jobs[iter->second].w_lines.push_back(line_start);
    }
    else
    {
      throw std::runtime_error("Invalid walk segment " + first_segment);
    }
  });

  // Sort the jobs to process largest ones first.
  std::sort(jobs.begin(), jobs.end());

  // Delete temporary structures before reporting time, as some structures are a bit complex.
  size_t num_components = components.size();
  components.clear(); node_to_job.clear(); graph.reset();
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Created " << jobs.size() << " jobs for " << num_components << " components in " << seconds << " seconds" << std::endl;
  }
  return jobs;
}

//------------------------------------------------------------------------------

std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<SequenceSource>>
gfa_to_gbwt(const std::string& gfa_filename, const GFAParsingParameters& parameters)
{
  // Metadata handling.
  MetadataBuilder metadata(parameters.path_name_regex, parameters.path_name_fields);

  // GFA parsing.
  GFAFile gfa_file(gfa_filename, parameters.show_progress);
  check_gfa_file(gfa_file, parameters);

  // Adjust batch size by GFA size and maximum path length.
  gbwt::size_type batch_size = determine_batch_size(gfa_file, parameters);

  // Parse segments.
  std::unique_ptr<SequenceSource> source;
  std::unique_ptr<EmptyGraph> graph;
  std::tie(source, graph) = parse_segments(gfa_file, parameters);

  // Parse links and create jobs.
  parse_links(gfa_file, *source, *graph, parameters);
  std::vector<ConstructionJob> jobs = determine_jobs(gfa_file, *source, graph, parameters);

  // Build the GBWT index.
  gbwt::Metadata final_metadata = parse_metadata(gfa_file, jobs, metadata, parameters);
  std::unique_ptr<gbwt::GBWT> gbwt_index = parse_paths(gfa_file, jobs, *source, parameters, batch_size);
  gbwt_index->addMetadata();
  gbwt_index->metadata = final_metadata;

  return std::make_pair(std::move(gbwt_index), std::move(source));
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
write_links(const GBWTGraph& graph, const SegmentCache& cache, std::ostream& out, const GFAExtractionParameters& parameters)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Writing links" << std::endl;
  }

  size_t threads = parameters.threads();
  std::vector<ManualTSVWriter> writers(threads, ManualTSVWriter(out));
  std::vector<size_t> links(threads, 0);

  if(graph.has_segment_names())
  {
    graph.for_each_link([&](const edge_t& edge, const std::string& from, const std::string& to) -> bool
    {
      size_t thread = omp_get_thread_num();
      ManualTSVWriter& writer = writers[thread];
      writer.put('L'); writer.newfield();
      writer.write(from); writer.newfield();
      writer.put((graph.get_is_reverse(edge.first) ? '-' : '+')); writer.newfield();
      writer.write(to); writer.newfield();
      writer.put((graph.get_is_reverse(edge.second) ? '-' : '+')); writer.newfield();
      writer.put('*'); writer.newline();
      links[thread]++;
      if(writer.full())
      {
        #pragma omp critical
        {
          writer.flush();
        }
      }
      return true;
    }, (threads > 1));
  }
  else
  {
    graph.for_each_edge([&](const edge_t& edge)
    {
      size_t thread = omp_get_thread_num();
      ManualTSVWriter& writer = writers[thread];
      writer.put('L'); writer.newfield();
      writer.write(cache.get(edge.first).first); writer.newfield();
      writer.put((graph.get_is_reverse(edge.first) ? '-' : '+')); writer.newfield();
      writer.write(cache.get(edge.second).first); writer.newfield();
      writer.put((graph.get_is_reverse(edge.second) ? '-' : '+')); writer.newfield();
      writer.put('*'); writer.newline();
      links[thread]++;
      if(writer.full())
      {
        #pragma omp critical
        {
          writer.flush();
        }
      }
    }, (threads > 1));
  }

  // Flush the writers and count the total number of links.
  size_t total_links = 0;
  for(size_t i = 0; i < threads; i++)
  {
    writers[i].flush();
    total_links += links[i];
  }

  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Wrote " << total_links << " links in " << seconds << " seconds" << std::endl;
  }
}

void
write_paths(const GBWTGraph& graph, const SegmentCache& cache, std::ostream& out, gbwt::size_type ref_sample, const GFAExtractionParameters& parameters)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Writing reference paths" << std::endl;
  }

  const gbwt::GBWT& index = *(graph.index);
  std::vector<gbwt::size_type> ref_paths = index.metadata.pathsForSample(ref_sample);
  std::vector<ManualTSVWriter> writers(parameters.threads(), ManualTSVWriter(out));

  #pragma omp parallel for schedule(dynamic, 1)
  for(gbwt::size_type path_id : ref_paths)
  {
    ManualTSVWriter& writer = writers[omp_get_thread_num()];
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
    #pragma omp critical
    {
      writer.flush();
    }
  }

  if(parameters.show_progress && !(ref_paths.empty()))
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Wrote " << ref_paths.size() << " paths in " << seconds << " seconds" << std::endl;
  }
}

void
write_walks(const GBWTGraph& graph, const SegmentCache& cache, std::ostream& out, gbwt::size_type ref_sample, const GFAExtractionParameters& parameters)
{
  double start = gbwt::readTimer();
  size_t walks = 0;
  if(parameters.show_progress)
  {
    std::cerr << "Writing walks" << std::endl;
  }

  const gbwt::GBWT& index = *(graph.index);
  std::vector<ManualTSVWriter> writers(parameters.threads(), ManualTSVWriter(out));

  #pragma omp parallel for schedule(dynamic, 1)
  for(gbwt::size_type path_id = 0; path_id < index.metadata.paths(); path_id++)
  {
    const gbwt::PathName& path_name = index.metadata.path(path_id);
    if(path_name.sample == ref_sample) { continue; }
    walks++;
    ManualTSVWriter& writer = writers[omp_get_thread_num()];
    gbwt::vector_type path = index.extract(gbwt::Path::encode(path_id, false));
    size_t length = 0;
    for(auto node : path) { length += graph.get_length(GBWTGraph::node_to_handle(node)); }
    writer.put('W'); writer.newfield();
    if(index.metadata.hasSampleNames()) { writer.write(index.metadata.sample(path_name.sample)); }
    else { writer.write(path_name.sample); }
    writer.newfield();
    writer.write(path_name.phase); writer.newfield();
    if(index.metadata.hasContigNames()) { writer.write(index.metadata.contig(path_name.contig)); }
    else { writer.write(path_name.contig); }
    writer.newfield();
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
    #pragma omp critical
    {
      writer.flush();
    }
  }

  if(parameters.show_progress && walks > 0)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Wrote " << walks << " walks in " << seconds << " seconds" << std::endl;
  }
}

void
write_all_paths(const GBWTGraph& graph, const SegmentCache& cache, std::ostream& out, const GFAExtractionParameters& parameters)
{
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Writing paths" << std::endl;
  }

  const gbwt::GBWT& index = *(graph.index);
  std::vector<ManualTSVWriter> writers(parameters.threads(), ManualTSVWriter(out));

  for(gbwt::size_type seq_id = 0; seq_id < index.sequences(); seq_id += 2)
  {
    ManualTSVWriter& writer = writers[omp_get_thread_num()];
    gbwt::size_type path_id = seq_id / 2;
    gbwt::vector_type path = index.extract(seq_id);
    writer.put('P'); writer.newfield();
    writer.write(path_id); writer.newfield();
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
    #pragma omp critical
    {
      writer.flush();
    }
  }

  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Wrote " << (index.sequences() / 2) << " paths in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

void
gbwt_to_gfa(const GBWTGraph& graph, std::ostream& out, const GFAExtractionParameters& parameters)
{
  bool sufficient_metadata = graph.index->hasMetadata() && graph.index->metadata.hasPathNames();

  // Cache segment names.
  double start = gbwt::readTimer();
  if(parameters.show_progress)
  {
    std::cerr << "Caching segments" << std::endl;
  }
  SegmentCache cache(graph);
  if(parameters.show_progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Cached " << cache.size() << " segments in " << seconds << " seconds" << std::endl;
  }

  // Cache and write the segments using a single writer.
  TSVWriter writer(out);
  writer.put('H'); writer.newfield();
  writer.write(std::string("VN:Z:1.0")); writer.newline();
  write_segments(graph, cache, writer, parameters.show_progress);
  writer.flush();

  // Write the links and paths using multiple threads.
  omp_set_num_threads(parameters.threads());
  write_links(graph, cache, out, parameters);
  if(sufficient_metadata)
  {
    gbwt::size_type ref_sample = graph.index->metadata.sample(REFERENCE_PATH_SAMPLE_NAME);
    write_paths(graph, cache, out, ref_sample, parameters);
    write_walks(graph, cache, out, ref_sample, parameters);
  }
  else { write_all_paths(graph, cache, out, parameters); }
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
