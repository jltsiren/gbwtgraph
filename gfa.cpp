#include <gbwtgraph/gfa.h>

#include <algorithm>
#include <functional>
#include <map>
#include <regex>
#include <string>
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
  int         fd;
  size_t      file_size;
  char*       ptr;

  struct field_type
  {
    const char* begin;
    const char* end;
    bool        has_next;

    field_type(const char* begin, const char* end, bool has_next) : begin(begin), end(end), has_next(has_next) {}

    size_t size() const { return this->end - this->begin; }
    bool empty() const { return (this->size() == 0); }

    std::string str() const { return std::string(this->begin, this->end); }
    std::string oriented_str() const { return std::string(this->begin, this->end - 1); }
    char front() const { return *(this->begin); }
    char back() const { return *(this->end - 1); }
    bool is_reverse() const { return (this->back() == '-'); }
  };

  GFAFile(const std::string& filename) :
    fd(-1), file_size(0), ptr(nullptr)
  {
    this->fd = ::open(filename.c_str(), O_RDONLY);
    if(this->fd < 0)
    {
      std::cerr << "GFAFile::GFAFile(): Cannot open file " << filename << std::endl;
      return;
    }

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
  }

  ~GFAFile()
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

  bool ok() const { return (this->fd >= 0 && this->file_size > 0 && this->ptr != nullptr); }
  size_t size() const { return this->file_size; }

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
    while(limit != this->end() && *limit != '\n' && *limit != '\t') { ++limit; }
    return field_type(iter, limit, (limit != this->end() && *limit == '\t'));
  }

  // Return the comma-separated subfield starting at the iterator.
  field_type get_subfield(const char* iter) const
  {
    const char* limit = iter;
    while(limit != this->end() && *limit != '\n' && *limit != '\t' && *limit != ',') { ++limit; }
    return field_type(iter, limit, (limit != this->end() && *limit == ','));
  }

  bool is_segment(const char* iter) const { return *iter == 'S'; }
  bool is_path(const char* iter) const { return *iter == 'P'; }

  /*
    Iterate over the file, calling segment() for all segments. Stops early if segment()
    returns false.
  */
  void for_each_segment(const std::function<bool(const std::string& name, const std::string& sequence)>& segment) const
  {
    if(!(this->ok())) { return; }

    // Iterate over the lines in the file.
    const char* iter = this->begin();
    while(iter != this->end())
    {
      if(this->is_segment(iter))
      {
        // Skip the record type field.
        field_type field = this->get_field(iter);
        if(!(field.has_next))
        {
          iter = this->next_line(field.end);
          continue;
        }

        // Segment name field.
        field = this->get_field(field.end + 1);
        if(!(field.has_next) || field.empty())
        {
          iter = this->next_line(field.end);
          continue;
        }
        std::string name = field.str();

        // Sequence field.
        field = this->get_field(field.end + 1);
        std::string sequence = field.str();
        if(!segment(name, sequence)) { return; }
        iter = this->next_line(field.end);
      }
      else
      {
        iter = this->next_line(iter);
      }
    }
  }

  /*
    Iterate over the file, calling path() for each path, path_segment() for
    each path segment, and finish_path() after parsing each path. Stops early
    if any call returns false.
  */
  void for_each_path(const std::function<bool(const std::string& name)>& path,
                     const std::function<bool(const std::string& name, bool is_reverse)>& path_segment,
                     const std::function<bool(const std::string& name)>& finish_path)
  {
    if(!(this->ok())) { return; }

    // Iterate over the lines in the file.
    const char* iter = this->begin();
    while(iter != this->end())
    {
      if(this->is_path(iter))
      {
        // Skip the record type field.
        field_type field = this->get_field(iter);
        if(!(field.has_next))
        {
          iter = this->next_line(field.end);
          continue;
        }

        // Path name field.
        field = this->get_field(field.end + 1);
        if(!(field.has_next) || field.empty())
        {
          iter = this->next_line(field.end);
          continue;
        }
        std::string path_name = field.str();
        if(!path(path_name)) { return; }

        // Segment names field.
        do
        {
          field = this->get_subfield(field.end + 1);
          if(field.size() >= 2)
          {
            std::string segment_name = field.oriented_str();
            if(!path_segment(segment_name, field.is_reverse())) { return; }
          }
        }
        while(field.has_next);

        if(!finish_path(path_name)) { return; }
        iter = this->next_line(field.end);
      }
      else
      {
        iter = this->next_line(iter);
      }
    }
  }
};

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

  explicit PathNameParser(const GFAParsingParameters& parameters) :
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

  bool parse(const std::string& name)
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

  bool empty() const { return this->path_names.empty(); }

  void setMetadata(gbwt::Metadata& metadata)
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
  if(parameters.show_progress)
  {
    std::cerr << "Opening GFA file " << gfa_filename << std::endl;
  }
  GFAFile gfa_file(gfa_filename);
  if(!(gfa_file.ok()))
  {
    std::cerr << "gfa_to_gbwt(): Cannot read file " << gfa_filename << std::endl;
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }

  // First pass: Segments.
  if(parameters.show_progress)
  {
    std::cerr << "Parsing segments" << std::endl;
  }
  SequenceSource* source = new SequenceSource();
  bool failed = false;
  gfa_file.for_each_segment([&](const std::string& name, const std::string& sequence) -> bool
  {
    nid_t node_id = 0;
    try { node_id = std::stoul(name); }
    catch(const std::invalid_argument&) {}
    if(node_id == 0)
    {
      std::cerr << "gfa_to_gbwt(): Invalid segment name " << name << std::endl;
      failed = true; return false;
    }
    if(source->has_node(node_id))
    {
      std::cerr << "gfa_to_gbwt(): Duplicate segment name " << name << std::endl;
      failed = true; return false;
    }
    source->add_node(node_id, sequence);
    return true;
  });
  if(failed)
  {
    std::cerr << "gfa_to_gbwt(): Invalid GFA file " << gfa_filename << std::endl;
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }
  if(parameters.show_progress)
  {
    std::cerr << "Parsed " << source->get_node_count() << " segments" << std::endl;
  }

  // GBWT construction.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type batch_size = std::min(static_cast<gbwt::size_type>(gfa_file.size()), parameters.batch_size);
  gbwt::GBWTBuilder builder(parameters.node_width, batch_size, parameters.sample_interval);
  builder.index.addMetadata();

  // Second pass: Metadata and path validation.
  if(parameters.show_progress)
  {
    std::cerr << "Parsing metadata" << std::endl;
  }
  gfa_file.for_each_path([&](const std::string& name) -> bool
  {
    if(!(parser.parse(name))) { failed = true; return false; }
    return true;
  }, [&](const std::string& name, bool) -> bool
  {
    nid_t node_id = 0;
    try { node_id = std::stoul(name); }
    catch(const std::invalid_argument&) {}
    if(!(source->has_node(node_id)))
    {
      std::cerr << "gfa_to_gbwt(): Invalid path segment name " << name << std::endl;
      failed = true; return false;
    }
    return true;
  }, [&](const std::string&) -> bool { return true; });
  if(failed)
  {
    std::cerr << "gfa_to_gbwt(): Invalid GFA file " << gfa_filename << std::endl;
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }
  if(parser.empty())
  {
    std::cerr << "gfa_to_gbwt(): No paths in GFA file " << gfa_filename << std::endl;
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }
  parser.setMetadata(builder.index.metadata);
  parser.clear();
  if(parameters.show_progress)
  {
    std::cerr << "Metadata: "; gbwt::operator<<(std::cerr, builder.index.metadata) << std::endl;
  }

  // Third pass: Paths.
  if(parameters.show_progress)
  {
    std::cerr << "Indexing paths" << std::endl;
  }
  gbwt::vector_type current_path;
  size_t path_count = 0;
  gfa_file.for_each_path([&](const std::string&) -> bool
  {
    return true;
  }, [&](const std::string& name, bool is_reverse) -> bool
  {
    nid_t node_id = std::stoul(name);
    current_path.push_back(gbwt::Node::encode(node_id, is_reverse));
    return true;
  }, [&](const std::string& name) -> bool
  {
    if(current_path.empty())
    {
      std::cerr << "gfa_to_gbwt(): Path " << name << " is empty" << std::endl;
      failed = true; return false;
    }
    builder.insert(current_path, true); current_path.clear(); path_count++;
    return true;
  });
  if(failed)
  {
    std::cerr << "gfa_to_gbwt(): Invalid GFA file " << gfa_filename << std::endl;
    return std::make_pair(std::unique_ptr<gbwt::GBWT>(nullptr), std::unique_ptr<SequenceSource>(nullptr));
  }
  builder.finish();
  if(parameters.show_progress)
  {
    std::cerr << "Indexed " << path_count << " paths" << std::endl;
  }

  return std::make_pair(std::unique_ptr<gbwt::GBWT>(new gbwt::GBWT(builder.index)), std::unique_ptr<SequenceSource>(source));
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
