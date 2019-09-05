#include <gbwtgraph/gfa.h>

#include <algorithm>
#include <functional>
#include <string>

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

namespace gbwtgraph
{

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

    this->ptr = static_cast<char*>(::mmap(nullptr, file_size, PROT_READ, MAP_FILE | MAP_SHARED, this->fd, 0));
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
    Iterate over the file, calling segment() for all segments and path_segment() for all
    path segments, assuming that the relevant names are non-empty. Argument 'first'
    indicates that this is the first segment of the path.
  */
  void for_each_line(const std::function<void(const std::string& name, const std::string& sequence)>& segment,
                     const std::function<void(const std::string& path_name, const std::string& segment_name, bool is_rev, bool first)>& path_segment)
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
        segment(name, sequence);
        iter = this->next_line(field.end);
      }
      else if(this->is_path(iter))
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

        // Segment names field.
        bool first_segment = true;
        do
        {
          field = this->get_subfield(field.end + 1);
          if(field.size() >= 2)
          {
            std::string segment_name = field.oriented_str();
            path_segment(path_name, segment_name, field.is_reverse(), first_segment);
            first_segment = false;
          }
        }
        while(field.has_next);
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

std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<SequenceSource>>
gfa_to_gbwt(const std::string& gfa_filename, gbwt::size_type node_width, gbwt::size_type batch_size, gbwt::size_type sample_interval)
{
  GFAFile gfa_file(gfa_filename);
  if(!(gfa_file.ok()))
  {
    std::cerr << "gfa_to_gbwt(): Cannot read file " << gfa_filename << std::endl;
    return std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<SequenceSource>>(nullptr, nullptr);
  }

  // Index the paths.
  gbwt::GBWTBuilder builder(node_width, batch_size, sample_interval);
  SequenceSource* source = new SequenceSource();
  gbwt::vector_type current_path;
  auto finish_path = [&]()
  {
    if(!(current_path.empty()))
    {
      builder.insert(current_path, true);
      current_path.clear();
    }
  };
  gfa_file.for_each_line([&](const std::string& name, const std::string& sequence)
  {
    source->add_node(std::stol(name), sequence);
  }, [&](const std::string&, const std::string& segment, bool is_rev, bool first)
  {
    if(first) { finish_path(); }
    current_path.push_back(gbwt::Node::encode(std::stol(segment), is_rev));
  });
  finish_path();
  builder.finish();

  // FIXME Metadata

  return std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<SequenceSource>>(new gbwt::GBWT(builder.index), source);
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
