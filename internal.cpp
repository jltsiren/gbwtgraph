#include <gbwtgraph/internal.h>

#include <algorithm>
#include <functional>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

BufferedHashSet::BufferedHashSet()
{
}

BufferedHashSet::~BufferedHashSet()
{
  // Wait for the insertion thread to finish.
  if(this->worker.joinable()) { this->worker.join(); }
}

void
BufferedHashSet::finish()
{
  // Flush the buffer if necessary.
  this->flush();

  // Wait for the insertion thread to finish.
  if(this->worker.joinable()) { this->worker.join(); }
}

void
insert_keys(std::unordered_set<std::string>& data, std::vector<std::string>& buffer)
{
  for(std::string& str : buffer) { data.insert(std::move(str)); }
}

void
BufferedHashSet::flush()
{
  // Wait for the insertion thread to finish.
  if(this->worker.joinable()) { this->worker.join(); }

  // Swap the input buffer and the internal buffer.
  this->internal_buffer.swap(this->input_buffer);
  this->input_buffer.clear();

  // Launch a new construction thread if necessary.
  if(this->internal_buffer.size() > 0)
  {
    this->worker = std::thread(insert_keys, std::ref(this->data), std::ref(this->internal_buffer));
  }
}

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

} // namespace gbwtgraph
