#include <gbwtgraph/internal.h>

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

} // namespace gbwtgraph