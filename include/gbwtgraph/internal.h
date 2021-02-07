#ifndef GBWTGRAPH_INTERNAL_H
#define GBWTGRAPH_INTERNAL_H

#include <gbwtgraph/utils.h>

#include <thread>
#include <unordered_set>
#include <vector>

/*
  internal.h: Internal support structures.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  A hash set with string keys. New keys are buffered and inserted in
  a background thread. The overall logic is similar to `GBWTBuilder`.
*/
struct BufferedHashSet
{
  BufferedHashSet();
  ~BufferedHashSet();

  void insert(view_type view)
  {
    if(this->input_buffer.size() >= BUFFER_SIZE) { this->flush(); }
    this->input_buffer.emplace_back(view.first, view.second);
  }

  // Insert all buffered keys into the hash table.
  void finish();

  // Buffer up to this many keys.
  constexpr static size_t BUFFER_SIZE = 4 * 1048576;

  std::unordered_set<std::string> data;

  // New keys are buffered here.
  std::vector<std::string> input_buffer;

  // The insertion thread uses this data.
  std::vector<std::string> internal_buffer;
  std::thread              worker;

  BufferedHashSet(const BufferedHashSet&) = delete;
  BufferedHashSet& operator=(const BufferedHashSet&) = delete;

private:
  void flush();
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_INTERNAL_H
