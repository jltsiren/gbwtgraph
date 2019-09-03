#include <gbwtgraph/utils.h>

#include <vector>

namespace gbwtgraph
{

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

} // namespace gbwtgraph

