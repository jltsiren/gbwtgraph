#ifndef GBWTGRAPH_IO_H
#define GBWTGRAPH_IO_H

#include <iostream>
#include <vector>

/*
  io.h: Internal I/O functions.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

namespace io {

constexpr static size_t BLOCK_SIZE = 4 * 1024 * 1024; // Elements per block.

// Serialize a simple element.
template<typename Element>
size_t
serialize(std::ostream& out, const Element& element)
{
  out.write(reinterpret_cast<const char*>(&element), sizeof(element));
  return sizeof(element);
}

// Load a simple element.
template<typename Element>
void
load(std::istream& in, Element& element)
{
  in.read(reinterpret_cast<char*>(&element), sizeof(element));
}

// Serialize the size of a container.
template<class Container>
size_t
serialize_size(std::ostream& out, const Container& c)
{
  size_t size = c.size();
  return serialize(out, size);
}

// Resize the container to the serialized size.
template<class Container>
void
load_size(std::istream& in, Container& c)
{
  size_t size = 0;
  load(in, size);
  c.resize(size);
}

// Serialize a vector of simple elements in blocks.
template<typename Element>
size_t
serialize_vector(std::ostream& out, const std::vector<Element>& v)
{
  size_t bytes = 0;

  bytes += serialize_size(out, v);

  // Data in blocks of BLOCK_SIZE elements.
  for(size_t i = 0; i < v.size(); i += BLOCK_SIZE)
  {
    size_t block_size = std::min(v.size() - i, BLOCK_SIZE);
    size_t byte_size = block_size * sizeof(Element);
    out.write(reinterpret_cast<const char*>(v.data() + i), byte_size);
    bytes += byte_size;
  }

  return bytes;
}

// Load a serialized vector of simple elements.
template<typename Element>
void
load_vector(std::istream& in, std::vector<Element>& v)
{
  load_size(in, v);

  // Data in blocks of BLOCK_SIZE elements.
  for(size_t i = 0; i < v.size(); i += BLOCK_SIZE)
  {
    size_t block_size = std::min(v.size() - i, BLOCK_SIZE);
    size_t byte_size = block_size * sizeof(Element);
    in.read(reinterpret_cast<char*>(v.data() + i), byte_size);
  }
}

// Serialize a hash table, replacing pointers with empty values.
// The hash table can be loaded with load_vector().
template<class CellType, class ValueType>
size_t
serialize_hash_table(std::ostream& out, const std::vector<CellType>& hash_table, const ValueType NO_VALUE)
{
  size_t bytes = 0;

  bytes += serialize_size(out, hash_table);

  // Data in blocks of BLOCK_SIZE elements. Replace pointers with NO_VALUE to ensure
  // that the file contents are deterministic.
  for(size_t i = 0; i < hash_table.size(); i += BLOCK_SIZE)
  {
    size_t block_size = std::min(hash_table.size() - i, BLOCK_SIZE);
    size_t byte_size = block_size * sizeof(CellType);
    std::vector<CellType> buffer(hash_table.begin() + i, hash_table.begin() + i + block_size);
    for(size_t j = 0; j < buffer.size(); j++)
    {
      if(buffer[j].first.is_pointer()) { buffer[j].second.value = NO_VALUE; }
    }
    out.write(reinterpret_cast<const char*>(buffer.data()), byte_size);
    bytes += byte_size;
  }

  return bytes;
}

} // namespace io

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_IO_H
