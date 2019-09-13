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
serialize(std::ostream& out, const Element& element, bool& ok)
{
  out.write(reinterpret_cast<const char*>(&element), sizeof(element));
  if(out.fail())
  {
    ok = false;
    return 0;
  }
  return sizeof(element);
}

// Load a simple element and return true if successful.
template<typename Element>
bool
load(std::istream& in, Element& element)
{
  in.read(reinterpret_cast<char*>(&element), sizeof(element));
  return (in.gcount() == sizeof(element));
}

// Serialize the size of a container.
template<class Container>
size_t
serialize_size(std::ostream& out, const Container& c, bool &ok)
{
  size_t size = c.size();
  return serialize(out, size, ok);
}

// Resize the container to the serialized size.
template<class Container>
bool
load_size(std::istream& in, Container& c)
{
  size_t size = 0;
  if(!load(in, size)) { return false; }
  c.resize(size);
  return true;
}

// Serialize a vector of simple elements in blocks.
template<typename Element>
size_t
serialize_vector(std::ostream& out, const std::vector<Element>& v, bool& ok)
{
  size_t bytes = 0;

  bytes += serialize_size(out, v, ok);

  // Data in blocks of BLOCK_SIZE elements.
  for(size_t i = 0; i < v.size(); i += BLOCK_SIZE)
  {
    size_t block_size = std::min(v.size() - i, BLOCK_SIZE);
    size_t byte_size = block_size * sizeof(Element);
    out.write(reinterpret_cast<const char*>(v.data() + i), byte_size);
    if(out.fail()) { ok = false; return bytes; }
    bytes += byte_size;
  }

  return bytes;
}

// Load a serialized vector of simple elements.
template<typename Element>
bool
load_vector(std::istream& in, std::vector<Element>& v)
{
  if(!load_size(in, v)) { return false; }

  // Data in blocks of BLOCK_SIZE elements.
  for(size_t i = 0; i < v.size(); i += BLOCK_SIZE)
  {
    size_t block_size = std::min(v.size() - i, BLOCK_SIZE);
    size_t byte_size = block_size * sizeof(Element);
    in.read(reinterpret_cast<char*>(v.data() + i), byte_size);
    if(static_cast<size_t>(in.gcount()) != byte_size) { return false; }
  }

  return true;
}

// Serialize a hash table, replacing pointers with empty values.
// The hash table can be loaded with load_vector().
template<class CellType, class ValueType>
size_t
serialize_hash_table(std::ostream& out, const std::vector<CellType>& hash_table,
                     const sdsl::bit_vector& is_pointer, const ValueType NO_VALUE, bool& ok)
{
  size_t bytes = 0;

  bytes += serialize_size(out, hash_table, ok);

  // Data in blocks of BLOCK_SIZE elements. Replace pointers with NO_VALUE to ensure
  // that the file contents are deterministic.
  for(size_t i = 0; i < hash_table.size(); i += BLOCK_SIZE)
  {
    size_t block_size = std::min(hash_table.size() - i, BLOCK_SIZE);
    size_t byte_size = block_size * sizeof(CellType);
    std::vector<CellType> buffer(hash_table.begin() + i, hash_table.begin() + i + block_size);
    for(size_t j = 0; j < buffer.size(); j++)
    {
      if(is_pointer[i + j]) { buffer[j].second.value = NO_VALUE; }
    }
    out.write(reinterpret_cast<const char*>(buffer.data()), byte_size);
    if(out.fail()) { ok = false; return bytes; }
    bytes += byte_size;
  }

  return bytes;
}

} // namespace io

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_IO_H
