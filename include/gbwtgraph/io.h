#ifndef GBWTGRAPH_IO_H
#define GBWTGRAPH_IO_H

#include <iostream>
#include <vector>

#include <gbwtgraph/minimizer.h>

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
template<class KeyType>
size_t
serialize_hash_table(std::ostream& out, const KmerIndex<KeyType>& index)
{
  using code_type = typename KmerIndex<KeyType>::code_type;
  using cell_type = typename KmerIndex<KeyType>::cell_type;
  size_t words_in_block = BLOCK_SIZE * index.cell_size();
  size_t cell_size = index.cell_size();

  size_t bytes = 0;
  bytes += serialize_size(out, index.hash_table);

  // Data in blocks of BLOCK_SIZE cells. Replace pointers with empty values
  // to ensure that the serialization is deterministic.
  for(size_t array_offset = 0; array_offset < index.hash_table.size(); array_offset += words_in_block)
  {
    size_t current_block_size = std::min(index.hash_table.size() - array_offset, words_in_block);
    size_t byte_size = current_block_size * sizeof(code_type);
    std::vector<code_type> buffer(index.hash_table.begin() + array_offset, index.hash_table.begin() + array_offset + current_block_size);
    for(size_t buffer_offset = 0; buffer_offset < buffer.size(); buffer_offset += cell_size)
    {
      cell_type cell = KmerIndex<KeyType>::hash_table_cell(buffer, buffer_offset);
      if(cell.first.is_pointer()) { index.clear_value(cell, false); }
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
