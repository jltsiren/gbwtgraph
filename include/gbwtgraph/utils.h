#ifndef GBWTGRAPH_UTILS_H
#define GBWTGRAPH_UTILS_H

#include <handlegraph/handle_graph.hpp>

#include <iostream>
#include <tuple>

/*
  utils.h: Common utilities and definitions.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Import type definitions from libhandlegraph.
using nid_t = handlegraph::nid_t;
using off_t = handlegraph::off_t;
using pos_t = handlegraph::pos_t;
using handle_t = handlegraph::handle_t;
using edge_t = handlegraph::edge_t;

using HandleGraph = handlegraph::HandleGraph;
using SerializableHandleGraph = handlegraph::SerializableHandleGraph;

//------------------------------------------------------------------------------

// Functions for pos_t manipulation.

inline pos_t
make_pos_t(nid_t id, bool is_rev, off_t off)
{
  return std::make_tuple(id, is_rev, off);
}

inline nid_t
id(const pos_t& pos)
{
  return std::get<0>(pos);
}

inline bool
is_rev(const pos_t& pos)
{
  return std::get<1>(pos);
}

inline off_t
offset(const pos_t& pos)
{
  return std::get<2>(pos);
}

inline nid_t&
get_id(pos_t& pos)
{
  return std::get<0>(pos);
}

inline bool&
get_is_rev(pos_t& pos)
{
  return std::get<1>(pos);
}

inline off_t&
get_offset(pos_t& pos)
{
  return std::get<2>(pos);
}

inline bool
is_empty(const pos_t& pos)
{
  return (id(pos) == 0);
}

// Get a pos_t for the same base in the other orientation.
inline pos_t
reverse_base_pos(const pos_t& pos, size_t node_length)
{
  return make_pos_t(id(pos), !is_rev(pos), (node_length - 1) - offset(pos));
}

inline std::ostream&
operator<<(std::ostream& out, const pos_t& pos)
{
  return out << id(pos) << (is_rev(pos) ? "-" : "+") << offset(pos);
}

//------------------------------------------------------------------------------

// Utility functions.


void reverse_complement_in_place(std::string& seq);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_UTILS_H
