#ifndef GBWTGRAPH_TESTS_SHARED_H
#define GBWTGRAPH_TESTS_SHARED_H

#include <vector>

#include <gbwt/dynamic_gbwt.h>

#include <gbwtgraph/minimizer.h>

/*
  shared.h: Utility functions and data definitions shared between the tests.
*/

namespace
{

//------------------------------------------------------------------------------

typedef std::pair<gbwtgraph::nid_t, std::string> node_type;
typedef std::pair<std::string, std::pair<gbwtgraph::nid_t, gbwtgraph::nid_t>> translation_type;

//------------------------------------------------------------------------------

inline gbwtgraph::view_type
get_view(const std::string& source)
{
  return gbwtgraph::view_type(source.data(), source.length());
}

//------------------------------------------------------------------------------

gbwt::vector_type alt_path
{
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(2, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(8, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
};

gbwt::vector_type short_path
{
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(7, false)),
  static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
};

// Build a GBWT with three paths including a duplicate.
inline gbwt::GBWT
build_gbwt_index()
{
  std::vector<gbwt::vector_type> paths
  {
    short_path, alt_path, short_path
  };

  // Determine node width in bits.
  gbwt::size_type node_width = 1, total_length = 0;
  for(auto& path : paths)
  {
    for(auto node : path)
    {
      node_width = std::max(node_width, gbwt::size_type(sdsl::bits::length(gbwt::Node::encode(node, true))));
    }
    total_length += 2 * (path.size() + 1);
  }

  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::GBWTBuilder builder(node_width, total_length);
  for(auto& path : paths) { builder.insert(path, true); }
  builder.finish();

  return gbwt::GBWT(builder.index);
}

inline gbwt::GBWT
build_gbwt_index_with_ref()
{
  std::vector<gbwt::vector_type> paths
  {
    short_path, alt_path, alt_path,
    short_path, alt_path, short_path
  };

  // Determine node width in bits.
  gbwt::size_type node_width = 1, total_length = 0;
  for(auto& path : paths)
  {
    for(auto node : path)
    {
      node_width = std::max(node_width, gbwt::size_type(sdsl::bits::length(gbwt::Node::encode(node, true))));
    }
    total_length += 2 * (path.size() + 1);
  }

  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::GBWTBuilder builder(node_width, total_length);
  for(auto& path : paths) { builder.insert(path, true); }
  builder.finish();

  return gbwt::GBWT(builder.index);
}

inline void
build_source(gbwtgraph::SequenceSource& source)
{
  source.add_node(1, "G");
  source.add_node(2, "A");
  source.add_node(3, "T");
  source.add_node(4, "GGG");
  source.add_node(5, "T");
  source.add_node(6, "A");
  source.add_node(7, "C");
  source.add_node(8, "A");
  source.add_node(9, "A");
}

//------------------------------------------------------------------------------

template<class KeyType>
typename gbwtgraph::MinimizerIndex<KeyType>::minimizer_type
get_minimizer(KeyType key, typename gbwtgraph::MinimizerIndex<KeyType>::offset_type offset = 0, bool orientation = false)
{
  return { key, key.hash(), offset, orientation };
}

template<class KeyType>
typename gbwtgraph::MinimizerIndex<KeyType>::minimizer_type
get_minimizer(std::string key, typename gbwtgraph::MinimizerIndex<KeyType>::offset_type offset = 0, bool orientation = false)
{
  return get_minimizer(KeyType::encode(key), offset, orientation);
}

//------------------------------------------------------------------------------

} // anonymous namespace

#endif // GBWTGRAPH_TESTS_SHARED_H
