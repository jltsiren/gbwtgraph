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

typedef gbwtgraph::view_type view_type;
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

// Build a GBWT with 6 paths, 3 duplicates of each.
// TODO: what's the "ref" here?
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

// Build a GBWT that contains 3 paths, two of which are named, reference paths in the metadata.
inline gbwt::GBWT
build_gbwt_index_with_named_paths()
{
    // Start with the no-metadata GBWT we usually use
    gbwt::GBWT built = build_gbwt_index();
    
    built.addMetadata();
    
    // Name the set of samples, including our special ref one
    built.metadata.setSamples({gbwtgraph::REFERENCE_PATH_SAMPLE_NAME, "Jouni Sirén"});
    
    // Name the set of contigs we are over.
    built.metadata.setContigs({"chr1", "chr2"});
    
    // We have a ref path on the first contig, and a ref path on the second
    // contig, and a sample path on the first contig.
    gbwt::PathName ref1;
    ref1.sample = 0;
    ref1.contig = 0;
    ref1.phase = 0;
    ref1.count = 0;
    built.metadata.addPath(ref1);
    gbwt::PathName ref2;
    ref2.sample = 0;
    ref2.contig = 1;
    ref2.phase = 0;
    ref2.count = 0;
    built.metadata.addPath(ref2);
    gbwt::PathName sample1;
    sample1.sample = 1;
    sample1.contig = 0;
    sample1.phase = 0;
    sample1.count = 0;
    built.metadata.addPath(sample1);
    
    // Record that we have 3 total haplotypes in the GBWT.
    built.metadata.setHaplotypes(3);
    
    return built;
}

// Dump a GBWT and its metadata.
inline void
dump_gbwt(const gbwt::GBWT& built) {
    // Check to make sure threads and metadata match.
    if (built.hasMetadata()) {
        for (size_t path_num = 0; path_num < built.metadata.paths(); path_num++) {
            auto extracted = built.extract(GBWTGraph::path_to_sequence(path_num, false));
            std::cerr << "GBWT stored forward path " << path_num << std::endl;
            auto path_name = built.metadata.path(path_num);
            std::cerr << "\tSample " << path_name.sample
                << " Contig " << path_name.contig
                << " Phase " << path_name.phase 
                << " Count " << path_name.count << std::endl;
            for (auto& gbwt_node : extracted) {
                std::cerr << "\t" << gbwt::Node::id(gbwt_node) << " " << gbwt::Node::is_reverse(gbwt_node) << std::endl;
            }
        }
    }
}

inline void
build_source(gbwtgraph::SequenceSource& source, bool with_translation = false)
{
  if(with_translation)
  {
    std::string seq = "GATGGGTACAA";
    source.translate_segment("s1", view_type(seq.data() + 0, 1), 3);
    source.translate_segment("s2", view_type(seq.data() + 1, 1), 3);
    source.translate_segment("s3", view_type(seq.data() + 2, 1), 3);
    source.translate_segment("s4", view_type(seq.data() + 3, 4), 3);
    source.translate_segment("s5", view_type(seq.data() + 7, 1), 3);
    source.translate_segment("s6", view_type(seq.data() + 8, 1), 3);
    source.translate_segment("s7", view_type(seq.data() + 9, 1), 3);
    source.translate_segment("s8", view_type(seq.data() + 10, 1), 3);
  }
  else
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
