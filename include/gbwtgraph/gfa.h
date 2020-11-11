#ifndef GBWTGRAPH_GFA_H
#define GBWTGRAPH_GFA_H

#include <memory>

#include <gbwt/dynamic_gbwt.h>

#include <gbwtgraph/utils.h>

/*
  gfa.h: Tools for building GBWTGraph from GFA.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

struct GFAParsingParameters
{
  gbwt::size_type node_width = gbwt::WORD_BITS;
  gbwt::size_type batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;
  gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;

  bool show_progress = false;

  /*
    path_name_regex is the regex used for parsing path names. Each submatch (part of the
    regex separated by parentheses) is a field. The fields are numbered according to
    preorder traversal from left to right, with 0 corresponding to the entire path name.

    path_name_fields[i] maps field i to a GBWT path name component. Possible values are:

      s or S  sample name
      c or C  contig name
      h       phase/haplotype identifier (0-based)
      H       phase/haplotype identifier (1-based)
      n       running count / fragment identifier (0-based)
      N       running count / fragment identifier (1-based)

    Any other character indicates that the field should not be used. If the string is too
    short, subsequent fields are not used. Each component may occur only once in the string.
  */
  std::string path_name_regex = ".*";
  std::string path_name_fields = "s";
};

//------------------------------------------------------------------------------

/*
  Build GBWT from GFA paths. This completely ignores link lines and makes the following
  assumptions:

    1. Links and paths have no overlaps between segments.
    2. There are no containments.
    3. Segment names are positive integers.

  The construction is done in several passes over a memory-mapped GFA file. The
  function returns the GBWT index and a sequence source for GBWTGraph construction.
*/
std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<SequenceSource>>
gfa_to_gbwt(const std::string& gfa_filename, const GFAParsingParameters& parameters = GFAParsingParameters());

extern const std::string GFA_EXTENSION; // ".gfa"

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_GFA_H
