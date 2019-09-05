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

/*
  Build GBWT from GFA paths. This completely ignores link lines and makes the following
  assumptions:

    1. Links and paths have no overlaps between segments.
    2. There are no containments.
    3. Segment names are integers.

  The construction is done in a single memory-mapped pass over the GFA file. The
  function returns the GBWT index and a sequence source for GBWTGraph construction.
*/
std::pair<std::unique_ptr<gbwt::GBWT>, std::unique_ptr<SequenceSource>>
gfa_to_gbwt(const std::string& gfa_filename,
            gbwt::size_type node_width = gbwt::WORD_BITS,
            gbwt::size_type batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE,
            gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_GFA_H
