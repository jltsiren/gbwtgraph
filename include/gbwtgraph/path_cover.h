#ifndef GBWTGRAPH_PATH_COVER_H
#define GBWTGRAPH_PATH_COVER_H

#include <gbwt/dynamic_gbwt.h>

#include <gbwtgraph/gbwtgraph.h>

/*
  path_cover.h: Build GBWT from a path cover of a HandleGraph.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

constexpr size_t PATH_COVER_DEFAULT_N = 16;
constexpr size_t PATH_COVER_DEFAULT_K = 4;
constexpr size_t PATH_COVER_MIN_K     = 2;

//------------------------------------------------------------------------------

/*
  Find a path cover of the graph with n paths per component and return a GBWT of the paths.
  The path cover is built greedily. Each time we extend a path, we choose the extension
  where the coverage of the k >= 2 node window is the lowest. Note that this is a maximum
  coverage algorithm that tries to maximize the number of windows covered by a fixed number
  of paths.

  This algorithm has been inspired by:

    Ghaffaari and Marschall: Fully-sensitive seed finding in sequence graphs using a
    hybrid index. Bioinformatics, 2019.

  Because the graph here may have cycles and orientation flips, some changes were necessary:

    - Instead of starting from a head node, we start from an arbitrary node with minimal
      coverage and extend in both directions.

    - We stop when the length of the path exceeds the size of the component.

    - The length of the window is in nodes instead of base pairs. We expect a sparse graph,
      where the nodes between variants are long.

    - When the path is shorter than k nodes, we consider the coverage of individual nodes
      instead of windows.

    - When determining window coverage, we consider the window equivalent to its reverse
      complement.
*/
gbwt::GBWT path_cover_gbwt(const HandleGraph& graph,
                           size_t n = PATH_COVER_DEFAULT_N, size_t k = PATH_COVER_DEFAULT_K,
                           gbwt::size_type batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE,
                           gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL,
                           bool show_progress = false);

//------------------------------------------------------------------------------

/*
  As above, but build the path cover using local haplotypes. Each window must be consistent
  with the haplotypes in the GBWT index. Instead of choosing the extension with the lowest
  coverage, this algorithm maximizes the D'Hondt ratio:

    true_coverage / (selected_coverage + 1).
*/

gbwt::GBWT local_haplotypes(const GBWTGraph& graph,
                            size_t n = PATH_COVER_DEFAULT_N, size_t k = PATH_COVER_DEFAULT_K,
                            gbwt::size_type batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE,
                            gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL,
                            bool show_progress = false);

//------------------------------------------------------------------------------

/*
  Return the weakly connected components in the graph. The components are sorted by
  the minimum node id, and the node ids in each component are also sorted.
*/
std::vector<std::vector<nid_t>> weakly_connected_components(const HandleGraph& graph);

/*
  Determine whether the given component is acyclic in a nice way. In such graphs,
  when we start from nodes with indegree 0 in forward orientation, we reach each node
  in a single orientation and find no cycles. Return the head nodes when the component
  passes the tests or an empty vector otherwise.
*/
std::vector<nid_t> is_nice_and_acyclic(const HandleGraph& graph, const std::vector<nid_t>& component);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_PATH_COVER_H
