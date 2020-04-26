#ifndef GBWTGRAPH_ALGORITHMS_H
#define GBWTGRAPH_ALGORITHMS_H

#include <gbwtgraph/gbwtgraph.h>

/*
  algorithms.h: GBWTGraph-specific graph algorithms. These should be faster than
  generic HandleGraph algorithms, as they can take advantage of the CachedGBWT
  layer.
*/

namespace gbwtgraph
{

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

#endif // GBWTGRAPH_ALGORITHMS_H
