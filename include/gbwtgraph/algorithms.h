#ifndef GBWTGRAPH_ALGORITHMS_H
#define GBWTGRAPH_ALGORITHMS_H

#include <unordered_set>

#include <gbwtgraph/gbwtgraph.h>

/*
  algorithms.h: Various graph algorithms.
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

  Ignores node ids that are not present in the graph.
*/
std::vector<nid_t> is_nice_and_acyclic(const HandleGraph& graph, const std::vector<nid_t>& component);

//------------------------------------------------------------------------------

/*
  Return a topological order of handles in the subgraph induced by the given node ids,
  or an empty vector if no such order exists.

  NOTE: Node ids that do not exist in the graph are ignored. In particular, some
  nodes of the original graph may be missing from the corresponding GBWTGraph if no
  path passes through them.

  If the subgraph is small, it may be a good idea to use CachedGBWTGraph instead of
  GBWTGraph.
*/
std::vector<handle_t> topological_order(const HandleGraph& graph, const std::unordered_set<nid_t>& subgraph);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_ALGORITHMS_H
