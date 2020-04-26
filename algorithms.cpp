#include <gbwtgraph/algorithms.h>

#include <limits>
#include <stack>
#include <unordered_map>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// A quick-and-dirty union-find data structure using path splitting and
// union by rank.
struct DisjointSets
{
  std::vector<size_t> parent;
  std::vector<std::uint8_t> rank; // Rank is at most ~log(size).
  nid_t offset;                   // Node i is stored at [i - offset].

  DisjointSets(size_t n, nid_t offset) :
    parent(n, 0), rank(n, 0), offset(offset)
  {
    for(size_t i = 0; i < this->size(); i++) { this->parent[i] = i; }
  }

  size_t size() const { return this->parent.size(); }

  void clear()
  {
    this->parent = std::vector<size_t>();
    this->rank = std::vector<std::uint8_t>();
  }

  size_t find(nid_t node)
  {
    size_t element = node - this->offset;
    while(this->parent[element] != element)
    {
      size_t next = this->parent[element];
      this->parent[element] = this->parent[next];
      element = next;
    }
    return element;
  }

  void set_union(nid_t node_a, nid_t node_b)
  {
    size_t a = this->find(node_a), b = this->find(node_b);
    if(a == b) { return; }
    if(this->rank[a] < this->rank[b]) { std::swap(a, b); }
    this->parent[b] = a;
    if(this->rank[b] == this->rank[a]) { this->rank[a]++; }
  }

  std::vector<std::vector<nid_t>> sets(const std::function<bool(nid_t)>& include_node)
  {
    std::vector<std::vector<nid_t>> result;
    std::unordered_map<size_t, size_t> root_to_set;
    for(nid_t node = this->offset; node < static_cast<nid_t>(this->offset + this->size()); node++)
    {
      if(!include_node(node)) { continue; }
      size_t root = this->find(node);
      if(root_to_set.find(root) == root_to_set.end())
      {
        root_to_set[root] = result.size();
        result.emplace_back();
      }
      result[root_to_set[root]].push_back(node);
    }
    return result;
  }
};

/*
  Return the weakly connected components in the graph.
*/
std::vector<std::vector<nid_t>>
weakly_connected_components(const HandleGraph& graph)
{
  nid_t min_id = graph.min_node_id(), max_id = graph.max_node_id();

  sdsl::bit_vector found(max_id + 1 - min_id, 0);
  DisjointSets components(max_id + 1 - min_id, min_id);
  graph.for_each_handle([&](const handle_t& handle)
  {
    nid_t start_id = graph.get_id(handle);
    if(found[start_id - min_id]) { return; }
    std::stack<handle_t> handles;
    handles.push(handle);
    while(!(handles.empty()))
    {
      handle_t h = handles.top(); handles.pop();
      nid_t id = graph.get_id(h);
      if(found[id - min_id]) { continue; }
      found[id - min_id] = true;
      auto handle_edge = [&](const handle_t& next)
      {
        components.set_union(id, graph.get_id(next));
        handles.push(next);
      };
      graph.follow_edges(h, false, handle_edge);
      graph.follow_edges(h, true, handle_edge);
    }
  }, false);

  return components.sets([&](nid_t node) -> bool { return graph.has_node(node); });
}

// FIXME non-destructive version should be faster
std::vector<nid_t>
is_nice_and_acyclic(const HandleGraph& graph, const std::vector<nid_t>& component)
{
  std::vector<nid_t> head_nodes;
  if(component.empty()) { return head_nodes; }

  constexpr size_t NOT_SEEN = std::numeric_limits<size_t>::max();
  std::unordered_map<nid_t, bool> orientations;
  std::unordered_map<nid_t, size_t> indegrees; // Nodes we have not visited from every incoming edge.
  std::stack<handle_t> active;

  // Find the head nodes.
  for(nid_t node : component)
  {
    handle_t handle = graph.get_handle(node, false);
    size_t indegree = graph.get_degree(handle, true);
    if(indegree == 0)
    {
      orientations[node] = false;
      head_nodes.push_back(node);
      active.push(handle);
    }
    else
    {
      indegrees[node] = NOT_SEEN; // Set the degree when we decide on the orientation.
    }
  }

  // Active nodes are the current head nodes. Process the successors, determine their
  // orientations, and decrement their indegrees. If the indegree becomes 0, activate
  // the node and remove it from the set of unfinished nodes.
  bool ok = true;
  while(!(active.empty()))
  {
    handle_t curr = active.top(); active.pop();
    graph.follow_edges(curr, false, [&](const handle_t& next)
    {
      nid_t next_id = graph.get_id(next);
      nid_t next_orientation = graph.get_is_reverse(next);
      auto orientation_iter = orientations.find(next_id);
      auto indegree_iter = indegrees.find(next_id);
      if(orientation_iter == orientations.end())
      {
        orientations[next_id] = next_orientation;
        indegree_iter->second = graph.get_degree(next, true);
      }
      else if(next_orientation != orientation_iter->second)
      {
        ok = false; return;
      }
      indegree_iter->second--;
      if(indegree_iter->second == 0)
      {
        active.push(next);
        indegrees.erase(indegree_iter);
      }
    });
    if(!ok) { break; }
  }
  if(!(indegrees.empty())) { ok = false; }

  if(!ok) { head_nodes.clear(); }
  return head_nodes;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
