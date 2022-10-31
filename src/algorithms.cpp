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

//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------

std::vector<nid_t>
is_nice_and_acyclic(const HandleGraph& graph, const std::vector<nid_t>& component)
{
  std::vector<nid_t> head_nodes;
  if(component.empty()) { return head_nodes; }

  constexpr size_t NOT_SEEN = std::numeric_limits<size_t>::max();
  std::unordered_map<nid_t, std::pair<size_t, bool>> nodes; // (remaining indegree, orientation)
  std::stack<handle_t> active;
  size_t found = 0; // Number of nodes that have become head nodes.

  // Find the head nodes.
  size_t missing_nodes = 0;
  for(nid_t node : component)
  {
    if(!(graph.has_node(node))) { missing_nodes++; continue; }
    handle_t handle = graph.get_handle(node, false);
    size_t indegree = graph.get_degree(handle, true);
    if(indegree == 0)
    {
      nodes[node] = std::make_pair(indegree, false);
      head_nodes.push_back(node);
      active.push(handle);
      found++;
    }
    else
    {
      nodes[node] = std::make_pair(NOT_SEEN, false);
    }
  }

  // Active nodes are the current head nodes. Process the successors, determine their
  // orientations, and decrement their indegrees. If the indegree becomes 0, activate
  // the node.
  bool ok = true;
  while(!(active.empty()))
  {
    handle_t curr = active.top(); active.pop();
    graph.follow_edges(curr, false, [&](const handle_t& next) -> bool
    {
      nid_t next_id = graph.get_id(next);
      bool next_orientation = graph.get_is_reverse(next);
      auto iter = nodes.find(next_id);
      if(iter->second.first == NOT_SEEN) // First visit to the node.
      {
        iter->second.first = graph.get_degree(next, true);
        iter->second.second = next_orientation;
      }
      else if(next_orientation != iter->second.second) // Already visited, wrong orientation.
      {
        ok = false; return false;
      }
      iter->second.first--;
      if(iter->second.first == 0)
      {
        active.push(next);
        found++;
      }
      return true;
    });
    if(!ok) { break; }
  }
  if(found != component.size() - missing_nodes) { ok = false; }

  if(!ok) { head_nodes.clear(); }
  return head_nodes;
}

//------------------------------------------------------------------------------

std::vector<handle_t>
topological_order(const HandleGraph& graph, const std::unordered_set<nid_t>& subgraph)
{
  std::vector<handle_t> result;
  result.reserve(2 * subgraph.size());
  if(subgraph.empty())
  {
    return result;
  }

  std::unordered_map<handle_t, size_t> indegrees;
  indegrees.reserve(2 * subgraph.size());
  std::stack<handle_t> active;

  // Add all handles to the map.
  size_t missing_nodes = 0;
  for(nid_t node : subgraph)
  {
    if(!(graph.has_node(node))) { missing_nodes++; continue; }
    indegrees[graph.get_handle(node, false)] = 0;
    indegrees[graph.get_handle(node, true)] = 0;
  }

  // Determine indegrees and activate head nodes.
  for(auto iter = indegrees.begin(); iter != indegrees.end(); ++iter)
  {
    graph.follow_edges(iter->first, true, [&](const handle_t& next) -> bool
    {
      if(indegrees.find(next) != indegrees.end()) { iter->second++; }
      return true;
    });
    if(iter->second == 0)
    {
      active.push(iter->first);
      result.push_back(iter->first);
    }
  }

  // Follow edges from active nodes and activate the nodes we have reached using
  // all incoming edges.
  while(!(active.empty()))
  {
    handle_t curr = active.top(); active.pop();
    graph.follow_edges(curr, false, [&](const handle_t& next) -> bool
    {
      auto iter = indegrees.find(next);
      if(iter == indegrees.end()) { return true; }
      iter->second--;
      if(iter->second == 0)
      {
        active.push(iter->first);
        result.push_back(iter->first);
      }
      return true;
    });
  }

  if(result.size() != 2 * (subgraph.size() - missing_nodes)) { result.clear(); }
  return result;
}

//------------------------------------------------------------------------------

void
ConstructionJobs::clear()
{
  this->nodes_per_job = {};
  this->node_to_job = {};
  this->components = 0;
}

ConstructionJobs
gbwt_construction_jobs(const HandleGraph& graph, size_t size_bound)
{
  ConstructionJobs jobs;

  std::vector<std::vector<nid_t>> components = weakly_connected_components(graph);
  jobs.components = components.size();

  size_t nodes = graph.get_node_count();
  jobs.node_to_job.reserve(nodes);

  for(size_t i = 0; i < components.size(); i++)
  {
    if(jobs.nodes_per_job.empty() || jobs.nodes_per_job.back() + components[i].size() > size_bound)
    {
      jobs.nodes_per_job.push_back(0);
    }
    jobs.nodes_per_job.back() += components[i].size();
    for(nid_t node_id : components[i]) { jobs.node_to_job[node_id] = jobs.size() - 1; }
  }

  return jobs;
}

std::vector<std::vector<handlegraph::net_handle_t>>
partition_chains(const handlegraph::SnarlDecomposition& snarls, const HandleGraph& graph, const ConstructionJobs& jobs)
{
  std::vector<std::vector<handlegraph::net_handle_t>> result(jobs.size());

  size_t unassigned = 0;
  snarls.for_each_child(snarls.get_root(), [&](const handlegraph::net_handle_t& chain)
  {
    bool assigned = false;
    snarls.for_each_child(chain, [&](const handlegraph::net_handle_t& child) -> bool
    {
      if(snarls.is_node(child))
      {
        nid_t node_id = graph.get_id(snarls.get_handle(child, &graph));
        size_t job_id = jobs(node_id);
        if(job_id < jobs.size())
        {
          result[job_id].push_back(chain);
          assigned = true;
        }
        return false;
      }
      else { return true; }
    });
    if(!assigned) { unassigned++; }
  });

  if(unassigned > 0)
  {
    std::cerr << "partition_chains(): Warning: Could not assign " << unassigned << " chains to jobs" << std::endl;
  }

  return result;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
