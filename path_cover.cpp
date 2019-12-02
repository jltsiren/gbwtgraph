#include <gbwtgraph/path_cover.h>

#include <algorithm>
#include <cassert>
#include <deque>
#include <limits>
#include <map>
#include <stack>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

std::vector<std::vector<nid_t>>
weakly_connected_components(const HandleGraph& graph)
{
  nid_t min_id = graph.min_node_id(), max_id = graph.max_node_id();

  std::vector<std::vector<nid_t>> components;
  sdsl::bit_vector found(max_id + 1 - min_id, 0);
  graph.for_each_handle([&](const handle_t& handle)
  {
    nid_t start_id = graph.get_id(handle);
    if(found[start_id - min_id]) { return; }
    components.emplace_back();
    std::stack<handle_t> handles;
    handles.push(handle);
    while(!(handles.empty()))
    {
      handle_t h = handles.top(); handles.pop();
      nid_t id = graph.get_id(h);
      if(found[id - min_id]) { continue; }
      found[id - min_id] = true;
      components.back().push_back(id);
      auto push = [&handles](const handle_t& next) { handles.push(next); };
      graph.follow_edges(h, false, push);
      graph.follow_edges(h, true, push);
    }
  }, false);

  return components;
}

// Similar to std::lower_bound().
template<class NodeCoverage>
size_t
find_first(std::vector<NodeCoverage>& array, nid_t id)
{
  size_t first = 0, mid = 0;
  size_t count = array.size();
  while(count > 0)
  {
    size_t step = count / 2;
    mid = first + step;
    if(array[mid].first < id) { first = mid + 1; count -= step + 1; }
    else { count = step; }
  }
  return first;
}

std::vector<handle_t>
reverse_complement(const HandleGraph& graph, std::vector<handle_t>& forward)
{
  std::vector<handle_t> result = forward;
  std::reverse(result.begin(), result.end());
  for(handle_t& handle : result) { handle = graph.flip(handle); }
  return result;
}

std::vector<handle_t>
forward_window(const HandleGraph& graph, const std::deque<handle_t>& path, const handle_t& successor, size_t k)
{
  std::vector<handle_t> forward;
  forward.reserve(k);
  forward.insert(forward.end(), path.end() - (k - 1), path.end());
  forward.push_back(successor);

  std::vector<handle_t> reverse = reverse_complement(graph, forward);
  return (forward < reverse ? forward : reverse);
}

std::vector<handle_t>
backward_window(const HandleGraph& graph, const std::deque<handle_t>& path, const handle_t& predecessor, size_t k)
{
  std::vector<handle_t> forward;
  forward.reserve(k);
  forward.push_back(predecessor);
  forward.insert(forward.end(), path.begin(), path.begin() + (k - 1));

  std::vector<handle_t> reverse = reverse_complement(graph, forward);
  return (forward < reverse ? forward : reverse);
}

template<class Coverage>
struct BestCoverage
{
  typedef typename Coverage::coverage_t coverage_t;

  coverage_t coverage;
  handle_t   handle;

  BestCoverage() : coverage(Coverage::worst_coverage()), handle() {}

  void update(const coverage_t& new_coverage, const handle_t& new_handle)
  {
    if(new_coverage < this->coverage)
    {
      this->coverage = new_coverage;
      this->handle = new_handle;
    }
  }
};

//------------------------------------------------------------------------------

/*
  The best candidate is the one with the lowest coverage so far.
*/

struct SimpleCoverage
{
  typedef HandleGraph graph_t;

  typedef size_t coverage_t;
  typedef std::pair<nid_t, coverage_t> node_coverage_t;

  static std::vector<node_coverage_t> init_node_coverage(const graph_t&, const std::vector<nid_t>& component)
  {
    std::vector<node_coverage_t> node_coverage;
    node_coverage.reserve(component.size());
    for(nid_t id : component) { node_coverage.emplace_back(id, static_cast<coverage_t>(0)); }
    return node_coverage;
  }

  static bool extend_forward(const graph_t& graph, std::deque<handle_t>& path, size_t k, std::vector<node_coverage_t>& node_coverage, std::map<std::vector<handle_t>, coverage_t>& path_coverage)
  {
    bool success = false;
    BestCoverage<SimpleCoverage> best;
    graph.follow_edges(path.back(), false, [&](const handle_t& next)
    {
      success = true;
      if(path.size() + 1 < k) // Node coverage.
      {
        size_t first = find_first(node_coverage, graph.get_id(next));
        best.update(node_coverage[first].second, next);
      }
      else
      {
        std::vector<handle_t> window = forward_window(graph, path, next, k);
        best.update(path_coverage[window], next);
      }
    });

    if(success)
    {
      if(path.size() + 1 >= k)
      {
        std::vector<handle_t> window = forward_window(graph, path, best.handle, k);
        increase_coverage(path_coverage[window]);
      }
      size_t first = find_first(node_coverage, graph.get_id(best.handle));
      increase_coverage(node_coverage[first]);
      path.push_back(best.handle);
    }

    return success;
  }

  static bool extend_backward(const graph_t& graph, std::deque<handle_t>& path, size_t k, std::vector<node_coverage_t>& node_coverage, std::map<std::vector<handle_t>, coverage_t>& path_coverage)
  {
    bool success = false;
    BestCoverage<SimpleCoverage> best;
    graph.follow_edges(path.front(), true, [&](const handle_t& prev)
    {
      success = true;
      if(path.size() + 1 < k) // Node coverage.
      {
        size_t first = find_first(node_coverage, graph.get_id(prev));
        best.update(node_coverage[first].second, prev);
      }
      else
      {
        std::vector<handle_t> window = backward_window(graph, path, prev, k);
        best.update(path_coverage[window], prev);
      }
    });

    if(success)
    {
      if(path.size() + 1 >= k)
      {
        std::vector<handle_t> window = backward_window(graph, path, best.handle, k);
        increase_coverage(path_coverage[window]);
      }
      size_t first = find_first(node_coverage, graph.get_id(best.handle));
      increase_coverage(node_coverage[first]);
      path.push_front(best.handle);
    }

    return success;
  }

  static void increase_coverage(coverage_t& coverage)
  {
    coverage++;
  }

  static void increase_coverage(node_coverage_t& node)
  {
    increase_coverage(node.second);
  }

  static coverage_t worst_coverage() { return std::numeric_limits<coverage_t>::max(); }
};

//------------------------------------------------------------------------------

/*
  The best candidate is the one with the highest ratio

    true_coverage / (selected_coverage + 1).
*/

struct LocalHaplotypes
{
  typedef GBWTGraph graph_t;

  struct coverage_t
  {
    size_t selected_coverage;
    size_t true_coverage;
    double score;

    // Give priority to high scores.
    bool operator<(const coverage_t& another) const
    {
      return (this->score > another.score);
    }

    coverage_t() : selected_coverage(0), true_coverage(0), score(0.0) {}

    coverage_t(size_t true_coverage) : selected_coverage(0), true_coverage(true_coverage), score(true_coverage) {}

    void compute_score() { this->score = this->true_coverage / (this->selected_coverage + 1.0); }
  };
  typedef std::pair<nid_t, coverage_t> node_coverage_t;

  static std::vector<node_coverage_t> init_node_coverage(const graph_t& graph, const std::vector<nid_t>& component)
  {
    std::vector<node_coverage_t> node_coverage;
    node_coverage.reserve(component.size());
    for(nid_t id : component)
    {
      coverage_t coverage;
      coverage.selected_coverage = 0;
      coverage.true_coverage = graph.index->nodeSize(gbwt::Node::encode(id, false));
      coverage.compute_score();
      node_coverage.emplace_back(id, coverage);
    }
    return node_coverage;
  }

  static bool extend_forward(const graph_t& graph, std::deque<handle_t>& path, size_t k, std::vector<node_coverage_t>& node_coverage, std::map<std::vector<handle_t>, coverage_t>& path_coverage)
  {
    bool success = false;
    BestCoverage<LocalHaplotypes> best;
    auto start = (path.size() + 1 < k ? path.begin() : path.end() - (k - 1));
    std::vector<handle_t> context(start, path.end());
    gbwt::BidirectionalState state = graph.bd_find(context);
    graph.follow_paths(state, false, [&](const gbwt::BidirectionalState& next) -> bool
    {
      success = true;
      handle_t handle = GBWTGraph::node_to_handle(next.forward.node);
      if(path.size() + 1 < k) // Node coverage.
      {
        size_t first = find_first(node_coverage, graph.get_id(handle));
        best.update(node_coverage[first].second, handle);
      }
      else
      {
        std::vector<handle_t> window = forward_window(graph, path, handle, k);
        // Insert empty coverage or find the existing coverage.
        auto result = path_coverage.insert({ window, coverage_t(next.size()) });
        best.update(result.first->second, handle);
      }
      return true;
    });

    if(success)
    {
      if(path.size() + 1 >= k)
      {
        std::vector<handle_t> window = forward_window(graph, path, best.handle, k);
        increase_coverage(path_coverage[window]);
      }
      size_t first = find_first(node_coverage, graph.get_id(best.handle));
      increase_coverage(node_coverage[first]);
      path.push_back(best.handle);
    }

    return success;
  }

  static bool extend_backward(const graph_t& graph, std::deque<handle_t>& path, size_t k, std::vector<node_coverage_t>& node_coverage, std::map<std::vector<handle_t>, coverage_t>& path_coverage)
  {
    bool success = false;
    BestCoverage<LocalHaplotypes> best;
    auto limit = (path.size() + 1 < k ? path.end() : path.begin() + (k - 1));
    std::vector<handle_t> context(path.begin(), limit);
    gbwt::BidirectionalState state = graph.bd_find(context);
    graph.follow_paths(state, true, [&](const gbwt::BidirectionalState& prev) -> bool
    {
      success = true;
      handle_t handle = GBWTGraph::node_to_handle(prev.backward.node);
      handle = graph.flip(handle); // Get the correct orientation.
      if(path.size() + 1 < k) // Node coverage.
      {
        size_t first = find_first(node_coverage, graph.get_id(handle));
        best.update(node_coverage[first].second, handle);
      }
      else
      {
        std::vector<handle_t> window = backward_window(graph, path, handle, k);
        // Insert empty coverage or find the existing coverage.
        auto result = path_coverage.insert({ window, coverage_t(prev.size()) });
        best.update(result.first->second, handle);
      }
      return true;
    });

    if(success)
    {
      if(path.size() + 1 >= k)
      {
        std::vector<handle_t> window = backward_window(graph, path, best.handle, k);
        increase_coverage(path_coverage[window]);
      }
      size_t first = find_first(node_coverage, graph.get_id(best.handle));
      increase_coverage(node_coverage[first]);
      path.push_front(best.handle);
    }

    return success;
  }

  static void increase_coverage(coverage_t& coverage)
  {
    coverage.selected_coverage++;
    coverage.compute_score();
  }

  static void increase_coverage(node_coverage_t& node)
  {
    increase_coverage(node.second);
  }

  static coverage_t worst_coverage() { return coverage_t(); }
};

//------------------------------------------------------------------------------

template<class Coverage>
gbwt::GBWT
generic_path_cover(const typename Coverage::graph_t& graph, size_t n, size_t k, gbwt::size_type batch_size, gbwt::size_type sample_interval, bool show_progress)
{
  typedef typename Coverage::coverage_t coverage_t;
  typedef typename Coverage::node_coverage_t node_coverage_t;

  // Sanity checks.
  size_t node_count = graph.get_node_count();
  if(node_count == 0 || n == 0) { return gbwt::GBWT(); }
  if(k < PATH_COVER_MIN_K)
  {
    std::cerr << "generic_path_cover(): Window length (" << k << ") must be at least " << PATH_COVER_MIN_K << std::endl;
    return gbwt::GBWT();
  }
  nid_t min_id = graph.min_node_id();
  if(min_id < 1)
  {
    std::cerr << "generic_path_cover(): Minimum node id (" << min_id << ") must be positive" << std::endl;
    return gbwt::GBWT();
  }
  nid_t max_id = graph.max_node_id();

  // Find weakly connected components, ignoring the directions of the edges.
  std::vector<std::vector<nid_t>> components = weakly_connected_components(graph);

  // GBWT construction parameters. Adjust the batch size down for small graphs.
  // We will also set basic metadata: n samples with each component as a separate contig.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type node_width = gbwt::bit_length(gbwt::Node::encode(max_id, true));
  batch_size = std::min(batch_size, static_cast<gbwt::size_type>(2 * n * (node_count + components.size())));
  gbwt::GBWTBuilder builder(node_width, batch_size, sample_interval);
  builder.index.addMetadata();

  // Handle each component separately.
  for(size_t contig = 0; contig < components.size(); contig++)
  {
    if(show_progress)
    {
      std::cerr << "Processing component " << (contig + 1) << " / " << components.size() << std::endl;
    }
    std::vector<nid_t>& component = components[contig];
    std::vector<node_coverage_t> node_coverage = Coverage::init_node_coverage(graph, component);
    component = std::vector<nid_t>(); // Save a little bit of memory.
    std::map<std::vector<handle_t>, coverage_t> path_coverage; // Path and its reverse complement are equivalent.

    // Generate n paths in the component.
    for(size_t i = 0; i < n; i++)
    {
      // Choose a starting node with minimum coverage and then sort the nodes by id.
      std::deque<handle_t> path;
      std::sort(node_coverage.begin(), node_coverage.end(), [](const node_coverage_t& a, const node_coverage_t& b) -> bool
      {
        return (a.second < b.second);
      });
      path.push_back(graph.get_handle(node_coverage.front().first, false));
      Coverage::increase_coverage(node_coverage.front());
      std::sort(node_coverage.begin(), node_coverage.end(), [](const node_coverage_t& a, const node_coverage_t& b) -> bool
      {
        return (a.first < b.first);
      });

      // Extend the path in both directions.
      bool success = true;
      while(success && path.size() < node_coverage.size())
      {
        success = false;
        success |= Coverage::extend_forward(graph, path, k, node_coverage, path_coverage);
        if(path.size() >= node_coverage.size()) { break; }
        success |= Coverage::extend_backward(graph, path, k, node_coverage, path_coverage);
      }

      // Insert the path and its name into the index.
      gbwt::vector_type buffer;
      buffer.reserve(path.size());
      for(handle_t handle : path)
      {
        buffer.push_back(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
      }
      builder.insert(buffer, true);
      builder.index.metadata.addPath(
      {
        static_cast<gbwt::PathName::path_name_type>(i),
        static_cast<gbwt::PathName::path_name_type>(contig),
        static_cast<gbwt::PathName::path_name_type>(0),
        static_cast<gbwt::PathName::path_name_type>(0)
      });
    }
  }

  // Finish the construction, add basic metadata, and return the GBWT.
  builder.finish();
  builder.index.metadata.setSamples(n);
  builder.index.metadata.setContigs(components.size());
  builder.index.metadata.setHaplotypes(n);
  if(show_progress)
  {
    gbwt::operator<<(std::cerr, builder.index.metadata) << std::endl;
  }
  return gbwt::GBWT(builder.index);
}

//------------------------------------------------------------------------------

gbwt::GBWT
path_cover_gbwt(const HandleGraph& graph, size_t n, size_t k, gbwt::size_type batch_size, gbwt::size_type sample_interval, bool show_progress)
{
  return generic_path_cover<SimpleCoverage>(graph, n, k, batch_size, sample_interval, show_progress);
}

gbwt::GBWT
local_haplotypes(const GBWTGraph& graph, size_t n, size_t k, gbwt::size_type batch_size, gbwt::size_type sample_interval, bool show_progress)
{
  return generic_path_cover<LocalHaplotypes>(graph, n, k, batch_size, sample_interval, show_progress);
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
