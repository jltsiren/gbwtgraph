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

typedef std::pair<nid_t, size_t> coverage_t;

std::vector<coverage_t>::iterator
find_first(std::vector<coverage_t>& array, nid_t id)
{
  auto iter = std::lower_bound(array.begin(), array.end(), coverage_t(id, 0));
  assert(iter != array.end());
  assert(iter->first == id);
  return iter;
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

//------------------------------------------------------------------------------

gbwt::GBWT
path_cover_gbwt(const HandleGraph& graph, size_t n, size_t k, gbwt::size_type batch_size, gbwt::size_type sample_interval)
{
  // Sanity checks.
  size_t node_count = graph.get_node_count();
  if(node_count == 0 || n == 0) { return gbwt::GBWT(); }
  if(k <= 1)
  {
    std::cerr << "path_cover_gbwt(): Window length (" << k << ") must be > 1" << std::endl;
    return gbwt::GBWT();
  }
  nid_t min_id = graph.min_node_id();
  if(min_id < 1)
  {
    std::cerr << "path_cover_gbwt(): Minimum node id (" << min_id << ") must be positive" << std::endl;
    return gbwt::GBWT();
  }
  nid_t max_id = graph.max_node_id();

  // Find weakly connected components, ignoring the directions of the edges.
  std::vector<std::vector<nid_t>> components;
  {
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
  }

  // GBWT construction parameters. Adjust the batch size down for small graphs.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type node_width = gbwt::bit_length(gbwt::Node::encode(max_id, true));
  batch_size = std::min(batch_size, static_cast<gbwt::size_type>(2 * n * (node_count + components.size())));
  gbwt::GBWTBuilder builder(node_width, batch_size, sample_interval);

  // Handle each component separately.
  for(std::vector<nid_t>& component : components)
  {
    std::vector<coverage_t> node_coverage;
    node_coverage.reserve(component.size());
    for(nid_t id : component) { node_coverage.emplace_back(id, 0); }
    component = std::vector<nid_t>(); // Save a little bit of memory.
    std::map<std::vector<handle_t>, size_t> path_coverage; // Path and its reverse complement are equivalent.

    // Generate n paths in the component.
    for(size_t i = 0; i < n; i++)
    {
      // Choose a starting node with minimum coverage and then sort the nodes by id.
      std::deque<handle_t> path;
      std::sort(node_coverage.begin(), node_coverage.end(), [](coverage_t a, coverage_t b) -> bool
      {
        return(a.second < b.second);
      });
      path.push_back(graph.get_handle(node_coverage.front().first, false));
      node_coverage.back().second++;
      std::sort(node_coverage.begin(), node_coverage.end(), [](coverage_t a, coverage_t b) -> bool
      {
        return(a.first < b.first);
      });

      // Extend the path in both directions.
      bool forward_success = true, backward_success = true;
      while((forward_success || backward_success) && path.size() < node_coverage.size())
      {
        size_t best_coverage;
        handle_t best;
        auto update_best = [&best_coverage, &best](size_t coverage, const handle_t& candidate)
        {
          if(coverage < best_coverage)
          {
            best_coverage = coverage;
            best = candidate;
          }
        };

        // Extend forward.
        forward_success = false;
        best_coverage = std::numeric_limits<size_t>::max();
        graph.follow_edges(path.back(), false, [&](const handle_t& next)
        {
          forward_success = true;
          if(path.size() + 1 < k) // Node coverage.
          {
            auto iter = find_first(node_coverage, graph.get_id(next));
            update_best(iter->second, next);
          }
          else
          {
            std::vector<handle_t> window = forward_window(graph, path, next, k);
            update_best(path_coverage[window], next);
          }
        });
        if(forward_success)
        {
          if(path.size() + 1 >= k)
          {
            std::vector<handle_t> window = forward_window(graph, path, best, k);
            path_coverage[window]++;
          }
          auto iter = find_first(node_coverage, graph.get_id(best));
          iter->second++;
          path.push_back(best);
          if(path.size() >= node_coverage.size()) { break; }
        }

        // Extend backward.
        backward_success = false;
        best_coverage = std::numeric_limits<size_t>::max();
        graph.follow_edges(path.front(), true, [&](const handle_t& prev)
        {
          backward_success = true;
          if(path.size() + 1 < k) // Node coverage.
          {
            auto iter = find_first(node_coverage, graph.get_id(prev));
            update_best(iter->second, prev);
          }
          else
          {
            std::vector<handle_t> window = backward_window(graph, path, prev, k);
            update_best(path_coverage[window], prev);
          }
        });
        if(backward_success)
        {
          if(path.size() + 1 >= k)
          {
            std::vector<handle_t> window = backward_window(graph, path, best, k);
            path_coverage[window]++;
          }
          auto iter = find_first(node_coverage, graph.get_id(best));
          iter->second++;
          path.push_front(best);
        }
      }

      // Insert the path into the index.
      gbwt::vector_type buffer;
      buffer.reserve(path.size());
      for(handle_t handle : path)
      {
        buffer.push_back(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
      }
      builder.insert(buffer, true);
    }
  }

  return gbwt::GBWT(builder.index);
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
