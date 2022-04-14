#include <gbwtgraph/path_cover.h>

#include <gbwtgraph/algorithms.h>

#include <algorithm>
#include <deque>
#include <limits>
#include <map>
#include <string>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

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
  if(path.size() + 1 < k) { k = path.size() + 1; } // Handle the short initial paths in DAGs.
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

  static std::vector<node_coverage_t> init_node_coverage(const graph_t& graph, const std::vector<nid_t>& component)
  {
    std::vector<node_coverage_t> node_coverage;
    node_coverage.reserve(component.size());
    for(nid_t id : component)
    {
      // We are actually interested in the intersection of this graph and the component.
      // For example, some nodes of the original graph may be missing from a GBWTGraph.
      if(!(graph.has_node(id))) { continue; }
      node_coverage.emplace_back(id, static_cast<coverage_t>(0));
    }
    return node_coverage;
  }

  static bool extend_forward(const graph_t& graph, std::deque<handle_t>& path, size_t k, std::vector<node_coverage_t>& node_coverage, std::map<std::vector<handle_t>, coverage_t>& path_coverage, bool acyclic)
  {
    bool success = false;
    BestCoverage<SimpleCoverage> best;
    graph.follow_edges(path.back(), false, [&](const handle_t& next)
    {
      success = true;
      if(!acyclic && path.size() + 1 < k) // Node coverage.
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
      if(acyclic || path.size() + 1 >= k)
      {
        std::vector<handle_t> window = forward_window(graph, path, best.handle, k);
        increase_coverage(path_coverage[window]);
      }
      if(!acyclic)
      {
        size_t first = find_first(node_coverage, graph.get_id(best.handle));
        increase_coverage(node_coverage[first]);
      }
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

  static std::string name() { return "SimpleCoverage"; }
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
      // We are actually interested in the intersection of this graph and the component.
      // For example, some nodes of the original graph may be missing from a GBWTGraph.
      // This also implies that the true coverage of selected nodes is nonzero.
      if(!(graph.has_node(id))) { continue; }
      coverage_t coverage;
      coverage.selected_coverage = 0;
      coverage.true_coverage = graph.index->nodeSize(gbwt::Node::encode(id, false));
      if(coverage.true_coverage == 0) { continue; }
      coverage.compute_score();
      node_coverage.emplace_back(id, coverage);
    }
    return node_coverage;
  }

  static bool extend_forward(const graph_t& graph, std::deque<handle_t>& path, size_t k, std::vector<node_coverage_t>& node_coverage, std::map<std::vector<handle_t>, coverage_t>& path_coverage, bool acyclic)
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
      if(!acyclic && path.size() + 1 < k) // Node coverage.
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
      if(acyclic || path.size() + 1 >= k)
      {
        std::vector<handle_t> window = forward_window(graph, path, best.handle, k);
        increase_coverage(path_coverage[window]);
      }
      if(!acyclic)
      {
        size_t first = find_first(node_coverage, graph.get_id(best.handle));
        increase_coverage(node_coverage[first]);
      }
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

  static std::string name() { return "LocalHaplotypes"; }
};

//------------------------------------------------------------------------------

bool
path_cover_sanity_checks(const HandleGraph& graph, size_t n, size_t k)
{
  if(graph.get_node_count() == 0 || n == 0) { return false; }
  if(k < PATH_COVER_MIN_K)
  {
    std::cerr << "path_cover_sanity_checks(): Window length (" << k << ") must be at least " << PATH_COVER_MIN_K << std::endl;
    return false;
  }
  nid_t min_id = graph.min_node_id();
  if(min_id < 1)
  {
    std::cerr << "path_cover_sanity_checks(): Minimum node id (" << min_id << ") must be positive" << std::endl;
    return false;
  }
  return true;
}

template<class Coverage>
bool
component_path_cover(const typename Coverage::graph_t& graph, gbwt::GBWTBuilder& builder, std::vector<std::vector<nid_t>>& components, size_t component_id, size_t n, size_t k, bool show_progress, size_t sample_id_offset, size_t contig_id)
{
  typedef typename Coverage::coverage_t coverage_t;
  typedef typename Coverage::node_coverage_t node_coverage_t;

  std::vector<nid_t>& component = components[component_id];
  size_t component_size = component.size();
  std::vector<nid_t> head_nodes = is_nice_and_acyclic(graph, component);
  bool acyclic = !(head_nodes.empty());
  if(show_progress)
  {
    std::cerr << Coverage::name() << ": Processing component " << (component_id + 1) << " / " << components.size();
    if(acyclic) { std::cerr << " (acyclic)"; }
    std::cerr << std::endl;
  }

  // Node coverage for the potential starting nodes.
  std::vector<node_coverage_t> node_coverage = Coverage::init_node_coverage(graph, (acyclic ? head_nodes : component));
  std::map<std::vector<handle_t>, coverage_t> path_coverage; // Path and its reverse complement are equivalent.

  // Node coverage will be empty if we cannot create this type of path cover for the component.
  // For example, if there are no haplotypes for LocalHaplotypes.
  if(node_coverage.empty())
  {
    if(show_progress)
    {
      std::cerr << Coverage::name() << ": Cannot find this type of path cover for the component" << std::endl;
    }
    return false;
  }

  // Now that we know that we can find a path cover, we can save a little bit of memory by deleting
  // the component.
  component = std::vector<nid_t>();

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

    // Extend the path forward if acyclic or in both directions otherwise.
    bool success = true;
    while(success && path.size() < component_size)
    {
      success = false;
      success |= Coverage::extend_forward(graph, path, k, node_coverage, path_coverage, acyclic);
      if(!acyclic && path.size() < component_size)
      {
        success |= Coverage::extend_backward(graph, path, k, node_coverage, path_coverage);
      }
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
      static_cast<gbwt::PathName::path_name_type>(sample_id_offset + i),
      static_cast<gbwt::PathName::path_name_type>(contig_id),
      static_cast<gbwt::PathName::path_name_type>(0),
      static_cast<gbwt::PathName::path_name_type>(0)
    });
  }

  return true;
}

void
finish_path_cover(gbwt::GBWTBuilder& builder, size_t n, const std::vector<std::string>& contig_names, size_t haplotypes, bool show_progress)
{
  builder.finish();
  
  // Record each added sample, with names if needed
  if(builder.index.metadata.hasSampleNames())
  {
    std::vector<std::string> sample_names;
    for(size_t i = 0; i < n; i++) { sample_names.push_back("path_cover_" + std::to_string(i)); }
    builder.index.metadata.addSamples(sample_names);
  }
  else
  {
    builder.index.metadata.setSamples(builder.index.metadata.samples() + n);
  }
  
  // Record each added contig, with names if needed
  if(builder.index.metadata.hasContigNames())
  {
    builder.index.metadata.addContigs(contig_names);
  }
  else
  {
    builder.index.metadata.setContigs(builder.index.metadata.contigs() + contig_names.size());
  }
  
  // Record the additional stored haplotypes
  builder.index.metadata.setHaplotypes(builder.index.metadata.haplotypes() + haplotypes);
  
  if(show_progress)
  {
    gbwt::operator<<(std::cerr, builder.index.metadata) << std::endl;
  }
}

void
finish_path_cover(gbwt::GBWTBuilder& builder, size_t n, size_t contigs, size_t haplotypes, bool show_progress)
{
  std::vector<std::string> contig_names(contigs);
  if(builder.index.metadata.hasContigNames())
  {
    // Need to fill in contig names
    for(size_t i = 0; i < contigs; i++) { contig_names[i] = ("path_cover_contig_" + std::to_string(i)); }
  }
  finish_path_cover(builder, n, contig_names, haplotypes, show_progress);
}

/*
  Get a contig number that each component touches, or std::numeric_limits<size_t>::max() if it doesn't touch anything.
  If path_filter is set, only considers paths that pass the filter.
*/
std::vector<size_t>
find_contigs_for_components(const PathHandleGraph* path_graph,
                            const gbwt::GBWTBuilder& builder,
                            const std::vector<std::vector<nid_t>>& components,
                            const std::function<bool(const path_handle_t&)>* path_filter = nullptr)
{
  // Each component is going to get a contig number. If a component shares a
  // node with a path, we need it to use the same contig number as that path.
  std::vector<size_t> component_contigs;
  
  for(auto& component : components)
  {
    // What contig if any did we find for this component
    size_t component_contig = std::numeric_limits<size_t>::max();
    for(auto& node_id : component)
    {
      // Look at every node in the component
      bool finished = path_graph->for_each_step_on_handle(path_graph->get_handle(node_id), [&](const step_handle_t& step)
      {
        // And then each step on a path on the node
        path_handle_t path = path_graph->get_path_handle_of_step(step);
        if (path_filter && !(*path_filter)(path)) {
          // Skip this path and keep going
          return true;
        }
        std::string path_name = path_graph->get_path_name(path);
        // And see if we have a contig for the path already
        size_t path_contig = builder.index.metadata.contig(path_name);
        if(path_contig != builder.index.metadata.contigs())
        {
          // Component touches a node that touches this path that is a contig already.
          // Component belongs on this contig.
          component_contig = path_contig;
          // Stop scanning the component
          return false;
        }
        // Otherwise continue scanning the component
        return true;
      });
      if(!finished)
      {
        // We're done early!
        break;
      }
    }
    
    // Save what we found, or that we found nothing.
    component_contigs.push_back(component_contig);
  }
  
  return component_contigs;
}


//------------------------------------------------------------------------------

void
store_named_paths(gbwt::GBWTBuilder& builder, const PathHandleGraph& graph, const std::function<bool(const path_handle_t&)>* path_filter)
{
  if(!builder.index.metadata.hasContigNames() && builder.index.metadata.contigs() > 0) {
    throw std::logic_error("Cannot add paths to an index with existing unnamed contigs");
  }
  if(!builder.index.metadata.hasSampleNames() && builder.index.metadata.samples() > 0) {
    throw std::logic_error("Cannot add paths to an index with existing unnamed samples");
  }
  if(!builder.index.metadata.hasPathNames() && builder.index.metadata.paths() > 0) {
    throw std::logic_error("Cannot add paths to an index with existing unnamed paths");
  }
  
  // TODO: Change this to use the path metadata API!
  
  size_t paths_to_add = 0;
  {
    // Work out what new contigs paths correspond to
    std::vector<std::string> new_contig_names;
    graph.for_each_path_handle([&](const path_handle_t& path)
    {
      if(graph.is_empty(path) || (path_filter && !(*path_filter)(path)))
      {
        // Skip this path
        return;
      }
      paths_to_add++;
      std::string path_name = graph.get_path_name(path);
      if(builder.index.metadata.contig(path_name) == builder.index.metadata.contigs())
      {
        // This is a new contig we will need to give a number to
        new_contig_names.push_back(path_name);
      }
    });
    
    if(paths_to_add == 0)
    {
      // We aren't actually adding anything
      return;
    }
    
    if(!new_contig_names.empty())
    {
      // Give numbers to these new contigs.
      builder.index.metadata.addContigs(new_contig_names);
    }
  }
  
  // What sample number should named generic path threads belong to?
  size_t generic_path_sample = builder.index.metadata.sample(gbwtgraph::NAMED_PATH_SAMPLE_PREFIX);
  if(generic_path_sample == builder.index.metadata.samples())
  {
    // No path sample is stored yet, so add one. It will end up at that index.
    builder.index.metadata.addSamples({gbwtgraph::NAMED_PATH_SAMPLE_PREFIX});
    // And allocate a haplotype for all of its contigs
    builder.index.metadata.setHaplotypes(builder.index.metadata.haplotypes() + 1);
  }
  
  // Now actually add all the paths
  graph.for_each_path_handle([&](const path_handle_t& path)
  {
    if(graph.is_empty(path) || (path_filter && !(*path_filter)(path)))
    {
      // Skip this path
      return;
    }
    std::string path_name = graph.get_path_name(path);
    gbwt::vector_type buffer;
    // TODO: This counting can be a whole path scan itself. Can we know and avoid it?
    buffer.reserve(graph.get_step_count(path));
    for(handle_t handle : graph.scan_path(path))
    {
      buffer.push_back(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
    }
    builder.insert(buffer, true); // Insert in both orientations.
    // Each generic path is the special sample (generic_path_sample), on a contig we need to look up.
    size_t contig_number = builder.index.metadata.contig(path_name);
    // Record that the path is in the special sample, on the contig we found or made.
    builder.index.metadata.addPath(generic_path_sample, contig_number, 0, 0);
  });
  
}

//------------------------------------------------------------------------------

gbwt::GBWT
path_cover_gbwt(const HandleGraph& graph,
                size_t n,
                size_t k,
                gbwt::size_type batch_size,
                gbwt::size_type sample_interval,
                bool include_named_paths,
                const std::function<bool(const path_handle_t&)>* path_filter,
                bool show_progress)
{
  // Sanity checks.
  if(!path_cover_sanity_checks(graph, n, k))
  {
    return gbwt::GBWT();
  }

  // Find weakly connected components, ignoring the direction of the edges.
  std::vector<std::vector<nid_t>> components = weakly_connected_components(graph);

  // GBWT construction parameters. Adjust the batch size down for small graphs.
  // We will also set basic metadata: n samples with each component as a separate contig.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type node_width = sdsl::bits::length(gbwt::Node::encode(graph.max_node_id(), true));
  gbwt::GBWTBuilder builder(node_width, batch_size, sample_interval);
  builder.index.addMetadata();
  
  // Each component is going to get a contig number. If a component shares a
  // node with a path, we need it to use the same contig number as that path.
  // If this is empty, we just assign fresh contigs to everything.
  std::vector<size_t> component_contigs;
  if(include_named_paths)
  {
    const PathHandleGraph* path_graph = dynamic_cast<const PathHandleGraph*>(&graph);
    if(path_graph && path_graph->get_path_count() > 0)
    {
      // Copy over all named paths
      store_named_paths(builder, *path_graph, path_filter);
      
      // Find the right contigs for all the components, by path name
      component_contigs = find_contigs_for_components(path_graph, builder, components, path_filter);
    }
  }
 
  // Handle each component separately.
  size_t base_sample = builder.index.metadata.samples();
  size_t next_contig = builder.index.metadata.contigs();
  for(size_t component = 0; component < components.size(); component++)
  {
    // Decide what contig this component belongs to
    size_t assigned_contig = component_contigs.empty() ? std::numeric_limits<size_t>::max() : component_contigs[component];
    size_t contig = assigned_contig == std::numeric_limits<size_t>::max() ? next_contig : assigned_contig;
    if(component_path_cover<SimpleCoverage>(graph, builder, components, component, n, k, show_progress, base_sample, contig))
    {
      if(assigned_contig == std::numeric_limits<size_t>::max())
      {
        // We used a new contig slot
        next_contig++;
      }
    }
  }

  // Finish the construction, add basic metadata, and return the GBWT.
  finish_path_cover(builder, n, next_contig - builder.index.metadata.contigs(), n, show_progress);
  return gbwt::GBWT(builder.index);
}

//------------------------------------------------------------------------------

gbwt::GBWT
local_haplotypes(const HandleGraph& graph,
                 const gbwt::GBWT& index,
                 size_t n,
                 size_t k,
                 gbwt::size_type batch_size,
                 gbwt::size_type sample_interval,
                 bool include_named_paths,
                 const std::function<bool(const path_handle_t&)>* path_filter,
                 bool show_progress)
{
  // Sanity checks.
  if(!path_cover_sanity_checks(graph, n, k))
  {
    return gbwt::GBWT();
  }

  // Find weakly connected components, ignoring the direction of the edges.
  std::vector<std::vector<nid_t>> components = weakly_connected_components(graph);

  // GBWT construction parameters. Adjust the batch size down for small graphs.
  // We will also set basic metadata: n samples with each component as a separate contig.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type node_width = sdsl::bits::length(gbwt::Node::encode(graph.max_node_id(), true));
  gbwt::GBWTBuilder builder(node_width, batch_size, sample_interval);
  builder.index.addMetadata();

  // If the graph we were given is a GBWTGraph using the same GBWT index, we can
  // use it directly for sampling local haplotypes. Otherwise we have to build a
  // temporary graph.
  const GBWTGraph* gbwt_graph = dynamic_cast<const GBWTGraph*>(&graph);
  if(gbwt_graph != nullptr)
  {
    if(gbwt_graph->index != &index) { gbwt_graph = nullptr; }
  }
  GBWTGraph created_gbwt_graph;
  if(gbwt_graph == nullptr)
  {
    if(show_progress)
    {
      std::cerr << "Building a temporary GBWTGraph" << std::endl;
    }
    created_gbwt_graph = GBWTGraph(index, graph);
    gbwt_graph = &created_gbwt_graph;
  }
  
  // Each component is going to get a contig number. If a component shares a
  // node with a path, we need it to use the same contig number as that path.
  // If this is empty, we just assign fresh contigs to everything.
  std::vector<size_t> component_contigs;
  if(include_named_paths)
  {
    const PathHandleGraph* path_graph = dynamic_cast<const PathHandleGraph*>(&graph);
    if(path_graph && path_graph->get_path_count() > 0)
    {
      // Copy over all named paths
      store_named_paths(builder, *path_graph, path_filter);
      
      // Find the right contigs for all the components, by path name
      component_contigs = find_contigs_for_components(path_graph, builder, components, path_filter);
    }
  }

  // Handle each component separately.
  size_t base_sample = builder.index.metadata.samples();
  size_t next_contig = builder.index.metadata.contigs();
  for(size_t component = 0; component < components.size(); component++)
  {
    // Decide what contig this component belongs to
    size_t assigned_contig = component_contigs.empty() ? std::numeric_limits<size_t>::max() : component_contigs[component];
    size_t contig = assigned_contig == std::numeric_limits<size_t>::max() ? next_contig : assigned_contig;
    // Revert to regular path cover if we cannot sample local haplotypes.
    if(component_path_cover<LocalHaplotypes>(*gbwt_graph, builder, components, component, n, k, show_progress, base_sample, contig) ||
       component_path_cover<SimpleCoverage>(graph, builder, components, component, n, k, show_progress, base_sample, contig))
    {
      if(assigned_contig == std::numeric_limits<size_t>::max())
      {
        // We used a new contig slot
        next_contig++;
      }
    }
  }

  // Finish the construction, add basic metadata, and return the GBWT.
  finish_path_cover(builder, n, next_contig - builder.index.metadata.contigs(), n, show_progress);
  return gbwt::GBWT(builder.index);
}

//------------------------------------------------------------------------------

size_t
augment_gbwt(const HandleGraph& graph,
             gbwt::DynamicGBWT& index,
             size_t n,
             size_t k,
             gbwt::size_type batch_size,
             gbwt::size_type sample_interval,
             bool show_progress)
{
  // Sanity checks.
  if(!path_cover_sanity_checks(graph, n, k)) { return 0; }
  if(!(index.bidirectional()))
  {
    std::cerr << "augment_gbwt(): The GBWT index must be bidirectional" << std::endl;
    return 0;
  }

  // Find weakly connected components, ignoring the direction of the edges.
  std::vector<std::vector<nid_t>> components = weakly_connected_components(graph);

  // GBWT construction parameters. Adjust the batch size down for small graphs.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type node_width = sdsl::bits::length(gbwt::Node::encode(graph.max_node_id(), true));
  gbwt::GBWTBuilder builder(node_width, batch_size, sample_interval);
  builder.swapIndex(index);

  // Handle each component separately, but only if there are no GBWT paths in it.
  std::vector<std::string> contig_names;
  size_t sample_id_offset = builder.index.metadata.samples();
  for(size_t contig = 0; contig < components.size(); contig++)
  {
    bool has_paths = false;
    for(nid_t nid : components[contig])
    {
      gbwt::node_type node = gbwt::Node::encode(nid, false);
      has_paths |= (builder.index.contains(node) && !(builder.index.empty(node)));
    }
    if(has_paths) { continue; }
    size_t contig_id = builder.index.metadata.contigs() + contig_names.size();
    if(component_path_cover<SimpleCoverage>(graph, builder, components, contig, n, k, show_progress, sample_id_offset, contig_id))
    {
      contig_names.emplace_back("component_" + std::to_string(contig));
    }
  }
  
  if(!contig_names.empty())
  {
    // Finish the construction, augment metadata, and return the number of augmented components.
    finish_path_cover(builder, n, contig_names, 0, show_progress);
  }
  builder.swapIndex(index);
  return contig_names.size();
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
