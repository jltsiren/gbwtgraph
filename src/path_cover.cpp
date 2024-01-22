#include <gbwtgraph/path_cover.h>

#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/internal.h>

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

// Create a path cover of the given component and insert it into the builder.
// This does not create metadata. Returns true if the path cover was successfully
// created.
template<class Coverage>
bool
component_path_cover(
  const typename Coverage::graph_t& graph,
  gbwt::GBWTBuilder& builder,
  const std::vector<nid_t>& component,
  size_t job_id, size_t component_id,
  const PathCoverParameters& parameters)
{
  typedef typename Coverage::coverage_t coverage_t;
  typedef typename Coverage::node_coverage_t node_coverage_t;

  std::vector<nid_t> head_nodes = is_nice_and_acyclic(graph, component);
  bool acyclic = !(head_nodes.empty());
  if(parameters.show_progress)
  {
    std::string msg =
      "Job " + std::to_string(job_id) +
      ", component " + std::to_string(component_id) +
      ": " + Coverage::name();
    if(acyclic) { msg += " (acyclic)"; }
    #pragma omp critical
    {
      std::cerr << msg << std::endl;
    }
  }

  // Node coverage for the potential starting nodes.
  std::vector<node_coverage_t> node_coverage = Coverage::init_node_coverage(graph, (acyclic ? head_nodes : component));
  std::map<std::vector<handle_t>, coverage_t> path_coverage; // Path and its reverse complement are equivalent.

  // Node coverage will be empty if we cannot create this type of path cover for the component.
  // For example, if there are no haplotypes for LocalHaplotypes.
  if(node_coverage.empty())
  {
    if(parameters.show_progress)
    {
      #pragma omp critical
      {
        std::cerr << Coverage::name() << ": Cannot find this type of path cover for component " << component_id << std::endl;
      }
    }
    return false;
  }

  // Generate n paths in the component.
  for(size_t i = 0; i < parameters.num_paths; i++)
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
    while(success && path.size() < component.size())
    {
      success = false;
      success |= Coverage::extend_forward(graph, path, parameters.context, node_coverage, path_coverage, acyclic);
      if(!acyclic && path.size() < component.size())
      {
        success |= Coverage::extend_backward(graph, path, parameters.context, node_coverage, path_coverage);
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
  }

  return true;
}

//------------------------------------------------------------------------------

void
store_named_paths(gbwt::GBWTBuilder& builder, const PathHandleGraph& graph, const std::function<bool(const path_handle_t&)>* path_filter)
{
  // What path senses should we copy across?
  std::unordered_set<PathSense> named_senses = {PathSense::GENERIC, PathSense::REFERENCE};

  store_paths(builder, graph, named_senses, path_filter);
}

void
store_paths(gbwt::GBWTBuilder& builder, const PathHandleGraph& graph, const std::unordered_set<PathSense>& senses, const std::function<bool(const path_handle_t&)>* path_filter)
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

  // Turn the metadata back into a MetadataBuilder.
  MetadataBuilder metadata_builder(builder.index.metadata);

  // We set this to true if we found any paths
  bool added_paths = false;

  // We can only use one sense per sample currently. We use this to work out what it should be.
  std::unordered_map<std::string, PathSense> sample_sense;

  // Work out what new contigs paths correspond to
  graph.for_each_path_matching(&senses, nullptr, nullptr, [&](const path_handle_t& path)
  {
    if(path_filter != nullptr && !(*path_filter)(path))
    {
      // The filter wants us to skip this path
      return;
    }
  
    // Get the sense
    PathSense sense = graph.get_sense(path);
    // Check for divergent senses for the sample
    std::string sample_name = graph.get_sample_name(path);

    auto it = sample_sense.find(sample_name);
    if(it != sample_sense.end())
    {
      if(it->second != sense)
      {
        // Complain to the user.
        std::cerr << "gbwtgraph::store_paths(): Warning: multiple senses of path in sample \""
                  << sample_name << "\"; only one will be stored";
      }
    }
    else
    {
      sample_sense.emplace_hint(it, sample_name, sense);
    }

    // Store the path metadata
    metadata_builder.add_path(
      sense,
      sample_name,
      graph.get_locus_name(path),
      graph.get_haplotype(path),
      graph.get_phase_block(path),
      graph.get_subrange(path)
    );

    // Add the path to the GBWT
    gbwt::vector_type buffer;
    // TODO: This counting can be a whole path scan itself. Can we know and avoid it?
    buffer.reserve(graph.get_step_count(path));
    for(handle_t handle : graph.scan_path(path))
    {
      buffer.push_back(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
    }
    builder.insert(buffer, true); // Insert in both orientations.

    added_paths = true;
  });

  if(!sample_sense.empty())
  {
    // Paths were stored

    // Commit back metadata changes
    builder.index.metadata = metadata_builder.get_metadata();

    // Commit path senses
    set_sample_path_senses(builder.index.tags, sample_sense);
  }
}

//------------------------------------------------------------------------------

// FIXME These should be exposed and tested.

// Assigns the paths we want to include to the construction jobs. Alse generates
// the metadata for them.
std::vector<std::vector<path_handle_t>>
assign_paths(
  const PathHandleGraph& graph,
  const ConstructionJobs& jobs,
  MetadataBuilder& metadata,
  const std::function<bool(const path_handle_t&)>* path_filter)
{
  std::vector<std::vector<path_handle_t>> result(jobs.size());
  std::unordered_set<PathSense> senses = { PathSense::GENERIC, PathSense::REFERENCE };

  graph.for_each_path_matching(&senses, nullptr, nullptr, [&](const path_handle_t& path)
  {
    // Check the path filter if we have one.
    if(path_filter != nullptr && !(*path_filter)(path)) { return; }

    // Find the job for this path.
    nid_t node = graph.get_id(graph.get_handle_of_step(graph.path_begin(path)));
    size_t job = jobs.job(node);
    if(job >= jobs.size()) { return; }

    result[job].push_back(path);
    metadata.add_path(
      graph.get_sense(path),
      graph.get_sample_name(path),
      graph.get_locus_name(path),
      graph.get_haplotype(path),
      graph.get_phase_block(path),
      graph.get_subrange(path),
      job
    );
  });

  return result;
}

// Inserts the selected paths into the GBWT builder.
void
insert_paths(
  const PathHandleGraph& graph,
  const std::vector<path_handle_t>& paths,
  gbwt::GBWTBuilder& builder,
  size_t job_id, bool show_progress)
{
  if(show_progress && paths.size() > 0)
  {
    #pragma omp critical
    {
      std::cerr << "Job " << job_id << ": Inserting " << paths.size() << " paths" << std::endl;
    }
  }
  for(const path_handle_t& path : paths)
  {
    gbwt::vector_type buffer;
    for(handle_t handle : graph.scan_path(path))
    {
      buffer.push_back(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
    }
    builder.insert(buffer, true); // Insert in both orientations.
  }
}

//------------------------------------------------------------------------------

gbwt::GBWT
path_cover_gbwt(
  const PathHandleGraph& graph,
  const PathCoverParameters& parameters,
  bool include_named_paths,
  const std::function<bool(const path_handle_t&)>* path_filter)
{
  // Sanity checks.
  if(!path_cover_sanity_checks(graph, parameters.num_paths, parameters.context))
  {
    return gbwt::GBWT();
  }

  // GBWT construction parameters.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type node_width = sdsl::bits::length(gbwt::Node::encode(graph.max_node_id(), true));
  size_t size_bound = graph.get_node_count() / std::max(size_t(1), parameters.approximate_num_jobs);

  // Determine GBWT construction jobs.
  ConstructionJobs jobs = gbwt_construction_jobs(graph, size_bound);
  MetadataBuilder metadata;

  // Assign the paths we want to include to construction jobs.
  std::vector<std::vector<path_handle_t>> paths_to_include(jobs.size());
  if(include_named_paths)
  {
    paths_to_include = assign_paths(graph, jobs, metadata, path_filter);
  }

  // Create metadata for the path cover.
  std::vector<std::string> contig_names = jobs.contig_names(graph);
  for(size_t component = 0; component < jobs.size(); component++)
  {
    size_t job = jobs.job_for_component(component);
    if(job >= paths_to_include.size()) { continue; }
    for(size_t i = 0; i < parameters.num_paths; i++)
    {
      metadata.add_haplotype("path_cover_" + std::to_string(i), contig_names[component], 0, 0, job);
    }
  }

  // Create the actual path cover.
  std::vector<gbwt::GBWT> partial_indexes(jobs.size());
  std::vector<std::vector<size_t>> components_per_job = jobs.components_per_job();
  int old_threads = omp_get_max_threads();
  omp_set_num_threads(parameters.parallel_jobs);
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t job = 0; job < jobs.size(); job++)
  {
    gbwt::GBWTBuilder builder(node_width, parameters.batch_size, parameters.sample_interval);

    insert_paths(graph, paths_to_include[job], builder, job, parameters.show_progress);
    for(size_t component : components_per_job[job])
    {
      component_path_cover<SimpleCoverage>(
        graph, builder, jobs.weakly_connected_components[component],
        job, component, parameters
      );
    }

    builder.finish();
    partial_indexes[job] = gbwt::GBWT(builder.index);
  }
  omp_set_num_threads(old_threads);

  // Merge the GBWTs and add metadata.
  if(parameters.show_progress)
  {
    std::cerr << "Merging " << partial_indexes.size() << " partial GBWTs" << std::endl;
  }
  gbwt::GBWT result(partial_indexes);
  result.addMetadata();
  result.metadata = metadata.get_metadata();

  return result;
}

//------------------------------------------------------------------------------

gbwt::GBWT
local_haplotypes(
  const PathHandleGraph& graph,
  const gbwt::GBWT& index,
  const PathCoverParameters& parameters,
  bool include_named_paths,
  const std::function<bool(const path_handle_t&)>* path_filter)
{
  // Sanity checks.
  if(!path_cover_sanity_checks(graph, parameters.num_paths, parameters.context))
  {
    return gbwt::GBWT();
  }

  // GBWT construction parameters.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type node_width = sdsl::bits::length(gbwt::Node::encode(graph.max_node_id(), true));
  size_t size_bound = graph.get_node_count() / std::max(size_t(1), parameters.approximate_num_jobs);

  // Determine GBWT construction jobs.
  ConstructionJobs jobs = gbwt_construction_jobs(graph, size_bound);
  MetadataBuilder metadata;

  // Assign the paths we want to include to construction jobs.
  std::vector<std::vector<path_handle_t>> paths_to_include(jobs.size());
  if(include_named_paths)
  {
    paths_to_include = assign_paths(graph, jobs, metadata, path_filter);
  }

  // Create metadata for the path cover.
  std::vector<std::string> contig_names = jobs.contig_names(graph);
  for(size_t component = 0; component < jobs.size(); component++)
  {
    size_t job = jobs.job_for_component(component);
    if(job >= paths_to_include.size()) { continue; }
    for(size_t i = 0; i < parameters.num_paths; i++)
    {
      metadata.add_haplotype("path_cover_" + std::to_string(i), contig_names[component], 0, 0, job);
    }
  }

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
    if(parameters.show_progress)
    {
      std::cerr << "Building a temporary GBWTGraph" << std::endl;
    }
    created_gbwt_graph = GBWTGraph(index, graph);
    gbwt_graph = &created_gbwt_graph;
  }

  // Create the actual path cover.
  std::vector<gbwt::GBWT> partial_indexes(jobs.size());
  std::vector<std::vector<size_t>> components_per_job = jobs.components_per_job();
  int old_threads = omp_get_max_threads();
  omp_set_num_threads(parameters.parallel_jobs);
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t job = 0; job < jobs.size(); job++)
  {
    gbwt::GBWTBuilder builder(node_width, parameters.batch_size, parameters.sample_interval);

    insert_paths(graph, paths_to_include[job], builder, job, parameters.show_progress);
    for(size_t component : components_per_job[job])
    {
      if(component_path_cover<LocalHaplotypes>(
        *gbwt_graph, builder, jobs.weakly_connected_components[component],
        job, component, parameters
      )) { continue; }
      component_path_cover<SimpleCoverage>(
        graph, builder, jobs.weakly_connected_components[component],
        job, component, parameters
      );
    }

    builder.finish();
    partial_indexes[job] = gbwt::GBWT(builder.index);
  }
  omp_set_num_threads(old_threads);

  // Merge the GBWTs and add metadata.
  if(parameters.show_progress)
  {
    std::cerr << "Merging " << partial_indexes.size() << " partial GBWTs" << std::endl;
  }
  gbwt::GBWT result(partial_indexes);
  result.addMetadata();
  result.metadata = metadata.get_metadata();

  return result;
}

//------------------------------------------------------------------------------

size_t
augment_gbwt(
  const HandleGraph& graph,
  gbwt::DynamicGBWT& index,
  const PathCoverParameters& parameters)
{
  // Sanity checks.
  if(!path_cover_sanity_checks(graph, parameters.num_paths, parameters.context)) { return 0; }
  if(!(index.bidirectional()))
  {
    std::cerr << "augment_gbwt(): The GBWT index must be bidirectional" << std::endl;
    return 0;
  }
  if(!(index.hasMetadata()))
  {
    std::cerr << "augment_gbwt(): The GBWT index must have metadata" << std::endl;
    return 0;
  }
  if(!(index.empty()))
  {
    if(!(index.metadata.hasSampleNames() && index.metadata.hasContigNames() && index.metadata.hasPathNames()))
    {
      std::cerr << "augment_gbwt(): The GBWT index must have sample, contig, and path names" << std::endl;
      return 0;
    }
  }

  // GBWT construction parameters.
  gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  gbwt::size_type node_width = sdsl::bits::length(gbwt::Node::encode(graph.max_node_id(), true));
  gbwt::GBWTBuilder builder(node_width, parameters.batch_size, parameters.sample_interval);
  builder.swapIndex(index);
  MetadataBuilder metadata(builder.index.metadata);

  // Find weakly connected components, ignoring the direction of the edges.
  std::vector<std::vector<nid_t>> components = weakly_connected_components(graph);

  // Handle each component separately, but only if there are no GBWT paths in it.
  size_t augmented_components = 0;
  for(size_t component = 0; component < components.size(); component++)
  {
    bool has_paths = false;
    for(nid_t nid : components[component])
    {
      gbwt::node_type node = gbwt::Node::encode(nid, false);
      if(builder.index.contains(node) && !(builder.index.empty(node)))
      {
        has_paths = true;
        break;
      }
    }
    if(has_paths) { continue; }
    component_path_cover<SimpleCoverage>(
      graph, builder, components[component],
      0, component, parameters
    );
    for(size_t i = 0; i < parameters.num_paths; i++)
    {
      metadata.add_haplotype("path_cover_" + std::to_string(i), "component_" + std::to_string(component), 0, 0, 0);
    }
    augmented_components++;
  }

  // Finish the construction and add the updated metadata.
  builder.finish();
  builder.index.metadata = metadata.get_metadata();
  builder.swapIndex(index);

  return augmented_components;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
