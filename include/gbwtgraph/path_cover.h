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

constexpr size_t LOCAL_HAPLOTYPES_DEFAULT_N = 64;
constexpr size_t PATH_COVER_DEFAULT_N       = 16;
constexpr size_t PATH_COVER_DEFAULT_K       = 4;
constexpr size_t PATH_COVER_MIN_K           = 2;
constexpr size_t PATH_COVER_DEFAULT_JOBS    = 32;

//------------------------------------------------------------------------------

// TODO: Remove this when vg no longer uses it.
/*
  Store the named paths from the given graph into the given GBWT builder.
  Generic named paths will go under a sample with the special
  REFERENCE_PATH_SAMPLE_NAME name. Creates new contigs with the appropriate
  names, or re-uses any found in the metadata.
  
  Skips empty paths.
  
  If the given filter function is set, and returns false for a path, that path
  is not added.

  NOTE: This function will be deprecated soon. Consider using `assign_paths` and
  `insert_paths` instead.
*/
void
store_named_paths(gbwt::GBWTBuilder& builder, const PathHandleGraph& graph, const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

// TODO: Remove this when vg no longer uses it.
/*
  Store paths from the given graph into the given GBWT builder. Generic named
  paths will go under a sample with the special REFERENCE_PATH_SAMPLE_NAME
  name. Creates new contigs with the appropriate names, or re-uses any found in
  the metadata.
  
  Only stores paths with the senses in the given set.
  
  Skips empty paths.
  
  If the given filter function is set, and returns false for a path, that path
  is not added.

  NOTE: This function will be deprecated soon. Consider using `assign_paths` and
  `insert_paths` instead.
*/
void
store_paths(gbwt::GBWTBuilder& builder, const PathHandleGraph& graph, const std::unordered_set<PathSense>& senses, const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

//------------------------------------------------------------------------------

struct PathCoverParameters
{
  // Number of paths per component.
  size_t num_paths = PATH_COVER_DEFAULT_N;

  // Length of the window in nodes.
  size_t context = PATH_COVER_DEFAULT_K;

  // Batch size for GBWT construction (in nodes).
  gbwt::size_type batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE;

  // Sample interval for GBWT construction (in nodes).
  gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;

  // Approximate number of construction jobs.
  size_t approximate_num_jobs = PATH_COVER_DEFAULT_JOBS;

  // Number of parallel GBWT construction jobs.
  size_t parallel_jobs = 1;

  // Show progress information.
  bool show_progress = false;
};

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

    - If the component is not a DAG, we start from an arbitrary node with minimal
      coverage and extend in both directions.

    - In a DAG, we start from the head node with the lowest coverage so far.

    - We stop when the length of the path reaches the size of the component.

    - The length of the window is in nodes instead of base pairs. We expect a sparse graph,
      where the nodes between variants are long.

    - If the component is not a DAG and the path is shorter than k - 1 nodes, we consider
      the coverage of individual nodes instead of windows.

    - When determining window coverage, we consider the window equivalent to its reverse
      complement.

  If include_named_paths is set, named paths from the graph will be stored. If a path_filter
  is supplied, only paths matching it will be stored.

  NOTE: The resulting GBWT does not have any samples assigned as reference samples. If
  reference paths were included, the REFERENCE_SAMPLE_LIST_GBWT_TAG tag should be set
  separately when appropriate.
*/
gbwt::GBWT path_cover_gbwt(
  const PathHandleGraph& graph,
  const PathCoverParameters& parameters,
  bool include_named_paths = false,
  const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

//------------------------------------------------------------------------------

/*
  As above, but build the path cover using local haplotypes. Each window must be consistent
  with the haplotypes in the GBWT index. Instead of choosing the extension with the lowest
  coverage, this algorithm maximizes the D'Hondt ratio:

    true_coverage / (selected_coverage + 1).

  In graph components without haplotypes in the GBWT index, this algorithm will revert to
  the regular path cover algorithm.
  
  If include_named_paths is set, named paths from the graph will be stored. If a path_filter
  is supplied, only paths matching it will be stored.

  NOTE: The resulting GBWT does not have any samples assigned as reference samples. If
  reference paths were included, the REFERENCE_SAMPLE_LIST_GBWT_TAG tag should be set
  separately when appropriate.
*/
gbwt::GBWT local_haplotypes(
  const PathHandleGraph& graph,
  const gbwt::GBWT& index,
  const PathCoverParameters& parameters,
  bool include_named_paths = false,
  const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

//------------------------------------------------------------------------------

/*
  Augment the given GBWT index with a path cover of the components that do not have any
  paths. This will add n new samples but no new haplotypes into the metadata. If the
  metadata contains sample/contig names, the path cover will use names path_cover_i and
  component_i. Returns the number of components that received a path cover.

  This algorithm is effectively a single construction job, as it takes an existing
  dynamic index. We further assume that most components are already covered by the GBWT
  and we only need to augment a few small ones.
*/
size_t augment_gbwt(
  const HandleGraph& graph,
  gbwt::DynamicGBWT& index,
  const PathCoverParameters& parameters);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_PATH_COVER_H
