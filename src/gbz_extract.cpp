#include <iostream>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <getopt.h>
#include <unistd.h>

#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/gbz.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

/*
  This tool extracts sequences from a GBZ graph, one sequence per selected path.
  The sequences are written in the same order as the paths in the graph and
  terminated by a newline character.

  The sequences are extracted either in the forward orientation or in both
  orientations. In the latter case, the reverse complement of each sequence
  follows the original sequence.
*/

const std::string tool_name = "GBZ Sequence Extractor";

struct Config
{
  Config(int argc, char** argv);

  std::string graph_name;
  std::string contig_name;

  // Extract sequences in both orientations.
  bool both_orientations = false;

  // TODO: Block size to assign multiple paths to a thread.

  bool verbose = false;
  size_t threads;
};

// Returns the path identifiers for the selected contig, or for all paths
// if the name is empty.
std::vector<gbwt::size_type> select_paths(const GBZ& graph, const std::string& contig_name, bool verbose);

// Extracts the nucleotide sequence corresponding to the given path id to the
// provided buffer. The buffer is cleared before the extraction.
void extract_path(const GBZ& graph, gbwt::size_type path_id, std::string& buffer);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  Config config(argc, argv);
  GBZ graph;

  double start = gbwt::readTimer();
  if(config.verbose)
  {
    std::cerr << "Loading the graph from " << config.graph_name << std::endl;
  }
  sdsl::simple_sds::load_from(graph, config.graph_name);
  std::vector<gbwt::size_type> paths = select_paths(graph, config.contig_name, config.verbose);

  // Start the extraction.
  if(config.verbose)
  {
    std::cerr << "Extracting " << paths.size() << " sequences";
    if(config.both_orientations) { std::cerr << " in both orientations"; }
    std::cerr << std::endl;
  }
  std::vector<std::thread> threads; threads.reserve(config.threads);
  std::vector<std::string> thread_buffers(config.threads);
  for(size_t to_read = 0; to_read < config.threads && to_read < paths.size(); to_read++)
  {
    threads.emplace_back(extract_path, std::ref(graph), paths[to_read], std::ref(thread_buffers[to_read]));
  }

  // Write the extracted sequences.
  size_t total_length = 0;
  for(size_t to_write = 0; to_write < paths.size(); to_write++)
  {
    // Get the buffer from the thread that was responsible for exracting this path.
    size_t thread_id = to_write % config.threads;
    std::string write_buffer;
    threads[thread_id].join();
    write_buffer.swap(thread_buffers[thread_id]);

    // Relaunch the thread if there are still more paths to extract.
    if(to_write + threads.size() < paths.size())
    {
      threads[thread_id] = std::thread(extract_path,
        std::ref(graph), paths[to_write + threads.size()], std::ref(thread_buffers[thread_id]));
    }

    // And write the extracted sequence.
    std::cout.write(write_buffer.data(), write_buffer.length());
    std::cout.put('\n');
    total_length += write_buffer.length();
    if(config.both_orientations)
    {
      reverse_complement_in_place(write_buffer);
      std::cout.write(write_buffer.data(), write_buffer.length());
      std::cout.put('\n');
      total_length += write_buffer.length();
    }
  }
  if(config.verbose)
  {
    std::cerr << "Total sequence length: " << total_length << " bp" << std::endl;
  }

  if(config.verbose)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Used " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------

size_t
default_threads()
{
  return std::thread::hardware_concurrency();
}

size_t
max_threads()
{
  return 2 * default_threads();
}

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: gbz_extract [options] graph.gbz > output" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Computational parameters:" << std::endl;
  std::cerr << "  -c, --contig NAME        select paths corresponding to contig NAME" << std::endl;
  std::cerr << "  -b, --both-orientations  extract sequences in both orientations" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Other options:" << std::endl;
  std::cerr << "  -t, --threads N          use N threads (default: " << default_threads() << ", max: " << max_threads() << ")" << std::endl;
  std::cerr << "  -v, --verbose            print progress information" << std::endl;
  std::cerr << "  -h, --help               print this help message" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

Config::Config(int argc, char** argv) :
  threads(default_threads())
{
  if(argc < 2) { printUsage(EXIT_SUCCESS); }

  // Data for `getopt_long()`.
  int c = 0, option_index = 0;
  option long_options[] =
  {
    { "contig", required_argument, 0, 'c' },
    { "both-orientations", no_argument, 0, 'b' },
    { "threads", required_argument, 0, 't' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 }
  };

  // Process options.
  while((c = getopt_long(argc, argv, "c:bt:vh", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'c':
      this->contig_name = optarg;
      break;
    case 'b':
      this->both_orientations = true;
      break;

    case 't':
      try { this->threads = std::stoul(optarg); }
      catch(std::exception& e)
      {
        std::cerr << "Cannot parse --threads " << optarg << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'v':
      this->verbose = true;
      break;
    case 'h':
      printUsage(EXIT_SUCCESS);
      break;

    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  // Positional arguments.
  if(optind +1 != argc) { printUsage(EXIT_FAILURE); }
  this->graph_name = argv[optind]; optind++;

  // Sanity checks.
  if(this->threads < 1 || this->threads > max_threads())
  {
    std::cerr << "Invalid number of threads: " << this->threads << " (must be between 1 and " << max_threads() << ")" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//------------------------------------------------------------------------------

size_t
component_for_path(const GBZ& graph, gbwt::size_type path_id, const std::unordered_map<nid_t, size_t>& node_to_component)
{
  gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
  nid_t node = gbwt::Node::id(graph.index.start(seq_id).first);
  auto iter = node_to_component.find(node);
  return (iter == node_to_component.end() ? std::numeric_limits<size_t>::max() : iter->second);
}

std::vector<gbwt::size_type>
select_paths(const GBZ& graph, const std::string& contig_name, bool verbose)
{
  std::vector<gbwt::size_type> result;
  if(contig_name.empty())
  {
    if(verbose)
    {
      std::cerr << "Selecting all " << (graph.index.sequences() / 2) << " paths" << std::endl;
    }
    result.reserve(graph.index.sequences() / 2);
    for(gbwt::size_type i = 0; i < graph.index.sequences() / 2; i++) { result.push_back(i); }
    return result;
  }

  // Find reference paths for the contig.
  if(verbose)
  {
    std::cerr << "Selecting reference paths for contig " << contig_name << std::endl;
  }
  const gbwt::Metadata& metadata = graph.index.metadata;
  if(!(graph.index.hasMetadata()) || !(metadata.hasContigNames()) || !(metadata.hasPathNames()))
  {
    std::cerr << "Error: The graph does not have metadata with contig names and path names" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  gbwt::size_type contig_id = metadata.contig(contig_name);
  if(contig_id >= metadata.contigs())
  {
    std::cerr << "Error: Contig " << contig_name << " not found in the graph" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::vector<gbwt::size_type> ref_paths = metadata.pathsForContig(contig_id);
  if(ref_paths.empty())
  {
    std::cerr << "Warning: Contig " << contig_name << " does not have any paths" << std::endl;
    return result;
  }
  if(verbose)
  {
    std::cerr << "Found " << ref_paths.size() << " reference paths" << std::endl;
  }

  // Find the graph components for the selected reference paths.
  if(verbose)
  {
    std::cerr << "Finding graph components for contig " << contig_name << std::endl;
  }
  auto components = weakly_connected_components(graph.graph);
  std::unordered_map<nid_t, size_t> node_to_component;
  for(size_t i = 0; i < components.size(); i++)
  {
    for(nid_t node : components[i])
    {
      node_to_component[node] = i;
    }
  }
  components = std::vector<std::vector<nid_t>>();
  std::unordered_set<size_t> selected_components;
  for(gbwt::size_type ref_path : ref_paths)
  {
    size_t component = component_for_path(graph, ref_path, node_to_component);
    if(component != std::numeric_limits<size_t>::max())
    {
      selected_components.insert(component);
    }
  }
  if(verbose)
  {
    std::cerr << "Found " << selected_components.size() << " components" << std::endl;
  }

  // Now select the paths.
  if(verbose)
  {
    std::cerr << "Selecting paths for contig " << contig_name << std::endl;
  }
  for(gbwt::size_type path_id = 0; path_id < metadata.paths(); path_id++)
  {
    size_t component = component_for_path(graph, path_id, node_to_component);
    if(selected_components.find(component) != selected_components.end())
    {
      result.push_back(path_id);
    }
  }
  if(verbose)
  {
    std::cerr << "Found " << result.size() << " paths" << std::endl;
  }

  return result;
}

//------------------------------------------------------------------------------

void
extract_path(const GBZ& graph, gbwt::size_type path_id, std::string& buffer)
{
  buffer.clear();
  gbwt::size_type seq_id = gbwt::Path::encode(path_id, false);
  gbwt::edge_type curr = graph.index.start(seq_id);
  while(curr.first != gbwt::ENDMARKER)
  {
    handle_t handle = GBWTGraph::node_to_handle(curr.first);
    view_type view = graph.graph.get_sequence_view(handle);
    buffer.append(view.first, view.second);
    curr = graph.index.LF(curr);
  }
}

//------------------------------------------------------------------------------
