#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

#include <getopt.h>
#include <unistd.h>

#include <omp.h>

#include <gbwtgraph/gbz.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

/*
  This tool finds all exact occurrences of the pattern in the haplotypes.
  It iterates over non-redundant haplotype windows. If a window contains the
  pattern, the tool outputs the following tab-separated fields:

  1. The path using GFA walk notation.
  2. The sequence as comma-separated node labels, excluding left/right flanks.
  3. Number of haplotypes containing the path.
  4. Length of the left flank.
  5. Length of the right flank.
  6. Reference coordinates for each node, if requested.

  Left (right) flank is the prefix (suffix) of the first (last) node not
  included in the window. If the path is a single node, the flanks are always
  empty.

  TODO: Make this a library function.
  TODO: How to avoid reporting the same window in both orientations if the pattern is palindromic?
*/

const std::string tool_name = "GBZ Haplotype Search";

struct Config
{
  Config(int argc, char** argv);

  std::string graph_name;
  std::string pattern;
  std::string reference_sample;

  bool progress;
  size_t threads;

  constexpr static size_t BUFFER_SIZE = 1024; // Lines.
};

// Returns a mapping GBWT node id -> (GBWT sequence id, sequence offset).
std::unordered_map<gbwt::node_type, std::pair<std::uint32_t, std::uint32_t>> reference_labels(const GBZ& graph, const Config& config);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  double start = gbwt::readTimer();
  Config config(argc, argv);
  omp_set_num_threads(config.threads);

  GBZ graph;
  if(config.progress)
  {
    std::cerr << "Loading the graph from " << config.graph_name << std::endl;
  }
  sdsl::simple_sds::load_from(graph, config.graph_name);

  std::unordered_map<gbwt::node_type, std::pair<std::uint32_t, std::uint32_t>> ref_labels = reference_labels(graph, config);
  auto create_label = [&](gbwt::size_type seq_id, gbwt::size_type offset, bool is_reverse) -> std::string
  {
    const gbwt::Metadata& metadata = graph.index.metadata;
    gbwt::FullPathName path_name = metadata.fullPath(seq_id);
    std::string result;
    result.push_back(is_reverse ? '<' : '>');
    result.append(path_name.sample_name);
    result.push_back('#');
    result.append(std::to_string(path_name.haplotype));
    result.push_back('#');
    result.append(path_name.contig_name);
    result.push_back('#');
    result.append(std::to_string(path_name.offset + offset));
    return result;
  };
  auto get_label = [&](handle_t handle) -> std::string
  {
    gbwt::node_type node = GBWTGraph::handle_to_node(handle);
    auto it = ref_labels.find(node);
    if(it != ref_labels.end())
    {
      return create_label(it->second.first, it->second.second, false);
    }
    it = ref_labels.find(gbwt::Node::reverse(node));
    if(it != ref_labels.end())
    {
      return create_label(it->second.first, it->second.second, true);
    }
    return ">-";
  };

  if(config.progress)
  {
    std::cerr << "Searching the haplotypes for " << config.pattern << std::endl;
  }
  std::vector<std::vector<std::string>> buffers(config.threads);
  size_t windows = 0;
  auto flush = [&](size_t thread_id)
  {
    std::vector<std::string>& buffer = buffers[thread_id];
    for(const std::string& line : buffer)
    {
      std::cout.write(line.data(), line.length());
      std::cout.put('\n');
    }
    windows += buffer.size();
    buffer.clear();
  };

  // FIXME: update
  for_each_nonredundant_window
  (
    graph.graph, config.pattern.length(),
    [&](const std::vector<handle_t>& path, const std::string& sequence)
    {
      if(sequence.find(config.pattern) == std::string::npos) { return; }

      // Path in the GFA walk notation.
      std::string line;
      for(handle_t handle : path)
      {
        gbwt::node_type node = GBWTGraph::handle_to_node(handle);
        line.push_back(gbwt::Node::is_reverse(node) ? '<' : '>');
        line += std::to_string(gbwt::Node::id(node));
      }

      // Path as a comma-separated list of oriented node labels.
      // Also determine the lengths of the left and right flanks.
      line.push_back('\t');
      size_t left_flank = 0, right_flank = 0;
      if(path.size() > 1)
      {
        size_t first_node_length = graph.graph.get_length(path.front());
        if(first_node_length >= config.pattern.length())
        {
          left_flank = first_node_length - config.pattern.length() + 1;
        }
        size_t offset = first_node_length - left_flank;
        line += std::string_view(sequence).substr(0, offset);
        for(size_t i = 1; i < path.size(); i++)
        {
          size_t node_length = graph.graph.get_length(path[i]);
          line.push_back(',');
          line += std::string_view(sequence).substr(offset, node_length);
          offset += node_length;
        }
        right_flank = offset - sequence.length();
      }
      else
      {
        line += sequence;
      }

      // Number of haplotypes containing the path.
      // TODO: for_each_haplotype_window should expose the search state.
      line.push_back('\t');
      gbwt::SearchState state = graph.index.find(GBWTGraph::handle_to_node(path.front()));
      for(size_t i = 1; i < path.size(); i++)
      {
        state = graph.index.extend(state, GBWTGraph::handle_to_node(path[i]));
      }
      line += std::to_string(state.size());

      // Length of the left and right flanks.
      line.push_back('\t');
      line += std::to_string(left_flank);
      line.push_back('\t');
      line += std::to_string(right_flank);

      // Reference coordinates for each node, if requested.
      // If the orientation of the final node is the opposite to the reference,
      // the coordinates are for the other end of the node. To determine the 
      // true coordinates, we need to use the previous field.
      if(!ref_labels.empty())
      {
        line.push_back('\t');
        for(handle_t handle : path)
        {
          line += get_label(handle);
        }
      }

      size_t thread_id = omp_get_thread_num();
      buffers[thread_id].push_back(std::move(line));
      if(buffers[thread_id].size() >= config.BUFFER_SIZE)
      {
        #pragma omp critical
        flush(thread_id);
      }
    }, true
  );
  for(size_t thread_id = 0; thread_id < config.threads; thread_id++)
  {
    flush(thread_id);
  }
  std::cout.flush();

  if(config.progress)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Found a match in " << windows << " windows in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------

size_t
default_threads()
{
  return std::thread::hardware_concurrency();
}

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: gbz_find [options] graph.gbz pattern > output" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -r, --reference NAME  output coordinates in this reference sample" << std::endl;
  std::cerr << "  -t, --threads N       use N threads (default and max: " << default_threads() << ")" << std::endl;
  std::cerr << "  -p, --progress        print progress information" << std::endl;
  std::cerr << "  -h, --help            print this help message" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

Config::Config(int argc, char** argv) :
  progress(false), threads(default_threads())
{
  if(argc < 2) { printUsage(EXIT_SUCCESS); }

  // Data for `getopt_long()`.
  int c = 0, option_index = 0;
  option long_options[] =
  {
    { "reference", required_argument, 0, 'r' },
    { "threads", required_argument, 0, 't' },
    { "progress", no_argument, 0, 'p' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 }
  };

  // Process options.
  while((c = getopt_long(argc, argv, "r:t:ph", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'r':
      this->reference_sample = optarg;
      break;
    case 't':
      try { this->threads = std::stoul(optarg); }
      catch(std::exception& e)
      {
        std::cerr << "Cannot parse --threads " << optarg << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'p':
      this->progress = true;
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
  if(optind + 2 != argc) { printUsage(EXIT_FAILURE); }
  this->graph_name = argv[optind]; optind++;
  this->pattern = argv[optind]; optind++;

  // Sanity checks.
  if(this->threads < 1 || this->threads > default_threads())
  {
    std::cerr << "Invalid number of threads: " << this->threads << " (must be between 1 and " << default_threads() << ")" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(this->pattern.empty())
  {
    std::cerr << "Pattern cannot be empty" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//------------------------------------------------------------------------------

std::unordered_map<gbwt::node_type, std::pair<std::uint32_t, std::uint32_t>>
reference_labels(const GBZ& graph, const Config& config)
{
  std::unordered_map<gbwt::node_type, std::pair<std::uint32_t, std::uint32_t>> result;
  if(config.reference_sample.empty()) { return result; }
  if(config.progress)
  {
    std::cerr << "Determining coordinates in reference sample " << config.reference_sample << std::endl;
  }

  const gbwt::Metadata& metadata = graph.index.metadata;
  gbwt::size_type sample = metadata.sample(config.reference_sample);
  std::vector<gbwt::size_type> paths = metadata.pathsForSample(sample);
  if(paths.empty())
  {
    if(config.progress)
    {
      std::cerr << "No paths found for the reference sample" << std::endl;
    }
    return result;
  }

  for(gbwt::size_type path_id : paths)
  {
    gbwt::vector_type path = graph.index.extract(gbwt::Path::encode(path_id, false));
    std::uint32_t offset = 0;
    for(gbwt::node_type node : path)
    {
      // TODO: Should we print a warning?
      if(result.find(node) != result.end()) { continue; }
      result[node] = std::pair<std::uint32_t, std::uint32_t>(path_id, offset);
      offset += graph.graph.get_length(GBWTGraph::node_to_handle(node));
    }
  }

  return result;
}

//------------------------------------------------------------------------------
