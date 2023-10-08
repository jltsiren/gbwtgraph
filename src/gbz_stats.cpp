#include <atomic>
#include <iostream>
#include <fstream>
#include <map>

#include <getopt.h>
#include <unistd.h>

#include <gbwtgraph/gbz.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

const std::string tool_name = "GBZ Statistics";

struct Config
{
  Config(int argc, char** argv);

  bool graph = false;
  bool gbwt = false;

  size_t window_length = 0;

  bool node_degrees = false;
  bool node_visits = false;

  bool record_bytes = false;
  bool record_runs = false;

  std::string filename;
};

void print_distribution(const std::map<size_t, size_t>& distribution, const std::string& header1, const std::string& header2);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  Config config(argc, argv);
  GBZ gbz;
  sdsl::simple_sds::load_from(gbz, config.filename);

  if(config.graph)
  {
    std::cout << "Nodes\t" << gbz.graph.get_node_count() << std::endl;
    std::cout << "Edges\t" << gbz.graph.get_edge_count() << std::endl;
    size_t total_length = 0;
    gbz.graph.for_each_handle([&](const handle_t& handle)
    {
      total_length += gbz.graph.get_length(handle);
    });
    std::cout << "Sequence\t" << total_length << std::endl;
  }

  if(config.gbwt)
  {
    gbwt::printStatistics(gbz.index, config.filename, std::cout);
  }

  if(config.window_length > 0)
  {
    std::atomic<size_t> fragments(0), total_length(0);
    for_each_haplotype_window(gbz.graph, config.window_length,
      [&](const std::vector<handle_t>&, const std::string& sequence) -> bool
      {
        fragments++; total_length += sequence.length();
        return true;
      }, true);
      std::cout << "Fragments\t" << fragments << std::endl;
      std::cout << "Total length\t" << total_length << std::endl;
  }

  if(config.node_degrees)
  {
    std::map<size_t, size_t> distribution;
    gbz.graph.for_each_handle([&](const handle_t& handle)
    {
      size_t degree = gbz.graph.get_degree(handle, false) + gbz.graph.get_degree(handle, true);
      distribution[degree]++;
    });
    print_distribution(distribution, "Degree", "Nodes");
  }

  if(config.node_visits)
  {
    std::map<size_t, size_t> distribution;
    for(gbwt::node_type node = gbz.index.firstNode(); node < gbz.index.sigma(); node += 2)
    {
      distribution[gbz.index.nodeSize(node)]++;
    }
    print_distribution(distribution, "Visits", "Nodes");
  }

  if(config.record_bytes)
  {
    std::map<size_t, size_t> distribution;
    for(gbwt::node_type node = gbz.index.firstNode(); node < gbz.index.sigma(); node++)
    {
      std::pair<gbwt::size_type, gbwt::size_type> range = gbz.index.bwt.getRange(gbz.index.toComp(node));
      distribution[range.second - range.first]++;
    }
    print_distribution(distribution, "Bytes", "Records");
  }

  if(config.record_runs)
  {
    std::map<size_t, size_t> distribution;
    for(gbwt::node_type node = gbz.index.firstNode(); node < gbz.index.sigma(); node++)
    {
      distribution[gbz.index.record(node).runs().first]++;
    }
    print_distribution(distribution, "Runs", "Records");
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: gbz_stats [options] graph.gbz" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Overall statistics:" << std::endl;
  std::cerr << "  -g, --graph         Graph statistics" << std::endl;
  std::cerr << "  -i, --gbwt          GBWT index statistics" << std::endl;
  std::cerr << "  -w, --windows N     N bp haplotype windows" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Nodes:" << std::endl;
  std::cerr << "  -d, --node-degrees  Node degree distribution" << std::endl;
  std::cerr << "  -v, --node-visits   Node visit distribution" << std::endl;
  std::cerr << std::endl;
  std::cerr << "GBWT records:" << std::endl;
  std::cerr << "  -b, --record-bytes  Record size distribution" << std::endl;
  std::cerr << "  -r, --record-runs   Run count distribution" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------

Config::Config(int argc, char** argv)
{
  if(argc < 2) { printUsage(EXIT_SUCCESS); }

  // Data for `getopt_long()`.
  int c = 0, option_index = 0;
  option long_options[] =
  {
    { "graph", no_argument, 0, 'g' },
    { "gbwt", no_argument, 0, 'i' },
    { "windows", required_argument, 0, 'w' },
    { "node-degrees", no_argument, 0, 'd' },
    { "node-visits", no_argument, 0, 'v' },
    { "record-bytes", no_argument, 0, 'b' },
    { "record-runs", no_argument, 0, 'r' },
  };

  // Process options.
  while((c = getopt_long(argc, argv, "giw:dvbr", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'g':
      this->graph = true;
      break;
    case 'i':
      this->gbwt = true;
      break;
    case 'w':
      try { this->window_length = std::stoul(optarg); }
      catch(std::exception& e)
      {
        std::cerr << "Cannot parse --window " << optarg << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;

    case 'd':
      this->node_degrees = true;
      break;
    case 'v':
      this->node_visits = true;
      break;

    case 'b':
      this->record_bytes = true;
      break;
    case 'r':
      this->record_runs = true;
      break;

    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  // Sanity checks.
  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  this->filename = argv[optind]; optind++;
}

//------------------------------------------------------------------------------

void
print_distribution(const std::map<size_t, size_t>& distribution, const std::string& header1, const std::string& header2)
{
  std::cout << header1 << "\t" << header2 << std::endl;
  for(auto iter = distribution.begin(); iter != distribution.end(); ++iter)
  {
    std::cout << iter->first << "\t" << iter->second << std::endl;
  }
}

//------------------------------------------------------------------------------
