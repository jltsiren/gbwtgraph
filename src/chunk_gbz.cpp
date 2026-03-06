#include <iostream>
#include <fstream>

#include <getopt.h>
#include <unistd.h>

#include <gbwtgraph/algorithms.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

const std::string tool_name = "Chunk GBZ";
const std::string default_prefix = "chunk";

struct Config
{
  Config(int argc, char** argv);

  std::string input;
  std::string output_prefix = default_prefix;

  ChunkParameters params;
};

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  double start = gbwt::readTimer();
  Config config(argc, argv);

  GBZ gbz;
  if(config.params.verbose)
  {
    std::cerr << "Loading the graph from " << config.input << std::endl;
  }
  sdsl::simple_sds::load_from(gbz, config.input);

  std::vector<std::pair<std::string, GBZ>> chunks = chunk_graph(gbz, config.params);

  if(config.params.verbose)
  {
    std::cerr << "Writing the output files" << std::endl;
  }
  for(size_t i = 0; i < chunks.size(); i++)
  {
    std::string filename = config.output_prefix + "_" + std::to_string(i) + "_" + chunks[i].first + ".gbz";
    sdsl::simple_sds::serialize_to(chunks[i].second, filename);
  }

  if(config.params.verbose)
  {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Used " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
  }
  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: chunk_gbz [options] graph.gbz" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Chunks the given graph into weakly connected components." << std::endl;
  std::cerr << "Output files are named <prefix>_<index>_<contig>.gbz." << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -c, --contig X  only extract components with this contig name" << std::endl;
  std::cerr << "  -j, --jobs N    run this many parallel jobs (default: 1)" << std::endl;
  std::cerr << "  -p, --prefix X  use this output prefix (default: " << default_prefix << ")" << std::endl;
  std::cerr << "  -v, --verbose   print progress information" << std::endl;
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
    { "contig", required_argument, nullptr, 'c' },
    { "jobs", required_argument, nullptr, 'j' },
    { "prefix", required_argument, nullptr, 'p' },
    { "verbose", no_argument, nullptr, 'v' },
    { 0, 0, 0, 0 }
  };

  // Process options.
  while((c = getopt_long(argc, argv, "c:j:p:v", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'c':
      this->params.contig_name = optarg;
      break;
    case 'j':
      try { this->params.parallel_jobs = std::stoul(optarg); }
      catch(std::exception& e)
      {
        std::cerr << "Cannot parse --jobs " << optarg << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'p':
      this->output_prefix = optarg;
      break;
    case 'v':
      this->params.verbose = true;
      break;

    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  // Sanity checks.
  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  this->input = argv[optind]; optind++;
}

//------------------------------------------------------------------------------
