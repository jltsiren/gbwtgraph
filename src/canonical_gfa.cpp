#include <getopt.h>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>

#include <handlegraph/algorithms/canonical_gfa.hpp>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

const std::string tool_name = "Canonical GFA";

struct Config
{
  Config(int argc, char** argv);

  bool integer_ids = false;
  bool handlegraph_algorithm = false;

  std::string filename;
};

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  Config config(argc, argv);

  // Load the graph.
  GBZ gbz;
  sdsl::simple_sds::load_from(gbz, config.filename);

  // Print the GFA.
  if(config.integer_ids && !config.handlegraph_algorithm)
  {
    gbwt_to_canonical_gfa(gbz.graph, std::cout);
  }
  else
  {
    handlegraph::algorithms::canonical_gfa(gbz.graph, std::cout, config.integer_ids);
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: canonical_gfa [options] graph.gbz" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -i, --integer-ids   order the nodes by integer ids" << std::endl;
  std::cerr << "  -H, --handlegraph   always use the libhandlegraph algorithm" << std::endl;
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
    { "integer-ids", no_argument, 0, 'i' },
    { "handlegraph", no_argument, 0, 'H' },
    { 0, 0, 0, 0 }
  };
  
  // Process options.
  while((c = getopt_long(argc, argv, "iH", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'i':
      this->integer_ids = true;
      break;
    case 'H':
      this->handlegraph_algorithm = true;
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
