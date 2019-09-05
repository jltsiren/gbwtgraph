#include <iostream>
#include <fstream>

#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

const std::string tool_name = "GFA to GBWTGraph";

void printUsage(int exit_code = EXIT_SUCCESS);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  // Parse command line options.
  int c = 0;
  while((c = getopt(argc, argv, "")) != -1)
  {
    switch(c)
    {
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  // Check command line options.
  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  std::string base_name = argv[optind]; optind++;

  // Initial output.
  Version::print(std::cout, tool_name);
  gbwt::printHeader("Base name"); std::cout << base_name << std::endl;
  std::cout << std::endl;

  double start = gbwt::readTimer();

  std::cout << "Parsing GFA and building GBWT" << std::endl;
  auto results = gfa_to_gbwt(base_name + GFA_EXTENSION);
  if(results.first.get() == nullptr || results.second.get() == nullptr)
  {
    std::cerr << "gfa2gbwt: Construction failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << std::endl;
  gbwt::printStatistics(*(results.first), base_name);

  std::cout << "Serializing GBWT" << std::endl;
  {
    std::string gbwt_name = base_name + gbwt::GBWT::EXTENSION;
    std::ofstream out(gbwt_name, std::ios_base::binary);
    if(!out)
    {
      std::cerr << "gfa2gbwt: Cannot open file " << gbwt_name << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    results.first->serialize(out);
    out.close();
  }
  std::cout << std::endl;

  std::cout << "Building GBWTGraph" << std::endl;
  GBWTGraph graph(*(results.first), *(results.second));
  std::cout << std::endl;

  std::cout << "Serializing GBWTGraph" << std::endl;
  {
    std::string graph_name = base_name + GBWTGraph::EXTENSION;
    std::ofstream out(graph_name, std::ios_base::binary);
    if(!out)
    {
      std::cerr << "gfa2gbwt: Cannot open file " << graph_name << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    graph.serialize(out);
    out.close();
  }
  std::cout << std::endl;

  double seconds = gbwt::readTimer() - start;

  std::cout << "GBWTGraph built in " << seconds << " seconds" << std::endl;
  std::cout << "Memory usage " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: gfa2gbwt base_name" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
