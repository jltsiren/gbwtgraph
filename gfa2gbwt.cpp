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

// FIXME new parameters: chopping limit
// FIXME use long options

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  GFAParsingParameters parameters;
  parameters.show_progress = true;

  // Parse command line options.
  int c = 0;
  while((c = getopt(argc, argv, "r:f:")) != -1)
  {
    switch(c)
    {
    case 'r':
      parameters.path_name_regex = optarg;
      break;
    case 'f':
      parameters.path_name_fields = optarg;
      break;
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
  std::cout << "Path name regex: " << parameters.path_name_regex << std::endl;
  std::cout << "Path name fields: " << parameters.path_name_fields << std::endl;
  auto result = gfa_to_gbwt(base_name + GFA_EXTENSION, parameters);
  if(result.first.get() == nullptr || result.second.get() == nullptr)
  {
    std::cerr << "gfa2gbwt: Construction failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "Serializing GBWT" << std::endl;
  {
    std::string gbwt_name = base_name + gbwt::GBWT::EXTENSION;
    std::ofstream out(gbwt_name, std::ios_base::binary);
    if(!out)
    {
      std::cerr << "gfa2gbwt: Cannot open file " << gbwt_name << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    result.first->serialize(out);
    out.close();
  }

  std::cout << "Building GBWTGraph" << std::endl;
  GBWTGraph graph(*(result.first), *(result.second));

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

  const std::unordered_map<std::string, std::pair<nid_t, nid_t>>& translation = result.second->segment_translation;
  if(!(translation.empty()))
  {
    std::cout << "Writing translation table" << std::endl;

    // FIXME source for extension
    std::string translation_name = base_name + ".trans";
    std::ofstream out(translation_name, std::ios_base::binary);
    for(auto iter = translation.begin(); iter != translation.end(); ++iter)
    {
      out << "T\t" << iter->first << "\t" << iter->second.first;
      for(nid_t i = iter->second.first + 1; i < iter->second.second; i++)
      {
        out << "," << i;
      }
      out << "\n";
    }
    out.close();
  }
  std::cout << std::endl;

  gbwt::printStatistics(*(result.first), base_name);

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

  std::cerr << "Usage: gfa2gbwt [options] base_name" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -r X   parse path names using regex X (default " << GFAParsingParameters::DEFAULT_REGEX << ")" << std::endl;
  std::cerr << "  -f X   map the submatches to fields X (default " << GFAParsingParameters::DEFAULT_FIELDS << ")" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Fields (case insensitive):" << std::endl;
  std::cerr << "  S      sample name" << std::endl;
  std::cerr << "  C      contig name" << std::endl;
  std::cerr << "  H      numerical haplotype identifier" << std::endl;
  std::cerr << "  F      numerical fragment identifier" << std::endl;
  std::cerr << "  other  ignore this field" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
