#include <iostream>
#include <fstream>
#include <string>

#include <getopt.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

enum input_type { input_gfa, input_gbz, input_gg };
enum output_type { output_gfa, output_gbz, output_gg };

struct Config
{
  Config(int argc, char** argv);

  GFAParsingParameters parameters;
  std::string base_name;

  input_type input = input_gfa;
  output_type output = output_gg;

  bool translation = false;
  bool show_progress = false;
};

const std::string tool_name = "GFA to GBWTGraph";

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  double start = gbwt::readTimer();
  Config config(argc, argv);

  // Initial output.
  if(config.show_progress)
  {
    Version::print(std::cerr, tool_name);
    gbwt::printHeader("Base name", std::cerr) << config.base_name << std::endl;
    std::cerr << std::endl;
  }

  // This is the data we are using.
  gbwt::GBWT index;
  GBWTGraph graph;
  std::unordered_map<std::string, std::pair<nid_t, nid_t>> translation;

  // TODO this should depend on input/output
  if(config.input == input_gfa)
  {
    if(config.show_progress)
    {
      std::cerr << "Parsing GFA and building GBWT" << std::endl;
      std::cerr << "Path name regex: " << config.parameters.path_name_regex << std::endl;
      std::cerr << "Path name fields: " << config.parameters.path_name_fields << std::endl;
    }
    auto result = gfa_to_gbwt(config.base_name + GFA_EXTENSION, config.parameters);
    if(result.first.get() == nullptr || result.second.get() == nullptr)
    {
      std::cerr << "gfa2gbwt: Construction failed" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    index.swap(*(result.first));
    if(config.show_progress)
    {
      std::cerr << "Building GBWTGraph" << std::endl;
    }
    graph = GBWTGraph(index, *(result.second));
    if(config.translation)
    {
      translation = result.second->segment_translation;
    }
  }

  if(config.output == output_gg)
  {
    if(config.show_progress)
    {
      std::cerr << "Serializing GBWT" << std::endl;
    }
    std::string gbwt_name = config.base_name + gbwt::GBWT::EXTENSION;
    std::ofstream gbwt_file(gbwt_name, std::ios_base::binary);
    if(!gbwt_file)
    {
      std::cerr << "gfa2gbwt: Cannot open file " << gbwt_name << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    index.serialize(gbwt_file);
    gbwt_file.close();

    if(config.show_progress)
    {
      std::cerr << "Serializing GBWTGraph" << std::endl;
    }
    std::string graph_name = config.base_name + GBWTGraph::EXTENSION;
    std::ofstream graph_file(graph_name, std::ios_base::binary);
    if(!graph_file)
    {
      std::cerr << "gfa2gbwt: Cannot open file " << graph_name << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    graph.serialize(graph_file);
    graph_file.close();
  }

  if(config.output == output_gbz)
  {
    if(config.show_progress)
    {
      std::cerr << "Compressing GBWTGraph" << std::endl;
    }
    std::string graph_name = config.base_name + GBWTGraph::COMPRESSED_EXTENSION;
    std::ofstream out(graph_name, std::ios_base::binary);
    if(!out)
    {
      std::cerr << "gfa2gbwt: Cannot open file " << graph_name << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    index.serialize(out);
    graph.compress(out);
    out.close();
  }

  // TODO we should extract this from GBWTGraph
  if(config.translation)
  {
    if(config.show_progress)
    {
      std::cerr << "Writing translation table" << std::endl;
    }
    std::string translation_name = config.base_name + SequenceSource::TRANSLATION_EXTENSION;
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

  if(config.show_progress)
  {
    std::cerr << std::endl;
    gbwt::printStatistics(index, config.base_name, std::cerr);
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Used " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    std::cerr << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: gfa2gbwt [mode] [options] basename" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Modes:" << std::endl;
  std::cerr << "  -b, --build-graph      read " << GFA_EXTENSION << ", write " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << " (default)" << std::endl;
  std::cerr << "  -c, --compress-gfa     read " << GFA_EXTENSION << ", write " << GBWTGraph::COMPRESSED_EXTENSION << std::endl;
  std::cerr << std::endl;
  std::cerr << "General options:" << std::endl;
  std::cerr << "  -p, --progress         show progress information" << std::endl;
  std::cerr << "  -t, --translation      write translation table into a " << SequenceSource::TRANSLATION_EXTENSION << " file" << std::endl;
  std::cerr << std::endl;
  std::cerr << "GFA parsing parameters:" << std::endl;
  std::cerr << "  -m, --max-node N       break > N bp segments into multiple nodes (default " << MAX_NODE_LENGTH << ")" << std::endl;
  std::cerr << "                         (minimizer index requires nodes of length <= 1024 bp)" << std::endl;
  std::cerr << "  -r, --path-regex STR   parse path names using regex STR (default " << GFAParsingParameters::DEFAULT_REGEX << ")" << std::endl;
  std::cerr << "  -f, --path-fields STR  map the submatches to fields STR (default " << GFAParsingParameters::DEFAULT_FIELDS << ")" << std::endl;
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

Config::Config(int argc, char** argv)
{
  if(argc < 2) { printUsage(EXIT_SUCCESS); }

  // Data for `getopt_long()`.
  int c = 0, option_index = 0;
  option long_options[] =
  {
    { "build-graph", no_argument, 0, 'b' },
    { "compress-gfa", no_argument, 0, 'c' },
    { "progress", no_argument, 0, 'p' },
    { "translation", no_argument, 0, 't' },
    { "max-node", required_argument, 0, 'm' },
    { "path-regex", required_argument, 0, 'r' },
    { "path-fields", required_argument, 0, 'f' },
  };

  // Process options.
  while((c = getopt_long(argc, argv, "bcptm:r:f:", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'b':
      this->input = input_gfa;
      this->output = output_gg;
      break;
    case 'c':
      this->input = input_gfa;
      this->output = output_gbz;
      break;

    case 'p':
      this->show_progress = true;
      this->parameters.show_progress = true;
      break;
    case 't':
      this->translation = true;
      break;

    case 'm':
      try { this->parameters.max_node_length = std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "gfa2gbwt: Invalid maximum node length: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'r':
      this->parameters.path_name_regex = optarg;
      break;
    case 'f':
      this->parameters.path_name_fields = optarg;
      break;

    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  // Sanity checks.
  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  this->base_name = argv[optind]; optind++;
}

//------------------------------------------------------------------------------
