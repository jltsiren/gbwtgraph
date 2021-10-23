#include <iostream>
#include <fstream>
#include <string>

#include <getopt.h>
#include <unistd.h>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>
#include <gbwtgraph/internal.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

enum input_type { input_gfa, input_gbz, input_graph };
enum output_type { output_gfa, output_gbz, output_graph };

struct Config
{
  Config(int argc, char** argv);

  GFAParsingParameters parameters;
  std::string basename;

  input_type input = input_gfa;
  output_type output = output_gbz;

  bool translation = false;
  bool show_progress = false;
  bool simple_sds_graph = false;
};

const std::string tool_name = "GFA to GBWTGraph";

void parse_gfa(GBZ& gbz, const Config& config);
void load_gbz(GBZ& gbz, const Config& config);
void load_graph(GBZ& gbz, const Config& config);

void write_gfa(const GBZ& gbz, const Config& config);
void write_gbz(const GBZ& gbz, const Config& config);
void write_graph(const GBZ& gbz, const Config& config);

void extract_translation(const GBZ& gbz, const Config& config);

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
    gbwt::printHeader("Base name", std::cerr) << config.basename << std::endl;
    std::cerr << std::endl;
  }

  // This is the data we are using.
  GBZ gbz;
  std::unordered_map<std::string, std::pair<nid_t, nid_t>> translation;

  // Handle the input.
  if(config.input == input_gfa)
  {
    parse_gfa(gbz, config);
  }
  else if(config.input == input_gbz)
  {
    load_gbz(gbz, config);
  }
  else if(config.input == input_graph)
  {
    load_graph(gbz, config);
  }

  // Handle the output.
  if(config.output == output_gfa)
  {
    write_gfa(gbz, config);
  }
  else if(config.output == output_gbz)
  {
    write_gbz(gbz, config);
  }
  else if(config.output == output_graph)
  {
    write_graph(gbz, config);
  }

  // Extract the translation.
  if(config.translation)
  {
    extract_translation(gbz, config);
  }

  if(config.show_progress)
  {
    std::cerr << std::endl;
    gbwt::printStatistics(gbz.index, config.basename, std::cerr);
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
  std::cerr << "  -b, --build-graph       read " << GFA_EXTENSION << ", write " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << std::endl;
  std::cerr << "  -e, --extract-gfa       read " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << ", write " << GFA_EXTENSION << std::endl;
  std::cerr << "  -c, --compress-gfa      read " << GFA_EXTENSION << ", write " << GBZ::EXTENSION << " (default)" << std::endl;
  std::cerr << "  -d, --decompress-gfa    read " << GBZ::EXTENSION << ", write " << GFA_EXTENSION << std::endl;
  std::cerr << "  -C, --compress-graph    read " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << ", write " << GBZ::EXTENSION << std::endl;
  std::cerr << "  -D, --decompress-graph  read " << GBZ::EXTENSION << ", write " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << std::endl;
  std::cerr << std::endl;
  std::cerr << "General options:" << std::endl;
  std::cerr << "  -p, --progress          show progress information" << std::endl;
  std::cerr << "  -t, --translation       write translation table into a " << SequenceSource::TRANSLATION_EXTENSION << " file" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Other options:" << std::endl;
  std::cerr << "  -s, --simple-sds-graph  serialize " << GBWTGraph::EXTENSION << " in simple-sds format instead of libhandlegraph format" << std::endl;
  std::cerr << "                          (this tool cannot read simple-sds graphs)" << std::endl;
  std::cerr << std::endl;
  std::cerr << "GFA parsing parameters:" << std::endl;
  std::cerr << "  -m, --max-node N        break > N bp segments into multiple nodes (default " << MAX_NODE_LENGTH << ")" << std::endl;
  std::cerr << "                          (minimizer index requires nodes of length <= 1024 bp)" << std::endl;
  std::cerr << "  -r, --path-regex STR    parse path names using regex STR (default " << GFAParsingParameters::DEFAULT_REGEX << ")" << std::endl;
  std::cerr << "  -f, --path-fields STR   map the submatches to fields STR (default " << GFAParsingParameters::DEFAULT_FIELDS << ")" << std::endl;
  std::cerr << "                          (the first submatch is the entire path name)" << std::endl;
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
    { "extract-gfa", no_argument, 0, 'e' },
    { "compress-gfa", no_argument, 0, 'c' },
    { "decompress-gfa", no_argument, 0, 'd' },
    { "compress-graph", no_argument, 0, 'C' },
    { "decompress-graph", no_argument, 0, 'D' },
    { "progress", no_argument, 0, 'p' },
    { "translation", no_argument, 0, 't' },
    { "simple-sds-graph", no_argument, 0, 's' },
    { "max-node", required_argument, 0, 'm' },
    { "path-regex", required_argument, 0, 'r' },
    { "path-fields", required_argument, 0, 'f' },
  };

  // Process options.
  while((c = getopt_long(argc, argv, "becdCDptsm:r:f:", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'b':
      this->input = input_gfa;
      this->output = output_graph;
      break;
    case 'e':
      this->input = input_graph;
      this->output = output_gfa;
      break;
    case 'c':
      this->input = input_gfa;
      this->output = output_gbz;
      break;
    case 'd':
      this->input = input_gbz;
      this->output = output_gfa;
      break;
    case 'C':
      this->input = input_graph;
      this->output = output_gbz;
      break;
    case 'D':
      this->input = input_gbz;
      this->output = output_graph;
      break;

    case 'p':
      this->show_progress = true;
      this->parameters.show_progress = true;
      break;
    case 't':
      this->translation = true;
      break;

    case 's':
      this->simple_sds_graph = true;
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
  this->basename = argv[optind]; optind++;
}

//------------------------------------------------------------------------------

void
parse_gfa(GBZ& gbz, const Config& config)
{
  std::string gfa_name = config.basename + GFA_EXTENSION;

  if(config.show_progress)
  {
    std::cerr << "Parsing GFA from " << gfa_name << " and building GBWT" << std::endl;
    std::cerr << "Path name regex: " << config.parameters.path_name_regex << std::endl;
    std::cerr << "Path name fields: " << config.parameters.path_name_fields << std::endl;
  }
  auto result = gfa_to_gbwt(gfa_name, config.parameters);
  if(result.first.get() == nullptr || result.second.get() == nullptr)
  {
    std::cerr << "gfa2gbwt: Construction failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(config.show_progress)
  {
    std::cerr << "Building GBWTGraph" << std::endl;
  }
  gbz = GBZ(result.first, result.second);
}

void
load_gbz(GBZ& gbz, const Config& config)
{
  std::string gbz_name = config.basename + GBZ::EXTENSION;
  if(config.show_progress)
  {
    std::cerr << "Decompressing GBWT and GBWTGraph from " << gbz_name << std::endl;
  }
  sdsl::simple_sds::load_from(gbz, gbz_name);
}

void
load_graph(GBZ& gbz, const Config& config)
{
  std::string gbwt_name = config.basename + gbwt::GBWT::EXTENSION;
  std::string graph_name = config.basename + GBWTGraph::EXTENSION;
  if(config.show_progress)
  {
    std::cerr << "Loading GBWT and GBWTGraph from " << gbwt_name << " and " << graph_name << std::endl;
  }
  gbz.load_from_files(gbwt_name, graph_name);
}

//------------------------------------------------------------------------------

void
write_gfa(const GBZ& gbz, const Config& config)
{
  std::string gfa_name = config.basename + GFA_EXTENSION;

  if(config.show_progress)
  {
    std::cerr << "Writing the GFA to " << gfa_name << std::endl;
  }
  std::ofstream out;
  out.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  out.open(gfa_name, std::ios_base::binary);
  gbwt_to_gfa(gbz.graph, out, config.show_progress);
  out.close();
}

void
write_gbz(const GBZ& gbz, const Config& config)
{
  std::string gbz_name = config.basename + GBZ::EXTENSION;
  if(config.show_progress)
  {
    std::cerr << "Compressing GBWT and GBWTGraph to " << gbz_name << std::endl;
  }
  sdsl::simple_sds::serialize_to(gbz, gbz_name);
}

void
write_graph(const GBZ& gbz, const Config& config)
{
  std::string gbwt_name = config.basename + gbwt::GBWT::EXTENSION;
  std::string graph_name = config.basename + GBWTGraph::EXTENSION;
  if(config.show_progress)
  {
    std::cerr << "Writing GBWT and GBWTGraph to " << gbwt_name << " and " << graph_name << std::endl;
  }
  gbz.serialize_to_files(gbwt_name, graph_name, config.simple_sds_graph);
}

//------------------------------------------------------------------------------

void
extract_translation(const GBZ& gbz, const Config& config)
{
  std::string translation_name = config.basename + SequenceSource::TRANSLATION_EXTENSION;
  if(config.show_progress)
  {
    std::cerr << "Writing the translation table to " << translation_name << std::endl;
  }

  std::ofstream out(translation_name, std::ios_base::binary);
  TSVWriter writer(out);
  gbz.graph.for_each_segment([&](const std::string& name, std::pair<nid_t, nid_t> nodes) -> bool
  {
    writer.put('T'); writer.newfield();
    writer.write(name); writer.newfield();
    writer.write(nodes.first);
    for(nid_t i = nodes.first + 1; i < nodes.second; i++)
    {
      writer.put(','); writer.write(i);
    }
    writer.newline();
    return true;
  });
}

//------------------------------------------------------------------------------
