#include <iostream>
#include <fstream>
#include <string>

#include <getopt.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>

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
  output_type output = output_graph;

  bool translation = false;
  bool show_progress = false;
};

const std::string tool_name = "GFA to GBWTGraph";

void parse_gfa(gbwt::GBWT& index, GBWTGraph& graph, const Config& config);
void load_gbz(gbwt::GBWT& index, GBWTGraph& graph, const Config& config);
void load_graph(gbwt::GBWT& index, GBWTGraph& graph, const Config& config);

void write_gfa(const gbwt::GBWT& index, const GBWTGraph& graph, const Config& config);
void write_gbz(const gbwt::GBWT& index, const GBWTGraph& graph, const Config& config);
void write_graph(const gbwt::GBWT& index, const GBWTGraph& graph, const Config& config);

void extract_translation(const GBWTGraph& graph, const Config& config);

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
  gbwt::GBWT index;
  GBWTGraph graph;
  std::unordered_map<std::string, std::pair<nid_t, nid_t>> translation;

  // Handle the input.
  if(config.input == input_gfa)
  {
    parse_gfa(index, graph, config);
  }
  else if(config.input == input_gbz)
  {
    load_gbz(index, graph, config);
  }
  else if(config.input == input_graph)
  {
    load_graph(index, graph, config);
  }

  // Handle the output.
  if(config.output == output_gfa)
  {
    write_gfa(index, graph, config);
  }
  else if(config.output == output_gbz)
  {
    write_gbz(index, graph, config);
  }
  else if(config.output == output_graph)
  {
    write_graph(index, graph, config);
  }

  // Extract the translation.
  if(config.translation)
  {
    extract_translation(graph, config);
  }

  if(config.show_progress)
  {
    std::cerr << std::endl;
    gbwt::printStatistics(index, config.basename, std::cerr);
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
  std::cerr << "  -b, --build-graph       read " << GFA_EXTENSION << ", write " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << " (default)" << std::endl;
  std::cerr << "  -e, --extract-gfa       read " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << ", write " << GFA_EXTENSION << std::endl;
  std::cerr << "  -c, --compress-gfa      read " << GFA_EXTENSION << ", write " << GBWTGraph::COMPRESSED_EXTENSION << std::endl;
  std::cerr << "  -d, --decompress-gfa    read " << GBWTGraph::COMPRESSED_EXTENSION << ", write " << GFA_EXTENSION << std::endl;
  std::cerr << "  -C, --compress-graph    read " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << ", write " << GBWTGraph::COMPRESSED_EXTENSION << std::endl;
  std::cerr << "  -D, --decompress-graph  read " << GBWTGraph::COMPRESSED_EXTENSION << ", write " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << std::endl;
  std::cerr << std::endl;
  std::cerr << "General options:" << std::endl;
  std::cerr << "  -p, --progress          show progress information" << std::endl;
  std::cerr << "  -t, --translation       write translation table into a " << SequenceSource::TRANSLATION_EXTENSION << " file" << std::endl;
  std::cerr << std::endl;
  std::cerr << "GFA parsing parameters:" << std::endl;
  std::cerr << "  -m, --max-node N        break > N bp segments into multiple nodes (default " << MAX_NODE_LENGTH << ")" << std::endl;
  std::cerr << "                          (minimizer index requires nodes of length <= 1024 bp)" << std::endl;
  std::cerr << "  -r, --path-regex STR    parse path names using regex STR (default " << GFAParsingParameters::DEFAULT_REGEX << ")" << std::endl;
  std::cerr << "  -f, --path-fields STR   map the submatches to fields STR (default " << GFAParsingParameters::DEFAULT_FIELDS << ")" << std::endl;
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
    { "max-node", required_argument, 0, 'm' },
    { "path-regex", required_argument, 0, 'r' },
    { "path-fields", required_argument, 0, 'f' },
  };

  // Process options.
  while((c = getopt_long(argc, argv, "becdCDptm:r:f:", long_options, &option_index)) != -1)
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
parse_gfa(gbwt::GBWT& index, GBWTGraph& graph, const Config& config)
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
  index.swap(*(result.first));

  if(config.show_progress)
  {
    std::cerr << "Building GBWTGraph" << std::endl;
  }
  graph = GBWTGraph(index, *(result.second));
}

void
load_gbz(gbwt::GBWT& index, GBWTGraph& graph, const Config& config)
{
  std::string gbz_name = config.basename + GBWTGraph::COMPRESSED_EXTENSION;

  if(config.show_progress)
  {
    std::cerr << "Decompressing GBWT and GBWTGraph from " << gbz_name << std::endl;
  }
  std::ifstream in(gbz_name, std::ios_base::binary);
  if(!in)
  {
    std::cerr << "gfa2gbwt: Cannot open file " << gbz_name << " for reading" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  index.load(in);
  graph.decompress(in, index);
  in.close();
}

void
load_graph(gbwt::GBWT& index, GBWTGraph& graph, const Config& config)
{
  std::string gbwt_name = config.basename + gbwt::GBWT::EXTENSION;
  std::string graph_name = config.basename + GBWTGraph::EXTENSION;

  if(config.show_progress)
  {
    std::cerr << "Loading GBWT from " << gbwt_name << std::endl;
  }
  if(!sdsl::load_from_file(index, gbwt_name))
  {
    std::cerr << "gfa2gbwt: Cannot load GBWT from " << gbwt_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(config.show_progress)
  {
    std::cerr << "Loading GBWTGraph from " << graph_name << std::endl;
  }
  graph.deserialize(graph_name);
  graph.set_gbwt(index);
}

//------------------------------------------------------------------------------

// TODO this should be a library function with tests
// TODO and progress information (number of each line type)
void
write_gfa(const gbwt::GBWT& index, const GBWTGraph& graph, const Config& config)
{
  std::string gfa_name = config.basename + GFA_EXTENSION;

  if(config.show_progress)
  {
    std::cerr << "Writing the GFA to " << gfa_name << std::endl;
  }
  std::ofstream out(gfa_name, std::ios_base::binary);

  // GFA header.
  // TODO This should use buffered writing.
  out << "H\tVN:Z:1.0\n";

  // S-lines.
  // TODO we should store views of all segment names for P-lines and W-lines
  if(graph.has_segment_names())
  {
    graph.for_each_segment([&](const std::string& name, std::pair<nid_t, nid_t> nodes) -> bool
    {
      out << "S\t" << name << "\t";
      for(nid_t id = nodes.first; id < nodes.second; id++)
      {
        out << graph.get_sequence(graph.get_handle(id, false));
      }
      out << "\n";
      return true;
    });
  }
  else
  {
    graph.for_each_handle([&](const handle_t& handle)
    {
      out << "S\t" << graph.get_id(handle) << "\t" << graph.get_sequence(handle) << "\n";
    });
  }

  // L-lines.
  if(graph.has_segment_names())
  {
    graph.for_each_link([&](const edge_t& edge, const std::string& from, const std::string& to) -> bool
    {
      out << "L\t"
          << from << (graph.get_is_reverse(edge.first) ? "\t-\t" : "\t+\t")
          << to << (graph.get_is_reverse(edge.second) ? "\t-\t" : "\t+\t")
          << "*\n";
      return true;
    });
  }
  else
  {
    graph.for_each_edge([&](const edge_t& edge)
    {
      out << "L\t"
          << graph.get_id(edge.first) << (graph.get_is_reverse(edge.first) ? "\t-\t" : "\t+\t")
          << graph.get_id(edge.second) << (graph.get_is_reverse(edge.second) ? "\t-\t" : "\t+\t")
          << "*\n";
    });
  }

  // P-lines.
  gbwt::size_type ref_sample = index.metadata.sample(REFERENCE_PATH_SAMPLE_NAME);
  std::vector<gbwt::size_type> ref_paths = index.metadata.pathsForSample(ref_sample);
  for(gbwt::size_type path_id : ref_paths)
  {
    gbwt::vector_type path = index.extract(gbwt::Path::encode(path_id, false));
    out << "P\t" << index.metadata.contig(index.metadata.path(path_id).contig) << "\t";
    size_t segments = path.size();
    if(graph.has_segment_names())
    {
      // We assume that the path is a valid concatenation of segments.
      size_t offset = 0;
      while(offset < path.size())
      {
        auto segment = graph.get_segment(GBWTGraph::node_to_handle(path[offset]));
        out << segment.first << (gbwt::Node::is_reverse(path[offset]) ? "-" : "+");
        offset += segment.second.second - segment.second.first;
        if(offset < path.size()) { out << ","; }
      }
    }
    else
    {
      for(size_t i = 0; i < path.size(); i++)
      {
        out << gbwt::Node::id(path[i]) << (gbwt::Node::is_reverse(path[i]) ? "-" : "+");
        if(i + 1 < path.size()) { out << ","; }
      }
    }
    out << "\t";
    for(size_t i = 1; i < segments; i++)
    {
      out << "*";
      if(i + 1 < segments) { out << ","; }
    }
    out << "\n";
  }

  // W-lines.
  for(gbwt::size_type path_id = 0; path_id < index.metadata.paths(); path_id++)
  {
    const gbwt::PathName& path_name = index.metadata.path(path_id);
    if(path_name.sample == ref_sample) { continue; }
    gbwt::vector_type path = index.extract(gbwt::Path::encode(path_id, false));
    size_t length = 0;
    for(auto node : path) { length += graph.get_length(GBWTGraph::node_to_handle(node)); }
    out << "W\t" << index.metadata.sample(path_name.sample)
        << "\t" << path_name.phase
        << "\t" << index.metadata.contig(path_name.contig)
        << "\t" << path_name.count << "\t" << (path_name.count + length) << "\t";
    if(graph.has_segment_names())
    {
      // We assume that the path is a valid concatenation of segments.
      size_t offset = 0;
      while(offset < path.size())
      {
        auto segment = graph.get_segment(GBWTGraph::node_to_handle(path[offset]));
        out << (gbwt::Node::is_reverse(path[offset]) ? "<" : ">") << segment.first;
        offset += segment.second.second - segment.second.first;
      }
    }
    else
    {
      for(auto node : path)
      {
        out << (gbwt::Node::is_reverse(node) ? "<" : ">") << gbwt::Node::id(node);
      }
    }
    out << "\n";
  }
}

void
write_gbz(const gbwt::GBWT& index, const GBWTGraph& graph, const Config& config)
{
  std::string gbz_name = config.basename + GBWTGraph::COMPRESSED_EXTENSION;

  if(config.show_progress)
  {
    std::cerr << "Compressing GBWT and GBWTGraph to " << gbz_name << std::endl;
  }
  std::ofstream out(gbz_name, std::ios_base::binary);
  if(!out)
  {
    std::cerr << "gfa2gbwt: Cannot open file " << gbz_name << " for writing" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  index.serialize(out);
  graph.compress(out);
  out.close();
}

void
write_graph(const gbwt::GBWT& index, const GBWTGraph& graph, const Config& config)
{
  std::string gbwt_name = config.basename + gbwt::GBWT::EXTENSION;
  std::string graph_name = config.basename + GBWTGraph::EXTENSION;

  if(config.show_progress)
  {
    std::cerr << "Writing GBWT to " << gbwt_name << std::endl;
  }
  if(!sdsl::store_to_file(index, gbwt_name))
  {
    std::cerr << "gfa2gbwt: Cannot write GBWT to " << gbwt_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(config.show_progress)
  {
    std::cerr << "Writing GBWTGraph to " << graph_name << std::endl;
  }
  graph.serialize(graph_name);
}

//------------------------------------------------------------------------------

void
extract_translation(const GBWTGraph& graph, const Config& config)
{
  std::string translation_name = config.basename + SequenceSource::TRANSLATION_EXTENSION;

  if(config.show_progress)
  {
    std::cerr << "Writing the translation table to " << translation_name << std::endl;
  }
  // TODO This should use buffered writing, like in GFA extraction.
  std::ofstream out(translation_name, std::ios_base::binary);
  graph.for_each_segment([&](const std::string& name, std::pair<nid_t, nid_t> nodes) -> bool
  {
    out << "T\t" << name << "\t" << nodes.first;
    for(nid_t i = nodes.first + 1; i < nodes.second; i++)
    {
      out << "," << i;
    }
    out << "\n";
    return true;
  });
  out.close();
}

//------------------------------------------------------------------------------
