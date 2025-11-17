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
enum output_type { output_gfa, output_gbz, output_graph, output_none };

struct Config
{
  Config(int argc, char** argv);

  GFAParsingParameters parameters;
  GFAExtractionParameters output_parameters;
  std::string basename;

  input_type input = input_gfa;
  output_type output = output_gbz;

  bool translation = false;
  bool show_progress = false;
  bool simple_sds_graph = false;
  bool bitvectors = false;
};

const std::string tool_name = "GFA to GBWTGraph";

void parse_gfa(GBZ& gbz, const Config& config); // May throw `std::runtime_error`.
void load_gbz(GBZ& gbz, const Config& config);
void load_graph(GBZ& gbz, const Config& config);

void write_gfa(const GBZ& gbz, const Config& config);
void write_gbz(const GBZ& gbz, const Config& config);
void write_graph(const GBZ& gbz, const Config& config);

void extract_translation(const GBZ& gbz, const Config& config);
void extract_bitvectors(const GBZ& gbz, const Config& config);

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
    if(config.input == input_gfa)
    {
      gbwt::printHeader("--approx-jobs", std::cerr) << config.parameters.approximate_num_jobs << std::endl;
      gbwt::printHeader("--parallel-jobs", std::cerr) << config.parameters.parallel_jobs << std::endl;
      gbwt::printHeader("--max-node", std::cerr) << config.parameters.max_node_length << std::endl;
      gbwt::printHeader("--path-regex", std::cerr) << config.parameters.path_name_formats.front().regex << std::endl;
      gbwt::printHeader("--path-fields", std::cerr) << config.parameters.path_name_formats.front().fields << std::endl;
      gbwt::printHeader("--path-sense", std::cerr) << (int)config.parameters.path_name_formats.front().sense << std::endl;
    }
    if(config.output == output_gfa)
    {
      gbwt::printHeader("--parallel-jobs", std::cerr) << config.output_parameters.num_threads << std::endl;
      gbwt::printHeader("--cache-records", std::cerr) << config.output_parameters.large_record_bytes << std::endl;
      gbwt::printHeader("--paths", std::cerr) << GFAExtractionParameters::mode_name(config.output_parameters.mode) << std::endl;
      if(!(config.output_parameters.use_translation))
      {
        std::cerr << "--no-translation" << std::endl;
      }
    }
    std::cerr << std::endl;
  }

  // This is the data we are using.
  GBZ gbz;
  std::unordered_map<std::string, std::pair<nid_t, nid_t>> translation;

  try
  {
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

    // Extract parts of the GBZ.
    if(config.translation)
    {
      extract_translation(gbz, config);
    }
    if(config.bitvectors)
    {
      extract_bitvectors(gbz, config);
    }
  }
  catch(const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(config.show_progress)
  {
    std::cerr << std::endl;
    if(config.input == input_gfa)
    {
      gbwt::printStatistics(gbz.index, config.basename, std::cerr);
    }
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
  std::cerr << "  -c, --compress-gfa      read " << GFA_EXTENSION << ", write " << GBZ::EXTENSION << " (default)" << std::endl;
  std::cerr << "  -d, --decompress-gfa    read " << GBZ::EXTENSION << ", write " << GFA_EXTENSION << std::endl;
  std::cerr << "  -b, --build-graph       read " << GFA_EXTENSION << ", write " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << std::endl;
  std::cerr << "  -e, --extract-gfa       read " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << ", write " << GFA_EXTENSION << std::endl;
  std::cerr << "  -C, --compress-graph    read " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << ", write " << GBZ::EXTENSION << std::endl;
  std::cerr << "  -D, --decompress-graph  read " << GBZ::EXTENSION << ", write " << gbwt::GBWT::EXTENSION << " and " << GBWTGraph::EXTENSION << std::endl;
  std::cerr << std::endl;
  std::cerr << "General options:" << std::endl;
  std::cerr << "  -p, --progress          show progress information" << std::endl;
  std::cerr << "  -t, --translation       write translation table into a " << SequenceSource::TRANSLATION_EXTENSION << " file" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Parallel options:" << std::endl;
  std::cerr << "  -j, --approx-jobs N     create approximately N GBWT construction jobs (default " << GFAParsingParameters::APPROXIMATE_NUM_JOBS << ")" << std::endl;
  std::cerr << "  -P, --parallel-jobs N   run N construction / extraction jobs in parallel (default 1)" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Output options:" << std::endl;
  std::cerr << "  -R, --cache-records N   cache > N-byte GBWT records for " << GFA_EXTENSION << " output (default " << GFAExtractionParameters::LARGE_RECORD_BYTES << ")" << std::endl;
  std::cerr << "  -s, --simple-sds-graph  serialize " << GBWTGraph::EXTENSION << " in simple-sds format instead of libhandlegraph format" << std::endl;
  std::cerr << "                          (this tool cannot read simple-sds graphs)" << std::endl;
  std::cerr << "      --paths STR         extract paths as STR (default, pan-sn, ref-only)" << std::endl;
  std::cerr << "      --pan-sn            extract paths as P-lines with PanSN names" << std::endl;
  std::cerr << "      --ref-only          extract only named paths as P-lines" << std::endl;
  std::cerr << "      --no-translation    do not use the node-to-segment translation" << std::endl;
  std::cerr << std::endl;
  std::cerr << "GFA parsing parameters:" << std::endl;
  std::cerr << "  -m, --max-node N        break > N bp segments into multiple nodes (default " << MAX_NODE_LENGTH << ")" << std::endl;
  std::cerr << "                          (minimizer index requires nodes of length <= 1024 bp)" << std::endl;
  std::cerr << "  -r, --path-regex STR    parse path names using regex STR (default " << GFAParsingParameters::DEFAULT_REGEX << ")" << std::endl;
  std::cerr << "  -f, --path-fields STR   map the submatches to fields STR (default " << GFAParsingParameters::DEFAULT_FIELDS << ")" << std::endl;
  std::cerr << "                          (the first submatch is the entire path name)" << std::endl;
  std::cerr << "      --path-sense INT    assign paths the sense INT (default " << (int)GFAParsingParameters::DEFAULT_SENSE << ")" << std::endl;
  std::cerr << "      --pan-sn            parse PanSN path names (sets --path-regex, --path-fields, and --path-sense)" << std::endl;
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

  constexpr int OPT_PATHS = 1000;
  constexpr int OPT_PAN_SN = 1001;
  constexpr int OPT_REF_ONLY = 1002;
  constexpr int OPT_NO_TRANSLATION = 1003;
  constexpr int OPT_PATH_SENSE = 1004;

  // Data for `getopt_long()`.
  int c = 0, option_index = 0;
  option long_options[] =
  {
    { "compress-gfa", no_argument, 0, 'c' },
    { "decompress-gfa", no_argument, 0, 'd' },
    { "build-graph", no_argument, 0, 'b' },
    { "extract-gfa", no_argument, 0, 'e' },
    { "compress-graph", no_argument, 0, 'C' },
    { "decompress-graph", no_argument, 0, 'D' },
    { "load-gbz", no_argument, 0, 'l' }, // Hidden.
    { "bitvectors", no_argument, 0, 'B' }, // Hidden.
    { "progress", no_argument, 0, 'p' },
    { "translation", no_argument, 0, 't' },
    { "approx-jobs", required_argument, 0, 'j' },
    { "parallel-jobs", required_argument, 0, 'P' },
    { "cache-records", required_argument, 0, 'R' },
    { "simple-sds-graph", no_argument, 0, 's' },
    { "paths", required_argument, 0, OPT_PATHS },
    { "pan-sn", no_argument, 0, OPT_PAN_SN },
    { "ref-only", no_argument, 0, OPT_REF_ONLY },
    { "no-translation", no_argument, 0, OPT_NO_TRANSLATION },
    { "max-node", required_argument, 0, 'm' },
    { "path-regex", required_argument, 0, 'r' },
    { "path-fields", required_argument, 0, 'f' },
    { "path-sense", required_argument, 0, OPT_PATH_SENSE },
    { 0, 0, 0, 0 }
  };

  // Process options.
  while((c = getopt_long(argc, argv, "cdbeCDlBptj:P:R:sm:r:f:", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'c':
      this->input = input_gfa;
      this->output = output_gbz;
      break;
    case 'd':
      this->input = input_gbz;
      this->output = output_gfa;
      break;
    case 'b':
      this->input = input_gfa;
      this->output = output_graph;
      break;
    case 'e':
      this->input = input_graph;
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

    case 'l':
      this->input = input_gbz;
      this->output = output_none;
      break;
    case 'B':
      this->bitvectors = true;
      break;

    case 'p':
      this->show_progress = true;
      this->parameters.show_progress = true;
      this->output_parameters.show_progress = true;
      break;
    case 't':
      this->translation = true;
      break;

    case 'j':
      try { this->parameters.approximate_num_jobs = std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "gfa2gbwt: Invalid number of jobs: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'P':
      try { this->parameters.parallel_jobs = std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "gfa2gbwt: Invalid number of parallel jobs: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      this->output_parameters.num_threads = this->parameters.parallel_jobs;
      break;

    case 'R':
      try { this->output_parameters.large_record_bytes = std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "gfa2gbwt: Invalid record caching threshold: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 's':
      this->simple_sds_graph = true;
      break;
    case OPT_PATHS:
      try { this->output_parameters.mode = GFAExtractionParameters::get_mode(optarg); }
      catch(const std::exception& e)
      {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case OPT_PAN_SN:
      this->output_parameters.mode = GFAExtractionParameters::mode_pan_sn;
      this->parameters.path_name_formats.clear();
      this->parameters.path_name_formats.emplace_back(
        GFAParsingParameters::PAN_SN_REGEX,
        GFAParsingParameters::PAN_SN_FIELDS,
        GFAParsingParameters::PAN_SN_SENSE
      );
      break;
    case OPT_REF_ONLY:
      this->output_parameters.mode = GFAExtractionParameters::mode_ref_only;
      break;
    case OPT_NO_TRANSLATION:
      this->output_parameters.use_translation = false;
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
      this->parameters.path_name_formats.front().regex = optarg;
      break;
    case 'f':
      this->parameters.path_name_formats.front().fields = optarg;
      break;
    case OPT_PATH_SENSE:
      // TODO: Use by-name sense parsing when available in libhandlegraph
      try { this->parameters.path_name_formats.front().sense = (PathSense) std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "gfa2gbwt: Invalid path sense: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (this->parameters.path_name_formats.front().sense > PathSense::HAPLOTYPE) {
        std::cerr << "gfa2gbwt: Invalid path sense: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
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
  }
  // This may throw an exception.
  auto result = gfa_to_gbwt(gfa_name, config.parameters);

  if(config.show_progress)
  {
    std::cerr << "Building GBWTGraph" << std::endl;
  }
  GraphName parent = result.second->graph_name();
  gbz = GBZ(result.first, result.second);

  if(config.show_progress)
  {
    std::cerr << "Computing graph name" << std::endl;
  }
  gbz.compute_pggname(&parent); // The heuristic will determine the relationship correctly.
  if(config.show_progress)
  {
    std::cerr << "Graph name: " << gbz.pggname() << std::endl;
  }
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
  gbwt_to_gfa(gbz, out, config.output_parameters);
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

sdsl::bit_vector
convert(const sdsl::sd_vector<>& source)
{
  sdsl::bit_vector result(source.size(), 0);
  for(auto iter = source.one_begin(); iter != source.one_end(); ++iter)
  {
    result[iter->second] = 1;
  }
  return result;
}

sdsl::bit_vector
convert(const sdsl::int_vector<0>& source, size_t n)
{
  // The sequences are stored in both orientations, so we only take every other value and divide by two.
  sdsl::bit_vector result(n / 2, 0);
  for(size_t i = 0; i < source.size(); i += 2)
  {
    result[source[i] / 2] = 1;
  }
  return result;
}

void
extract_bitvectors(const GBZ& gbz, const Config& config)
{
  std::string bv_gbwt = config.basename + ".bv_gbwt";
  std::string bv_seq = config.basename + ".bv_seq";

  if(config.show_progress)
  {
    std::cerr << "Writing bitvectors to " << bv_gbwt << " and " << bv_seq << std::endl;
  }
  sdsl::bit_vector bv = convert(gbz.index.bwt.index);
  sdsl::simple_sds::serialize_to(bv, bv_gbwt);
  bv = convert(gbz.graph.sequences.index, gbz.graph.sequences.length() + 1);
  sdsl::simple_sds::serialize_to(bv, bv_seq);
}

//------------------------------------------------------------------------------
