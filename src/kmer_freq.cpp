#include <iostream>
#include <fstream>
#include <map>

#include <getopt.h>
#include <unistd.h>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/index.h>

using namespace gbwtgraph;

//------------------------------------------------------------------------------

const std::string tool_name = "Kmer Frequencies";

struct Config
{
  Config(int argc, char** argv);

  constexpr static size_t DEFAULT_K = 29;

  size_t k = DEFAULT_K;

  bool distribution = false;
  size_t threshold = 0;
  size_t threads = 1;
  size_t hash_table_size = KmerIndex<Key64>::INITIAL_CAPACITY;

  bool space_efficient = false;

  bool verbose = true; // TODO: Do we need an option to silence this?

  std::string filename;
};

void print_distribution(const std::map<size_t, size_t>& distribution, const std::string& header1, const std::string& header2);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  typedef KmerIndex<Key64> index_type;

  double start = gbwt::readTimer();
  Config config(argc, argv);
  omp_set_num_threads(config.threads);

  // Load the graph.
  GBZ gbz;
  if(config.verbose)
  {
    std::cerr << "Loading GBZ from " << config.filename << std::endl;
  }
  sdsl::simple_sds::load_from(gbz, config.filename);

  // (number of hits, number of kmers)
  std::map<size_t, size_t> distribution;
  auto update_distribution = [&](const index_type& index)
  {
    index.for_each_kmer([&](index_type::key_type, index_type::multi_value_type values)
    {
      distribution[values.second]++;
    });
  };

  // Build the kmer indexes and extract the kmer frequency distribution.
  double checkpoint = gbwt::readTimer();
  size_t total_kmers = 0, total_unique = 0, total_values = 0;
  if(config.space_efficient)
  {
    for(char base : { 'A', 'C', 'G', 'T' })
    {
      if(config.verbose)
      {
        std::cerr << "Building a " << config.k << "-mer index for middle base " << base << std::endl;
      }
      index_type index(config.hash_table_size);
      build_kmer_index(gbz.graph, index, config.k, [&](Key64 key) -> bool { return (key.access(config.k, config.k / 2) == base); });
      if(config.verbose)
      {
        std::cerr << index.size() << " kmers (" << index.unique_keys() << " unique) with " << index.number_of_values() << " hits" << std::endl;
      }
      total_kmers += index.size(); total_unique += index.unique_keys(); total_values += index.number_of_values();
      update_distribution(index);
    }
  }
  else
  {
    if(config.verbose)
    {
      std::cerr << "Building all " << config.k << "-mer indexes in parallel" << std::endl;
    }
    std::array<index_type, 4> indexes
    {
      index_type(config.hash_table_size), index_type(config.hash_table_size), index_type(config.hash_table_size), index_type(config.hash_table_size)
    };
    build_kmer_indexes(gbz.graph, indexes, config.k);
    for(auto& index : indexes)
    {
      total_kmers += index.size(); total_unique += index.unique_keys(); total_values += index.number_of_values();
      update_distribution(index);
    }
  }
  if(config.verbose)
  {
    double seconds = gbwt::readTimer() - checkpoint;
    std::cerr << "Built the kmer indexes in " << seconds << " seconds" << std::endl;
    std::cerr << total_kmers << " kmers (" << total_unique << " unique) with " << total_values << " hits" << std::endl;
  }

  if(config.threshold > 0)
  {
    size_t count = 0;
    for(auto iter = distribution.begin(); iter != distribution.end(); ++iter)
    {
      if(iter->first > config.threshold) { count += iter->second; }
    }
    std::cout << count << std::endl;
  }

  if(config.distribution)
  {
    print_distribution(distribution, "Frequency", "Kmers");
  }

  // Deleting the hash table may take some time, so we delete it before returning final usage.
  if(config.verbose)
  {
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

  std::cerr << "Usage: kmer_freq [options] graph.gbz" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -d, --distribution   print kmer frequency distribution" << std::endl;
  std::cerr << "  -f, --frequent N     count the number of kmers with frequency > N" << std::endl;
  std::cerr << "  -k, --kmer-length N  count N-mers (default: " << Config::DEFAULT_K << ")" << std::endl;
  std::cerr << "  -s, --save-memory    use the slower space-efficient algorithm" << std::endl;
  std::cerr << "  -t, --threads N      use N parallel threads" << std::endl;
  std::cerr << "  -H, --hash-table N   initialize the hash table with 2^N cells" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------

Config::Config(int argc, char** argv)
{
  if(argc < 2) { printUsage(EXIT_SUCCESS); }

  size_t max_threads = omp_get_max_threads();
  size_t min_width = 10;
  size_t max_width = 36;

  // Data for `getopt_long()`.
  int c = 0, option_index = 0;
  option long_options[] =
  {
    { "distribution", no_argument, 0, 'd' },
    { "frequent", required_argument, 0, 'f' },
    { "kmer-length", required_argument, 0, 'k' },
    { "save-memory", no_argument, 0, 's' },
    { "threads", required_argument, 0, 't' },
    { "hash-tables", required_argument, 0, 'H' },
    { 0, 0, 0, 0 }
  };

  // Process options.
  while((c = getopt_long(argc, argv, "df:k:st:H:", long_options, &option_index)) != -1)
  {
    switch(c)
    {
    case 'd':
      this->distribution = true;
      break;
    case 'f':
      try { this->threshold = std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "kmer_freq: Invalid frequency threshold: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'k':
      try { this->k = std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "kmer_freq: Invalid kmer length: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if(this->k < 2 || this->k > Key64::KMER_MAX_LENGTH)
      {
        std::cerr << "kmer_freq: Invalid kmer length: " << this->k << " (must be from 2 to " << Key64::KMER_MAX_LENGTH << ")" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 's':
      this->space_efficient = true;
      break;
    case 't':
      try { this->threads = std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "kmer_freq: Invalid number of threads: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if(this->threads == 0 || this->threads > max_threads)
      {
        std::cerr << "kmer_freq: Invalid number of threads: " << this->threads << " (must be from 1 to " << max_threads << ")" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'H':
      size_t width;
      try { width = std::stoul(optarg); }
      catch(const std::invalid_argument&)
      {
        std::cerr << "kmer_freq: Invalid hash table width: " << optarg << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if(width < min_width || width > max_width)
      {
        std::cerr << "kmer_freq: Invalid hash table width: " << width << " (must be from " << min_width << " to " << max_width << ")" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      this->hash_table_size = size_t(1) << width;
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
