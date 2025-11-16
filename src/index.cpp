#include <gbwtgraph/index.h>
#include <gbwtgraph/internal.h>

#include <mutex>

#include <omp.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

const std::string PATH_NAME_FIELDS_TAG = "path_name_fields";

//------------------------------------------------------------------------------

template<typename KeyType>
void index_haplotypes
(
  const GBZ& gbz, MinimizerIndex<KeyType>& index,
  const std::function<const KmerEncoding::code_type*(const pos_t&)>& get_payload
)
{
  using minimizer_type = typename MinimizerIndex<KeyType>::minimizer_type;
  using code_type = KmerEncoding::code_type;
  using value_type = KmerEncoding::value_type;

  struct alignas(CACHE_LINE_SIZE) Cache
  {
    std::vector<std::pair<minimizer_type, pos_t>> data;
  };

  // Take the graph reference for convenience.
  // Also copy graph name and relationships from GBZ tags to MinimizerIndex tags.
  const GBWTGraph& graph = gbz.graph;
  GraphName graph_name = gbz.graph_name();
  graph_name.set_tags(index.get_tags());

  // Minimizer caching. We only generate the payloads after we have removed duplicate positions.
  int threads = omp_get_max_threads();
  constexpr size_t MINIMIZER_CACHE_SIZE = 1024;
  std::vector<Cache> cache(threads);
  auto flush_cache = [&](int thread_id)
  {
    auto& current_cache = cache[thread_id].data;
    gbwt::removeDuplicates(current_cache, false);
    std::vector<code_type> payloads;
    size_t payload_size = index.payload_size();
    payloads.reserve(current_cache.size() * payload_size);
    for (auto& minimizer : current_cache)
    {
      const code_type* payload = get_payload(minimizer.second);
      for(size_t i = 0; i < payload_size; i++) { payloads.push_back(payload[i]); }
    }
    #pragma omp critical (minimizer_index)
    {
      for(size_t i = 0, payload_offset = 0; i < current_cache.size(); i++, payload_offset += payload_size)
      {
        value_type value(Position(current_cache[i].second), &payloads[payload_offset]);
        index.insert(current_cache[i].first, value);
      }
    }
    current_cache.clear();
  };

  // Minimizer finding.
  auto find_minimizers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    std::vector<minimizer_type> minimizers = index.minimizers(seq); // Calls syncmers() when appropriate.
    auto iter = traversal.begin();
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();
    for(minimizer_type& minimizer : minimizers)
    {
      if(minimizer.empty()) { continue; }

      // Find the node covering minimizer starting position.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= minimizer.offset)
      {
        node_start += node_length;
        ++iter;
        node_length = graph.get_length(*iter);
      }
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start };
      if(minimizer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      if(!Position::valid_offset(pos))
      {
        #pragma omp critical (cerr)
        {
          std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large" << std::endl;
        }
        std::exit(EXIT_FAILURE);
      }
      cache[thread_id].data.emplace_back(minimizer, pos);
    }
    if(cache[thread_id].data.size() >= MINIMIZER_CACHE_SIZE) { flush_cache(thread_id); }
  };

  // Traverse all haplotype windows
  for_each_haplotype_window(graph, index.window_bp(), find_minimizers, (threads > 1));
  for(int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }
}

// Instantiate the template, as we did not define the function in the header.

template void index_haplotypes<Key64>
(
  const GBZ& gbz, MinimizerIndex<Key64>& index,
  const std::function<const KmerEncoding::code_type*(const pos_t&)>& get_payload
);

template void index_haplotypes<Key128>
(
  const GBZ& gbz, MinimizerIndex<Key128>& index,
  const std::function<const KmerEncoding::code_type*(const pos_t&)>& get_payload
);

//------------------------------------------------------------------------------

template<typename KeyType>
void index_haplotypes_with_paths
(
  const GBZ& gbz, MinimizerIndex<KeyType>& index,
  const std::function<const KmerEncoding::code_type*(const pos_t&)>& get_payload
)
{
  if(index.payload_size() == 0)
  {
    std::string msg = "index_haplotypes_with_paths(): Payload size must be greater than zero";
    throw std::runtime_error(msg);
  }

  using minimizer_type = typename MinimizerIndex<KeyType>::minimizer_type;
  using code_type = KmerEncoding::code_type;
  using value_type = KmerEncoding::value_type;

  // This version needs the path corresponding to the minimizer.
  struct alignas(CACHE_LINE_SIZE) Cache
  {
      std::vector<std::tuple<minimizer_type, pos_t, std::vector<gbwt::node_type>>> data;
  };

  // Take the graph reference for convenience.
  // Also copy graph name and relationships from GBZ tags to MinimizerIndex tags.
  const GBWTGraph& graph = gbz.graph;
  GraphName graph_name = gbz.graph_name();
  graph_name.set_tags(index.get_tags());

  int threads = omp_get_max_threads();
  constexpr size_t MINIMIZER_CACHE_SIZE = 1024;
  std::vector<Cache> cache(threads);
  const gbwt::Metadata& metadata = graph.index->metadata;
  PathIdMap path_ids_map(metadata);
  index.set_tag(PATH_NAME_FIELDS_TAG, PathIdMap::key_type_str(path_ids_map.key_type()));

  // Minimizer caching. We only generate the payloads after we have removed duplicate positions.
  auto flush_cache = [&](int thread_id)
  {
    // Haplotypes for search states we have already seen.
    std::unordered_map<gbwt::SearchState, code_type, SearchStateHasher> state_cache;
    auto cached_gbwt = graph.get_cache();
    auto& current_cache = cache[thread_id].data;
    gbwt::removeDuplicates(current_cache, false);

    // Determine the payloads.
    std::vector<code_type> payloads;
    size_t payload_size = index.payload_size();
    payloads.reserve(current_cache.size() * payload_size);
    for(auto& minimizer : current_cache)
    {
      // Fill in the payload except for the last word.
      const code_type* payload = get_payload(std::get<1>(minimizer));
      for(size_t i = 0; i < payload_size - 1; i++) { payloads.push_back(payload[i]); }

      // Determine the haplotypes that contain the minimizer.
      auto& traversed_nodes = std::get<2>(minimizer);
      gbwt::SearchState state = cached_gbwt.find(traversed_nodes.begin(), traversed_nodes.end());
      auto iter = state_cache.find(state);
      if(iter != state_cache.end())
      {
        payloads.push_back(iter->second);
      }
      else
      {
        code_type haps = 0;
        for(gbwt::size_type seq_id : cached_gbwt.locate(state))
        {
          const gbwt::PathName& path = metadata.path(gbwt::Path::id(seq_id));
          auto id = path_ids_map.id(path);
          haps |= (code_type(1) << id);
        }
        payloads.push_back(haps);
        state_cache[state] = haps;
      }
    }

    // And now flush the cache.
    #pragma omp critical (minimizer_index)
    {
      for(size_t i = 0, payload_offset = 0; i < current_cache.size(); i++, payload_offset += payload_size)
      {
        value_type value(Position(std::get<1>(current_cache[i])), &payloads[payload_offset]);
        index.insert(std::get<0>(current_cache[i]), value);
      }
    }
    current_cache.clear();
  };

  // Minimizer finding logic.
  auto find_minimizers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    std::vector<minimizer_type> minimizers = index.minimizers(seq); // Calls syncmers() when appropriate.
    auto iter = traversal.begin();
    size_t path_offset = 0;
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();

    for(minimizer_type& minimizer : minimizers)
    {
      if (minimizer.empty()) { continue; }

      // Find the node covering minimizer starting position and the path corresponding to the minimizer.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= minimizer.offset)
      {
        node_start += node_length;
        ++iter; path_offset++;
        node_length = graph.get_length(*iter);
      }
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start };
      std::vector<gbwt::node_type> traversed_nodes = extract_kmer_path(graph, traversal, path_offset, offset(pos), index.k(), minimizer.is_reverse);
      if(minimizer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      if(!Position::valid_offset(pos))
      {
        #pragma omp critical (cerr)
        {
          std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large" << std::endl;
        }
        std::exit(EXIT_FAILURE);
      }

      cache[thread_id].data.emplace_back(minimizer, pos, traversed_nodes);
    }

    if(cache[thread_id].data.size() >= MINIMIZER_CACHE_SIZE) { flush_cache(thread_id); }
  };

  // Traverse all haplotype windows
  for_each_haplotype_window(graph, index.window_bp(), find_minimizers, (threads > 1));
  for(int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }
}

// Instantiate the template, as we did not define the function in the header.

template void index_haplotypes_with_paths<Key64>
(
  const GBZ& gbz, MinimizerIndex<Key64>& index,
  const std::function<const KmerEncoding::code_type*(const pos_t&)>& get_payload
);

template void index_haplotypes_with_paths<Key128>
(
  const GBZ& gbz, MinimizerIndex<Key128>& index,
  const std::function<const KmerEncoding::code_type*(const pos_t&)>& get_payload
);

//------------------------------------------------------------------------------

template<class KeyType>
std::vector<Kmer<KeyType>>
canonical_kmers(std::string::const_iterator begin, std::string::const_iterator end, size_t k)
{
  std::vector<Kmer<KeyType>> result;
  if(k == 0 || k > KeyType::KMER_MAX_LENGTH)
  {
    std::cerr << "canonical_kmers(): k must be between 1 and " << KeyType::KMER_MAX_LENGTH << std::endl;
    return result;
  }

  size_t valid_chars = 0, offset = 0;
  KeyType forward_key, reverse_key;
  std::string::const_iterator iter = begin;
  while(iter != end)
  {
    forward_key.forward(k, *iter, valid_chars);
    reverse_key.reverse(k, *iter);
    if(valid_chars >= k)
    {
      if(forward_key < reverse_key)
      {
        result.push_back({ forward_key, forward_key.hash(), offset_type(offset - (k - 1)), false });
      }
      else
      {
        result.push_back({ reverse_key, reverse_key.hash(), offset_type(offset), true });
      }
    }
    ++iter; offset++;
  }

  std::sort(result.begin(), result.end());
  return result;
}

// Instantiate the template, as we did not define the function in the header.

template std::vector<Kmer<Key64>>
canonical_kmers<Key64>(std::string::const_iterator begin, std::string::const_iterator end, size_t k);

template std::vector<Kmer<Key128>>
canonical_kmers<Key128>(std::string::const_iterator begin, std::string::const_iterator end, size_t k);

//------------------------------------------------------------------------------

template<class KeyType>
void build_kmer_index(const GBWTGraph& graph, KmerIndex<KeyType>& index, size_t k, const std::function<bool(KeyType)>& include)
{
  using kmer_type = Kmer<KeyType>;
  using code_type = KmerEncoding::code_type;
  using value_type = KmerEncoding::value_type;
  struct alignas(CACHE_LINE_SIZE) Cache
  {
    std::vector<std::pair<kmer_type, Position>> data;
  };

  // Kmer caching.
  int threads = omp_get_max_threads();
  std::vector<Cache> cache(threads);
  constexpr size_t KMER_CACHE_SIZE = 1024;
  std::vector<code_type> null_payload(index.payload_size(), 0);
  auto flush_cache = [&](int thread_id)
  {
    auto& current_cache = cache[thread_id].data;
    gbwt::removeDuplicates(current_cache, false);
    #pragma omp critical (kmer_index)
    {
      for(auto& kmer : current_cache)
      {
        value_type value(Position(kmer.second), null_payload.data());
        index.insert(kmer.first.key, value, kmer.first.hash);
      }
    }
    current_cache.clear();
  };

  // Kmer finding.
  auto find_kmers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    std::vector<kmer_type> kmers = canonical_kmers<KeyType>(seq, k);
    auto iter = traversal.begin();
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();
    for(auto kmer : kmers)
    {
      if(kmer.empty()) { continue; }

      // Find the node covering kmer starting position.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= kmer.offset)
      {
        node_start += node_length;
        ++iter;
        node_length = graph.get_length(*iter);
      }
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), kmer.offset - node_start };
      if(kmer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      if(include(kmer.key)) { cache[thread_id].data.emplace_back(kmer, Position(pos)); }
    }
    if(cache[thread_id].data.size() >= KMER_CACHE_SIZE) { flush_cache(thread_id); }
  };

  // Count all kmers.
  for_each_haplotype_window(graph, k, find_kmers, (threads > 1));
  for(int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }
}

// Instantiate the template, as we did not define the function in the header.

template void
build_kmer_index<Key64>(const GBWTGraph& graph, KmerIndex<Key64>& index, size_t k, const std::function<bool(Key64)>& include);

template void
build_kmer_index<Key128>(const GBWTGraph& graph, KmerIndex<Key128>& index, size_t k, const std::function<bool(Key128)>& include);

//------------------------------------------------------------------------------

template<class KeyType>
void build_kmer_indexes(const GBWTGraph& graph, std::array<KmerIndex<KeyType>, 4>& indexes, size_t k)
{
  size_t payload_size = indexes[0].payload_size();
  for(size_t i = 1; i < indexes.size(); i++)
  {
    if(indexes[i].payload_size() != payload_size)
    {
      throw std::runtime_error("build_kmer_indexes(): All indexes must have the same payload size");
    }
  }

  using kmer_type = Kmer<KeyType>;
  using code_type = KmerEncoding::code_type;
  using value_type = KmerEncoding::value_type;
  struct alignas(CACHE_LINE_SIZE) Cache
  {
    std::array<std::vector<std::pair<kmer_type, Position>>, 4> data;
  };

  // Kmer caching.
  int threads = omp_get_max_threads();
  std::array<std::mutex, 4> mutexes;
  std::vector<Cache> cache(threads);
  constexpr size_t KMER_CACHE_SIZE = 1024;
  std::vector<code_type> null_payload(payload_size, 0);
  auto flush_cache = [&](int thread_id, size_t index)
  {
    auto& current_cache = cache[thread_id].data[index];
    gbwt::removeDuplicates(current_cache, false);
    {
      std::lock_guard<std::mutex> lock(mutexes[index]);
      for(auto& kmer : current_cache)
      {
        value_type value(Position(kmer.second), null_payload.data());
        indexes[index].insert(kmer.first.key, value, kmer.first.hash);
      }
    }
    current_cache.clear();
  };

  // Kmer finding.
  auto find_kmers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    std::vector<kmer_type> kmers = canonical_kmers<KeyType>(seq, k);
    auto iter = traversal.begin();
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();
    for(auto kmer : kmers)
    {
      if(kmer.empty()) { continue; }

      // Find the node covering kmer starting position.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= kmer.offset)
      {
        node_start += node_length;
        ++iter;
        node_length = graph.get_length(*iter);
      }
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), kmer.offset - node_start };
      if(kmer.is_reverse) { pos = reverse_base_pos(pos, node_length); }

      // Insert the kmer into the right index.
      size_t index = kmer.key.access_raw(k, k / 2);
      if(index < indexes.size()) { cache[thread_id].data[index].emplace_back(kmer, Position(pos)); }
    }
    for(size_t i = 0; i < indexes.size(); i++)
    {
      if(cache[thread_id].data[i].size() >= KMER_CACHE_SIZE) { flush_cache(thread_id, i); }
    }
  };

  // Count all kmers.
  for_each_haplotype_window(graph, k, find_kmers, (threads > 1));
  for(int thread_id = 0; thread_id < threads; thread_id++)
  {
    for(size_t index = 0; index < indexes.size(); index++) { flush_cache(thread_id, index); }
  }
}

// Instantiate the template, as we did not define the function in the header.

template void
build_kmer_indexes<Key64>(const GBWTGraph& graph, std::array<KmerIndex<Key64>, 4>& indexes, size_t k);

template void
build_kmer_indexes<Key128>(const GBWTGraph& graph, std::array<KmerIndex<Key128>, 4>& indexes, size_t k);

//------------------------------------------------------------------------------

template<class KeyType>
std::vector<KeyType>
frequent_kmers(const GBWTGraph& graph, size_t k, size_t threshold, bool space_efficient, size_t hash_table_size)
{
  using index_type = KmerIndex<KeyType>;
  using multi_value_type = KmerEncoding::multi_value_type;

  std::vector<KeyType> result;
  auto select_frequent = [&](const index_type& index)
  {
    index.for_each_kmer([&](KeyType key, multi_value_type values)
    {
      if(threshold == 0 || values.second > threshold)
      {
        result.push_back(key);
        result.back().clear_pointer();
      }
    });
  };

  if(space_efficient)
  {
    for(char base : { 'A', 'C', 'G', 'T' })
    {
      index_type index(hash_table_size, 0);
      std::function<bool(KeyType)> include = [&](KeyType key) -> bool { return (key.access(k, k / 2) == base); };
      build_kmer_index(graph, index, k, include);
      select_frequent(index);
    }
  }
  else
  {
    std::array<index_type, 4> indexes
    {
      index_type(hash_table_size, 0),
      index_type(hash_table_size, 0),
      index_type(hash_table_size, 0),
      index_type(hash_table_size, 0),
    };
    build_kmer_indexes(graph, indexes, k);
    for(auto& index : indexes) { select_frequent(index); }
  }

  gbwt::parallelQuickSort(result.begin(), result.end());
  return result;
}

// Instantiate the template, as we did not define the function in the header.

template std::vector<Key64>
frequent_kmers<Key64>(const GBWTGraph& graph, size_t k, size_t threshold, bool space_efficient, size_t hash_table_size);

template std::vector<Key128>
frequent_kmers<Key128>(const GBWTGraph& graph, size_t k, size_t threshold, bool space_efficient, size_t hash_table_size);

//------------------------------------------------------------------------------

} // namespace gbwtgraph
