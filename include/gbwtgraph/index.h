#ifndef GBWTGRAPH_CONSTRUCTION_H
#define GBWTGRAPH_CONSTRUCTION_H

#include <array>
#include <cstdlib>
#include <functional>
#include <mutex>

#include <omp.h>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/minimizer.h>
#include <gbwtgraph/internal/path_names_index.hpp>
/*
  index.h: Minimizer index construction from GBWTGraph.
*/

namespace gbwtgraph
{
using detail::SearchStateKey;
using detail::PathIDMap;
//------------------------------------------------------------------------------

// TODO: These algorithms are basically the same. Is there a clean way of
// merging the implementations?

/*
  Index the haplotypes in the graph. Insert the minimizers into the provided
  index. Function argument get_payload is used to generate the payload for each
  position stored in the index. The number of threads can be set through OpenMP.

  We do a lot of redundant work by traversing both orientations and finding
  almost the same minimizers in each orientation. If we consider only the
  windows starting in forward (reverse) orientation, we may miss windows that
  cross from a reverse node to a forward node (from a forward node to a reverse
  node).
*/
template<typename KeyType, typename PayloadType>
void index_haplotypes(const GBWTGraph& graph,
                      MinimizerIndex<KeyType, PositionPayload<PayloadType>>& index,
                      const std::function<PayloadType(const pos_t&)>& get_payload)
{
    using minimizer_type = typename MinimizerIndex<KeyType, PositionPayload<PayloadType>>::minimizer_type;
    //typedef typename MinimizerIndex<KeyType, PositionPayload<PayloadType>::minimizer_type minimizer_type;

    struct alignas(CACHE_LINE_SIZE) Cache
    {
        std::vector<std::tuple<minimizer_type, pos_t, std::vector<gbwt::node_type>>> data;
    };

    int threads = omp_get_max_threads();
    constexpr size_t MINIMIZER_CACHE_SIZE = 1024;

    std::vector<Cache> cache(threads);
    auto meta = graph.index->metadata;
    auto path_ids_map = PathIDMap(meta);

    // Flush cache for one thread
    auto flush_cache = [&](int thread_id) {
        std::unordered_map<SearchStateKey, uint64_t> state_cache;
        auto cached_gbwt = graph.get_cache();
        auto& current_cache = cache[thread_id].data;
        gbwt::removeDuplicates(current_cache, false);

        std::vector<PayloadType> payloads;
        payloads.reserve(current_cache.size());

        for (auto& minimizer : current_cache)
        {
            auto data = get_payload(std::get<1>(minimizer));
            auto& traversed_nodes = std::get<2>(minimizer);

            auto state = cached_gbwt.find(traversed_nodes.begin(), traversed_nodes.end());
            SearchStateKey key{state.node, state.range};
            uint64_t haps = 0;

            if (auto it = state_cache.find(key); it != state_cache.end())
            {
                haps = it->second;
            }
            else
            {
                for (auto h : cached_gbwt.locate(state))
                {
                    auto path = meta.path(h / 2);
                    auto id = path_ids_map.id(path);
                    haps |= (1UL << id);
                }
                state_cache[key] = haps;
            }
            data.paths = haps;
            payloads.push_back(data);
        }

        #pragma omp critical (minimizer_index)
        {
            for (size_t i = 0; i < current_cache.size(); ++i)
            {
                index.insert(std::get<0>(current_cache[i]),
                             {Position::encode(std::get<1>(current_cache[i])), payloads[i]});
            }
        }

        current_cache.clear();
    };

    // Minimizer finding logic
    auto find_minimizers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
    {
        std::vector<minimizer_type> minimizers = index.minimizers(seq);
        auto iter = traversal.begin();
        size_t node_start = 0;
        int thread_id = omp_get_thread_num();

        for (minimizer_type& minimizer : minimizers)
        {
            if (minimizer.empty()) continue;

            size_t node_length = graph.get_length(*iter);
            while (node_start + node_length <= minimizer.offset)
            {
                node_start += node_length;
                ++iter;
                node_length = graph.get_length(*iter);
            }

            pos_t pos{graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start};
            if (minimizer.is_reverse)
                pos = reverse_base_pos(pos, node_length);

            if (!Position::valid_offset(pos))
            {
                #pragma omp critical (cerr)
                std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large" << std::endl;
                std::exit(EXIT_FAILURE);
            }

            if constexpr (std::is_same_v<PayloadType, Payload>)
            {
                // Simple payload: no path traversal
                cache[thread_id].data.emplace_back(minimizer, pos, std::vector<gbwt::node_type>{});
            }
            else
            {
                // Extended payload: collect traversed nodes
                auto node_iter = iter;
                size_t total_length = node_length;
                std::vector<gbwt::node_type> traversed_nodes{graph.handle_to_node(*node_iter)};

                const auto extend = [&](auto step_fwd) {
                    while (total_length < index.k() && step_fwd())
                    {
                        total_length += graph.get_length(*node_iter);
                        traversed_nodes.push_back(graph.handle_to_node(*node_iter));
                    }
                };

                if (minimizer.is_reverse)
                {
                    extend([&] {
                        if (node_iter == traversal.begin()) return false;
                        --node_iter;
                        return true;
                    });
                    std::reverse(traversed_nodes.begin(), traversed_nodes.end());
                }
                else
                {
                    extend([&] {
                        if (std::next(node_iter) == traversal.end()) return false;
                        ++node_iter;
                        return true;
                    });
                }

                cache[thread_id].data.emplace_back(minimizer, pos, traversed_nodes);
            }
        }

        if (cache[thread_id].data.size() >= MINIMIZER_CACHE_SIZE)
        {
            flush_cache(thread_id);
        }
    };

    // Traverse all haplotype windows
    for_each_haplotype_window(graph, index.window_bp(), find_minimizers, (threads > 1));

    for (int thread_id = 0; thread_id < threads; ++thread_id)
    {
        flush_cache(thread_id);
    }
}

//------------------------------------------------------------------------------

/*
  Index the haplotypes in the graph. Insert the minimizers into the provided
  index. This version is used for minimizer indexes without payloads. The number
  of threads can be set through OpenMP.
*/
template<class KeyType>
void
index_haplotypes(const GBWTGraph& graph, MinimizerIndex<KeyType, Position>& index)
{
  typedef typename MinimizerIndex<KeyType, Position>::minimizer_type minimizer_type;
  struct alignas(CACHE_LINE_SIZE) Cache
  {
    std::vector<std::pair<minimizer_type, Position>> data;
  };

  int threads = omp_get_max_threads();

  // Minimizer caching. We only generate the payloads after we have removed duplicate positions.
  std::vector<Cache> cache(threads);
  constexpr size_t MINIMIZER_CACHE_SIZE = 1024;
  auto flush_cache = [&](int thread_id)
  {
    auto& current_cache = cache[thread_id].data;
    gbwt::removeDuplicates(current_cache, false);
    #pragma omp critical (minimizer_index)
    {
      for(auto& minimizer : current_cache)
      {
        index.insert(minimizer.first, minimizer.second);
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
      cache[thread_id].data.emplace_back(minimizer, Position::encode(pos));
    }
    if(cache[thread_id].data.size() >= MINIMIZER_CACHE_SIZE) { flush_cache(thread_id); }
  };

  for_each_haplotype_window(graph, index.window_bp(), find_minimizers, (threads > 1));
  for(int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }
}

//------------------------------------------------------------------------------

/*
  Returns all canonical kmers in the string specified by the iterators. A
  canonical kmer is the smallest of a kmer and its reverse complement. The
  return value is a vector of kmers sorted by their offsets.
*/
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

/*
  Returns all canonical kmers in the string. A canonical kmer is the smaller
  of a kmer and its reverse complement. The return value is a vector of kmers
  sorted by their offsets.
*/
template<class KeyType>
std::vector<Kmer<KeyType>>
canonical_kmers(const std::string& seq, size_t k)
{
  return canonical_kmers<KeyType>(seq.begin(), seq.end(), k);
}

//------------------------------------------------------------------------------

/*
  Index the haplotypes in the graph. Insert the kmers into the provided index
  if `include(key)` returns `true`. The number of threads can be set through
  OpenMP.
*/
template<class KeyType, class Predicate>
void
build_kmer_index(const GBWTGraph& graph, KmerIndex<KeyType, Position>& index, size_t k, const Predicate& include)
{
  typedef KeyType key_type;
  typedef Kmer<key_type> kmer_type;
  struct alignas(CACHE_LINE_SIZE) Cache
  {
    std::vector<std::pair<kmer_type, Position>> data;
  };

  // Kmer caching.
  int threads = omp_get_max_threads();
  std::vector<Cache> cache(threads);
  constexpr size_t KMER_CACHE_SIZE = 1024;
  auto flush_cache = [&](int thread_id)
  {
    auto& current_cache = cache[thread_id].data;
    gbwt::removeDuplicates(current_cache, false);
    #pragma omp critical (kmer_index)
    {
      for(auto& kmer : current_cache)
      {
        index.insert(kmer.first.key, kmer.second, kmer.first.hash);
      }
    }
    current_cache.clear();
  };

  // Kmer finding.
  auto find_kmers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    std::vector<kmer_type> kmers = canonical_kmers<key_type>(seq, k);
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
      if(include(kmer.key)) { cache[thread_id].data.emplace_back(kmer, Position::encode(pos)); }
    }
    if(cache[thread_id].data.size() >= KMER_CACHE_SIZE) { flush_cache(thread_id); }
  };

  // Count all kmers.
  for_each_haplotype_window(graph, k, find_kmers, (threads > 1));
  for(int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }
}

//------------------------------------------------------------------------------

/*
  Index the haplotypes in the graph. Insert the kmers into the provided indexes
  according to the middle base. The number of threads can be set through OpenMP.
*/
template<class KeyType>
void
build_kmer_indexes(const GBWTGraph& graph, std::array<KmerIndex<KeyType, Position>, 4>& indexes, size_t k)
{
  typedef KeyType key_type;
  typedef Kmer<key_type> kmer_type;
  struct alignas(CACHE_LINE_SIZE) Cache
  {
    std::array<std::vector<std::pair<kmer_type, Position>>, 4> data;
  };

  // Kmer caching.
  int threads = omp_get_max_threads();
  std::array<std::mutex, 4> mutexes;
  std::vector<Cache> cache(threads);
  constexpr size_t KMER_CACHE_SIZE = 1024;
  auto flush_cache = [&](int thread_id, size_t index)
  {
    auto& current_cache = cache[thread_id].data[index];
    gbwt::removeDuplicates(current_cache, false);
    {
      std::lock_guard<std::mutex> lock(mutexes[index]);
      for(auto& kmer : current_cache)
      {
        indexes[index].insert(kmer.first.key, kmer.second, kmer.first.hash);
      }
    }
    current_cache.clear();
  };

  // Kmer finding.
  auto find_kmers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    std::vector<kmer_type> kmers = canonical_kmers<key_type>(seq, k);
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
      if(index < indexes.size()) { cache[thread_id].data[index].emplace_back(kmer, Position::encode(pos)); }
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

//------------------------------------------------------------------------------

/*
  Returns all kmers in the haplotypes with more than `threshold` hits in the
  graph in sorted order. Considers a kmer and its reverse complement the same.
  The number of threads can be set through OpenMP.  If the number of kmers can
  be estimated in advance, providing a hash table size can save time and
  memory.

  There are two versions of this algorithm:

  * The fast version (default) uses four kmer indexes with the given hash
    table size and inserts each kmer into the index determined by the
    middle base. This version uses more memory and parallelizes better.

  * The space-efficient version does four passes over the graph. In each
    pass, it only considers the kmers with a specific middle base. This
    versions uses less memory, is slower, and does not parallelize as
    well.
*/
template<class KeyType>
std::vector<KeyType>
frequent_kmers(const GBWTGraph& graph, size_t k, size_t threshold, bool space_efficient, size_t hash_table_size = KmerIndex<KeyType, Position>::INITIAL_CAPACITY)
{
  typedef KmerIndex<KeyType, Position> index_type;
  typedef typename index_type::cell_type cell_type;

  std::vector<KeyType> result;
  auto select_frequent = [&](const index_type& index)
  {
    index.for_each_kmer([&](const cell_type& cell)
    {
      if(threshold == 0 || (cell.first.is_pointer() && cell.second.pointer->size() > threshold))
      {
        result.push_back(cell.first);
        result.back().clear_pointer();
      }
    });
  };

  if(space_efficient)
  {
    for(char base : { 'A', 'C', 'G', 'T' })
    {
      index_type index(hash_table_size);
      build_kmer_index(graph, index, k, [&](KeyType key) -> bool { return (key.access(k, k / 2) == base); });
      select_frequent(index);
    }
  }
  else
  {
    std::array<index_type, 4> indexes
    {
      index_type(hash_table_size), index_type(hash_table_size), index_type(hash_table_size), index_type(hash_table_size),
    };
    build_kmer_indexes(graph, indexes, k);
    for(auto& index : indexes) { select_frequent(index); }
  }

  gbwt::parallelQuickSort(result.begin(), result.end());
  return result;
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_CONSTRUCTION_H
