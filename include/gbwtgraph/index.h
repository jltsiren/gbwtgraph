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

/*
  Index the haplotypes in the graph. Insert the minimizers into the provided
  index. Function argument get_payload is used to generate the payload for each
  position stored in the index. The payload is provided as a pointer to
  index.payload_size() words of payload. The pointer must remain valid for the
  duration of the call. If payload size is zero, the pointer may be null.

  The number of threads can be set through OpenMP.

  We do a lot of redundant work by traversing both orientations and finding
  almost the same minimizers in each orientation. If we consider only the
  windows starting in forward (reverse) orientation, we may miss windows that
  cross from a reverse node to a forward node (from a forward node to a reverse
  node).
*/
template<typename KeyType>
void index_haplotypes
(
  const GBWTGraph& graph, MinimizerIndex<KeyType>& index,
  const std::function<const KmerEncoding::code_type*(const pos_t&)>& get_payload
);

// FIXME: better documentation
/*
  Index the haplotypes in the graph. This version requires that
  index.payload_size() > 0. It gets the first index.payload_size() - 1 words of
  payload from the get_payload function and uses the last word to store a
  representation of the haplotypes (paths) that the minimizer occurs in.

  Throws std::runtime_error if index.payload_size() == 0.
*/
template<typename KeyType>
void index_haplotypes_with_paths
(
  const GBWTGraph& graph, MinimizerIndex<KeyType>& index,
  const std::function<const KmerEncoding::code_type*(const pos_t&)>& get_payload
);

//------------------------------------------------------------------------------

/*
  Returns all canonical kmers in the string specified by the iterators. A
  canonical kmer is the smallest of a kmer and its reverse complement. The
  return value is a vector of kmers sorted by their offsets.
*/
template<class KeyType>
std::vector<Kmer<KeyType>>
canonical_kmers(std::string::const_iterator begin, std::string::const_iterator end, size_t k);

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
  if include(key) returns true. The number of threads can be set through
  OpenMP.

  This is mostly intended for indexes without payload. If payload size is
  nonzero, the payload will be empty.
*/
template<class KeyType>
void build_kmer_index(const GBWTGraph& graph, KmerIndex<KeyType>& index, size_t k, const std::function<bool(KeyType)>& include);

/*
  Index the haplotypes in the graph. Insert the kmers into the provided indexes
  according to the middle base. The number of threads can be set through OpenMP.

  This is mostly intended for indexes without payload. If payload size is
  nonzero, the payload will be empty.

  Throws std::runtime_error if the indexes do not have the same payload size.
*/
template<class KeyType>
void build_kmer_indexes(const GBWTGraph& graph, std::array<KmerIndex<KeyType>, 4>& indexes, size_t k);

/*
  Returns all kmers in the haplotypes with more than threshold hits in the
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
frequent_kmers(const GBWTGraph& graph, size_t k, size_t threshold, bool space_efficient, size_t hash_table_size = KmerIndex<KeyType>::INITIAL_CAPACITY);

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_CONSTRUCTION_H
