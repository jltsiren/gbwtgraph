#include <gtest/gtest.h>

#include <map>
#include <set>
#include <vector>

#include <gbwtgraph/index.h>
#include <gbwtgraph/internal.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

using KeyTypes = ::testing::Types<Key64, Key128>;

const std::vector<size_t> payload_sizes { 0, 1, 2, 3 };

//------------------------------------------------------------------------------

template<class KeyType>
class KmerExtraction : public ::testing::Test
{

};

TYPED_TEST_CASE(KmerExtraction, KeyTypes);

TYPED_TEST(KmerExtraction, CanonicalKmers)
{
  typedef TypeParam key_type;

  std::string sequence = "GATTACAGATTA";
  std::vector<Kmer<key_type>> truth
  {
    { key_type::encode("ATC"), 0, 2, true },
    { key_type::encode("AAT"), 0, 3, true },
    { key_type::encode("TAA"), 0, 4, true },
    { key_type::encode("GTA"), 0, 5, true },
    { key_type::encode("ACA"), 0, 4, false },
    { key_type::encode("CAG"), 0, 5, false },
    { key_type::encode("AGA"), 0, 6, false },
    { key_type::encode("ATC"), 0, 9, true },
    { key_type::encode("AAT"), 0, 10, true },
    { key_type::encode("TAA"), 0, 11, true },
  };
  std::sort(truth.begin(), truth.end());

  auto kmers = canonical_kmers<TypeParam>(sequence, 3);
  ASSERT_EQ(kmers.size(), truth.size()) << "Invalid number of kmers";
  for(size_t i = 0; i < truth.size(); i++)
  {
    EXPECT_EQ(kmers[i], truth[i]) << "Kmer " << i << ": expected " << truth[i] << ", got " << kmers[i];
  }
}

TYPED_TEST(KmerExtraction, InvalidCharacters)
{
  typedef TypeParam key_type;

  std::string sequence = "GATTAxAGATTA";
  std::vector<Kmer<key_type>> truth
  {
    { key_type::encode("ATC"), 0, 2, true },
    { key_type::encode("AAT"), 0, 3, true },
    { key_type::encode("TAA"), 0, 4, true },
    { key_type::encode("AGA"), 0, 6, false },
    { key_type::encode("ATC"), 0, 9, true },
    { key_type::encode("AAT"), 0, 10, true },
    { key_type::encode("TAA"), 0, 11, true },
  };
  std::sort(truth.begin(), truth.end());

  auto kmers = canonical_kmers<TypeParam>(sequence, 3);
  ASSERT_EQ(kmers.size(), truth.size()) << "Invalid number of kmers";
  for(size_t i = 0; i < truth.size(); i++)
  {
    EXPECT_EQ(kmers[i], truth[i]) << "Kmer " << i << ": expected " << truth[i] << ", got " << kmers[i];
  }
}

//------------------------------------------------------------------------------

template<class KeyType>
class IndexConstruction : public ::testing::Test
{
public:
  gbwt::GBWT index;
  SequenceSource source;
  GBWTGraph graph;

  void SetUp() override
  {
    this->index = build_gbwt_index();
    build_source(this->source);
    this->graph = GBWTGraph(this->index, this->source);
  }

  static std::vector<Kmer<KeyType>> get_kmers(const KmerIndex<KeyType>&, const std::string& str, size_t k)
  {
    return canonical_kmers<KeyType>(str, k);
  }

  static std::vector<Kmer<KeyType>> get_kmers(const MinimizerIndex<KeyType>& index, const std::string& str, size_t)
  {
    return index.minimizers(str);
  }

  template<class IndexType>
  void insert_values(const IndexType& index, const gbwt::vector_type& path, std::map<KeyType, std::set<owned_value_type>>& result, size_t k) const
  {
    // Convert the path to a string and find the minimizers.
    std::string str = path_to_string(this->graph, path);
    auto kmers = get_kmers(index, str, k);

    // Insert the minimizers into the result.
    auto iter = path.begin();
    size_t node_start = 0;
    for(auto kmer : kmers)
    {
      if(kmer.empty()) { continue; }
      handle_t handle = GBWTGraph::node_to_handle(*iter);
      size_t node_length = this->graph.get_length(handle);
      while(node_start + node_length <= kmer.offset)
      {
        node_start += node_length;
        ++iter;
        handle = GBWTGraph::node_to_handle(*iter);
        node_length = this->graph.get_length(handle);
      }
      pos_t pos { this->graph.get_id(handle), this->graph.get_is_reverse(handle), kmer.offset - node_start };
      if(kmer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      owned_value_type value = create_value(pos, index.payload_size(), hash(pos));
      result[kmer.key].insert(value);
    }
  }

  // This version stores path information in the last word of payload.
  template<class IndexType>
  void insert_values(const IndexType& index, const std::vector<gbwt::vector_type>& paths, std::map<KeyType, std::set<owned_value_type>>& result, size_t k) const
  {
    ASSERT_GE(index.payload_size(), 1) << "Index must have payload to store path information";

    for(size_t path_id = 0; path_id < paths.size(); path_id++)
    {
      // Convert the path to a string and find the minimizers.
      const gbwt::vector_type& path = paths[path_id];
      std::string str = path_to_string(this->graph, path);
      auto kmers = get_kmers(index, str, k);

      // Insert the minimizers into the result.
      auto iter = path.begin();
      size_t node_start = 0;
      for(auto kmer : kmers)
      {
        if(kmer.empty()) { continue; }
        handle_t handle = GBWTGraph::node_to_handle(*iter);
        size_t node_length = this->graph.get_length(handle);
        while(node_start + node_length <= kmer.offset)
        {
          node_start += node_length;
          ++iter;
          handle = GBWTGraph::node_to_handle(*iter);
          node_length = this->graph.get_length(handle);
        }
        pos_t pos { this->graph.get_id(handle), this->graph.get_is_reverse(handle), kmer.offset - node_start };
        if(kmer.is_reverse) { pos = reverse_base_pos(pos, node_length); }

        // This value will be inserted/updated in the set.
        owned_value_type value = create_value(pos, index.payload_size(), hash(pos));
        value.second.back() = 1 << path_id;

        auto found = result.find(kmer.key);
        if(found == result.end())
        {
          // Case 1: first occurrence of the kmer.
          result[kmer.key].insert(value);
          continue;
        }
        bool same_hit = false;
        owned_value_type old_value;
        for(auto& existing_value : found->second)
        {
          if(existing_value.first == value.first)
          {
            // Case 2: same position already exists; merge the sets of paths.
            value.second.back() |= existing_value.second.back();
            same_hit = true;
            old_value = existing_value;
            break;
          }
        }
        if(same_hit)
        {
          // Update the existing value with the new path information.
          found->second.erase(old_value);
          found->second.insert(value);
        }
        else
        {
          // Case 3: new position for the kmer.
          result[kmer.key].insert(value);
        }
      }
    }
  }

  std::map<KeyType, std::set<owned_value_type>> filter_values(const std::map<KeyType, std::set<owned_value_type>>& source, char middle_base) const
  {
    // NOTE: This assumes 3-mers.
    std::map<KeyType, std::set<owned_value_type>> result;
    for(auto iter = source.begin(); iter != source.end(); ++iter)
    {
      if(iter->first.access(3, 1) == middle_base) { result[iter->first] = iter->second; }
    }
    return result;
  }

  static KmerEncoding::multi_value_type find_values(const KmerIndex<KeyType>& index, KeyType key)
  {
    return index.find(key);
  }

  static KmerEncoding::multi_value_type find_values(const MinimizerIndex<KeyType>& index, KeyType key)
  {
    return index.find(get_minimizer<KeyType>(key));
  }

  template<class IndexType>
  void check_index(const IndexType& index, const std::map<KeyType, std::set<owned_value_type>>& correct_values, bool with_paths) const
  {
    using code_type = KmerEncoding::code_type;
    std::string payload_msg = " with payload size " + std::to_string(index.payload_size());

    size_t values = 0;
    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      values += iter->second.size();
    }
    ASSERT_EQ(index.size(), correct_values.size()) << "Wrong number of keys" << payload_msg;
    ASSERT_EQ(index.number_of_values(), values) << "Wrong number of values" << payload_msg;

    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      KmerEncoding::multi_value_type values = find_values(index, iter->first);
      EXPECT_TRUE(same_values(values, iter->second, index.payload_size(), with_paths)) << "Wrong values for key " << iter->first << payload_msg;
    }
  }
};

TYPED_TEST_CASE(IndexConstruction, KeyTypes);

TYPED_TEST(IndexConstruction, WithoutPayload)
{
  // Determine the correct minimizer occurrences.
  MinimizerIndex<TypeParam> index(3, 2, 0);
  std::map<TypeParam, std::set<owned_value_type>> correct_values;
  this->insert_values(index, alt_path, correct_values, index.k());
  this->insert_values(index, short_path, correct_values, index.k());

  // Check that we managed to index them.
  index_haplotypes(this->graph, index, [](const pos_t&) { return nullptr; });
  this->check_index(index, correct_values, false);
}

TYPED_TEST(IndexConstruction, WithPayload)
{
  using key_type = TypeParam;
  using index_type = MinimizerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    // Determine the correct minimizer occurrences.
    index_type index(3, 2, payload_size);
    std::map<TypeParam, std::set<owned_value_type>> correct_values;
    this->insert_values(index, alt_path, correct_values, index.k());
    this->insert_values(index, short_path, correct_values, index.k());

    // Extract the payloads.
    std::map<pos_t, std::vector<std::uint64_t>> payloads;
    for(auto& pair : correct_values)
    {
      for(auto& value : pair.second)
      {
        payloads[value.first.decode()] = value.second;
      }
    }

    // Check that we managed to index them.
    index_haplotypes(this->graph, index, [&](const pos_t& pos) { return payloads[pos].data(); });
    this->check_index(index, correct_values, false);
  }
}

TYPED_TEST(IndexConstruction, WithPaths)
{
  using key_type = TypeParam;
  using index_type = MinimizerIndex<key_type>;

  // We need path metadata for indexing with paths.
  ASSERT_EQ(this->index.sequences(), 6) << "Metadata generation code does not match test index";
  this->index.addMetadata();
  this->index.metadata.setSamples(3);
  this->index.metadata.setHaplotypes(3);
  this->index.metadata.setContigs(1);
  for(size_t path_id = 0; path_id < 3; path_id++)
  {
    this->index.metadata.addPath(path_id, 0, 0, 0);
  }

  for(size_t payload_size : payload_sizes)
  {
    if(payload_size == 0)
    {
      index_type index(3, 2, payload_size);
      EXPECT_THROW
      (
        index_haplotypes_with_paths(this->graph, index, [](const pos_t&) { return nullptr; }),
        std::runtime_error
      ) << "Expected exception when indexing with paths and zero payload size";
      continue;
    }

    // Determine the correct minimizer occurrences.
    index_type index(3, 2, payload_size);
    std::map<TypeParam, std::set<owned_value_type>> correct_values;
    std::vector<gbwt::vector_type> paths { short_path, alt_path, short_path };
    this->insert_values(index, paths, correct_values, index.k());

    // Extract the payloads.
    std::map<pos_t, std::vector<std::uint64_t>> payloads;
    for(auto& pair : correct_values)
    {
      for(auto& value : pair.second)
      {
        // We may have different payloads for different extensions of the same position.
        // But because we only use the shared prefix (of length payload_size() - 1),
        // we can overwrite the existing value.
        payloads[value.first.decode()] = value.second;
      }
    }

    // Check that we managed to index them.
    index_haplotypes_with_paths(this->graph, index, [&](const pos_t& pos) { return payloads[pos].data(); });
    this->check_index(index, correct_values, true);
  }
}

TYPED_TEST(IndexConstruction, CanonicalKmers)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  // Determine the correct canonical kmers.
  index_type index(0);
  std::map<TypeParam, std::set<owned_value_type>> correct_values;
  this->insert_values(index, alt_path, correct_values, 3);
  this->insert_values(index, short_path, correct_values, 3);

  // Check that we managed to index them.
  std::function<bool(key_type)> include = [](key_type) -> bool { return true; };
  build_kmer_index(this->graph, index, 3, include);
  this->check_index(index, correct_values, false);
}

TYPED_TEST(IndexConstruction, CanonicalKmersByMiddleBase)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  // Determine the correct canonical kmers with C as the middle base.
  index_type index(0);
  std::map<TypeParam, std::set<owned_value_type>> all_values;
  this->insert_values(index, alt_path, all_values, 3);
  this->insert_values(index, short_path, all_values, 3);
  std::map<TypeParam, std::set<owned_value_type>> correct_values = this->filter_values(all_values, 'C');

  // Check that we managed to index them.
  std::function<bool(key_type)> include = [&](key_type key) -> bool { return (key.access(3, 1) == 'C'); };
  build_kmer_index(this->graph, index, 3, include);
  this->check_index(index, correct_values, false);
}

TYPED_TEST(IndexConstruction, MultipleKmerIndexes)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  // Determine the correct canonical kmers.
  std::array<index_type, 4> indexes { index_type(0), index_type(0), index_type(0), index_type(0) };
  std::map<TypeParam, std::set<owned_value_type>> all_values;
  this->insert_values(indexes[0], alt_path, all_values, 3);
  this->insert_values(indexes[0], short_path, all_values, 3);

  // Build the indexes.
  build_kmer_indexes(this->graph, indexes, 3);

  // Check that the kmers were partitioned correctly.
  std::string bases = "ACGT";
  for(size_t i = 0; i < bases.length(); i++)
  {
    std::map<TypeParam, std::set<owned_value_type>> correct_values = this->filter_values(all_values, bases[i]);
    this->check_index(indexes[i], correct_values, false);
  }
}

//------------------------------------------------------------------------------

template<class KeyType>
class KmerCounting : public ::testing::Test
{
public:
  gbwt::GBWT index;
  SequenceSource source;
  GBWTGraph graph;

  void SetUp() override
  {
    this->index = build_gbwt_index();
    build_source(this->source);
    this->graph = GBWTGraph(this->index, this->source);
  }
};

TYPED_TEST_CASE(KmerCounting, KeyTypes);

TYPED_TEST(KmerCounting, FrequentKmers)
{
  // alt_path:   G A GGG T A - A A
  // short_path: G - GGG T A C - A
  // CCC: two reverse occurrences
  // GTA: one forward and one reverse occurrence

  size_t k = 3;
  std::vector<TypeParam> correct;
  correct.push_back(TypeParam::encode("CCC"));
  correct.push_back(TypeParam::encode("GTA"));

  auto kmers = frequent_kmers<TypeParam>(this->graph, k, 1, false);
  ASSERT_EQ(kmers.size(), correct.size()) << "Invalid number of frequent kmers";
  for(size_t i = 0; i < kmers.size(); i++)
  {
    EXPECT_EQ(kmers[i], correct[i]) << "Invalid kmer " << i << ": expected " << correct[i].decode(k) << ", got " << kmers[i].decode(k);
  }

  auto space_efficient = frequent_kmers<TypeParam>(this->graph, k, 1, true);
  ASSERT_EQ(space_efficient, kmers) << "Invalid frequent kmers using the space-efficient algorithm";
}

//------------------------------------------------------------------------------

class PathIdTest : public ::testing::Test
{
  static void shuffle(std::vector<gbwt::PathName>& paths, size_t seed = 0xACDCABBACAFEBEEF)
  {
    if(paths.empty()) { return; }
    std::mt19937_64 rng(seed);
    for(size_t i = paths.size() - 1; i > 0; i--)
    {
      size_t j = rng() % (i + 1);
      std::swap(paths[i], paths[j]);
    }
  }

  gbwt::Metadata generate_test_case(size_t samples, size_t haplotypes, size_t contigs, size_t fragments) const
  {
    gbwt::Metadata metadata;
    metadata.setSamples(samples);
    metadata.setHaplotypes(samples * haplotypes);
    metadata.setContigs(contigs);

    std::vector<gbwt::PathName> paths;
    for(size_t s = 0; s < samples; s++)
    {
      for(size_t c = 0; c < contigs; c++)
      {
        for(size_t h = 0; h < haplotypes; h++)
        {
          for(size_t f = 0; f < fragments; f++)
          {
            gbwt::PathName path(s, c, h, f);
            paths.push_back(path);
          }
        }
      }
    }
    shuffle(paths); // Shuffle to test that order does not matter.
    for(const gbwt::PathName& path : paths) { metadata.addPath(path); }

    return metadata;
  }

  std::unordered_map<gbwt::PathName, size_t, PathNameHasher> generate_truth(size_t samples, size_t haplotypes, size_t contigs, size_t fragments) const
  {
    std::unordered_map<gbwt::PathName, size_t, PathNameHasher> truth;

    bool use_contigs = (samples * haplotypes * contigs <= PathIdMap::MAX_HAPLOTYPES);
    bool use_haplotypes = (samples * haplotypes <= PathIdMap::MAX_HAPLOTYPES) & !use_contigs;
    bool use_samples = (samples <= PathIdMap::MAX_HAPLOTYPES) & !use_haplotypes & !use_contigs;

    // Note that the iteration order does not match the order of the fields.
    // Haplotypes are more significant than contigs.
    size_t next = 0;
    for(size_t s = 0; s < samples; s++)
    {
      for(size_t h = 0; h < haplotypes; h++)
      {
        for(size_t c = 0; c < contigs; c++)
        {
          for(size_t f = 0; f < fragments; f++)
          {
            gbwt::PathName path(s, c, h, f);
            truth[path] = next;
          }
          if(use_contigs) { next++; }
        }
        if(use_haplotypes) { next++; }
      }
      if(use_samples) { next++; }
    }

    return truth;
  }

public:
  void path_id_map_test(size_t samples, size_t haplotypes, size_t contigs, size_t fragments) const
  {
    std::string test_case = " with (" + std::to_string(samples) + " samples, "
      + std::to_string(haplotypes) + " haplotypes/sample, "
      + std::to_string(contigs) + " contigs, "
      + std::to_string(fragments) + " fragments/sequence)";

    gbwt::Metadata metadata = this->generate_test_case(samples, haplotypes, contigs, fragments);
    PathIdMap path_id_map(metadata);
    auto truth = this->generate_truth(samples, haplotypes, contigs, fragments);

    for(gbwt::size_type i = 0; i < metadata.paths(); i++)
    {
      const gbwt::PathName& path = metadata.path(i);
      size_t expected_id = truth[path];
      size_t actual_id = path_id_map.id(path);
      EXPECT_LT(actual_id, PathIdMap::MAX_HAPLOTYPES) << "Too large path id for path " << path_name_to_string(path) << test_case;
      EXPECT_EQ(actual_id, expected_id) << "Wrong path id for path " << path_name_to_string(path) << test_case;
    }

    gbwt::PathName missing(samples, 0, 0, 0);
    size_t missing_id = path_id_map.id(missing);
    EXPECT_EQ(missing_id, size_t(0)) << "Non-zero path id for missing path " << path_name_to_string(missing) << test_case;
  }

  std::vector<handle_t> get_path(const gbwt::GBWT& index, gbwt::size_type seq_id) const
  {
    gbwt::vector_type gbwt_path = index.extract(seq_id);
    std::vector<handle_t> path;
    for(size_t i = 0; i < gbwt_path.size(); i++)
    {
      path.push_back(GBWTGraph::node_to_handle(gbwt_path[i]));
    }
    return path;
  }

  std::vector<size_t> path_offset_by_seq_offset(const GBWTGraph& graph, const std::vector<handle_t>& path) const
  {
    std::vector<size_t> result;
    for(size_t i = 0; i < path.size(); i++)
    {
      handle_t handle = path[i];
      size_t node_length = graph.get_length(handle);
      for(size_t j = 0; j < node_length; j++)
      {
        result.push_back(i);
      }
    }
    return result;
  }

  std::vector<gbwt::node_type> true_kmer_path(const std::vector<handle_t>& path, const std::vector<size_t>& expanded, size_t seq_offset, size_t k, bool is_reverse) const
  {
    std::vector<gbwt::node_type> result;
    if(is_reverse)
    {
      size_t prev = path.size();
      for(size_t i = 0; i < k && i <= seq_offset; i++)
      {
        if(expanded[seq_offset - i] == prev) { continue; }
        prev = expanded[seq_offset - i];
        result.push_back(gbwt::Node::reverse(GBWTGraph::handle_to_node(path[prev])));
      }
    }
    else
    {
      size_t prev = path.size();
      for(size_t i = 0; i < k && seq_offset + i < expanded.size(); i++)
      {
        if(expanded[seq_offset + i] == prev) { continue; }
        prev = expanded[seq_offset + i];
        result.push_back(GBWTGraph::handle_to_node(path[prev]));
      }
    }
    return result;
  }
};

TEST_F(PathIdTest, PathIdMap)
{
  this->path_id_map_test(0, 0, 0, 0);  // Empty metadata.
  this->path_id_map_test(1, 1, 1, 1);  // Single path.
  this->path_id_map_test(2, 2, 2, 2);  // Multiple paths.
  this->path_id_map_test(4, 4, 4, 1);  // Should still use (sample, haplotype, contig).
  this->path_id_map_test(5, 1, 13, 1); // Fall back to (sample, haplotype).
  this->path_id_map_test(32, 2, 2, 1); // Should still use (sample, haplotype).
  this->path_id_map_test(13, 5, 1, 1); // Fall back to (sample).
  this->path_id_map_test(230, 2, 24, 3); // Large test case, where everything should map to 0.
}

TEST_F(PathIdTest, ExtractKmerPath)
{
  gbwt::GBWT index = build_gbwt_index();
  SequenceSource source;
  build_source(source);
  GBWTGraph graph(index, source);

  size_t k = 3;
  for(gbwt::size_type seq_id = 0; seq_id < index.sequences(); seq_id++)
  {
    std::vector<handle_t> path = this->get_path(index, seq_id);
    std::vector<size_t> expanded = this->path_offset_by_seq_offset(graph, path);
    size_t seq_offset = 0;
    for(size_t path_offset = 0; path_offset < path.size(); path_offset++)
    {
      size_t node_length = graph.get_length(path[path_offset]);
      for(size_t node_offset = 0; node_offset < node_length; node_offset++, seq_offset++)
      {

        for(bool is_reverse : { false, true })
        {
          std::string test_case = " for sequence " + std::to_string(seq_id)
            + ", path offset " + std::to_string(path_offset)
            + ", node offset " + std::to_string(node_offset)
            + (is_reverse ? " (reverse)" : " (forward)");
          std::vector<gbwt::node_type> kmer_path = extract_kmer_path(graph, path, path_offset, node_offset, k, is_reverse);
          std::vector<gbwt::node_type> truth = this->true_kmer_path(path, expanded, seq_offset, k, is_reverse);
          EXPECT_EQ(kmer_path, truth) << "Wrong kmer path" << test_case;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

} // namespace
