#include <gtest/gtest.h>

#include <map>
#include <set>
#include <vector>

#include <gbwtgraph/index.h>

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
  void check_index(const IndexType& index, const std::map<KeyType, std::set<owned_value_type>>& correct_values) const
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
      EXPECT_TRUE(same_values(values, iter->second, index.payload_size())) << "Wrong values for key " << iter->first << payload_msg;
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
  this->check_index(index, correct_values);
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
    this->check_index(index, correct_values);
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
  this->check_index(index, correct_values);
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
  this->check_index(index, correct_values);
}

TYPED_TEST(IndexConstruction, MultipleKmerIndexes)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  // Determine the correct canonical kmers.
  std::array<index_type, 4> indexes{ index_type(0), index_type(0), index_type(0), index_type(0) };
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
    this->check_index(indexes[i], correct_values);
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

} // namespace
