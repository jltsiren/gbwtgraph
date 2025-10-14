#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <gbwtgraph/minimizer.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

using KeyTypes = ::testing::Types<Key64, Key128>;

const std::vector<size_t> payload_sizes { 0, 1, 2, 3 };

//------------------------------------------------------------------------------

// FIXME: add a test with different payload sizes

template<class KeyType>
class KmerIndexManipulation : public ::testing::Test
{
};

TYPED_TEST_CASE(KmerIndexManipulation, KeyTypes);

TYPED_TEST(KmerIndexManipulation, EmptyIndex)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(payload_size);
    std::string msg = "New index with payload size " + std::to_string(payload_size);
    EXPECT_EQ(index.size(), 0) << msg << " contains keys";
    EXPECT_TRUE(index.empty()) << msg << " is not empty";
    EXPECT_EQ(index.number_of_values(), 0) << msg << " has values";
    EXPECT_EQ(index.unique_keys(), 0) << msg << " has unique keys";
    EXPECT_EQ(index.payload_size(), payload_size) << msg << " has wrong payload size";

    index_type copy(index);
    EXPECT_EQ(index, copy) << "A copy of a new index with payload size " << payload_size << " is not identical to the original";
  }
}

TYPED_TEST(KmerIndexManipulation, Contents)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(payload_size);
    index_type copy(index);

    // Different contents.
    {
      pos_t pos = make_pos_t(1, false, 3);
      auto value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key_type(1), value);
      EXPECT_NE(index, copy) << "Empty index is identical to nonempty index with payload size " << payload_size;
    }

    // Same key, different value.
    {
      pos_t pos = make_pos_t(2, false, 3);
      auto value = create_value(pos, payload_size, hash(pos));
      insert_value(copy, key_type(1), value);
      EXPECT_NE(index, copy) << "Indexes with different values are identical with payload size " << payload_size;
    }

    // Same contents.
    {
      copy = index;
      EXPECT_EQ(index, copy) << "A copy of a nonempty index is not identical to the original with payload size " << payload_size;
    }
  }
}

TYPED_TEST(KmerIndexManipulation, Swap)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type first(payload_size), second(payload_size);
    pos_t first_pos = make_pos_t(1, false, 3);
    auto first_value = create_value(first_pos, payload_size, hash(first_pos));
    pos_t second_pos = make_pos_t(2, false, 3);
    auto second_value = create_value(second_pos, payload_size, hash(second_pos));
    insert_value(first, key_type(1), first_value);
    insert_value(second, key_type(2), second_value);

    index_type first_copy(first), second_copy(second);
    first.swap(second);
    EXPECT_NE(first, first_copy) << "Swapping did not change the first index with payload size " << payload_size;
    EXPECT_EQ(first, second_copy) << "The first index was not swapped correctly with payload size " << payload_size;
    EXPECT_EQ(second, first_copy) << "The second index was not swapped correctly with payload size " << payload_size;
    EXPECT_NE(second, second_copy) << "Swapping did not change the second index with payload size " << payload_size;
  }
}

//------------------------------------------------------------------------------

template<class KeyType>
class CorrectKmers : public ::testing::Test
{
public:
  using key_type = KeyType;
  using index_type = KmerIndex<KeyType>;
  using code_type = KmerEncoding::code_type;
  using multi_value_type = KmerEncoding::multi_value_type;
  using result_type = std::map<key_type, std::set<owned_value_type>>;

  size_t total_keys;

  virtual void SetUp()
  {
    this->total_keys = 16;
  }

  void check_kmer_index
  (
    const index_type& index,
    const result_type& correct_values,
    size_t keys, size_t values, size_t unique
  )
  {
    std::string payload_msg = " with payload size " + std::to_string(index.payload_size());

    ASSERT_EQ(index.size(), keys) << "Wrong number of keys" << payload_msg;
    ASSERT_EQ(index.number_of_values(), values) << "Wrong number of values" << payload_msg;
    EXPECT_EQ(index.unique_keys(), unique) << "Wrong number of unique keys" << payload_msg;

    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      size_t count = index.count(iter->first);
      EXPECT_EQ(count, iter->second.size()) << "Wrong number of occurrences for key " << iter->first << payload_msg;
      if(count != iter->second.size()) { continue; }

      multi_value_type values = index.find(iter->first);
      EXPECT_TRUE(same_values(values, iter->second, index.payload_size())) << "Wrong values for key " << iter->first << payload_msg;
    }
  }
};

TYPED_TEST_CASE(CorrectKmers, KeyTypes);

TYPED_TEST(CorrectKmers, UniqueKeys)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(payload_size);
    size_t keys = 0, values = 0, unique = 0;
    typename TestFixture::result_type correct_values;

    for(size_t i = 1; i <= this->total_keys; i++)
    {
      key_type key(i);
      pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
      correct_values[key].insert(value);
      keys++; values++; unique++;
    }
    this->check_kmer_index(index, correct_values, keys, values, unique);
  }
}

TYPED_TEST(CorrectKmers, MissingKeys)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(payload_size);
    for(size_t i = 1; i <= this->total_keys; i++)
    {
      key_type key(i);
      pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
    }

    std::string payload_msg = " with payload size " + std::to_string(index.payload_size());
    for(size_t i = this->total_keys + 1; i <= 2 * this->total_keys; i++)
    {
      key_type key(i);
      EXPECT_EQ(index.count(key), size_t(0)) << "Non-zero occurrences for key " << i << payload_msg;
      EXPECT_EQ(index.find(key).second, size_t(0)) << "Non-empty occurrences for key " << i << payload_msg;
    }
  }
}

TYPED_TEST(CorrectKmers, EmptyKeysValues)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(payload_size);
    std::string payload_msg = " with payload size " + std::to_string(index.payload_size());

    key_type no_key = key_type::no_key();
    insert_value(index, no_key, create_value(make_pos_t(1, false, 0), payload_size, 0));
    EXPECT_EQ(index.count(no_key), size_t(0)) << "Non-zero occurrences for empty key" << payload_msg;
    EXPECT_EQ(index.find(no_key).second, size_t(0)) << "Non-empty value for empty key" << payload_msg;

    key_type key(this->total_keys + 1);
    insert_value(index, key, create_value(Position::no_pos().decode(), payload_size, 0));
    EXPECT_EQ(index.count(key), size_t(0)) << "Non-zero occurrences after inserting empty value" << payload_msg;
    EXPECT_EQ(index.find(key).second, size_t(0)) << "Non-empty value after inserting empty value" << payload_msg;
  }
}

TYPED_TEST(CorrectKmers, MultipleOccurrences)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(payload_size);
    size_t keys = 0, values = 0, unique = 0;
    typename TestFixture::result_type correct_values;

    for(size_t i = 1; i <= this->total_keys; i++)
    {
      key_type key(i);
      pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
      correct_values[key].insert(value);
      keys++; values++; unique++;
    }
    for(size_t i = 1; i <= this->total_keys; i += 2)
    {
      key_type key(i);
      pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
      correct_values[key].insert(value);
      values++; unique--;
    }
    for(size_t i = 1; i <= this->total_keys; i += 4)
    {
      key_type key(i);
      pos_t pos = make_pos_t(i + 2, i & 1, (i + 2) & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
      correct_values[key].insert(value);
      values++;
    }
    this->check_kmer_index(index, correct_values, keys, values, unique);
  }
}

TYPED_TEST(CorrectKmers, DuplicateValues)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(payload_size);
    size_t keys = 0, values = 0, unique = 0;
    typename TestFixture::result_type correct_values;

    for(size_t i = 1; i <= this->total_keys; i++)
    {
      key_type key(i);
      pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
      correct_values[key].insert(value);
      keys++; values++; unique++;
    }
    for(size_t i = 1; i <= this->total_keys; i += 2)
    {
      key_type key(i);
      pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
      correct_values[key].insert(value);
      values++; unique--;
    }
    for(size_t i = 1; i <= this->total_keys; i += 4)
    {
      // Also check that inserting duplicates does not change the payload.
      key_type key(i);
      pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos) + 1);
      insert_value(index, key, value);
    }
    this->check_kmer_index(index, correct_values, keys, values, unique);
  }
}

TYPED_TEST(CorrectKmers, Rehashing)
{
  using key_type = TypeParam;
  using index_type = KmerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(payload_size);
    std::string payload_msg = " with payload size " + std::to_string(index.payload_size());
    size_t keys = 0, values = 0, unique = 0;
    typename TestFixture::result_type correct_values;
    size_t threshold = index.capacity();

    for(size_t i = 1; i <= threshold; i++)
    {
      key_type key(i);
      pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
      correct_values[key].insert(value);
      keys++; values++; unique++;
    }
    ASSERT_EQ(index.capacity(), threshold) << "Index capacity changed at threshold" << payload_msg;

    {
      size_t i = threshold + 1;
      key_type key(i);
      pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(index, key, value);
      correct_values[key].insert(value);
      keys++; values++; unique++;
    }
    EXPECT_GT(index.capacity(), threshold) << "Index capacity did not increase after threshold" << payload_msg;

    this->check_kmer_index(index, correct_values, keys, values, unique);
  }
}

//------------------------------------------------------------------------------

template<class KeyType>
class ObjectManipulation : public ::testing::Test
{
};

TYPED_TEST_CASE(ObjectManipulation, KeyTypes);

TYPED_TEST(ObjectManipulation, EmptyIndex)
{
  using key_type = TypeParam;
  using index_type = MinimizerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type default_index(payload_size);
    index_type default_copy(default_index);
    index_type alt_index(15, 6, payload_size);
    index_type alt_copy(alt_index);

    std::string payload_msg = " with payload size " + std::to_string(payload_size);
    EXPECT_EQ(default_index, default_copy) << "A copy of the default index is not identical to the original" << payload_msg;
    EXPECT_EQ(alt_index, alt_copy) << "A copy of a parameterized index is not identical to the original" << payload_msg;
    EXPECT_NE(default_index, alt_index) << "Default and parameterized indexes are identical" << payload_msg;
  }
}

TYPED_TEST(ObjectManipulation, SyncmerIndex)
{
  using key_type = TypeParam;
  using index_type = MinimizerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type default_minimizer(payload_size);
    index_type default_syncmer(payload_size, true);
    index_type short_minimizer(15, 6, payload_size);
    index_type short_syncmer(15, 6, payload_size, true);

    std::string payload_msg = " with payload size " + std::to_string(payload_size);
    EXPECT_FALSE(default_minimizer.uses_syncmers()) << "Default minimizer index uses syncmers" << payload_msg;
    EXPECT_TRUE(default_syncmer.uses_syncmers()) << "Default syncmer index uses minimizers" << payload_msg;
    EXPECT_FALSE(short_minimizer.uses_syncmers()) << "Parameterized minimizer index uses syncmers" << payload_msg;
    EXPECT_TRUE(short_syncmer.uses_syncmers()) << "Parameterized syncmer index uses minimizers" << payload_msg;
  }
}

TYPED_TEST(ObjectManipulation, Contents)
{
  using key_type = TypeParam;
  using index_type = MinimizerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type default_index(payload_size);
    index_type default_copy(default_index);
    std::string payload_msg = " with payload size " + std::to_string(payload_size);

    // Different contents.
    {
      auto minimizer = get_minimizer(key_type(1));
      pos_t pos = make_pos_t(1, false, 3);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(default_index, minimizer, value);
      EXPECT_NE(default_index, default_copy) << "Empty index is identical to nonempty index" << payload_msg;
    }

    // Same key, different value.
    {
      auto minimizer = get_minimizer(key_type(1));
      pos_t pos = make_pos_t(2, false, 3);
      owned_value_type value = create_value(pos, payload_size, hash(pos));
      insert_value(default_copy, minimizer, value);
      EXPECT_NE(default_index, default_copy) << "Indexes with different values are identical" << payload_msg;
    }

    // Same contents.
    default_copy = default_index;
    EXPECT_EQ(default_index, default_copy) << "A copy of a nonempty index is not identical to the original" << payload_msg;
  }
}

TYPED_TEST(ObjectManipulation, Swap)
{
  using key_type = TypeParam;
  using index_type = MinimizerIndex<key_type>;

  for(size_t payload_size : payload_sizes)
  {
    index_type first(payload_size), second(payload_size);
    auto first_minimizer = get_minimizer(key_type(1));
    pos_t first_pos = make_pos_t(1, false, 3);
    owned_value_type first_value = create_value(first_pos, payload_size, hash(first_pos));
    insert_value(first, first_minimizer, first_value);
    auto second_minimizer = get_minimizer(key_type(2));
    pos_t second_pos = make_pos_t(2, false, 3);
    owned_value_type second_value = create_value(second_pos, payload_size, hash(second_pos));
    insert_value(second, second_minimizer, second_value);

    index_type first_copy(first), second_copy(second);
    first.swap(second);

    std::string payload_msg = " with payload size " + std::to_string(payload_size);
    EXPECT_NE(first, first_copy) << "Swapping did not change the first index" << payload_msg;
    EXPECT_EQ(first, second_copy) << "The first index was not swapped correctly" << payload_msg;
    EXPECT_EQ(second, first_copy) << "The second index was not swapped correctly" << payload_msg;
    EXPECT_NE(second, second_copy) << "Swapping did not change the second index" << payload_msg;
  }
}

//------------------------------------------------------------------------------

template<class KeyType>
class Serialization : public ::testing::Test
{
};

TYPED_TEST_CASE(Serialization, KeyTypes);

TYPED_TEST(Serialization, Serialize)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(15, 6, payload_size);
    auto first_minimizer = get_minimizer(key_type(1));
    pos_t first_pos = make_pos_t(1, false, 3);
    owned_value_type first_value = create_value(first_pos, payload_size, hash(first_pos));
    insert_value(index, first_minimizer, first_value);
    auto second_minimizer = get_minimizer(key_type(2));
    insert_value(index, second_minimizer, first_value);
    pos_t second_pos = make_pos_t(2, false, 3);
    owned_value_type second_value = create_value(second_pos, payload_size, hash(second_pos));
    insert_value(index, second_minimizer, second_value);

    std::string filename = gbwt::TempFile::getName("minimizer");
    std::ofstream out(filename, std::ios_base::binary);
    index.serialize(out);
    out.close();

    index_type copy(payload_size + 1); // Payload size should be overwritten.
    std::ifstream in(filename, std::ios_base::binary);
    copy.deserialize(in);
    in.close();
    gbwt::TempFile::remove(filename);

    EXPECT_EQ(index, copy) << "Loaded index is not identical to the original with payload size " << payload_size;
  }
}

TYPED_TEST(Serialization, WeightedMinimizers)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;

  for(size_t payload_size : payload_sizes)
  {
    index_type index(15, 6, payload_size);
    std::string payload_msg = " with payload size " + std::to_string(payload_size);
    ASSERT_FALSE(index.uses_weighted_minimizers()) << "Weighted minimizers are in use by default" << payload_msg;
    index.add_frequent_kmers({ key_type::encode("GATTACACATGATTA"), key_type::encode("TATTAGATTACATTA") }, 3);
    ASSERT_TRUE(index.uses_weighted_minimizers()) << "Weighted minimizers could not be enabled" << payload_msg;

    auto first_minimizer = get_minimizer(key_type(1));
    pos_t first_pos = make_pos_t(1, false, 3);
    owned_value_type first_value = create_value(first_pos, payload_size, hash(first_pos));
    insert_value(index, first_minimizer, first_value);
    auto second_minimizer = get_minimizer(key_type(2));
    insert_value(index, second_minimizer, first_value);
    pos_t second_pos = make_pos_t(2, false, 3);
    owned_value_type second_value = create_value(second_pos, payload_size, hash(second_pos));
    insert_value(index, second_minimizer, second_value);

    std::string filename = gbwt::TempFile::getName("minimizer");
    std::ofstream out(filename, std::ios_base::binary);
    index.serialize(out);
    out.close();

    index_type copy(payload_size + 1); // Payload size should be overwritten.
    std::ifstream in(filename, std::ios_base::binary);
    copy.deserialize(in);
    in.close();
    gbwt::TempFile::remove(filename);

    EXPECT_EQ(copy, index) << "Loaded index is not identical to the original" << payload_msg;
  }
}

//------------------------------------------------------------------------------

template<class KeyType>
class KeyEncodeDecode : public ::testing::Test
{
public:
  std::string get_string(size_t k) const
  {
    std::string result;
    for(size_t i = 0; i < k; i += 7) { result += "GATTACA"; }
    return result.substr(0, k);
  }
};

TYPED_TEST_CASE(KeyEncodeDecode, KeyTypes);

TYPED_TEST(KeyEncodeDecode, SimpleSequence)
{
  for(size_t k = 1; k <= TypeParam::KMER_MAX_LENGTH; k++)
  {
    for(auto base : {'A', 'C', 'G', 'T'})
    {
      std::string kmer(k, base);
      TypeParam encoded = TypeParam::encode(kmer);
      std::string decoded = encoded.decode(k);
      EXPECT_EQ(decoded, kmer) << "Decoded kmer is not identical to original";
    }
  }
}

TYPED_TEST(KeyEncodeDecode, ComplexSequence)
{
  for(size_t k = 1; k <= TypeParam::KMER_MAX_LENGTH; k++)
  {
    std::string kmer = this->get_string(k);
    TypeParam encoded = TypeParam::encode(kmer);
    std::string decoded = encoded.decode(k);
    EXPECT_EQ(decoded, kmer) << "Decoded kmer is not identical to original";
  }
}

TYPED_TEST(KeyEncodeDecode, RandomAccess)
{
  for(size_t k = 1; k <= TypeParam::KMER_MAX_LENGTH; k++)
  {
    std::string kmer = this->get_string(k);
    TypeParam encoded = TypeParam::encode(kmer);
    for(size_t i = 0; i < k; i++)
    {
      ASSERT_EQ(encoded.access(k, i), kmer[i]) << "Invalid access(" << k << ", " << i << ")";
    }
  }
}

TYPED_TEST(KeyEncodeDecode, ReverseComplement)
{
  for(size_t k = 1; k <= TypeParam::KMER_MAX_LENGTH; k++)
  {
    std::string kmer = this->get_string(k);
    TypeParam encoded = TypeParam::encode(kmer);
    TypeParam rc = encoded.reverse_complement(k);
    std::string decoded = rc.decode(k);
    std::string truth = reverse_complement(kmer);
    EXPECT_EQ(decoded, truth) << "Incorrect reverse complement";
  }
}

//------------------------------------------------------------------------------

/*
  Order of 3-mers using Key64:
  AAT < TGT < TTG < TAT < ATA < TCG < ATT < ACA < GAA < ACT < TAC < CGA < CAA < GTA < TTC < AGT

  Merged with reverse complements:
  AAT < TGT < TTG < TAT < TCG < GAA < ACT < TAC
  ATT   ACA   CAA   ATA   CGA   TTC   AGT   GTA

  Order of 3-mers using Key128:
  ATT < TCG < ATA < TAT < TTG < TGT < AAT < AGT < TTC < GTA < CAA < TAC < CGA < ACT < ACA < GAA

  Merged with reverse complements:
  ATT < TCG < ATA < TTG < TGT < AGT < TTC < GTA
  AAT   CGA   TAT   CAA   ACA   ACT   GAA   TAC
*/

template<class KeyType>
class MinimizerExtraction : public ::testing::Test
{
public:
  std::string str, rev, repetitive;

  virtual void SetUp()
  {
    this->str = "CGAATACAATACT";
    this->rev = reverse_complement(this->str);
    this->repetitive = "TATATA";
  }
};

TYPED_TEST_CASE(MinimizerExtraction, KeyTypes);

TYPED_TEST(MinimizerExtraction, KeyEncoding)
{
  typedef TypeParam key_type;

  size_t k = key_type::KMER_MAX_LENGTH - 1;
  std::string bases = "ACGT";

  key_type forward_key, reverse_key;
  size_t valid_chars = 0;
  std::string correct;
  for(size_t i = 0; i < 2 * k + 1; i++)
  {
    char c = bases[i % bases.length()];
    forward_key.forward(k, c, valid_chars);
    reverse_key.reverse(k, c);
    correct += c;
  }
  correct = correct.substr(correct.length() - k);
  std::string reverse = reverse_complement(correct);

  EXPECT_EQ(forward_key.decode(k), correct) << "Invalid forward key " << forward_key;
  EXPECT_EQ(reverse_key.decode(k), reverse) << "Invalid reverse key " << reverse_key;
}

TYPED_TEST(MinimizerExtraction, AllMinimizers)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 2, 0);
  std::vector<minimizer_type> correct;
  if(key_type::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<key_type>("TCG", 2, true),
      get_minimizer<key_type>("ATA", 3, false),
      get_minimizer<key_type>("ATT", 4, true),
      get_minimizer<key_type>("TGT", 7, true),
      get_minimizer<key_type>("TTG", 8, true),
      get_minimizer<key_type>("ATA", 8, false),
      get_minimizer<key_type>("ATT", 9, true),
      get_minimizer<key_type>("AGT", 12, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<key_type>("TCG", 2, true),
      get_minimizer<key_type>("AAT", 2, false),
      get_minimizer<key_type>("TAT", 5, true),
      get_minimizer<key_type>("TGT", 7, true),
      get_minimizer<key_type>("AAT", 7, false),
      get_minimizer<key_type>("TAT", 10, true),
      get_minimizer<key_type>("ACT", 10, false)
    };
  }

  std::vector<minimizer_type> result = index.minimizers(this->str.begin(), this->str.end());
  ASSERT_EQ(result, correct) << "Did not find the correct minimizers";
  std::vector<minimizer_type> fallback = index.syncmers(this->str.begin(), this->str.end());
  EXPECT_EQ(fallback, correct) << "Did not find the correct minimizers using syncmers()";
}

TYPED_TEST(MinimizerExtraction, WeightedMinimizers)
{
  // This is the same test as above, but we want to avoid kmer ATA. With Key64,
  // we get the lowest-priority kmer TAC as the minimizer in both windows. With
  // Key128, we get its reverse complement GTA instead.
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 2, 0);
  index.add_frequent_kmers({ key_type::encode("ATA") }, 3);
  std::vector<minimizer_type> correct;
  if(key_type::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<key_type>("TCG", 2, true),
      get_minimizer<key_type>("ATT", 4, true),
      get_minimizer<key_type>("GTA", 6, true),
      get_minimizer<key_type>("TGT", 7, true),
      get_minimizer<key_type>("TTG", 8, true),
      get_minimizer<key_type>("ATT", 9, true),
      get_minimizer<key_type>("GTA", 11, true),
      get_minimizer<key_type>("AGT", 12, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<key_type>("TCG", 2, true),
      get_minimizer<key_type>("AAT", 2, false),
      get_minimizer<key_type>("TAC", 4, false),
      get_minimizer<key_type>("TGT", 7, true),
      get_minimizer<key_type>("AAT", 7, false),
      get_minimizer<key_type>("TAC", 9, false),
      get_minimizer<key_type>("ACT", 10, false)
    };
  }

  std::vector<minimizer_type> result = index.minimizers(this->str.begin(), this->str.end());
  ASSERT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, ClosedSyncmers)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(5, 3, 0, true);
  std::vector<minimizer_type> correct;
  if(key_type::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<key_type>("ATACA", 3, false),
      get_minimizer<key_type>("ATTCG", 4, true),
      get_minimizer<key_type>("GTATT", 6, true),
      get_minimizer<key_type>("TTGTA", 8, true),
      get_minimizer<key_type>("ATACT", 8, false),
      get_minimizer<key_type>("ATTGT", 9, true),
      get_minimizer<key_type>("GTATT", 11, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<key_type>("AATAC", 2, false),
      get_minimizer<key_type>("ATACA", 3, false),
      get_minimizer<key_type>("ATTCG", 4, true),
      get_minimizer<key_type>("ACAAT", 5, false),
      get_minimizer<key_type>("AATAC", 7, false),
      get_minimizer<key_type>("AGTAT", 12, true)
    };
  }

  std::vector<minimizer_type> result = index.syncmers(this->str.begin(), this->str.end());
  ASSERT_EQ(result, correct) << "Did not find the correct syncmers";
  std::vector<minimizer_type> fallback = index.minimizers(this->str.begin(), this->str.end());
  EXPECT_EQ(fallback, correct) << "Did not find the correct syncmers using minimizers()";
}

TYPED_TEST(MinimizerExtraction, AllMinimizersWithRegions)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 2, 0);
  std::vector<std::tuple<minimizer_type, size_t, size_t>> correct;
  if(key_type::KEY_BITS == 128)
  {
    // 0000000000111
    // 0123456789012
    // CGAATACAATACT
    // CGA-
    // *--*
    //    ATA+
    //    *--*
    //   AAT-
    //  *---*
    //      ACA-
    //     *--*
    //       CAA-
    //      *--*
    //         ATA+
    //         *--*
    //        AAT-
    //       *---*
    //           ACT-
    //          *--*
    // Each window length will be at least k + w - 1 = 4 bases.
   
    correct =
    {
      std::make_tuple(get_minimizer<key_type>("TCG", 2, true), 0, 4),
      std::make_tuple(get_minimizer<key_type>("ATA", 3, false), 3, 4),
      std::make_tuple(get_minimizer<key_type>("ATT", 4, true), 1, 5),
      std::make_tuple(get_minimizer<key_type>("TGT", 7, true), 4, 4),
      std::make_tuple(get_minimizer<key_type>("TTG", 8, true), 5, 4),
      std::make_tuple(get_minimizer<key_type>("ATA", 8, false), 8, 4),
      std::make_tuple(get_minimizer<key_type>("ATT", 9, true), 6, 5),
      std::make_tuple(get_minimizer<key_type>("AGT", 12, true), 9, 4)
    };
  }
  else
  {
    // 0000000000111
    // 0123456789012
    // CGAATACAATACT
    // CGA-
    // *--*
    //   AAT+ 
    //  *---*
    //    ATA-
    //    *--*
    //      ACA-
    //     *---* 
    //        AAT+
    //       *---*
    //         ATA-
    //         *--*
    //           ACT+
    //          *--*
    
    // Note that no minimizers actually get replaced in this example; they all take their second possible window.

    correct =
    {
      std::make_tuple(get_minimizer<key_type>("TCG", 2, true), 0, 4),
      std::make_tuple(get_minimizer<key_type>("AAT", 2, false), 1, 5),
      std::make_tuple(get_minimizer<key_type>("TAT", 5, true), 3, 4),
      std::make_tuple(get_minimizer<key_type>("TGT", 7, true), 4, 5),
      std::make_tuple(get_minimizer<key_type>("AAT", 7, false), 6, 5),
      std::make_tuple(get_minimizer<key_type>("TAT", 10, true), 8, 4),
      std::make_tuple(get_minimizer<key_type>("ACT", 10, false), 9, 4)
    };
  }
  std::vector<std::tuple<minimizer_type, size_t, size_t>> result = index.minimizer_regions(this->str.begin(), this->str.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TEST(MinimizerExtraction, HardMinimizersWithRegion)
{
  using TypeParam = Key128;
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  // Here's a case I caught not working correctly.
  index_type index(29, 11, 0);
  std::string seq("ACCAGTTTTTTACACAAGCTGCTCTTTCCCTCAATTGTTCATTTGTCTCCTTGTCCAGGT");
  
  // No replacements happen
  
  // 000000000011111111112222222222333333333344444444445555555555
  // 012345678901234567890123456789012345678901234567890123456789
  // ACCAGTTTTTTACACAAGCTGCTCTTTCCCTCAATTGTTCATTTGTCTCCTTGTCCAGGT
  //       TTTTTACACAAGCTGCTCTTTCCCTCAAT-                        (#1)
  // *-------------------------------------------*                
  //             CACAAGCTGCTCTTTCCCTCAATTGTTCA+                  (#2)
  //        *------------------------------------------*         
  //                  GCTGCTCTTTCCCTCAATTGTTCATTTGT-             (#3)
  //              *-----------------------------------------*    
  //                         TTTCCCTCAATTGTTCATTTGTCTCCTTG+      (#4)
  //                   *----------------------------------------*
  // These then get sorted by minimizer offset, at the end for reverse strand minimizers.
  
  std::vector<std::tuple<minimizer_type, size_t, size_t>> correct;
  correct =
  {
    std::make_tuple(get_minimizer<key_type>("CACAAGCTGCTCTTTCCCTCAATTGTTCA", 12, false), 7, 44), // #2
    std::make_tuple(get_minimizer<key_type>("TTTCCCTCAATTGTTCATTTGTCTCCTTG", 24, false), 18, 42), // #4
    std::make_tuple(get_minimizer<key_type>("ATTGAGGGAAAGAGCAGCTTGTGTAAAAA", 34, true), 0, 45), // #1
    std::make_tuple(get_minimizer<key_type>("ACAAATGAACAATTGAGGGAAAGAGCAGC", 45, true), 13, 43) // #3
  };
  std::vector<std::tuple<minimizer_type, size_t, size_t>> result = index.minimizer_regions(seq.begin(), seq.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, WindowLength)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 3, 0);
  std::vector<minimizer_type> correct;
  if(key_type::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<key_type>("ATA", 3, false),
      get_minimizer<key_type>("ATT", 4, true),
      get_minimizer<key_type>("TTG", 8, true),
      get_minimizer<key_type>("ATA", 8, false),
      get_minimizer<key_type>("ATT", 9, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<key_type>("AAT", 2, false),
      get_minimizer<key_type>("TGT", 7, true),
      get_minimizer<key_type>("AAT", 7, false),
      get_minimizer<key_type>("TAT", 10, true)
    };
  }
  std::vector<minimizer_type> result = index.minimizers(this->str.begin(), this->str.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

// We have multiple occurrences of the same minimizer in the initial window and after it.
// When we find the first occurrence after the window, one of the original occurrences is
// still in the buffer.
TYPED_TEST(MinimizerExtraction, AllOccurrences)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 3, 0);
  std::vector<minimizer_type> correct;
  if(key_type::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<key_type>("ATA", 1, false),
      get_minimizer<key_type>("ATA", 2, true),
      get_minimizer<key_type>("ATA", 3, false),
      get_minimizer<key_type>("ATA", 4, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<key_type>("TAT", 0, false),
      get_minimizer<key_type>("TAT", 2, false),
      get_minimizer<key_type>("TAT", 3, true),
      get_minimizer<key_type>("TAT", 5, true)
    };
  }
  std::vector<minimizer_type> result = index.minimizers(this->repetitive.begin(), this->repetitive.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, WeirdSyncmers)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  // The string is the reverse complement of itself. The middle smers AAT and ATT
  // are the smallest ones using both key types.
  std::string weird = "CAATTG";
  index_type index(5, 3, 0, true);
  std::vector<minimizer_type> correct;
  if(key_type::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<key_type>("AATTG", 1, false),
      get_minimizer<key_type>("AATTG", 4, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<key_type>("CAATT", 0, false),
      get_minimizer<key_type>("CAATT", 5, true)
    };
  }
  std::vector<minimizer_type> result = index.syncmers(weird);
  EXPECT_EQ(result, correct) << "Did not find the correct syncmers";
}

TYPED_TEST(MinimizerExtraction, InvalidMinimizerCharacters)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  std::string weird = "CGAATAxAATACT";
  index_type index(3, 2, 0);
  std::vector<minimizer_type> correct;
  if(key_type::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<key_type>("TCG", 2, true),
      get_minimizer<key_type>("ATA", 3, false),
      get_minimizer<key_type>("ATT", 4, true),
      get_minimizer<key_type>("ATA", 8, false),
      get_minimizer<key_type>("ATT", 9, true),
      get_minimizer<key_type>("AGT", 12, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<key_type>("TCG", 2, true),
      get_minimizer<key_type>("AAT", 2, false),
      get_minimizer<key_type>("TAT", 5, true),
      get_minimizer<key_type>("AAT", 7, false),
      get_minimizer<key_type>("TAT", 10, true),
      get_minimizer<key_type>("ACT", 10, false)
    };
  }
  std::vector<minimizer_type> result = index.minimizers(weird.begin(), weird.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, InvalidSyncmerCharacters)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  std::string weird = "CGAATAxAATACT";
  index_type index(5, 3, 0, true);
  std::vector<minimizer_type> correct;
  if(key_type::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<key_type>("ATTCG", 4, true),
      get_minimizer<key_type>("ATACT", 8, false),
      get_minimizer<key_type>("GTATT", 11, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<key_type>("ATTCG", 4, true),
      get_minimizer<key_type>("AATAC", 7, false),
      get_minimizer<key_type>("AGTAT", 12, true)
    };
  }
  std::vector<minimizer_type> result = index.minimizers(weird.begin(), weird.end());
  EXPECT_EQ(result, correct) << "Did not find the correct syncmers";
}

TYPED_TEST(MinimizerExtraction, BothOrientations)
{
  typedef TypeParam key_type;
  typedef MinimizerIndex<key_type> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 2, 0);
  std::vector<minimizer_type> forward_minimizers = index.minimizers(this->str.begin(), this->str.end());
  std::vector<minimizer_type> reverse_minimizers = index.minimizers(this->rev.begin(), this->rev.end());
  ASSERT_EQ(forward_minimizers.size(), reverse_minimizers.size()) << "Different number of minimizers in forward and reverse orientations";
  for(size_t i = 0; i < forward_minimizers.size(); i++)
  {
    minimizer_type& f = forward_minimizers[i];
    minimizer_type& r = reverse_minimizers[forward_minimizers.size() - 1 - i];
    EXPECT_EQ(f.key, r.key) << "Wrong key for minimizer " << i;
    EXPECT_EQ(f.offset, this->str.length() - 1 - r.offset) << "Wrong offset for minimizer " << i;
    EXPECT_NE(f.is_reverse, r.is_reverse) << "Wrong orientation for minimizer " << i;
  }
}

//------------------------------------------------------------------------------

class HitsInSubgraphTest : public ::testing::Test
{
public:
  using value_type = KmerEncoding::value_type;
  using multi_value_type = KmerEncoding::multi_value_type;
  using owned_multi_value_type = std::vector<KmerEncoding::code_type>;
  using result_type = std::vector<owned_value_type>;  

  void check_results
  (
    const MinimizerIndex<Key64>& index,
    const std::unordered_set<nid_t>& subgraph,
    const owned_multi_value_type& hits,
    const result_type& expected_result,
    const std::string& test_case
  ) const
  {
    std::vector<nid_t> sorted_subgraph(subgraph.begin(), subgraph.end());
    std::sort(sorted_subgraph.begin(), sorted_subgraph.end());

    result_type result;
    auto report_hit = [&](value_type value) {
      result.push_back(create_value(value, index.payload_size()));
    };
    
    multi_value_type h(hits.data(), hits.size() / index.value_size());
    hits_in_subgraph(index, h, subgraph, report_hit);
    ASSERT_EQ(result.size(), expected_result.size()) << test_case << ": Incorrect number of results with the naive algorithm";
    ASSERT_EQ(result, expected_result) << test_case << ": Incorrect results with the naive algorithm";

    result.clear();
    hits_in_subgraph(index, h, sorted_subgraph, report_hit);
    ASSERT_EQ(result.size(), expected_result.size()) << test_case << ": Incorrect number of results with exponential search";
    ASSERT_EQ(result, expected_result) << test_case << ": Incorrect results with exponential search";
  }

  void encode_hit(owned_multi_value_type& hits, const owned_value_type& value, size_t payload_size) const
  {
    size_t offset = hits.size();
    hits.resize(offset + KmerIndex<Key64>::POS_SIZE + payload_size);
    value.first.write(hits.data() + offset);
    offset += KmerIndex<Key64>::POS_SIZE;
    for(size_t j = 0; j < payload_size; j++) { hits[offset + j] = value.second[j]; }
  }

  std::tuple<std::unordered_set<nid_t>, owned_multi_value_type, result_type> create_test_case
  (
    const MinimizerIndex<Key64>& index,
    size_t universe_size, size_t begin, size_t end,
    double interval_prob, double outlier_prob, double hit_prob, size_t random_seed
  ) const
  {
    std::unordered_set<nid_t> subgraph;
    owned_multi_value_type hits;
    result_type expected_result;

    std::mt19937_64 rng(random_seed);
    for(size_t i = 0; i < universe_size; i++)
    {
      // Is node i in the subgraph?
      double random_value = rng() / static_cast<double>(std::numeric_limits<std::mt19937_64::result_type>::max());
      bool in_subgraph = false;
      if(i >= begin && i < end)
      {
        if(random_value <= interval_prob)
        {
          subgraph.insert(i); in_subgraph = true;
        }
      }
      else if(random_value <= outlier_prob)
      {
        subgraph.insert(i); in_subgraph = true;
      }

      // Are there hits for node i?
      random_value = rng() / static_cast<double>(std::numeric_limits<std::mt19937_64::result_type>::max());
      size_t hit_count = 0;
      while(random_value <= hit_prob)
      {
        pos_t pos = make_pos_t(i, hit_count & 1, hit_count & Position::OFF_MASK);
        owned_value_type value = create_value(pos, index.payload_size(), hash(pos));
        this->encode_hit(hits, value, index.payload_size());
        if(in_subgraph) { expected_result.push_back(value); }
        hit_count++;
        random_value = rng() / static_cast<double>(std::numeric_limits<std::mt19937_64::result_type>::max());
      }
    }

    return std::make_tuple(subgraph, hits, expected_result);
  }
};

TEST_F(HitsInSubgraphTest, EmptySets)
{
  for(size_t payload_size : payload_sizes)
  {
    std::string payload_msg = " with payload size " + std::to_string(payload_size);
    MinimizerIndex<Key64> index(payload_size);
    std::unordered_set<nid_t> subgraph;
    owned_multi_value_type hits;
    result_type expected_result;

    std::string test_case = "Empty subgraph and hits" + payload_msg;
    this->check_results(index, subgraph, hits, expected_result, test_case);

    subgraph.insert(42);
    test_case = "Nonempty subgraph, empty hits" + payload_msg;
    this->check_results(index, subgraph, hits, expected_result, test_case);

    subgraph.clear();
    pos_t pos = make_pos_t(42, false, 0);
    owned_value_type value = create_value(pos, index.payload_size(), hash(pos));
    this->encode_hit(hits, value, index.payload_size());
    test_case = "Empty subgraph, nonempty hits" + payload_msg;
    this->check_results(index, subgraph, hits, expected_result, test_case);
  }
}

TEST_F(HitsInSubgraphTest, SmallSets)
{
  constexpr size_t UNIVERSE_SIZE = 1024;
  constexpr size_t INTERVALS = 10;
  constexpr size_t INTERVAL_START = 90;
  constexpr size_t INTERVAL_LENGTH = 30;
  constexpr double INTERVAL_PROB = 0.9;
  constexpr double OUTLIER_PROB = 0.01;
  constexpr double HIT_PROB = 0.1;

  for(size_t payload_size : payload_sizes)
  {
    MinimizerIndex<Key64> index(payload_size);
    for(size_t i = 1; i <= INTERVALS; i++)
    {
      std::unordered_set<nid_t> subgraph;
      owned_multi_value_type hits;
      result_type expected_result;
      size_t start = i * INTERVAL_START;
      size_t random_seed = i * 0xDEADBEEF;
      std::tie(subgraph, hits, expected_result) = this->create_test_case
      (
        index,
        UNIVERSE_SIZE, start, start + INTERVAL_LENGTH,
        INTERVAL_PROB, OUTLIER_PROB, HIT_PROB, random_seed
      );
      std::string test_case = "Payload size " + std::to_string(payload_size) + ", set " + std::to_string(i);
      this->check_results(index, subgraph, hits, expected_result, test_case);
    }
  }
}

TEST_F(HitsInSubgraphTest, LargeSets)
{
  constexpr size_t UNIVERSE_SIZE = 1048576;
  constexpr size_t INTERVALS = 3;
  constexpr size_t INTERVAL_START = 256 * 1024;
  constexpr size_t INTERVAL_LENGTH = 120;
  constexpr double INTERVAL_PROB = 0.9;
  constexpr double OUTLIER_PROB = 0.0001;
  constexpr double HIT_PROB = 0.05;

  for(size_t payload_size : payload_sizes)
  {
    MinimizerIndex<Key64> index(payload_size);
    for(size_t i = 1; i <= INTERVALS; i++)
    {
      std::unordered_set<nid_t> subgraph;
      owned_multi_value_type hits;
      result_type expected_result;
      size_t start = i * INTERVAL_START;
      size_t random_seed = i * 0xDEADBEEF;
      std::tie(subgraph, hits, expected_result) = this->create_test_case
      (
        index,
        UNIVERSE_SIZE, start, start + INTERVAL_LENGTH,
        INTERVAL_PROB, OUTLIER_PROB, HIT_PROB, random_seed
      );
      std::string test_case = "Payload size " + std::to_string(payload_size) + ", set " + std::to_string(i);
      this->check_results(index, subgraph, hits, expected_result, test_case);
    }
  }
}

//------------------------------------------------------------------------------

} // namespace
