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

std::pair<Position, std::vector<std::uint64_t>>
create_value(pos_t pos, size_t payload_size, std::uint64_t payload)
{
  return std::make_pair(Position(pos), std::vector<std::uint64_t>(payload_size, payload));
}

template<class KeyType>
void insert_value(KmerIndex<KeyType>& index, KeyType key, const std::pair<Position, std::vector<std::uint64_t>>& value)
{
  index.insert(key, std::make_pair(value.first, value.second.data()));
}

const std::vector<size_t> payload_sizes { 0, 1, 2, 3 };

//------------------------------------------------------------------------------

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
      auto value = create_value(make_pos_t(1, false, 3), payload_size, hash(1, false, 3));
      insert_value(index, key_type(1), value);
      EXPECT_NE(index, copy) << "Empty index is identical to nonempty index with payload size " << payload_size;
    }

    // Same key, different value.
    {
      auto value = create_value(make_pos_t(2, false, 3), payload_size, hash(2, false, 3));
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
    auto first_value = create_value(make_pos_t(1, false, 3), payload_size, hash(1, false, 3));
    auto second_value = create_value(make_pos_t(2, false, 3), payload_size, hash(2, false, 3));
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

// FIXME: from here

template<class KeyType>
class CorrectKmers : public ::testing::Test
{
public:
  using key_type = KeyType;
  using index_type = KmerIndex<KeyType>;
  using value_type = typename index_type::value_type; // FIXME: is this what we want?
  using result_type = std::map<key_type, std::set<value_type>>;

  size_t total_keys;

  virtual void SetUp()
  {
    this->total_keys = 16;
  }

  // FIXME: payload size in error messages?
  void check_kmer_index_index(const index_type& index,
                             const result_type& correct_values,
                             size_t keys, size_t values, size_t unique)
  {
    ASSERT_EQ(index.size(), keys) << "Wrong number of keys";
    ASSERT_EQ(index.number_of_values(), values) << "Wrong number of values";
    EXPECT_EQ(index.unique_keys(), unique) << "Wrong number of unique keys";

    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      std::vector<value_type> correct(iter->second.begin(), iter->second.end());

      size_t count = index.count(iter->first);
      EXPECT_EQ(count, correct.size()) << "Wrong number of occurrences for key " << iter->first;
      if(count != correct.size()) { continue; }

      std::pair<const value_type*, size_t> occs = index.find(iter->first);
      std::vector<value_type> found(occs.first, occs.first + occs.second);
      EXPECT_EQ(found, correct) << "Wrong positions for key " << iter->first;
    }
  }
};

TYPED_TEST_CASE(CorrectKmers, KeyTypes);

TYPED_TEST(CorrectKmers, UniqueKeys)
{
  typedef TypeParam index_type;
  typedef typename index_type::key_type key_type;
  typedef typename index_type::value_type value_type;

  index_type index;
  size_t keys = 0, values = 0, unique = 0;
  typename TestFixture::result_type correct_values;

  for(size_t i = 1; i <= this->total_keys; i++)
  {
    key_type key(i);
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
    correct_values[key].insert(value);
    keys++; values++; unique++;
  }
  this->check_kmer_index_index(index, correct_values, keys, values, unique);
}

TYPED_TEST(CorrectKmers, MissingKeys)
{
  typedef TypeParam index_type;
  typedef typename index_type::key_type key_type;
  typedef typename index_type::value_type value_type;

  index_type index;
  for(size_t i = 1; i <= this->total_keys; i++)
  {
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    index.insert(key_type(i), create_value<value_type>(pos, payload));
  }
  for(size_t i = this->total_keys + 1; i <= 2 * this->total_keys; i++)
  {
    key_type key(i);
    EXPECT_EQ(index.count(key), size_t(0)) << "Non-zero occurrences for key " << i;
    std::pair<const value_type*, size_t> correct(nullptr, 0);
    EXPECT_EQ(index.find(key), correct) << "Non-empty occurrences for key " << i;
  }
}

TYPED_TEST(CorrectKmers, EmptyKeysValues)
{
  typedef TypeParam index_type;
  typedef typename index_type::key_type key_type;
  typedef typename index_type::value_type value_type;

  index_type index;
  std::pair<const value_type*, size_t> correct(nullptr, 0);

  index.insert(key_type::no_key(), create_value<value_type>(make_pos_t(1, false, 0), Payload::create(0)));
  EXPECT_EQ(index.count(key_type::no_key()), size_t(0)) << "Non-zero occurrences for empty key";
  EXPECT_EQ(index.find(key_type::no_key()), correct) << "Non-empty value for empty key";

  key_type key(this->total_keys + 1);
  index.insert(key, create_value<value_type>(Position::no_value().decode(), Payload::create(0)));
  EXPECT_EQ(index.count(key), size_t(0)) << "Non-zero occurrences after inserting empty value";
  EXPECT_EQ(index.find(key), correct) << "Non-empty value after inserting empty value";
}

TYPED_TEST(CorrectKmers, MultipleOccurrences)
{
  typedef TypeParam index_type;
  typedef typename index_type::key_type key_type;
  typedef typename index_type::value_type value_type;

  index_type index;
  size_t keys = 0, values = 0, unique = 0;
  typename TestFixture::result_type correct_values;

  for(size_t i = 1; i <= this->total_keys; i++)
  {
    key_type key(i);
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
    correct_values[key].insert(value);
    keys++; values++; unique++;
  }
  for(size_t i = 1; i <= this->total_keys; i += 2)
  {
    key_type key(i);
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, (i + 1) & Position::OFF_MASK));
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
    correct_values[key].insert(value);
    values++; unique--;
  }
  for(size_t i = 1; i <= this->total_keys; i += 4)
  {
    key_type key(i);
    pos_t pos = make_pos_t(i + 2, i & 1, (i + 2) & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, (i + 2) & Position::OFF_MASK));
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
    correct_values[key].insert(value);
    values++;
  }
  this->check_kmer_index_index(index, correct_values, keys, values, unique);
}

TYPED_TEST(CorrectKmers, DuplicateValues)
{
  typedef TypeParam index_type;
  typedef typename index_type::key_type key_type;
  typedef typename index_type::value_type value_type;

  index_type index;
  size_t keys = 0, values = 0, unique = 0;
  typename TestFixture::result_type correct_values;

  for(size_t i = 1; i <= this->total_keys; i++)
  {
    key_type key(i);
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
    correct_values[key].insert(value);
    keys++; values++; unique++;
  }
  for(size_t i = 1; i <= this->total_keys; i += 2)
  {
    key_type key(i);
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, (i + 1) & Position::OFF_MASK));
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
    correct_values[key].insert(value);
    values++; unique--;
  }
  for(size_t i = 1; i <= this->total_keys; i += 4)
  {
    // Also check that inserting duplicates does not change the payload.
    key_type key(i);
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, (i + 1) & Position::OFF_MASK) + 1);
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
  }
  this->check_kmer_index_index(index, correct_values, keys, values, unique);
}

TYPED_TEST(CorrectKmers, Rehashing)
{
  typedef TypeParam index_type;
  typedef typename index_type::key_type key_type;
  typedef typename index_type::value_type value_type;

  index_type index;
  size_t keys = 0, values = 0, unique = 0;
  typename TestFixture::result_type correct_values;
  size_t threshold = index.capacity();

  for(size_t i = 1; i <= threshold; i++)
  {
    key_type key(i);
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
    correct_values[key].insert(value);
    keys++; values++; unique++;
  }
  ASSERT_EQ(index.capacity(), threshold) << "Index capacity changed at threshold";

  {
    size_t i = threshold + 1;
    key_type key(i);
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    value_type value = create_value<value_type>(pos, payload);
    index.insert(key, value);
    correct_values[key].insert(value);
    keys++; values++; unique++;
  }
  EXPECT_GT(index.capacity(), threshold) << "Index capacity did not increase after threshold";

  this->check_kmer_index_index(index, correct_values, keys, values, unique);
}

//------------------------------------------------------------------------------

template<class IndexType>
class ObjectManipulation : public ::testing::Test
{
};

TYPED_TEST_CASE(ObjectManipulation, MinimizerIndexes);

TYPED_TEST(ObjectManipulation, EmptyIndex)
{
  TypeParam default_index;
  TypeParam default_copy(default_index);
  TypeParam alt_index(15, 6);
  TypeParam alt_copy(alt_index);
  EXPECT_EQ(default_index, default_copy) << "A copy of the default index is not identical to the original";
  EXPECT_EQ(alt_index, alt_copy) << "A copy of a parameterized index is not identical to the original";
  EXPECT_NE(default_index, alt_index) << "Default and parameterized indexes are identical";
}

TYPED_TEST(ObjectManipulation, SyncmerIndex)
{
  TypeParam default_minimizer;
  TypeParam default_syncmer(true);
  TypeParam short_minimizer(15, 6);
  TypeParam short_syncmer(15, 6, true);
  EXPECT_FALSE(default_minimizer.uses_syncmers()) << "Default minimizer index uses syncmers";
  EXPECT_TRUE(default_syncmer.uses_syncmers()) << "Default syncmer index uses minimizers";
  EXPECT_FALSE(short_minimizer.uses_syncmers()) << "Parameterized minimizer index uses syncmers";
  EXPECT_TRUE(short_syncmer.uses_syncmers()) << "Parameterized syncmer index uses minimizers";
}

TYPED_TEST(ObjectManipulation, Contents)
{
  typedef TypeParam index_type;
  typedef typename index_type::key_type key_type;
  typedef typename index_type::value_type value_type;

  index_type default_index;
  index_type default_copy(default_index);

  // Different contents.
  default_index.insert(get_minimizer<key_type>(1), create_value<value_type>(make_pos_t(1, false, 3), Payload::create(hash(1, false, 3))));
  EXPECT_NE(default_index, default_copy) << "Empty index is identical to nonempty index";

  // Same key, different value.
  default_copy.insert(get_minimizer<key_type>(1), create_value<value_type>(make_pos_t(2, false, 3), Payload::create(hash(2, false, 3))));
  EXPECT_NE(default_index, default_copy) << "Indexes with different values are identical";

  // Same contents.
  default_copy = default_index;
  EXPECT_EQ(default_index, default_copy) << "A copy of a nonempty index is not identical to the original";
}

TYPED_TEST(ObjectManipulation, Swap)
{
  typedef TypeParam index_type;
  typedef typename index_type::key_type key_type;
  typedef typename index_type::value_type value_type;

  index_type first, second;
  first.insert(get_minimizer<key_type>(1), create_value<value_type>(make_pos_t(1, false, 3), Payload::create(hash(1, false, 3))));
  second.insert(get_minimizer<key_type>(2), create_value<value_type>(make_pos_t(2, false, 3), Payload::create(hash(2, false, 3))));

  index_type first_copy(first), second_copy(second);
  first.swap(second);
  EXPECT_NE(first, first_copy) << "Swapping did not change the first index";
  EXPECT_EQ(first, second_copy) << "The first index was not swapped correctly";
  EXPECT_EQ(second, first_copy) << "The second index was not swapped correctly";
  EXPECT_NE(second, second_copy) << "Swapping did not change the second index";
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
  typedef PositionPayload<Payload> value_type;
  typedef MinimizerIndex<key_type, value_type> index_type;

  index_type index(15, 6);
  index.insert(get_minimizer<key_type>(1), create_value<value_type>(make_pos_t(1, false, 3), Payload::create(hash(1, false, 3))));
  index.insert(get_minimizer<key_type>(2), create_value<value_type>(make_pos_t(1, false, 3), Payload::create(hash(1, false, 3))));
  index.insert(get_minimizer<key_type>(2), create_value<value_type>(make_pos_t(2, false, 3), Payload::create(hash(2, false, 3))));

  std::string filename = gbwt::TempFile::getName("minimizer");
  std::ofstream out(filename, std::ios_base::binary);
  index.serialize(out);
  out.close();

  index_type copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.deserialize(in);
  in.close();
  gbwt::TempFile::remove(filename);

  EXPECT_EQ(index, copy) << "Loaded index is not identical to the original";
}

TYPED_TEST(Serialization, WeightedMinimizers)
{
  typedef TypeParam key_type;
  typedef PositionPayload<Payload> value_type;
  typedef MinimizerIndex<key_type, value_type> index_type;

  index_type index(15, 6);
  ASSERT_FALSE(index.uses_weighted_minimizers()) << "Weighted minimizers are in use by default";
  index.add_frequent_kmers({ key_type::encode("GATTACACATGATTA"), key_type::encode("TATTAGATTACATTA") }, 3);
  ASSERT_TRUE(index.uses_weighted_minimizers()) << "Weighted minimizers could not be enabled";

  index.insert(get_minimizer<key_type>(1), create_value<value_type>(make_pos_t(1, false, 3), Payload::create(hash(1, false, 3))));
  index.insert(get_minimizer<key_type>(2), create_value<value_type>(make_pos_t(1, false, 3), Payload::create(hash(1, false, 3))));
  index.insert(get_minimizer<key_type>(2), create_value<value_type>(make_pos_t(2, false, 3), Payload::create(hash(2, false, 3))));

  std::string filename = gbwt::TempFile::getName("minimizer");
  std::ofstream out(filename, std::ios_base::binary);
  index.serialize(out);
  out.close();

  index_type copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.deserialize(in);
  in.close();
  gbwt::TempFile::remove(filename);

  EXPECT_EQ(copy, index) << "Loaded index is not identical to the original";
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 2);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 2);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(5, 3, true);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 2);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  // Here's a case I caught not working correctly.
  index_type index(29, 11);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 3);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 3);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  // The string is the reverse complement of itself. The middle smers AAT and ATT
  // are the smallest ones using both key types.
  std::string weird = "CAATTG";
  index_type index(5, 3, true);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  std::string weird = "CGAATAxAATACT";
  index_type index(3, 2);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  std::string weird = "CGAATAxAATACT";
  index_type index(5, 3, true);
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
  typedef MinimizerIndex<key_type, Position> index_type;
  typedef Kmer<key_type> minimizer_type;

  index_type index(3, 2);
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
  typedef std::vector<std::pair<pos_t, Payload>> result_type;

  

  void check_results(const std::unordered_set<nid_t>& subgraph, const std::vector<PositionPayload<Payload>>& hits,
                     const result_type& expected_result, const std::string& test_case)
  {
    std::vector<nid_t> sorted_subgraph(subgraph.begin(), subgraph.end());
    std::sort(sorted_subgraph.begin(), sorted_subgraph.end());

    result_type result;
    using Report = std::function<void(pos_t, Payload)>;
    Report cb1 = [&](pos_t pos, Payload payload){
      result.emplace_back(pos, payload);
    };
    
    hits_in_subgraph(hits.size(), hits.data(), subgraph, cb1);
    ASSERT_EQ(result, expected_result) << test_case << ": Incorrect results with the naive algorithm";

    result.clear();
    hits_in_subgraph(hits.size(), hits.data(), sorted_subgraph, cb1);
    ASSERT_EQ(result, expected_result) << test_case << ": Incorrect results with exponential search";
  }

  std::tuple<std::unordered_set<nid_t>, std::vector<PositionPayload<Payload>>, result_type>
  create_test_case(size_t universe_size, size_t begin, size_t end,
                   double interval_prob, double outlier_prob, double hit_prob, size_t random_seed) const
  {
    std::unordered_set<nid_t> subgraph;
    std::vector<PositionPayload<Payload>> hits;
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
        Payload payload = Payload::create(hit_count);
        hits.push_back({ Position::encode(pos), payload });
        if(in_subgraph) { expected_result.emplace_back(pos, payload); }
        hit_count++;
        random_value = rng() / (double)std::numeric_limits<size_t>::max();
      }
    }

    return std::make_tuple(subgraph, hits, expected_result);
  }
};

TEST_F(HitsInSubgraphTest, EmptySets)
{
  std::unordered_set<nid_t> subgraph;
  std::vector<PositionPayload<Payload>> hits;
  result_type expected_result;

  this->check_results(subgraph, hits, expected_result, "Empty subgraph and hits");

  subgraph.insert(42);
  this->check_results(subgraph, hits, expected_result, "Empty hits");

  subgraph.clear();
  hits.push_back({ Position::create(42), Payload::create(42) });
  this->check_results(subgraph, hits, expected_result, "Empty subgraph");
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

  for(size_t i = 1; i <= INTERVALS; i++)
  {
    std::unordered_set<nid_t> subgraph;
    std::vector<PositionPayload<Payload>> hits;
    result_type expected_result;
    size_t start = i * INTERVAL_START;
    size_t random_seed = i * 0xDEADBEEF;
    std::tie(subgraph, hits, expected_result) =
      this->create_test_case(UNIVERSE_SIZE, start, start + INTERVAL_LENGTH, INTERVAL_PROB, OUTLIER_PROB, HIT_PROB, random_seed);
    this->check_results(subgraph, hits, expected_result, "Set " + std::to_string(i));
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

  for(size_t i = 1; i <= INTERVALS; i++)
  {
    std::unordered_set<nid_t> subgraph;
    std::vector<PositionPayload<Payload>> hits;
    result_type expected_result;
    size_t start = i * INTERVAL_START;
    size_t random_seed = i * 0xDEADBEEF;
    std::tie(subgraph, hits, expected_result) =
      this->create_test_case(UNIVERSE_SIZE, start, start + INTERVAL_LENGTH, INTERVAL_PROB, OUTLIER_PROB, HIT_PROB, random_seed);
    this->check_results(subgraph, hits, expected_result, "Set " + std::to_string(i));
  }
}

//------------------------------------------------------------------------------

} // namespace
