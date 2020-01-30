#include <gtest/gtest.h>

#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <gbwtgraph/minimizer.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

using KeyTypes = ::testing::Types<Key64, Key128>;

//------------------------------------------------------------------------------

template<class KeyType>
class ObjectManipulation : public ::testing::Test
{
};

TYPED_TEST_CASE(ObjectManipulation, KeyTypes);

TYPED_TEST(ObjectManipulation, EmptyIndex)
{
  MinimizerIndex<TypeParam> default_index;
  MinimizerIndex<TypeParam> default_copy(default_index);
  MinimizerIndex<TypeParam> alt_index(15, 6);
  MinimizerIndex<TypeParam> alt_copy(alt_index);
  EXPECT_EQ(default_index, default_copy) << "A copy of the default index is not identical to the original";
  EXPECT_EQ(alt_index, alt_copy) << "A copy of a parametrized index is not identical to the original";
  EXPECT_NE(default_index, alt_index) << "Default and parametrized indexes are identical";
}

TYPED_TEST(ObjectManipulation, Contents)
{
  MinimizerIndex<TypeParam> default_index;
  MinimizerIndex<TypeParam> default_copy(default_index);

  // Different contents.
  default_index.insert(get_minimizer<TypeParam>(1), make_pos_t(1, false, 3));
  EXPECT_NE(default_index, default_copy) << "Empty index is identical to nonempty index";

  // Same key, different value.
  default_copy.insert(get_minimizer<TypeParam>(1), make_pos_t(2, false, 3));
  EXPECT_NE(default_index, default_copy) << "Indexes with different values are identical";

  // Same contents.
  default_copy = default_index;
  EXPECT_EQ(default_index, default_copy) << "A copy of a nonempty index is not identical to the original";
}

TYPED_TEST(ObjectManipulation, Swap)
{
  MinimizerIndex<TypeParam> first, second;
  first.insert(get_minimizer<TypeParam>(1), make_pos_t(1, false, 3));
  second.insert(get_minimizer<TypeParam>(2), make_pos_t(2, false, 3));

  MinimizerIndex<TypeParam> first_copy(first), second_copy(second);
  first.swap(second);
  EXPECT_NE(first, first_copy) << "Swapping did not change the first index";
  EXPECT_EQ(first, second_copy) << "The first index was not swapped correctly";
  EXPECT_EQ(second, first_copy) << "The second index was not swapped correctly";
  EXPECT_NE(second, second_copy) << "Swapping did not change the second index";
}

TYPED_TEST(ObjectManipulation, Serialization)
{
  MinimizerIndex<TypeParam> index(15, 6);
  index.insert(get_minimizer<TypeParam>(1), make_pos_t(1, false, 3));
  index.insert(get_minimizer<TypeParam>(2), make_pos_t(1, false, 3));
  index.insert(get_minimizer<TypeParam>(2), make_pos_t(2, false, 3));

  std::string filename = gbwt::TempFile::getName("minimizer");
  std::ofstream out(filename, std::ios_base::binary);
  index.serialize(out);
  out.close();

  MinimizerIndex<TypeParam> copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.deserialize(in);
  in.close();
  gbwt::TempFile::remove(filename);

  EXPECT_EQ(index, copy) << "Loaded index is not identical to the original";
}

//------------------------------------------------------------------------------

/*
  Order of 3-mers using Key64:
  AAT < TGT < TTG < TAT < ATA < TCG < ATT < ACA < GAA < ACT < TAC < CGA < CAA < GTA < TTC < AGT

  Order of 3-mers using Key128:
  ATT < TCG < ATA < TAT < TTG < TGT < AAT < AGT < TTC < GTA < CAA < TAC < CGA < ACT < ACA < GAA
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

  void check_weighted_minimizers(const std::string& query, size_t k, size_t w)
  {
    MinimizerIndex<KeyType> index(k, w);
    std::vector<typename MinimizerIndex<KeyType>::minimizer_type> result = index.minimizers(query);
    std::vector<std::pair<typename MinimizerIndex<KeyType>::minimizer_type, size_t>> weighted = index.weighted_minimizers(query);

    std::stringstream ss;
    ss << "(" << k << ", " << w << ")-minimizers in " << query;
    std::string test_description = ss.str();
    size_t correct_weight = query.length() + 2 - k - w;

    ASSERT_EQ(weighted.size(), result.size()) << "Wrong number of weighted " << test_description;
    size_t total_weight = 0;
    bool same_minimizers = true;
    for(size_t i = 0; i < result.size(); i++)
    {
      if(weighted[i].first != result[i]) { same_minimizers = false; }
      total_weight += weighted[i].second;
    }
    EXPECT_TRUE(same_minimizers) << "Incorrect weighted " << test_description;
    EXPECT_EQ(total_weight, correct_weight) << "Incorrect total weight for " << test_description;
  }
};

TYPED_TEST_CASE(MinimizerExtraction, KeyTypes);

TYPED_TEST(MinimizerExtraction, KeyEncoding)
{
  size_t k = TypeParam::KMER_MAX_LENGTH - 1;
  std::string bases = "ACGT";

  TypeParam forward_key, reverse_key;
  size_t valid_chars = 0;
  std::stringstream result;
  for(size_t i = 0; i < 2 * k + 1; i++)
  {
    char c = bases[i % bases.length()];
    forward_key.forward(k, c, valid_chars);
    reverse_key.reverse(k, c);
    result << c;
  }
  std::string correct = result.str();
  correct = correct.substr(correct.length() - k);
  std::string reverse = reverse_complement(correct);

  EXPECT_EQ(forward_key.decode(k), correct) << "Invalid forward key " << forward_key;
  EXPECT_EQ(reverse_key.decode(k), reverse) << "Invalid reverse key " << reverse_key;
}

TYPED_TEST(MinimizerExtraction, LeftmostOccurrence)
{
  MinimizerIndex<TypeParam> index(3, 2);
  typename MinimizerIndex<TypeParam>::minimizer_type correct = (TypeParam::KEY_BITS == 128 ?
    get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 3 * 1, 4, true) : // ATT
    get_minimizer<TypeParam>(0 * 16 + 0 * 4 + 3 * 1, 2));       // AAT
  typename MinimizerIndex<TypeParam>::minimizer_type result = index.minimizer(this->str.begin(), this->str.end());
  EXPECT_EQ(result, correct) << "The leftmost minimizer was not found";
}

TYPED_TEST(MinimizerExtraction, AllMinimizers)
{
  MinimizerIndex<TypeParam> index(3, 2);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>(3 * 16 + 1 * 4 + 2 * 1, 2, true),   // TCG
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 3, false),  // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 3 * 1, 4, true),   // ATT
      get_minimizer<TypeParam>(3 * 16 + 2 * 4 + 3 * 1, 7, true),   // TGT
      get_minimizer<TypeParam>(3 * 16 + 3 * 4 + 2 * 1, 8, true),   // TTG
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 8, false),  // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 3 * 1, 9, true),   // ATT
      get_minimizer<TypeParam>(0 * 16 + 2 * 4 + 3 * 1, 12, true)   // AGT
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>(3 * 16 + 1 * 4 + 2 * 1, 2, true),   // TCG
      get_minimizer<TypeParam>(0 * 16 + 0 * 4 + 3 * 1, 2, false),  // AAT
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 5, true),   // TAT
      get_minimizer<TypeParam>(3 * 16 + 2 * 4 + 3 * 1, 7, true),   // TGT
      get_minimizer<TypeParam>(0 * 16 + 0 * 4 + 3 * 1, 7, false),  // AAT
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 10, true),  // TAT
      get_minimizer<TypeParam>(0 * 16 + 1 * 4 + 3 * 1, 10, false)  // ACT
    };
  }
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.minimizers(this->str.begin(), this->str.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, WindowLength)
{
  MinimizerIndex<TypeParam> index(3, 3);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 3, false),  // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 3 * 1, 4, true),   // ATT
      get_minimizer<TypeParam>(3 * 16 + 3 * 4 + 2 * 1, 8, true),   // TTG
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 8, false),  // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 3 * 1, 9, true)    // ATT
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>(0 * 16 + 0 * 4 + 3 * 1, 2, false),  // AAT
      get_minimizer<TypeParam>(3 * 16 + 2 * 4 + 3 * 1, 7, true),   // TGT
      get_minimizer<TypeParam>(0 * 16 + 0 * 4 + 3 * 1, 7, false),  // AAT
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 10, true)   // TAT
    };
  }
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.minimizers(this->str.begin(), this->str.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

// We have multiple occurrences of the same minimizer in the initial window and after it.
// When we find the first occurrence after the window, one of the original occurrences is
// still in the buffer.
TYPED_TEST(MinimizerExtraction, AllOccurrences)
{
  MinimizerIndex<TypeParam> index(3, 3);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 1, false), // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 2, true),  // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 3, false), // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 4, true)   // ATA
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 0, false),  // TAT
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 2, false),  // TAT
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 3, true),   // TAT
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 5, true)    // TAT
    };
  }
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.minimizers(this->repetitive.begin(), this->repetitive.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, WeightedMinimizers)
{
  this->check_weighted_minimizers(this->str, 3, 2);
  this->check_weighted_minimizers(this->repetitive, 3, 3);
}

TYPED_TEST(MinimizerExtraction, InvalidCharacters)
{
  std::string weird = "CGAATAxAATACT";
  MinimizerIndex<TypeParam> index(3, 2);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>(3 * 16 + 1 * 4 + 2 * 1, 2, true),   // TCG
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 3, false),  // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 3 * 1, 4, true),   // ATT
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 0 * 1, 8, false),  // ATA
      get_minimizer<TypeParam>(0 * 16 + 3 * 4 + 3 * 1, 9, true),   // ATT
      get_minimizer<TypeParam>(0 * 16 + 2 * 4 + 3 * 1, 12, true)   // AGT
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>(3 * 16 + 1 * 4 + 2 * 1, 2, true),   // TCG
      get_minimizer<TypeParam>(0 * 16 + 0 * 4 + 3 * 1, 2, false),  // AAT
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 5, true),   // TAT
      get_minimizer<TypeParam>(0 * 16 + 0 * 4 + 3 * 1, 7, false),  // AAT
      get_minimizer<TypeParam>(3 * 16 + 0 * 4 + 3 * 1, 10, true),  // TAT
      get_minimizer<TypeParam>(0 * 16 + 1 * 4 + 3 * 1, 10, false)  // ACT
    };
  }
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.minimizers(weird.begin(), weird.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, BothOrientations)
{
  MinimizerIndex<TypeParam> index(3, 2);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> forward_minimizers = index.minimizers(this->str.begin(), this->str.end());
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> reverse_minimizers = index.minimizers(this->rev.begin(), this->rev.end());
  ASSERT_EQ(forward_minimizers.size(), reverse_minimizers.size()) << "Different number of minimizers in forward and reverse orientations";
  for(size_t i = 0; i < forward_minimizers.size(); i++)
  {
    typename MinimizerIndex<TypeParam>::minimizer_type& f = forward_minimizers[i];
    typename MinimizerIndex<TypeParam>::minimizer_type& r = reverse_minimizers[forward_minimizers.size() - 1 - i];
    EXPECT_EQ(f.key, r.key) << "Wrong key for minimizer " << i;
    EXPECT_EQ(f.offset, this->str.length() - 1 - r.offset) << "Wrong offset for minimizer " << i;
    EXPECT_NE(f.is_reverse, r.is_reverse) << "Wrong orientation for minimizer " << i;
  }
}

//------------------------------------------------------------------------------

template<class KeyType>
class CorrectKmers : public ::testing::Test
{
public:
  typedef std::map<KeyType, std::set<pos_t>> result_type;

  size_t total_keys;

  virtual void SetUp()
  {
    this->total_keys = 16;
  }

  void check_minimizer_index(const MinimizerIndex<KeyType>& index,
                             const result_type& correct_values,
                             size_t keys, size_t values, size_t unique)
  {
    ASSERT_EQ(index.size(), keys) << "Wrong number of keys";
    ASSERT_EQ(index.values(), values) << "Wrong number of values";
    EXPECT_EQ(index.unique_keys(), unique) << "Wrong number of unique keys";

    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      typename MinimizerIndex<KeyType>::minimizer_type minimizer = get_minimizer<KeyType>(iter->first);
      std::vector<pos_t> correct(iter->second.begin(), iter->second.end());

      size_t count = index.count(minimizer);
      EXPECT_EQ(count, correct.size()) << "Wrong number of occurrences for key " << iter->first;
      if(count != correct.size()) { continue; }

      std::vector<pos_t> occs = index.find(minimizer);
      EXPECT_EQ(occs, correct) << "Wrong positions for key " << iter->first;

      std::pair<size_t, const typename MinimizerIndex<KeyType>::code_type*> raw_occs = index.count_and_find(minimizer);
      bool ok = true;
      // Calling anything related to the test framework may invalidate the pointers.
      if(raw_occs.first == correct.size())
      {
        for(size_t i = 0; i < raw_occs.first; i++)
        {
          if(MinimizerIndex<KeyType>::decode(raw_occs.second[i]) != correct[i])
          {
            ok = false; break;
          }
        }
      }
      EXPECT_EQ(raw_occs.first, correct.size()) << "Wrong number of raw occurrences for key " << iter->first;
      EXPECT_TRUE(ok) << "Wrong raw positions for key " << iter->first;
    }
  }
};

TYPED_TEST_CASE(CorrectKmers, KeyTypes);

TYPED_TEST(CorrectKmers, UniqueKeys)
{
  MinimizerIndex<TypeParam> index;
  size_t keys = 0, values = 0, unique = 0;
  typename TestFixture::result_type correct_values;

  for(size_t i = 1; i <= this->total_keys; i++)
  {
    pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
    correct_values[i].insert(pos);
    keys++; values++; unique++;
  }
  this->check_minimizer_index(index, correct_values, keys, values, unique);
}

TYPED_TEST(CorrectKmers, MissingKeys)
{
  MinimizerIndex<TypeParam> index;
  using code_type = typename MinimizerIndex<TypeParam>::code_type;
  for(size_t i = 1; i <= this->total_keys; i++)
  {
    index.insert(get_minimizer<TypeParam>(i), make_pos_t(i, i & 1, i & MinimizerIndex<TypeParam>::OFF_MASK));
  }
  for(size_t i = this->total_keys + 1; i <= 2 * this->total_keys; i++)
  {
    typename MinimizerIndex<TypeParam>::minimizer_type minimizer = get_minimizer<TypeParam>(i);
    EXPECT_EQ(index.count(minimizer), static_cast<size_t>(0)) << "Non-zero occurrences for key " << i;
    EXPECT_TRUE(index.find(minimizer).empty()) << "Non-empty value for key " << i;
    std::pair<size_t, const code_type*> correct_raw(0, nullptr);
    EXPECT_EQ(index.count_and_find(minimizer), correct_raw) << "Non-empty raw occurrences for key " << i;
  }
}

TYPED_TEST(CorrectKmers, EmptyKeysValues)
{
  MinimizerIndex<TypeParam> index;
  using code_type = typename MinimizerIndex<TypeParam>::code_type;
  std::pair<size_t, const code_type*> correct_raw(0, nullptr);

  typename MinimizerIndex<TypeParam>::minimizer_type empty_key = get_minimizer<TypeParam>(MinimizerIndex<TypeParam>::key_type::no_key());
  index.insert(empty_key, make_pos_t(1, false, 0));
  EXPECT_EQ(index.count(empty_key), static_cast<size_t>(0)) << "Non-zero occurrences for empty key";
  EXPECT_TRUE(index.find(empty_key).empty()) << "Non-empty value for empty key";
  EXPECT_EQ(index.count_and_find(empty_key), correct_raw) << "Non-empty raw occurrences for empty key";

  typename MinimizerIndex<TypeParam>::minimizer_type key = get_minimizer<TypeParam>(this->total_keys + 1);
  index.insert(key, MinimizerIndex<TypeParam>::decode(MinimizerIndex<TypeParam>::NO_VALUE));
  EXPECT_EQ(index.count(key), static_cast<size_t>(0)) << "Non-zero occurrences after inserting empty value";
  EXPECT_TRUE(index.find(key).empty()) << "Non-empty value after inserting empty value";
  EXPECT_EQ(index.count_and_find(key), correct_raw) << "Non-empty raw occurrences after inserting empty value";
}

TYPED_TEST(CorrectKmers, MultipleOccurrences)
{
  MinimizerIndex<TypeParam> index;
  size_t keys = 0, values = 0, unique = 0;
  typename TestFixture::result_type correct_values;

  for(size_t i = 1; i <= this->total_keys; i++)
  {
    pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
    correct_values[i].insert(pos);
    keys++; values++; unique++;
  }
  for(size_t i = 1; i <= this->total_keys; i += 2)
  {
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
    correct_values[i].insert(pos);
    values++; unique--;
  }
  for(size_t i = 1; i <= this->total_keys; i += 4)
  {
    pos_t pos = make_pos_t(i + 2, i & 1, (i + 2) & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
    correct_values[i].insert(pos);
    values++;
  }
  this->check_minimizer_index(index, correct_values, keys, values, unique);
}

TYPED_TEST(CorrectKmers, DuplicateValues)
{
  MinimizerIndex<TypeParam> index;
  size_t keys = 0, values = 0, unique = 0;
  typename TestFixture::result_type correct_values;

  for(size_t i = 1; i <= this->total_keys; i++)
  {
    pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
    correct_values[i].insert(pos);
    keys++; values++; unique++;
  }
  for(size_t i = 1; i <= this->total_keys; i += 2)
  {
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
    correct_values[i].insert(pos);
    values++; unique--;
  }
  for(size_t i = 1; i <= this->total_keys; i += 4)
  {
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
  }
  this->check_minimizer_index(index, correct_values, keys, values, unique);
}

TYPED_TEST(CorrectKmers, Rehashing)
{
  MinimizerIndex<TypeParam> index;
  size_t keys = 0, values = 0, unique = 0;
  typename TestFixture::result_type correct_values;
  size_t threshold = index.max_keys();

  for(size_t i = 1; i <= threshold; i++)
  {
    pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
    correct_values[i].insert(pos);
    keys++; values++; unique++;
  }
  ASSERT_EQ(index.max_keys(), threshold) << "Index capacity changed at threshold";

  {
    size_t i = threshold + 1;
    pos_t pos = make_pos_t(i, i & 1, i & MinimizerIndex<TypeParam>::OFF_MASK);
    index.insert(get_minimizer<TypeParam>(i), pos);
    correct_values[i].insert(pos);
    keys++; values++; unique++;
  }
  EXPECT_GT(index.max_keys(), threshold) << "Index capacity not increased after threshold";

  this->check_minimizer_index(index, correct_values, keys, values, unique);
}

//------------------------------------------------------------------------------

} // namespace
