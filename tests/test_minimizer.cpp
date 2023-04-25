#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <sstream>
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
  EXPECT_EQ(alt_index, alt_copy) << "A copy of a parameterized index is not identical to the original";
  EXPECT_NE(default_index, alt_index) << "Default and parameterized indexes are identical";
}

TYPED_TEST(ObjectManipulation, SyncmerIndex)
{
  MinimizerIndex<TypeParam> default_minimizer;
  MinimizerIndex<TypeParam> default_syncmer(true);
  MinimizerIndex<TypeParam> short_minimizer(15, 6);
  MinimizerIndex<TypeParam> short_syncmer(15, 6, true);
  EXPECT_FALSE(default_minimizer.uses_syncmers()) << "Default minimizer index uses syncmers";
  EXPECT_TRUE(default_syncmer.uses_syncmers()) << "Default syncmer index uses minimizers";
  EXPECT_FALSE(short_minimizer.uses_syncmers()) << "Parameterized minimizer index uses syncmers";
  EXPECT_TRUE(short_syncmer.uses_syncmers()) << "Parameterized syncmer index uses minimizers";
}

TYPED_TEST(ObjectManipulation, Contents)
{
  MinimizerIndex<TypeParam> default_index;
  MinimizerIndex<TypeParam> default_copy(default_index);

  // Different contents.
  default_index.insert(get_minimizer<TypeParam>(1), make_pos_t(1, false, 3), Payload::create(hash(1, false, 3)));
  EXPECT_NE(default_index, default_copy) << "Empty index is identical to nonempty index";

  // Same key, different value.
  default_copy.insert(get_minimizer<TypeParam>(1), make_pos_t(2, false, 3), Payload::create(hash(2, false, 3)));
  EXPECT_NE(default_index, default_copy) << "Indexes with different values are identical";

  // Same contents.
  default_copy = default_index;
  EXPECT_EQ(default_index, default_copy) << "A copy of a nonempty index is not identical to the original";
}

TYPED_TEST(ObjectManipulation, Swap)
{
  MinimizerIndex<TypeParam> first, second;
  first.insert(get_minimizer<TypeParam>(1), make_pos_t(1, false, 3), Payload::create(hash(1, false, 3)));
  second.insert(get_minimizer<TypeParam>(2), make_pos_t(2, false, 3), Payload::create(hash(2, false, 3)));

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
  index.insert(get_minimizer<TypeParam>(1), make_pos_t(1, false, 3), Payload::create(hash(1, false, 3)));
  index.insert(get_minimizer<TypeParam>(2), make_pos_t(1, false, 3), Payload::create(hash(1, false, 3)));
  index.insert(get_minimizer<TypeParam>(2), make_pos_t(2, false, 3), Payload::create(hash(2, false, 3)));

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

template<class KeyType>
class KeyEncodeDecode : public ::testing::Test
{
};

TYPED_TEST_CASE(KeyEncodeDecode, KeyTypes);

TYPED_TEST(KeyEncodeDecode, SimpleSequence)
{
  for (size_t k = 1; k <= TypeParam::KMER_MAX_LENGTH; k++) {
    for (auto base : {'A', 'C', 'G', 'T'}) {
      std::stringstream gen;
      for (size_t i = 0; i < k; i++) {
        gen << base;
      }
      std::string kmer = gen.str();
      
      TypeParam encoded = TypeParam::encode(kmer);
      std::string decoded = encoded.decode(k);
      
      EXPECT_EQ(decoded, kmer) << "Decoded kmer is not identical to original";
    }
  }
}

TYPED_TEST(KeyEncodeDecode, ComplexSequence)
{
  for (size_t k = 1; k <= TypeParam::KMER_MAX_LENGTH; k++) {
    std::stringstream gen;
    for (size_t i = 0; i < k; i+= 7) {
      gen << "GATTACA";
    }
    std::string kmer = gen.str().substr(0, k);
    
    TypeParam encoded = TypeParam::encode(kmer);
    std::string decoded = encoded.decode(k);
    
    EXPECT_EQ(decoded, kmer) << "Decoded kmer is not identical to original";
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

TYPED_TEST(MinimizerExtraction, AllMinimizers)
{
  MinimizerIndex<TypeParam> index(3, 2);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>("TCG", 2, true),
      get_minimizer<TypeParam>("ATA", 3, false),
      get_minimizer<TypeParam>("ATT", 4, true),
      get_minimizer<TypeParam>("TGT", 7, true),
      get_minimizer<TypeParam>("TTG", 8, true),
      get_minimizer<TypeParam>("ATA", 8, false),
      get_minimizer<TypeParam>("ATT", 9, true),
      get_minimizer<TypeParam>("AGT", 12, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>("TCG", 2, true),
      get_minimizer<TypeParam>("AAT", 2, false),
      get_minimizer<TypeParam>("TAT", 5, true),
      get_minimizer<TypeParam>("TGT", 7, true),
      get_minimizer<TypeParam>("AAT", 7, false),
      get_minimizer<TypeParam>("TAT", 10, true),
      get_minimizer<TypeParam>("ACT", 10, false)
    };
  }

  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.minimizers(this->str.begin(), this->str.end());
  ASSERT_EQ(result, correct) << "Did not find the correct minimizers";
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> fallback = index.syncmers(this->str.begin(), this->str.end());
  EXPECT_EQ(fallback, correct) << "Did not find the correct minimizers using syncmers()";
}

TYPED_TEST(MinimizerExtraction, ClosedSyncmers)
{
  MinimizerIndex<TypeParam> index(5, 3, true);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>("ATACA", 3, false),
      get_minimizer<TypeParam>("ATTCG", 4, true),
      get_minimizer<TypeParam>("GTATT", 6, true),
      get_minimizer<TypeParam>("TTGTA", 8, true),
      get_minimizer<TypeParam>("ATACT", 8, false),
      get_minimizer<TypeParam>("ATTGT", 9, true),
      get_minimizer<TypeParam>("GTATT", 11, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>("AATAC", 2, false),
      get_minimizer<TypeParam>("ATACA", 3, false),
      get_minimizer<TypeParam>("ATTCG", 4, true),
      get_minimizer<TypeParam>("ACAAT", 5, false),
      get_minimizer<TypeParam>("AATAC", 7, false),
      get_minimizer<TypeParam>("AGTAT", 12, true)
    };
  }

  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.syncmers(this->str.begin(), this->str.end());
  ASSERT_EQ(result, correct) << "Did not find the correct syncmers";
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> fallback = index.minimizers(this->str.begin(), this->str.end());
  EXPECT_EQ(fallback, correct) << "Did not find the correct syncmers using minimizers()";
}

TYPED_TEST(MinimizerExtraction, AllMinimizersWithRegions)
{
  MinimizerIndex<TypeParam> index(3, 2);
  std::vector<std::tuple<typename MinimizerIndex<TypeParam>::minimizer_type, size_t, size_t>> correct;
  if(TypeParam::KEY_BITS == 128)
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
      std::make_tuple(get_minimizer<TypeParam>("TCG", 2, true), 0, 4),
      std::make_tuple(get_minimizer<TypeParam>("ATA", 3, false), 3, 4),
      std::make_tuple(get_minimizer<TypeParam>("ATT", 4, true), 1, 5),
      std::make_tuple(get_minimizer<TypeParam>("TGT", 7, true), 4, 4),
      std::make_tuple(get_minimizer<TypeParam>("TTG", 8, true), 5, 4),
      std::make_tuple(get_minimizer<TypeParam>("ATA", 8, false), 8, 4),
      std::make_tuple(get_minimizer<TypeParam>("ATT", 9, true), 6, 5),
      std::make_tuple(get_minimizer<TypeParam>("AGT", 12, true), 9, 4)
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
      std::make_tuple(get_minimizer<TypeParam>("TCG", 2, true), 0, 4),
      std::make_tuple(get_minimizer<TypeParam>("AAT", 2, false), 1, 5),
      std::make_tuple(get_minimizer<TypeParam>("TAT", 5, true), 3, 4),
      std::make_tuple(get_minimizer<TypeParam>("TGT", 7, true), 4, 5),
      std::make_tuple(get_minimizer<TypeParam>("AAT", 7, false), 6, 5),
      std::make_tuple(get_minimizer<TypeParam>("TAT", 10, true), 8, 4),
      std::make_tuple(get_minimizer<TypeParam>("ACT", 10, false), 9, 4)
    };
  }
  std::vector<std::tuple<typename MinimizerIndex<TypeParam>::minimizer_type, size_t, size_t>> result = index.minimizer_regions(this->str.begin(), this->str.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TEST(MinimizerExtraction, HardMinimizersWithRegion)
{
  using TypeParam = Key128;

  // Here's a case I caught not working correctly.
  MinimizerIndex<TypeParam> index(29, 11);
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
  
  std::vector<std::tuple<typename MinimizerIndex<TypeParam>::minimizer_type, size_t, size_t>> correct;
  correct =
  {
    std::make_tuple(get_minimizer<TypeParam>("CACAAGCTGCTCTTTCCCTCAATTGTTCA", 12, false), 7, 44), // #2
    std::make_tuple(get_minimizer<TypeParam>("TTTCCCTCAATTGTTCATTTGTCTCCTTG", 24, false), 18, 42), // #4
    std::make_tuple(get_minimizer<TypeParam>("ATTGAGGGAAAGAGCAGCTTGTGTAAAAA", 34, true), 0, 45), // #1
    std::make_tuple(get_minimizer<TypeParam>("ACAAATGAACAATTGAGGGAAAGAGCAGC", 45, true), 13, 43) // #3
  };
  std::vector<std::tuple<typename MinimizerIndex<TypeParam>::minimizer_type, size_t, size_t>> result = index.minimizer_regions(seq.begin(), seq.end());
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
      get_minimizer<TypeParam>("ATA", 3, false),
      get_minimizer<TypeParam>("ATT", 4, true),
      get_minimizer<TypeParam>("TTG", 8, true),
      get_minimizer<TypeParam>("ATA", 8, false),
      get_minimizer<TypeParam>("ATT", 9, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>("AAT", 2, false),
      get_minimizer<TypeParam>("TGT", 7, true),
      get_minimizer<TypeParam>("AAT", 7, false),
      get_minimizer<TypeParam>("TAT", 10, true)
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
      get_minimizer<TypeParam>("ATA", 1, false),
      get_minimizer<TypeParam>("ATA", 2, true),
      get_minimizer<TypeParam>("ATA", 3, false),
      get_minimizer<TypeParam>("ATA", 4, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>("TAT", 0, false),
      get_minimizer<TypeParam>("TAT", 2, false),
      get_minimizer<TypeParam>("TAT", 3, true),
      get_minimizer<TypeParam>("TAT", 5, true)
    };
  }
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.minimizers(this->repetitive.begin(), this->repetitive.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, WeirdSyncmers)
{
  // The string is the reverse complement of itself. The middle smers AAT and ATT
  // are the smallest ones using both key types.
  std::string weird = "CAATTG";
  MinimizerIndex<TypeParam> index(5, 3, true);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>("AATTG", 1, false),
      get_minimizer<TypeParam>("AATTG", 4, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>("CAATT", 0, false),
      get_minimizer<TypeParam>("CAATT", 5, true)
    };
  }
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.syncmers(weird);
  EXPECT_EQ(result, correct) << "Did not find the correct syncmers";
}

TYPED_TEST(MinimizerExtraction, InvalidMinimizerCharacters)
{
  std::string weird = "CGAATAxAATACT";
  MinimizerIndex<TypeParam> index(3, 2);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>("TCG", 2, true),
      get_minimizer<TypeParam>("ATA", 3, false),
      get_minimizer<TypeParam>("ATT", 4, true),
      get_minimizer<TypeParam>("ATA", 8, false),
      get_minimizer<TypeParam>("ATT", 9, true),
      get_minimizer<TypeParam>("AGT", 12, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>("TCG", 2, true),
      get_minimizer<TypeParam>("AAT", 2, false),
      get_minimizer<TypeParam>("TAT", 5, true),
      get_minimizer<TypeParam>("AAT", 7, false),
      get_minimizer<TypeParam>("TAT", 10, true),
      get_minimizer<TypeParam>("ACT", 10, false)
    };
  }
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.minimizers(weird.begin(), weird.end());
  EXPECT_EQ(result, correct) << "Did not find the correct minimizers";
}

TYPED_TEST(MinimizerExtraction, InvalidSyncmerCharacters)
{
  std::string weird = "CGAATAxAATACT";
  MinimizerIndex<TypeParam> index(5, 3, true);
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> correct;
  if(TypeParam::KEY_BITS == 128)
  {
    correct =
    {
      get_minimizer<TypeParam>("ATTCG", 4, true),
      get_minimizer<TypeParam>("ATACT", 8, false),
      get_minimizer<TypeParam>("GTATT", 11, true)
    };
  }
  else
  {
    correct =
    {
      get_minimizer<TypeParam>("ATTCG", 4, true),
      get_minimizer<TypeParam>("AATAC", 7, false),
      get_minimizer<TypeParam>("AGTAT", 12, true)
    };
  }
  std::vector<typename MinimizerIndex<TypeParam>::minimizer_type> result = index.minimizers(weird.begin(), weird.end());
  EXPECT_EQ(result, correct) << "Did not find the correct syncmers";
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
  typedef std::map<KeyType, std::set<std::pair<pos_t, Payload>>> result_type;

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
      std::vector<std::pair<pos_t, Payload>> correct(iter->second.begin(), iter->second.end());

      size_t count = index.count(minimizer);
      EXPECT_EQ(count, correct.size()) << "Wrong number of occurrences for key " << iter->first;
      if(count != correct.size()) { continue; }

      std::vector<std::pair<pos_t, Payload>> occs = index.find(minimizer);
      EXPECT_EQ(occs, correct) << "Wrong positions for key " << iter->first;

      std::pair<size_t, const PositionPayload*> raw_occs = index.count_and_find(minimizer);
      bool ok = true;
      // Calling anything related to the test framework may invalidate the pointers.
      if(raw_occs.first == correct.size())
      {
        for(size_t i = 0; i < raw_occs.first; i++)
        {
          if(raw_occs.second[i].position.decode() != correct[i].first ||
             raw_occs.second[i].payload != correct[i].second)
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
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
    correct_values[i].insert(std::make_pair(pos, payload));
    keys++; values++; unique++;
  }
  this->check_minimizer_index(index, correct_values, keys, values, unique);
}

TYPED_TEST(CorrectKmers, MissingKeys)
{
  MinimizerIndex<TypeParam> index;
  for(size_t i = 1; i <= this->total_keys; i++)
  {
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
  }
  for(size_t i = this->total_keys + 1; i <= 2 * this->total_keys; i++)
  {
    typename MinimizerIndex<TypeParam>::minimizer_type minimizer = get_minimizer<TypeParam>(i);
    EXPECT_EQ(index.count(minimizer), static_cast<size_t>(0)) << "Non-zero occurrences for key " << i;
    EXPECT_TRUE(index.find(minimizer).empty()) << "Non-empty value for key " << i;
    std::pair<size_t, const PositionPayload*> correct_raw(0, nullptr);
    EXPECT_EQ(index.count_and_find(minimizer), correct_raw) << "Non-empty raw occurrences for key " << i;
  }
}

TYPED_TEST(CorrectKmers, EmptyKeysValues)
{
  MinimizerIndex<TypeParam> index;
  std::pair<size_t, const PositionPayload*> correct_raw(0, nullptr);

  typename MinimizerIndex<TypeParam>::minimizer_type empty_key = get_minimizer<TypeParam>(MinimizerIndex<TypeParam>::key_type::no_key());
  index.insert(empty_key, make_pos_t(1, false, 0));
  EXPECT_EQ(index.count(empty_key), static_cast<size_t>(0)) << "Non-zero occurrences for empty key";
  EXPECT_TRUE(index.find(empty_key).empty()) << "Non-empty value for empty key";
  EXPECT_EQ(index.count_and_find(empty_key), correct_raw) << "Non-empty raw occurrences for empty key";

  typename MinimizerIndex<TypeParam>::minimizer_type key = get_minimizer<TypeParam>(this->total_keys + 1);
  index.insert(key, Position::no_value().decode());
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
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
    correct_values[i].insert(std::make_pair(pos, payload));
    keys++; values++; unique++;
  }
  for(size_t i = 1; i <= this->total_keys; i += 2)
  {
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, (i + 1) & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
    correct_values[i].insert(std::make_pair(pos, payload));
    values++; unique--;
  }
  for(size_t i = 1; i <= this->total_keys; i += 4)
  {
    pos_t pos = make_pos_t(i + 2, i & 1, (i + 2) & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, (i + 2) & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
    correct_values[i].insert(std::make_pair(pos, payload));
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
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
    correct_values[i].insert(std::make_pair(pos, payload));
    keys++; values++; unique++;
  }
  for(size_t i = 1; i <= this->total_keys; i += 2)
  {
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, (i + 1) & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
    correct_values[i].insert(std::make_pair(pos, payload));
    values++; unique--;
  }
  for(size_t i = 1; i <= this->total_keys; i += 4)
  {
    // Also check that inserting duplicates does not change the payload.
    pos_t pos = make_pos_t(i + 1, i & 1, (i + 1) & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, (i + 1) & Position::OFF_MASK) + 1);
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
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
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
    correct_values[i].insert(std::make_pair(pos, payload));
    keys++; values++; unique++;
  }
  ASSERT_EQ(index.max_keys(), threshold) << "Index capacity changed at threshold";

  {
    size_t i = threshold + 1;
    pos_t pos = make_pos_t(i, i & 1, i & Position::OFF_MASK);
    Payload payload = Payload::create(hash(i, i & 1, i & Position::OFF_MASK));
    index.insert(get_minimizer<TypeParam>(i), pos, payload);
    correct_values[i].insert(std::make_pair(pos, payload));
    keys++; values++; unique++;
  }
  EXPECT_GT(index.max_keys(), threshold) << "Index capacity not increased after threshold";

  this->check_minimizer_index(index, correct_values, keys, values, unique);
}

//------------------------------------------------------------------------------

class HitsInSubgraphTest : public ::testing::Test
{
public:
  typedef std::vector<std::pair<pos_t, Payload>> result_type;

  void check_results(const std::unordered_set<nid_t>& subgraph, const std::vector<PositionPayload>& hits,
                     const result_type& expected_result, const std::string& test_case)
  {
    std::vector<nid_t> sorted_subgraph(subgraph.begin(), subgraph.end());
    std::sort(sorted_subgraph.begin(), sorted_subgraph.end());

    result_type result;
    hits_in_subgraph(hits.size(), hits.data(), subgraph, [&](pos_t pos, Payload payload)
    {
      result.emplace_back(pos, payload);
    });
    ASSERT_EQ(result, expected_result) << test_case << ": Incorrect results with the naive algorithm";

    result.clear();
    hits_in_subgraph(hits.size(), hits.data(), sorted_subgraph, [&](pos_t pos, Payload payload)
    {
      result.emplace_back(pos, payload);
    });
    ASSERT_EQ(result, expected_result) << test_case << ": Incorrect results with exponential search";
  }

  std::tuple<std::unordered_set<nid_t>, std::vector<PositionPayload>, result_type>
  create_test_case(size_t universe_size, size_t begin, size_t end,
                   double interval_prob, double outlier_prob, double hit_prob, size_t random_seed) const
  {
    std::unordered_set<nid_t> subgraph;
    std::vector<PositionPayload> hits;
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
  std::vector<PositionPayload> hits;
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
    std::vector<PositionPayload> hits;
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
    std::vector<PositionPayload> hits;
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
