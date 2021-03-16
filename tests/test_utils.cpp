#include <gtest/gtest.h>

#include <gbwtgraph/utils.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class SourceTest : public ::testing::Test
{
public:
  void check_nodes(const SequenceSource& source, const std::vector<node_type>& truth) const
  {
    ASSERT_EQ(source.get_node_count(), truth.size()) << "Incorrect number of nodes";
    for(const node_type& node : truth)
    {
      ASSERT_TRUE(source.has_node(node.first)) << "Node id " << node.first << " is missing from the sequence source";
      EXPECT_EQ(source.get_length(node.first), node.second.length()) << "Invalid sequence length for node " << node.first;
      EXPECT_EQ(source.get_sequence(node.first), node.second) << "Invalid sequence for node " << node.first;
      view_type view = source.get_sequence_view(node.first);
      EXPECT_EQ(std::string(view.first, view.second), node.second) << "Invalid sequence view for node " << node.first;
    }
  }

  void check_translation(const SequenceSource& source, const std::vector<translation_type>& truth) const
  {
    ASSERT_EQ(source.uses_translation(), !(truth.empty())) << "Segment translation is not used as expected";
    if(!(source.uses_translation())) { return; }
    ASSERT_EQ(source.segment_translation.size(), truth.size()) << "Invalid number of segments";
    for(const translation_type& translation : truth)
    {
      EXPECT_EQ(source.get_translation(translation.first), translation.second) << "Invalid translation for " << translation.first;
    }
  }

  void check_inverse_translation(const SequenceSource& source) const
  {
    auto result = source.invert_translation();
    ASSERT_EQ(result.first.size(), source.segment_translation.size()) << "Invalid number of segments in inverse translation";
    ASSERT_EQ(result.second.size(), size_t(source.next_id)) << "Invalid number of nodes in inverse translation";
    ASSERT_EQ(result.first.size(), result.second.ones()) << "Inconsistent number of segments in inverse translation";
    auto iter = result.second.one_begin();
    for(size_t i = 0; i < result.first.size(); i++)
    {
      std::string segment = result.first.str(i);
      std::pair<nid_t, nid_t> translation = source.get_translation(segment);
      EXPECT_EQ(nid_t(iter->second), translation.first) << "Invalid start for segment " << segment;
      ++iter;
      EXPECT_EQ(nid_t(iter->second), translation.second) << "Invalid limit for segment " << segment;
    }
  }
};

TEST_F(SourceTest, EmptySource)
{
  SequenceSource source;
  std::vector<node_type> nodes;
  std::vector<translation_type> translation;

  this->check_nodes(source, nodes);
  this->check_translation(source, translation);
  this->check_inverse_translation(source);
}

TEST_F(SourceTest, AddNodes)
{
  SequenceSource source;
  build_source(source);

  std::vector<node_type> nodes =
  {
    { nid_t(1), "G" },
    { nid_t(2), "A" },
    { nid_t(3), "T" },
    { nid_t(4), "GGG" },
    { nid_t(5), "T" },
    { nid_t(6), "A" },
    { nid_t(7), "C" },
    { nid_t(8), "A" },
    { nid_t(9), "A" },
  };
  std::vector<translation_type> translation;

  this->check_nodes(source, nodes);
  this->check_translation(source, translation);
  this->check_inverse_translation(source);
}

TEST_F(SourceTest, TranslateSegments)
{
  SequenceSource source;
  std::vector<std::pair<std::string, std::string>> segments =
  {
    { "s1", "G" },
    { "s2", "A" },
    { "s3", "T" },
    { "s4", "GGGT" },
    { "s6", "A" },
    { "s7", "C" },
    { "s8", "A" },
    { "s9", "A" },
  };
  for(const std::pair<std::string, std::string>& segment : segments)
  {
    source.translate_segment(segment.first, get_view(segment.second), 3);
  }

  std::vector<node_type> nodes =
  {
    { nid_t(1), "G" },
    { nid_t(2), "A" },
    { nid_t(3), "T" },
    { nid_t(4), "GGG" },
    { nid_t(5), "T" },
    { nid_t(6), "A" },
    { nid_t(7), "C" },
    { nid_t(8), "A" },
    { nid_t(9), "A" },
  };
  std::vector<translation_type> translation =
  {
    { "s1", { 1, 2 } },
    { "s2", { 2, 3 } },
    { "s3", { 3, 4 } },
    { "s4", { 4, 6 } },
    { "s6", { 6, 7 } },
    { "s7", { 7, 8 } },
    { "s8", { 8, 9 } },
    { "s9", { 9, 10 } },
  };

  this->check_nodes(source, nodes);
  this->check_translation(source, translation);  
  this->check_inverse_translation(source);
}

//------------------------------------------------------------------------------

class StringArrayTest : public ::testing::Test
{
public:
  void check_array(const StringArray& array, const std::vector<std::string>& truth) const
  {
    ASSERT_EQ(array.size(), truth.size()) << "Incorrect array size";
    ASSERT_EQ(array.empty(), truth.empty()) << "Incorrect emptiness";
    size_t total_length = 0;
    for(const std::string& s : truth) { total_length += s.length(); }
    ASSERT_EQ(array.length(), total_length) << "Incorrect total length of the strings";

    for(size_t i = 0; i < array.size(); i++)
    {
      std::string correct = truth[i];
      ASSERT_EQ(array.str(i), correct) << "Incorrect string " << i;
      EXPECT_EQ(array.length(i), correct.length()) << "Incorrect length for string " << i;
      view_type view = array.view(i);
      std::string from_view(view.first, view.second);
      EXPECT_EQ(from_view, correct) << "Incorrect view of string " << i;
    }
  }
};

TEST_F(StringArrayTest, DefaultEmptyArray)
{
  std::vector<std::string> truth;
  StringArray array;
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, FromEmptyArray)
{
  std::vector<std::string> truth;
  StringArray array(truth);
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, NonEmptyArray)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray array(truth);
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, SerializeEmpty)
{
  std::vector<std::string> truth;
  StringArray original(truth);

  std::string filename = gbwt::TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.serialize(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.deserialize(in);
  in.close();
  ASSERT_EQ(copy, original) << "Serialization changed the empty array";

  gbwt::TempFile::remove(filename);
}

TEST_F(StringArrayTest, CompressEmpty)
{
  std::vector<std::string> truth;
  StringArray original(truth);

  std::string filename = gbwt::TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.compress(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.decompress(in);
  in.close();
  ASSERT_EQ(copy, original) << "Compression changed the empty array";

  gbwt::TempFile::remove(filename);
}

TEST_F(StringArrayTest, SerializeNonEmpty)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray original(truth);

  std::string filename = gbwt::TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.serialize(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.deserialize(in);
  in.close();
  ASSERT_EQ(copy, original) << "Serialization changed the non-empty array";

  gbwt::TempFile::remove(filename);
}

TEST_F(StringArrayTest, CompressNonEmpty)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray original(truth);

  std::string filename = gbwt::TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.compress(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.decompress(in);
  in.close();
  ASSERT_EQ(copy, original) << "Compression changed the non-empty array";

  gbwt::TempFile::remove(filename);
}

TEST_F(StringArrayTest, CompressWithEmptyString)
{
  // Here we test that the compression still works when there is an empty
  // string in the middle and the sd_vector used for the offsets contains
  // duplicate values.
  std::vector<std::string> truth
  {
    "first",
    "second",
    "",
    "fourth"
  };
  StringArray original(truth);

  std::string filename = gbwt::TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.compress(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.decompress(in);
  in.close();
  ASSERT_EQ(copy, original) << "Compression changed the array with an empty string in the middle";

  gbwt::TempFile::remove(filename);
}

//------------------------------------------------------------------------------

} // namespace
