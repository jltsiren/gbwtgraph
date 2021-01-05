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
      handle_t handle = source.get_handle(node.first, false);
      EXPECT_EQ(source.get_length(handle), node.second.length()) << "Invalid sequence length for node " << node.first;
      EXPECT_EQ(source.get_sequence(handle), node.second) << "Invalid sequence for node " << node.first;
      view_type view = source.get_sequence_view(handle);
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
};

TEST_F(SourceTest, EmptySource)
{
  SequenceSource source;
  std::vector<node_type> nodes;
  std::vector<translation_type> translation;

  this->check_nodes(source, nodes);
  this->check_translation(source, translation);
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
}

//------------------------------------------------------------------------------

} // namespace
