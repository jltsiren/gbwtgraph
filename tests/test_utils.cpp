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
    if(source.uses_translation())
    {
      ASSERT_EQ(source.segment_translation.size(), truth.size()) << "Invalid number of segments";
      for(const translation_type& translation : truth)
      {
        EXPECT_EQ(source.get_translation(translation.first), translation.second) << "Invalid translation for " << translation.first;
        EXPECT_EQ(source.force_translate(translation.first), translation.second) << "Invalid forced translation for " << translation.first;
      }
    }
    else
    {
      for(auto iter = source.nodes.begin(); iter != source.nodes.end(); ++iter)
      {
        std::string segment = std::to_string(iter->first);
        std::pair<nid_t, nid_t> translation(iter->first, iter->first + 1);
        EXPECT_EQ(source.force_translate(segment), translation) << "Invalid forced translation for " << segment;
      }
    }
  }

  void check_inverse_translation(const SequenceSource& source, const std::function<bool(std::pair<nid_t, nid_t>)>& is_present) const
  {
    auto result = source.invert_translation(is_present);
    ASSERT_EQ(result.first.size(), source.segment_translation.size()) << "Invalid number of segments in inverse translation";
    ASSERT_EQ(result.second.size(), size_t(source.next_id)) << "Invalid number of nodes in inverse translation";
    ASSERT_EQ(result.first.size(), result.second.ones()) << "Inconsistent number of segments in inverse translation";
    auto iter = result.second.one_begin();
    nid_t start = iter->second;
    for(size_t i = 0; i < result.first.size(); i++)
    {
      ++iter;
      nid_t limit = iter->second;
      std::string segment = result.first.str(i);
      if(is_present(std::make_pair(start, limit)))
      {
        std::pair<nid_t, nid_t> translation = source.get_translation(segment);
        EXPECT_EQ(start, translation.first) << "Invalid start for segment " << segment;
        EXPECT_EQ(limit, translation.second) << "Invalid limit for segment " << segment;
      }
      else
      {
        EXPECT_TRUE(segment.empty()) << "Got a name for an unused segment from " << start << " to " << limit;
      }
      start = limit;
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
  this->check_inverse_translation(source, [](std::pair<nid_t, nid_t>) -> bool { return true; });
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
  this->check_inverse_translation(source, [](std::pair<nid_t, nid_t>) -> bool { return true; });
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
  this->check_inverse_translation(source, [](std::pair<nid_t, nid_t>) -> bool { return true; });

  // Check an inverse translation without the names for s4 and s6.
  this->check_inverse_translation(source, [](std::pair<nid_t, nid_t> nodes) -> bool
  {
    return (nodes.first < 4 || nodes.first > 6);
  });
}

//------------------------------------------------------------------------------

} // namespace
