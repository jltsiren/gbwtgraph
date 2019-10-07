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

class IndexConstruction : public ::testing::Test
{
public:
  typedef std::map<DefaultMinimizerIndex::key_type, std::set<pos_t>> result_type;

  gbwt::GBWT index;
  SequenceSource source;
  GBWTGraph graph;
  DefaultMinimizerIndex mi;

  IndexConstruction() :
    mi(3, 2)
  {
  }

  void SetUp() override
  {
    this->index = build_gbwt_index();
    build_source(this->source);
    this->graph = GBWTGraph(this->index, this->source);
  }

  void insert_values(const gbwt::vector_type& path, result_type& result) const
  {
    // Convert the path to a string and find the minimizers.
    std::string str;
    for(gbwt::node_type node : path)
    {
      str += this->graph.get_sequence(GBWTGraph::node_to_handle(node));
    }
    std::vector<DefaultMinimizerIndex::minimizer_type> minimizers = this->mi.minimizers(str);

    // Insert the minimizers into the result.
    auto iter = path.begin();
    size_t node_start = 0;
    for(DefaultMinimizerIndex::minimizer_type minimizer : minimizers)
    {
      if(minimizer.empty()) { continue; }
      handle_t handle = GBWTGraph::node_to_handle(*iter);
      size_t node_length = this->graph.get_length(handle);
      while(node_start + node_length <= minimizer.offset)
      {
        node_start += node_length;
        ++iter;
        handle = GBWTGraph::node_to_handle(*iter);
        node_length = this->graph.get_length(handle);
      }
      pos_t pos { this->graph.get_id(handle), this->graph.get_is_reverse(handle), minimizer.offset - node_start };
      if(minimizer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      result[minimizer.key].insert(pos);
    }
  }

  void check_minimizer_index(const result_type& correct_values)
  {
    size_t values = 0;
    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      values += iter->second.size();
    }
    ASSERT_EQ(this->mi.size(), correct_values.size()) << "Wrong number of keys";
    ASSERT_EQ(this->mi.values(), values) << "Wrong number of values";

    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      std::vector<pos_t> result = this->mi.find(get_minimizer(iter->first));
      std::vector<pos_t> correct(iter->second.begin(), iter->second.end());
      EXPECT_EQ(result, correct) << "Wrong positions for key " << iter->first;
    }
  }
};

TEST_F(IndexConstruction, DefaultMinimizerIndex)
{
  // Determine the correct minimizer occurrences.
  std::map<DefaultMinimizerIndex::key_type, std::set<pos_t>> correct_values;
  this->insert_values(alt_path, correct_values);
  this->insert_values(short_path, correct_values);

  // Check that we managed to index them.
  index_haplotypes(this->graph, this->mi);
  this->check_minimizer_index(correct_values);
}

//------------------------------------------------------------------------------

} // namespace
