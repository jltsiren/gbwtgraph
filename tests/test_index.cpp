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

template<class ValueType>
ValueType
create_value(pos_t pos, Payload payload) = delete;

template<>
Position
create_value<Position>(pos_t pos, Payload)
{
  return Position::encode(pos);
}

template<>
PositionPayload
create_value<PositionPayload>(pos_t pos, Payload payload)
{
  return { Position::encode(pos), payload };
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

  template<class ValueType>
  void insert_values(const MinimizerIndex<KeyType, ValueType>& index, const gbwt::vector_type& path, std::map<KeyType, std::set<ValueType>>& result) const
  {
    typedef typename MinimizerIndex<KeyType, ValueType>::minimizer_type minimizer_type;

    // Convert the path to a string and find the minimizers.
    std::string str;
    for(gbwt::node_type node : path)
    {
      str += this->graph.get_sequence(GBWTGraph::node_to_handle(node));
    }
    std::vector<minimizer_type> minimizers = index.minimizers(str);

    // Insert the minimizers into the result.
    auto iter = path.begin();
    size_t node_start = 0;
    for(auto minimizer : minimizers)
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
      result[minimizer.key].emplace(create_value<ValueType>(pos, Payload::create(hash(pos))));
    }
  }

  template<class ValueType>
  void check_minimizer_index(const MinimizerIndex<KeyType, ValueType>& index, const std::map<KeyType, std::set<ValueType>>& correct_values)
  {
    size_t values = 0;
    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      values += iter->second.size();
    }
    ASSERT_EQ(index.size(), correct_values.size()) << "Wrong number of keys";
    ASSERT_EQ(index.number_of_values(), values) << "Wrong number of values";

    for(auto iter = correct_values.begin(); iter != correct_values.end(); ++iter)
    {
      auto values = index.find(get_minimizer<KeyType>(iter->first));
      std::vector<ValueType> result(values.first, values.first + values.second);
      std::vector<ValueType> correct(iter->second.begin(), iter->second.end());
      EXPECT_EQ(result, correct) << "Wrong positions for key " << iter->first;
    }
  }
};

TYPED_TEST_CASE(IndexConstruction, KeyTypes);

TYPED_TEST(IndexConstruction, WithoutPayload)
{
  // Determine the correct minimizer occurrences.
  MinimizerIndex<TypeParam, Position> index(3, 2);
  std::map<TypeParam, std::set<Position>> correct_values;
  this->insert_values(index, alt_path, correct_values);
  this->insert_values(index, short_path, correct_values);

  // Check that we managed to index them.
  index_haplotypes(this->graph, index);
  this->check_minimizer_index(index, correct_values);
}

TYPED_TEST(IndexConstruction, WithPayload)
{
  // Determine the correct minimizer occurrences.
  MinimizerIndex<TypeParam, PositionPayload> index(3, 2);
  std::map<TypeParam, std::set<PositionPayload>> correct_values;
  this->insert_values(index, alt_path, correct_values);
  this->insert_values(index, short_path, correct_values);

  // Check that we managed to index them.
  index_haplotypes(this->graph, index, [](const pos_t& pos) -> Payload
  {
    return Payload::create(hash(pos));
  });
  this->check_minimizer_index(index, correct_values);
}

//------------------------------------------------------------------------------

} // namespace
