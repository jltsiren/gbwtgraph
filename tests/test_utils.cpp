#include <gtest/gtest.h>

#include <gbwtgraph/utils.h>
#include <gbwtgraph/gfa.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

//------------------------------------------------------------------------------

class DigestTest : public ::testing::Test
{
public:
  std::vector<std::pair<std::string, std::string>> test_vectors
  {
    { "", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" },
    { "abc", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad" },
    { "message digest", "f7846f55cf23e14eebeab5b4e1550cad5b509e3348fbc4efa3a1413d393cb650" },
    { "abcdefghijklmnopqrstuvwxyz", "71c480df93d6ae2f1efad1447c66c9525e316218cf51fc8d9ed832f2daf18b73" },
    { "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
      "db4bfcbd4da0cd85a60c3c37d3fbd8805c77f15fc6b1fdfe614ee0a7c8fdb4c0" },
    { "12345678901234567890123456789012345678901234567890123456789012345678901234567890",
      "f371bc4a311f2b009eef952dd83ca80e2b60026c8e935592d0f9c308453c813e" }
  };
};

TEST_F(DigestTest, ByCharacter)
{
  for(const auto& vector : this->test_vectors)
  {
    DigestStream digest_stream;
    std::string input = vector.first;
    for(char c : input)
    {
      digest_stream.put(c);
    }
    std::string digest = digest_stream.finish();
    EXPECT_EQ(digest, vector.second) << "Wrong digest for input \"" << input << "\"";
  }
}

TEST_F(DigestTest, ByBlock)
{
  for(const auto& vector : this->test_vectors)
  {
    DigestStream digest_stream;
    std::string input = vector.first;
    digest_stream.write(input.data(), input.length());
    std::string digest = digest_stream.finish();
    EXPECT_EQ(digest, vector.second) << "Wrong digest for input \"" << input << "\"";
  }
}

//------------------------------------------------------------------------------

class GraphNameTest : public ::testing::Test
{
public:
  // NOTE: Set to true to print descriptions of graph relationships during tests.
  bool print_descriptions = false;

  const std::string GRAPH = "A";
  const std::vector<std::pair<std::string, std::string>> SUBGRAPH
  {
    { "A", "B" },
    { "C", "D" },
    { "D", "E" }
  };
  const std::vector<std::pair<std::string, std::string>> TRANSLATION
  {
    { "B", "C" },
    { "C", "F" }
  };

  GraphName build_manual() const
  {
    GraphName gn(this->GRAPH);
    for(const auto& pair : this->SUBGRAPH)
    {
      gn.add_subgraph(view_type(pair.first), view_type(pair.second));
    }
    for(const auto& pair : this->TRANSLATION)
    {
      gn.add_translation(view_type(pair.first), view_type(pair.second));
    }
    return gn;
  }

  // The first object contains all tags and the second only those relevant to graph names.
  std::pair<gbwt::Tags, gbwt::Tags> build_tags() const
  {
    gbwt::Tags all_tags, graph_name_tags;
    all_tags.set(Version::SOURCE_KEY, Version::SOURCE_VALUE);

    all_tags.set(GraphName::GBZ_NAME_TAG, this->GRAPH);
    graph_name_tags.set(GraphName::GBZ_NAME_TAG, this->GRAPH);

    std::string subgraph_value;
    for(size_t i = 0; i < this->SUBGRAPH.size(); i++)
    {
      if(i > 0) { subgraph_value += ";"; }
      subgraph_value += this->SUBGRAPH[i].first + "," + this->SUBGRAPH[i].second;
    }
    all_tags.set(GraphName::GBZ_SUBGRAPH_TAG, subgraph_value);
    graph_name_tags.set(GraphName::GBZ_SUBGRAPH_TAG, subgraph_value);

    std::string translation_value;
    for(size_t i = 0; i < this->TRANSLATION.size(); i++)
    {
      if(i > 0) { translation_value += ";"; }
      translation_value += this->TRANSLATION[i].first + "," + this->TRANSLATION[i].second;
    }
    all_tags.set(GraphName::GBZ_TRANSLATION_TAG, translation_value);
    graph_name_tags.set(GraphName::GBZ_TRANSLATION_TAG, translation_value);

    return { all_tags, graph_name_tags };
  }

  // The first vector contains all header lines and the second only those relevant to graph names.
  std::pair<std::vector<std::string>, std::vector<std::string>> build_gfa_headers() const
  {
    std::vector<std::string> all_headers;
    all_headers.push_back("H\tVN:Z:1.1");
    all_headers.push_back("H\tNM:Z:" + this->GRAPH);

    for(size_t i = 0; i < this->SUBGRAPH.size(); i++)
    {
      std::string subgraph_value = "H\tSG:Z:" + this->SUBGRAPH[i].first + "," + this->SUBGRAPH[i].second;;
      all_headers.push_back(subgraph_value);
    }

    for(size_t i = 0; i < this->TRANSLATION.size(); i++)
    {
      std::string translation_value = "H\tTL:Z:" + this->TRANSLATION[i].first + "," + this->TRANSLATION[i].second;;
      all_headers.push_back(translation_value);
    }

    std::vector<std::string> graph_name_headers(all_headers.begin() + 1, all_headers.end());
    return { all_headers, graph_name_headers };
  }

  // The first vector contains all header lines and the second only those relevant to graph names.
  std::pair<std::vector<std::string>, std::vector<std::string>> build_gaf_headers() const
  {
    std::vector<std::string> all_headers;
    all_headers.push_back("@HD\tVN:Z:1.0");
    all_headers.push_back("@RN\t" + this->GRAPH);

    for(size_t i = 0; i < this->SUBGRAPH.size(); i++)
    {
      std::string subgraph_value = "@SG\t" + this->SUBGRAPH[i].first + "\t" + this->SUBGRAPH[i].second;;
      all_headers.push_back(subgraph_value);
    }

    for(size_t i = 0; i < this->TRANSLATION.size(); i++)
    {
      std::string translation_value = "@TL\t" + this->TRANSLATION[i].first + "\t" + this->TRANSLATION[i].second;;
      all_headers.push_back(translation_value);
    }

    std::vector<std::string> graph_name_headers(all_headers.begin() + 1, all_headers.end());
    return { all_headers, graph_name_headers };
  }

  size_t expected_lines(size_t steps, bool no_path) const
  {
    size_t result = 2; // from_desc + to_desc
    result += steps; // one line per relationship
    result++; // "With:" line
    result += steps + 1; // one line per graph
    if(no_path) { result++; } // final graph
    return result;
  }

  size_t count_lines(const std::string& description) const
  {
    return std::count(description.begin(), description.end(), '\n');
  }

  void describe_relationships
  (
    const GraphName& from, const GraphName& to,
    const std::string& from_desc, const std::string& to_desc,
    size_t steps, bool no_path
  )
  {
    std::string desc = from.describe_relationship(to, from_desc, to_desc);
    size_t expected_lines = this->expected_lines(steps, no_path);
    size_t actual_lines = this->count_lines(desc);
    EXPECT_EQ(actual_lines, expected_lines)
      << "Wrong number of lines in description of relationships between "
      << from.name() << " and " << to.name();
    if(this->print_descriptions)
    {
      std::cout << desc << std::endl;
    }
  }
};

TEST_F(GraphNameTest, Empty)
{
  gbwt::Tags empty_tags;
  std::vector<std::string> empty_headers;

  GraphName manual;
  EXPECT_EQ(manual.name(), "") << "Default constructor sets non-empty name";
  EXPECT_FALSE(manual.has_name()) << "Default constructor sets a name";

  GraphName from_tags(empty_tags);
  EXPECT_EQ(from_tags.name(), "") << "Non-empty name from empty tags";
  EXPECT_FALSE(from_tags.has_name()) << "Importing empty tags sets a name";
  EXPECT_FALSE(from_tags.same(manual)) << "Empty names refer to the same graph (from tags)";
  gbwt::Tags to_tags;
  manual.set_tags(to_tags);
  EXPECT_EQ(to_tags, empty_tags) << "Tags set from default GraphName are not empty";
  this->describe_relationships(from_tags, manual, "tags", "manually built", 0, true);

  GraphName from_headers(empty_headers);
  EXPECT_EQ(from_headers.name(), "") << "Non-empty name from empty headers";
  EXPECT_FALSE(from_headers.has_name()) << "Importing empty headers sets a name";
  EXPECT_FALSE(from_headers.same(manual)) << "Empty names refer to the same graph (from headers)";
  std::vector<std::string> to_headers = manual.gfa_header_lines();
  EXPECT_EQ(to_headers, empty_headers) << "GFA headers from default GraphName are not empty";
  to_headers = manual.gaf_header_lines();
  EXPECT_EQ(to_headers, empty_headers) << "GAF headers from default GraphName are not empty";
};

TEST_F(GraphNameTest, Tags)
{
  gbwt::Tags all_tags;
  gbwt::Tags graph_name_tags;
  std::tie(all_tags, graph_name_tags) = this->build_tags();
  GraphName from_tags(all_tags);
  EXPECT_EQ(from_tags.name(), this->GRAPH) << "Wrong graph name from tags";
  EXPECT_TRUE(from_tags.has_name()) << "Importing non-empty tags does not set a name";

  GraphName manual = this->build_manual();
  EXPECT_TRUE(from_tags.same(manual)) << "GraphName from tags is not same as manually built";
  EXPECT_EQ(from_tags, manual) << "GraphName from tags is not equal to manually built";

  gbwt::Tags to_tags;
  manual.set_tags(to_tags);
  EXPECT_EQ(to_tags, graph_name_tags) << "Tags set from manually built GraphName do not match original tags";

  GraphName empty;
  empty.set_tags(to_tags);
  EXPECT_FALSE(to_tags.contains(GraphName::GBZ_NAME_TAG)) << "Name tag not removed from tags set from empty GraphName";
  EXPECT_FALSE(to_tags.contains(GraphName::GBZ_SUBGRAPH_TAG)) << "Subgraph tag not removed from tags set from empty GraphName";
  EXPECT_FALSE(to_tags.contains(GraphName::GBZ_TRANSLATION_TAG)) << "Translation tag not removed from tags set from empty GraphName";
}

TEST_F(GraphNameTest, GFAHeaders)
{
  std::vector<std::string> all_headers;
  std::vector<std::string> graph_name_headers;
  std::tie(all_headers, graph_name_headers) = this->build_gfa_headers();
  GraphName from_headers(all_headers);
  EXPECT_EQ(from_headers.name(), this->GRAPH) << "Wrong graph name from GFA headers";
  EXPECT_TRUE(from_headers.has_name()) << "Importing non-empty GFA headers does not set a name";

  GraphName manual = this->build_manual();
  EXPECT_TRUE(from_headers.same(manual)) << "GraphName from GFA headers is not same as manually built";
  EXPECT_EQ(from_headers, manual) << "GraphName from GFA headers is not equal to manually built";

  std::vector<std::string> to_headers = manual.gfa_header_lines();
  EXPECT_EQ(to_headers, graph_name_headers) << "GFA headers from manually built GraphName do not match original headers";

  std::vector<std::string> weird_headers
  {
    "H\tVN:Z:1.1",
    "H\tNM:Z:A:B", // We allow colons in string values.
    "H\txy:Z:C:D", // Including unknown typed fields.
  };
  GraphName from_weird(weird_headers);
  EXPECT_EQ(from_weird.name(), "A:B") << "Wrong graph name from weird GFA headers";
}

TEST_F(GraphNameTest, GAFHeaders)
{
  std::vector<std::string> all_headers;
  std::vector<std::string> graph_name_headers;
  std::tie(all_headers, graph_name_headers) = this->build_gaf_headers();
  GraphName from_headers(all_headers);
  EXPECT_EQ(from_headers.name(), this->GRAPH) << "Wrong graph name from GAF headers";
  EXPECT_TRUE(from_headers.has_name()) << "Importing non-empty GAF headers does not set a name";

  GraphName manual = this->build_manual();
  EXPECT_TRUE(from_headers.same(manual)) << "GraphName from GAF headers is not same as manually built";
  EXPECT_EQ(from_headers, manual) << "GraphName from GAF headers is not equal to manually built";

  std::vector<std::string> to_headers = manual.gaf_header_lines();
  EXPECT_EQ(to_headers, graph_name_headers) << "GAF headers from manually built GraphName do not match original headers";
}

TEST_F(GraphNameTest, SubgraphPath)
{
  GraphName a = this->build_manual();
  GraphName a_empty(this->GRAPH);

  // Same graph.
  EXPECT_TRUE(a.subgraph_of(a_empty)) << "Graph is not subgraph of itself (relationships in subgraph)";
  EXPECT_TRUE(a_empty.subgraph_of(a)) << "Graph is not subgraph of itself (relationships in supergraph)";
  EXPECT_TRUE(a_empty.subgraph_of(a_empty)) << "Graph is not subgraph of itself (no relationships)";
  this->describe_relationships(a, a_empty, "the original graph", "the same graph", 0, false);

  // Single step.
  {
    GraphName b_empty("B");
    GraphName b = b_empty;
    b.add_relationships(a);
    EXPECT_TRUE(a.subgraph_of(b_empty)) << "A is not subgraph of B (relationships in subgraph)";
    EXPECT_TRUE(a_empty.subgraph_of(b)) << "A is not subgraph of B (relationships in supergraph)";
    EXPECT_FALSE(b.subgraph_of(a)) << "B is subgraph of A";
    this->describe_relationships(a, b, "subgraph (first)", "supergraph (second)", 1, false);
    this->describe_relationships(b, a, "supergraph (first)", "subgraph (second)", 1, false);
  }

  // No path.
  GraphName c_empty("C");
  GraphName c = c_empty;
  c.add_relationships(a);
  EXPECT_FALSE(a.subgraph_of(c)) << "A is subgraph of C";
  EXPECT_FALSE(c.subgraph_of(a)) << "C is subgraph of A";
  // Skip description, as there is a path with a translation.

  // Multiple steps.
  {
    GraphName e_empty("E");
    GraphName e = e_empty;
    e.add_relationships(a);
    EXPECT_TRUE(c.subgraph_of(e_empty)) << "C is not subgraph of E (relationships in subgraph)";
    EXPECT_TRUE(c_empty.subgraph_of(e)) << "C is not subgraph of E (relationships in supergraph)";
    EXPECT_FALSE(e.subgraph_of(c)) << "E is subgraph of C";
    this->describe_relationships(c, e, "subgraph", "supergraph", 2, false);
  }

  // Relationships split between graphs.
  {
    GraphName from("from");
    from.add_subgraph(view_type("from"), view_type("middle"));
    GraphName to("to");
    to.add_subgraph(view_type("middle"), view_type("to"));
    EXPECT_TRUE(from.subgraph_of(to)) << "from is not subgraph of to (relationships split)";
    EXPECT_FALSE(to.subgraph_of(from)) << "to is subgraph of from";
    this->describe_relationships(from, to, "subgraph", "supergraph", 2, false);
  }
}

TEST_F(GraphNameTest, TranslationPath)
{
  GraphName a = this->build_manual();
  GraphName a_empty(this->GRAPH);

  // Same graph.
  EXPECT_TRUE(a.translates_to(a_empty)) << "Graph does not translate to itself (relationships in translation)";
  EXPECT_TRUE(a_empty.translates_to(a)) << "Graph does not translate to itself (no relationships)";
  EXPECT_TRUE(a_empty.translates_to(a_empty)) << "Graph does not translate to itself (no relationships)";
  this->describe_relationships(a, a_empty, "the original graph", "the same graph", 0, false);

  // Single subgraph step.
  GraphName b_empty("B");
  GraphName b = b_empty;
  b.add_relationships(a);
  EXPECT_TRUE(a.subgraph_of(b_empty)) << "A is not subgraph of B (relationships in subgraph)";
  EXPECT_TRUE(a_empty.subgraph_of(b)) << "A is not subgraph of B (relationships in supergraph)";
  EXPECT_FALSE(b.subgraph_of(a)) << "B is subgraph of A";
  this->describe_relationships(a, b, "subgraph", "supergraph", 1, false);

  // Single translation step.
  GraphName c_empty("C");
  GraphName c = c_empty;
  c.add_relationships(a);
  EXPECT_TRUE(b.translates_to(c_empty)) << "B does not translate to C (relationships in source graph)";
  EXPECT_TRUE(b_empty.translates_to(c)) << "B does not translate to C (relationships in target graph)";
  EXPECT_FALSE(c.translates_to(b)) << "C translates to B";
  this->describe_relationships(b, c, "source graph", "target graph", 1, false);

  // Mixed steps.
  EXPECT_TRUE(a.translates_to(c_empty)) << "A does not translate to C (relationships in source graph)";
  EXPECT_TRUE(a_empty.translates_to(c)) << "A does not translate to C (relationships in target graph)";
  EXPECT_FALSE(c.translates_to(a)) << "C translates to A";
  this->describe_relationships(a, c, "source graph", "target graph", 2, false);

  // Multiple subgraph steps.
  {
    GraphName e_empty("E");
    GraphName e = e_empty;
    e.add_relationships(a);
    EXPECT_TRUE(c.translates_to(e_empty)) << "C does not translate to E (relationships in source graph)";
    EXPECT_TRUE(c_empty.translates_to(e)) << "C does not translate to E (relationships in target graph)";
    EXPECT_FALSE(e.translates_to(c)) << "E translates to C";
    this->describe_relationships(c, e, "source graph", "target graph", 2, false);
  }

  // Multiple translation steps.
  {
    GraphName f_empty("F");
    GraphName f = f_empty;
    f.add_relationships(a);
    EXPECT_TRUE(b.translates_to(f_empty)) << "B does not translate to F (relationships in source graph)";
    EXPECT_TRUE(b_empty.translates_to(f)) << "B does not translate to F (relationships in target graph)";
    EXPECT_FALSE(f.translates_to(b)) << "F translates to B";
    this->describe_relationships(b, f, "source graph", "target graph",  2, false);
  }

  // Relationships split between graphs.
  {
    GraphName from("from");
    from.add_translation(view_type("from"), view_type("middle"));
    GraphName to("to");
    to.add_translation(view_type("middle"), view_type("to"));
    EXPECT_TRUE(from.translates_to(to)) << "from does not translate to to (relationships split)";
    EXPECT_FALSE(to.translates_to(from)) << "to translates to from";
    this->describe_relationships(from, to, "source graph", "target graph", 2, false);
  }
}

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

class MetadataBuilderTest : public ::testing::Test
{
public:
  void check_empty(const MetadataBuilder& builder) const
  {
    gbwt::Metadata metadata = builder.get_metadata();
    ASSERT_EQ(metadata.samples(), gbwt::size_type(0)) << "Invalid number of samples";
    ASSERT_EQ(metadata.contigs(), gbwt::size_type(0)) << "Invalid number of contigs";
    ASSERT_EQ(metadata.paths(), gbwt::size_type(0)) << "Invalid number of paths";
  }

  void create_example(
    std::vector<std::string>& samples,
    std::vector<std::string>& contigs,
    std::vector<gbwt::FullPathName>& paths,
    bool generic_reference) const
  {
    std::string reference_sample = (generic_reference ? REFERENCE_PATH_SAMPLE_NAME : "GRCh38");
    size_t reference_haplotype = (generic_reference ? GBWTGraph::NO_PHASE : 0);
    samples.push_back(reference_sample);
    samples.push_back("HG002");
    samples.push_back("HG003");

    contigs.push_back("chr1");
    contigs.push_back("chr2");

    paths.push_back({ reference_sample, "chr1", reference_haplotype, 0 });
    paths.push_back({ reference_sample, "chr2", reference_haplotype, 0 });
    paths.push_back({ "HG002", "chr1", 1, 0 });
    paths.push_back({ "HG002", "chr1", 2, 0 });
    paths.push_back({ "HG002", "chr2", 1, 0 });
    paths.push_back({ "HG002", "chr2", 2, 0 });
    paths.push_back({ "HG003", "chr1", 1, 0 });
    paths.push_back({ "HG003", "chr1", 2, 0 });
    paths.push_back({ "HG003", "chr2", 1, 0 });
    paths.push_back({ "HG003", "chr2", 2, 0 });
  }

  void add_hg004(
    std::vector<std::string>& samples,
    std::vector<gbwt::FullPathName>& paths) const
  {
    samples.push_back("HG004");
    paths.push_back({ "HG004", "chr1", 1, 0 });
    paths.push_back({ "HG004", "chr1", 2, 0 });
    paths.push_back({ "HG004", "chr2", 1, 0 });
    paths.push_back({ "HG004", "chr2", 2, 0 });
  }

  size_t get_job(const gbwt::FullPathName& path) const
  {
    if(path.contig_name == "chr1") { return 0; }
    if(path.contig_name == "chr2") { return 1; }
    return 0;
  }

  void add_haplotypes(MetadataBuilder& builder, const std::vector<gbwt::FullPathName>& paths, size_t from, bool assign_job)
  {
    for(size_t i = from; i < paths.size(); i++)
    {
      const gbwt::FullPathName& path = paths[i];
      size_t job = (assign_job ? get_job(path) : 0);
      if(path.sample_name == REFERENCE_PATH_SAMPLE_NAME)
      {
        builder.add_generic_path(path.contig_name, job);
      }
      else
      {
        builder.add_haplotype(path.sample_name, path.contig_name, path.haplotype, path.offset, job);
      }
    }
  }

  void add_walks(MetadataBuilder& builder, const std::vector<gbwt::FullPathName>& paths)
  {
    for(const gbwt::FullPathName& path : paths)
    {
      if(path.sample_name == REFERENCE_PATH_SAMPLE_NAME)
      {
        builder.add_generic_path(path.contig_name);
      }
      else
      {
        std::string haplotype = std::to_string(path.haplotype);
        std::string start = std::to_string(path.offset);
        builder.add_walk(path.sample_name, haplotype, path.contig_name, start);
      }
    }
  }

  void add_walks_no_interval(MetadataBuilder& builder, const std::vector<gbwt::FullPathName>& paths)
  {
    std::string no_interval = "*";
    for(const gbwt::FullPathName& path : paths)
    {
      if(path.sample_name == REFERENCE_PATH_SAMPLE_NAME)
      {
        builder.add_generic_path(path.contig_name);
      }
      else
      {
        std::string haplotype = std::to_string(path.haplotype);
        builder.add_walk(path.sample_name, haplotype, path.contig_name, no_interval);
      }
    }
  }

  void add_named_paths(MetadataBuilder& builder, const std::vector<gbwt::FullPathName>& paths)
  {
    for(const gbwt::FullPathName& path : paths)
    {
      std::string name;
      if(path.sample_name == REFERENCE_PATH_SAMPLE_NAME)
      {
        name = path.contig_name;
      }
      else
      {
        name = path.sample_name + "#" + std::to_string(path.haplotype) + "#" + path.contig_name;
      }
      builder.add_path(name);
    }
  }

  void check_metadata(
    const gbwt::Metadata& metadata,
    const std::vector<std::string>& samples,
    const std::vector<std::string>& contigs,
    const std::vector<gbwt::FullPathName>& paths) const
  {
    ASSERT_EQ(metadata.samples(), samples.size()) << "Invalid number of samples";
    for(size_t i = 0; i < samples.size(); i++)
    {
      EXPECT_EQ(metadata.sample(i), samples[i]) << "Invalid sample name at index " << i;
    }

    ASSERT_EQ(metadata.contigs(), contigs.size()) << "Invalid number of contigs";
    for(size_t i = 0; i < contigs.size(); i++)
    {
      EXPECT_EQ(metadata.contig(i), contigs[i]) << "Invalid contig name at index " << i;
    }

    ASSERT_EQ(metadata.paths(), paths.size()) << "Invalid number of paths";
    for(size_t i = 0; i < paths.size(); i++)
    {
      gbwt::PathName path = metadata.path(i);
      EXPECT_EQ(metadata.sample(path.sample), paths[i].sample_name) << "Invalid sample name for path " << i;
      EXPECT_EQ(metadata.contig(path.contig), paths[i].contig_name) << "Invalid contig name for path " << i;
      EXPECT_EQ(path.phase, paths[i].haplotype) << "Invalid haplotype for path " << i;
      EXPECT_EQ(path.count, paths[i].offset) << "Invalid offset for path " << i;
    }
  }
};

TEST_F(MetadataBuilderTest, Empty)
{
  this->check_empty(MetadataBuilder());
}

TEST_F(MetadataBuilderTest, GenericPathsAndHaplotypes)
{
  std::vector<std::string> samples, contigs;
  std::vector<gbwt::FullPathName> paths;
  this->create_example(samples, contigs, paths, true);

  MetadataBuilder builder;
  this->add_haplotypes(builder, paths, 0, false);
  this->check_metadata(builder.get_metadata(), samples, contigs, paths);
}

TEST_F(MetadataBuilderTest, GFAPathsAndWalks)
{
  std::vector<std::string> samples, contigs;
  std::vector<gbwt::FullPathName> paths;
  this->create_example(samples, contigs, paths, true);

  MetadataBuilder builder;
  this->add_walks(builder, paths);
  this->check_metadata(builder.get_metadata(), samples, contigs, paths);
}

TEST_F(MetadataBuilderTest, GFAWalksNoInterval)
{
  std::vector<std::string> samples, contigs;
  std::vector<gbwt::FullPathName> paths;
  this->create_example(samples, contigs, paths, true);

  MetadataBuilder builder;
  this->add_walks_no_interval(builder, paths);
  this->check_metadata(builder.get_metadata(), samples, contigs, paths);
}

TEST_F(MetadataBuilderTest, PanSN)
{
  std::vector<std::string> samples, contigs;
  std::vector<gbwt::FullPathName> paths;
  this->create_example(samples, contigs, paths, false);

  MetadataBuilder builder(
    GFAParsingParameters::PAN_SN_REGEX,
    GFAParsingParameters::PAN_SN_FIELDS,
    GFAParsingParameters::PAN_SN_SENSE
  );
  this->add_haplotypes(builder, paths, 0, false);
  this->check_metadata(builder.get_metadata(), samples, contigs, paths);
}

TEST_F(MetadataBuilderTest, Clear)
{
  std::vector<std::string> samples, contigs;
  std::vector<gbwt::FullPathName> paths;
  this->create_example(samples, contigs, paths, true);

  MetadataBuilder builder;
  this->add_haplotypes(builder, paths, 0, false);
  builder.clear();
  this->check_empty(builder);
}

TEST_F(MetadataBuilderTest, MultipleFormats)
{
  std::vector<std::string> samples, contigs;
  std::vector<gbwt::FullPathName> paths;
  this->create_example(samples, contigs, paths, true);

  MetadataBuilder builder;
  builder.add_path_name_format(
    GFAParsingParameters::PAN_SN_REGEX,
    GFAParsingParameters::PAN_SN_FIELDS,
    GFAParsingParameters::PAN_SN_SENSE
  );
  builder.add_path_name_format(".*", "C", PathSense::GENERIC);

  this->add_named_paths(builder, paths);
  this->check_metadata(builder.get_metadata(), samples, contigs, paths);
}

TEST_F(MetadataBuilderTest, FromMetadata)
{
  std::vector<std::string> samples, contigs;
  std::vector<gbwt::FullPathName> paths;
  this->create_example(samples, contigs, paths, true);
  size_t old_paths = paths.size();

  MetadataBuilder builder;
  this->add_haplotypes(builder, paths, 0, false);
  gbwt::Metadata metadata = builder.get_metadata();

  MetadataBuilder new_builder(metadata);
  this->add_hg004(samples, paths);
  this->add_haplotypes(new_builder, paths, old_paths, false);
  this->check_metadata(new_builder.get_metadata(), samples, contigs, paths);
}

TEST_F(MetadataBuilderTest, MultipleJobs)
{
  std::vector<std::string> samples, contigs;
  std::vector<gbwt::FullPathName> paths;
  this->create_example(samples, contigs, paths, true);

  MetadataBuilder builder;
  this->add_haplotypes(builder, paths, 0, true);

  std::vector<gbwt::FullPathName> reordered_paths;
  for(size_t job = 0; job < contigs.size(); job++)
  {
    for(size_t i = 0; i < paths.size(); i++)
    {
      if(this->get_job(paths[i]) == job)
      {
        reordered_paths.push_back(paths[i]);
      }
    }
  }
  this->check_metadata(builder.get_metadata(), samples, contigs, reordered_paths);
}

//------------------------------------------------------------------------------

} // namespace
