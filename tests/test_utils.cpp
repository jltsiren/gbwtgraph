#include <gtest/gtest.h>

#include <cstring>

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
      gn.add_subgraph(pair.first, pair.second);
    }
    for(const auto& pair : this->TRANSLATION)
    {
      gn.add_translation(pair.first, pair.second);
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
    from.add_subgraph("from", "middle");
    GraphName to("to");
    to.add_subgraph("middle", "to");
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
    from.add_translation("from", "middle");
    GraphName to("to");
    to.add_translation("middle", "to");
    EXPECT_TRUE(from.translates_to(to)) << "from does not translate to to (relationships split)";
    EXPECT_FALSE(to.translates_to(from)) << "to translates to from";
    this->describe_relationships(from, to, "source graph", "target graph", 2, false);
  }
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
    std::string reference_sample = (generic_reference ? GENERIC_PATH_SAMPLE_NAME : "GRCh38");
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
      if(path.sample_name == GENERIC_PATH_SAMPLE_NAME)
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
      if(path.sample_name == GENERIC_PATH_SAMPLE_NAME)
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
      if(path.sample_name == GENERIC_PATH_SAMPLE_NAME)
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
      if(path.sample_name == GENERIC_PATH_SAMPLE_NAME)
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

class DedupStringArrayTest : public ::testing::Test
{
public:
  constexpr static size_t LINEARIZED = std::numeric_limits<size_t>::max();

  // An original string: either an independent string or a substring of a superstring.
  struct Source
  {
    std::string seq; // Used if `index == LINEARIZED`.
    size_t index, offset, length;
  };

  std::vector<std::string> superstrings { "GATTACACATTAG", "CCCGGGAAATTT" };

  std::vector<Source> sources
  {
    { "AC", LINEARIZED, 0, 0 },
    { "", 0, 0, 5 }, // GATTA
    { "", 0, 4, 4 }, // ACAC (overlaps the previous one)
    { "", LINEARIZED, 0, 0 }, // An empty string in the linearized prefix.
    { "TTTT", LINEARIZED, 0, 0 },
    { "", 1, 3, 6 }, // GGGAAA
    { "", 0, 9, 4 }, // TTAG
    { "", 1, 5, 0 }, // An empty string pointing to a superstring.
    { "", 0, 0, 5 }, // GATTA again; superstring offset 0 after the last linearized string.
    { "", 1, 0, 12 }, // An entire superstring.
  };

  std::string source_seq(const Source& source) const
  {
    if(source.index == LINEARIZED) { return source.seq; }
    return this->superstrings[source.index].substr(source.offset, source.length);
  }

  DedupStringArray::StringSource string_source(size_t i) const
  {
    const Source& source = this->sources[i];
    if(source.index == LINEARIZED)
    {
      return { std::string_view(source.seq), this->superstrings.size() };
    }
    std::string_view superstring(this->superstrings[source.index]);
    return { superstring.substr(source.offset, source.length), source.index };
  }

  DedupStringArray build(bool reverse_complements) const
  {
    return DedupStringArray
    (
      this->sources.size(), reverse_complements, this->superstrings,
      [this](size_t i) -> DedupStringArray::StringSource { return this->string_source(i); }
    );
  }

  // Expected strings in the array.
  std::vector<std::string> truth(bool reverse_complements) const
  {
    std::vector<std::string> result;
    for(const Source& source : this->sources)
    {
      std::string seq = this->source_seq(source);
      result.push_back(seq);
      if(reverse_complements) { result.push_back(reverse_complement(seq)); }
    }
    return result;
  }

  // Expected linearization status for the strings in the array.
  // Empty strings are always in the linearized prefix.
  std::vector<bool> expected_linearized(bool reverse_complements) const
  {
    std::vector<bool> result;
    for(const Source& source : this->sources)
    {
      bool linearized = (source.index == LINEARIZED || source.length == 0);
      result.push_back(linearized);
      if(reverse_complements) { result.push_back(linearized); }
    }
    return result;
  }

  std::string superstring_suffix(bool reverse_complements) const
  {
    std::string result;
    for(const std::string& superstring : this->superstrings) { result += superstring; }
    if(reverse_complements) { result += reverse_complement(result); }
    return result;
  }

  void check_array
  (
    const DedupStringArray& array,
    const std::vector<std::string>& truth, const std::vector<bool>& linearized,
    bool reverse_complements, bool linear,
    const std::string& test_case
  ) const
  {
    ASSERT_EQ(array.size(), truth.size()) << test_case << ": Wrong number of strings";
    ASSERT_EQ(linearized.size(), truth.size()) << test_case << ": Invalid test data";
    EXPECT_EQ(array.empty(), truth.empty()) << test_case << ": Wrong emptiness";
    EXPECT_EQ(array.contains_reverse_complements(), reverse_complements) << test_case << ": Wrong reverse complement flag";
    EXPECT_EQ(array.linear(), linear) << test_case << ": Wrong linearity";

    size_t total_length = 0;
    std::string prefix;
    for(size_t i = 0; i < truth.size(); i++)
    {
      EXPECT_EQ(array.length(i), truth[i].length()) << test_case << ": Wrong length for string " << i;
      EXPECT_EQ(array.view(i), truth[i]) << test_case << ": Wrong view for string " << i;
      EXPECT_EQ(array.str(i), truth[i]) << test_case << ": Wrong string " << i;
      EXPECT_EQ(array.linearized(i), linearized[i]) << test_case << ": Wrong linearization status for string " << i;
      total_length += truth[i].length();
      if(linearized[i]) { prefix += truth[i]; }
    }
    EXPECT_EQ(array.length(), total_length) << test_case << ": Wrong total length";

    ASSERT_EQ(array.offsets.size(), array.size() + 1) << test_case << ": Wrong number of offsets";
    EXPECT_EQ(array.data.size(), array.data_array_length) << test_case << ": Wrong data array length";
    EXPECT_EQ(array.offsets[array.size()], array.data.size()) << test_case << ": Wrong sentinel";
    ASSERT_EQ(array.linearized_length, prefix.length()) << test_case << ": Wrong linearized length";
    ASSERT_LE(array.linearized_length, array.data.size()) << test_case << ": Linearized prefix exceeds data";
    EXPECT_EQ(std::string(array.data.data(), array.linearized_length), prefix) << test_case << ": Wrong linearized prefix";
  }

  void check_serialization(const DedupStringArray& original, const std::string& test_case) const
  {
    std::stringstream ss;
    original.simple_sds_serialize(ss);
    size_t bytes = ss.str().size();
    EXPECT_EQ(original.simple_sds_size() * sizeof(std::uint64_t), bytes) << test_case << ": Wrong serialized size";

    DedupStringArray copy;
    ASSERT_NO_THROW(copy.simple_sds_load(ss)) << test_case << ": Loading failed";
    EXPECT_TRUE(copy.equivalent(original)) << test_case << ": Loaded array is not equivalent";
    EXPECT_TRUE(copy.identical(original)) << test_case << ": Loaded array is not identical";
  }

  void expect_invalid(const DedupStringArray& array, const std::string& test_case) const
  {
    std::stringstream ss;
    array.simple_sds_serialize(ss);
    DedupStringArray copy;
    EXPECT_THROW(copy.simple_sds_load(ss), sdsl::simple_sds::InvalidData) << test_case << ": Invalid data was accepted";
  }

  // Serializes a valid array, replaces header field `field` with `value`, and
  // expects that loading fails.
  // Header fields: 0 = total length, 1 = linearized length, 2 = data array length.
  void expect_invalid_header(const DedupStringArray& array, size_t field, size_t value, const std::string& test_case) const
  {
    std::stringstream ss;
    array.simple_sds_serialize(ss);
    std::string serialized = ss.str();
    ASSERT_GE(serialized.size(), (field + 1) * sizeof(std::uint64_t)) << test_case << ": Serialized data is too short";
    std::uint64_t word = value;
    std::memcpy(serialized.data() + field * sizeof(std::uint64_t), &word, sizeof(word));

    std::istringstream in(serialized);
    DedupStringArray copy;
    EXPECT_THROW(copy.simple_sds_load(in), sdsl::simple_sds::InvalidData) << test_case << ": Invalid data was accepted";
  }

  void superstring_test(bool reverse_complements) const
  {
    std::string test_case = (reverse_complements ? "With reverse complements" : "Without reverse complements");
    DedupStringArray array = this->build(reverse_complements);
    std::vector<std::string> truth = this->truth(reverse_complements);
    this->check_array(array, truth, this->expected_linearized(reverse_complements), reverse_complements, false, test_case);

    // The data array should end with the superstrings.
    std::string suffix = this->superstring_suffix(reverse_complements);
    ASSERT_EQ(array.data.size(), array.linearized_length + suffix.length()) << test_case << ": Wrong data array length";
    EXPECT_EQ(std::string(array.data.data() + array.linearized_length, suffix.length()), suffix) << test_case << ": Wrong superstrings";

    gbwt::StringArray expected(truth);
    EXPECT_EQ(array.to_string_array(), expected) << test_case << ": Wrong conversion to StringArray";
    this->check_serialization(array, test_case);

    DedupStringArray linearized = array;
    linearized.linearize();
    std::vector<bool> all_linearized(truth.size(), true);
    this->check_array(linearized, truth, all_linearized, reverse_complements, true, test_case + " (linearized)");
    EXPECT_EQ(linearized.data_array_length, linearized.length()) << test_case << ": Wrong data array length after linearization";
    EXPECT_TRUE(linearized.equivalent(array)) << test_case << ": Linearized array is not equivalent";
    EXPECT_EQ(linearized.to_string_array(), expected) << test_case << ": Wrong conversion to StringArray after linearization";
    this->check_serialization(linearized, test_case + " (linearized)");
  }
};

TEST_F(DedupStringArrayTest, Empty)
{
  DedupStringArray array;
  this->check_array(array, {}, {}, false, true, "Empty");
  EXPECT_EQ(array.to_string_array(), gbwt::StringArray()) << "Wrong conversion to StringArray";
  this->check_serialization(array, "Empty");

  DedupStringArray linearized = array;
  linearized.linearize();
  EXPECT_TRUE(linearized.identical(array)) << "Linearization changed an empty array";

  // Empty arrays built in other ways.
  for(bool rc : { false, true })
  {
    std::string test_case = std::string("Empty") + (rc ? " with reverse complements" : "");
    DedupStringArray from_string_array(gbwt::StringArray(), rc);
    this->check_array(from_string_array, {}, {}, rc, true, test_case + " (from StringArray)");
    this->check_serialization(from_string_array, test_case + " (from StringArray)");

    DedupStringArray from_sources(0, rc, {}, [](size_t) -> DedupStringArray::StringSource { return { std::string_view(), 0 }; });
    this->check_array(from_sources, {}, {}, rc, true, test_case + " (from sources)");
    this->check_serialization(from_sources, test_case + " (from sources)");
  }
}

TEST_F(DedupStringArrayTest, FromStringArray)
{
  for(bool rc : { false, true })
  {
    std::string test_case = (rc ? "With reverse complements" : "Without reverse complements");
    std::vector<std::string> truth = this->truth(rc);
    std::vector<bool> all_linearized(truth.size(), true);
    gbwt::StringArray source(truth);

    DedupStringArray array(source, rc);
    this->check_array(array, truth, all_linearized, rc, true, test_case);
    EXPECT_EQ(array.data_array_length, array.length()) << test_case << ": Wrong data array length";
    EXPECT_EQ(array.to_string_array(), source) << test_case << ": Wrong conversion to StringArray";
    this->check_serialization(array, test_case);

    gbwt::StringArray to_move = source;
    DedupStringArray moved(std::move(to_move), rc);
    EXPECT_TRUE(moved.identical(array)) << test_case << ": Move constructor produced a different array";
    EXPECT_TRUE(to_move.empty()) << test_case << ": Source was not reset";
    EXPECT_EQ(to_move.length(), size_t(0)) << test_case << ": Source still contains data";
  }
}

TEST_F(DedupStringArrayTest, Superstrings)
{
  this->superstring_test(false);
}

TEST_F(DedupStringArrayTest, SuperstringsWithReverseComplements)
{
  this->superstring_test(true);
}

TEST_F(DedupStringArrayTest, EmptyStringWithoutView)
{
  // An empty string referring to a superstring with a default-constructed view
  // that does not point to the superstring.
  for(bool rc : { false, true })
  {
    std::string test_case = (rc ? "With reverse complements" : "Without reverse complements");
    this->sources = { { "AC", LINEARIZED, 0, 0 }, { "", 0, 0, 0 }, { "", 0, 0, 5 } };
    DedupStringArray array
    (
      this->sources.size(), rc, this->superstrings,
      [this](size_t i) -> DedupStringArray::StringSource
      {
        if(i == 1) { return { std::string_view(), 0 }; }
        return this->string_source(i);
      }
    );
    this->check_array(array, this->truth(rc), this->expected_linearized(rc), rc, false, test_case);
    this->check_serialization(array, test_case);
  }
}

TEST_F(DedupStringArrayTest, InvalidSources)
{
  std::string other = "GATTACA";
  std::string_view first(this->superstrings[0]), second(this->superstrings[1]);
  std::vector<std::pair<DedupStringArray::StringSource, std::string>> invalid
  {
    { { std::string_view(other), 0 }, "View outside the superstrings" },
    { { second.substr(0, 4), 0 }, "View in the wrong superstring" },
    { { std::string_view(first.data(), first.size() + 1), 0 }, "View extends past the end of the superstring" },
  };

  for(bool rc : { false, true })
  {
    for(const auto& test : invalid)
    {
      std::string test_case = test.second + (rc ? " (with reverse complements)" : "");
      auto sequence = [&](size_t i) -> DedupStringArray::StringSource
      {
        if(i == 1) { return test.first; }
        return this->string_source(i);
      };
      EXPECT_THROW(DedupStringArray(3, rc, this->superstrings, sequence), std::invalid_argument) << test_case << ": Invalid source was accepted";
    }
  }
}

TEST_F(DedupStringArrayTest, UnusedSuperstrings)
{
  for(bool rc : { false, true })
  {
    std::string test_case = (rc ? "With reverse complements" : "Without reverse complements");
    this->sources = { { "AC", LINEARIZED, 0, 0 }, { "", LINEARIZED, 0, 0 }, { "GATTA", LINEARIZED, 0, 0 } };
    DedupStringArray array = this->build(rc);
    std::vector<std::string> truth = this->truth(rc);
    std::vector<bool> all_linearized(truth.size(), true);

    // All strings are in the linearized prefix, but the data array still
    // contains the superstrings. Hence the array is not linear.
    this->check_array(array, truth, all_linearized, rc, false, test_case);
    EXPECT_EQ(array.to_string_array(), gbwt::StringArray(truth)) << test_case << ": Wrong conversion to StringArray";
    this->check_serialization(array, test_case);

    array.linearize();
    this->check_array(array, truth, all_linearized, rc, true, test_case + " (linearized)");
    EXPECT_EQ(array.data_array_length, array.length()) << test_case << ": Superstrings were not removed";
    EXPECT_EQ(array.to_string_array(), gbwt::StringArray(truth)) << test_case << ": Wrong conversion to StringArray after linearization";
  }
}

TEST_F(DedupStringArrayTest, Swap)
{
  DedupStringArray first = this->build(true);
  DedupStringArray second(gbwt::StringArray(this->truth(false)), false);
  DedupStringArray first_copy = first, second_copy = second;

  first.swap(second);
  EXPECT_TRUE(first.identical(second_copy)) << "First array is not the original second array";
  EXPECT_TRUE(second.identical(first_copy)) << "Second array is not the original first array";
  EXPECT_EQ(first.data_array_length, second_copy.data_array_length) << "Data array length was not swapped";
  EXPECT_EQ(second.data_array_length, first_copy.data_array_length) << "Data array length was not swapped";
}

TEST_F(DedupStringArrayTest, InvalidHeader)
{
  for(bool rc : { false, true })
  {
    std::string test_case = (rc ? "With reverse complements" : "Without reverse complements");
    DedupStringArray array = this->build(rc);

    // Sanity check: the header is where we expect it to be.
    this->expect_invalid_header(array, 0, array.total_length + 2, test_case + ": wrong total length");
    this->expect_invalid_header(array, 1, array.total_length + 2, test_case + ": linearized length exceeds total length");
    this->expect_invalid_header(array, 2, array.linearized_length - 2, test_case + ": linearized length exceeds data array length");
    this->expect_invalid_header(array, 2, array.data_array_length + 2, test_case + ": wrong data array length");
  }

  // With reverse complements, all lengths must be even.
  DedupStringArray array = this->build(true);
  this->expect_invalid_header(array, 0, array.total_length + 1, "Odd total length");
  this->expect_invalid_header(array, 1, array.linearized_length + 1, "Odd linearized length");
  this->expect_invalid_header(array, 2, array.data_array_length + 1, "Odd data array length");
}

TEST_F(DedupStringArrayTest, InvalidIndex)
{
  // Index of a source string that is not in the linearized prefix.
  size_t source_index = 1;

  for(bool rc : { false, true })
  {
    std::string test_case = (rc ? "With reverse complements" : "Without reverse complements");
    DedupStringArray original = this->build(rc);
    size_t i = (rc ? 2 * source_index : source_index);
    ASSERT_FALSE(original.linearized(i)) << test_case << ": Invalid test setup";

    DedupStringArray array = original;
    array.offsets[i] = array.linearized_length - 1;
    this->expect_invalid(array, test_case + ": string partially in the linearized prefix");

    array = original;
    array.offsets[i] = (rc ? array.linearized_length + (array.data_array_length - array.linearized_length) / 2 : array.data_array_length);
    this->expect_invalid(array, test_case + ": string outside the data array");
  }

  // Without reverse complements: a linearized string in the wrong position.
  {
    DedupStringArray array = this->build(false);
    ASSERT_TRUE(array.linearized(0) && array.linearized(4)) << "Invalid test setup";
    array.offsets[4] = 0;
    this->expect_invalid(array, "Linearized string in the wrong position");
  }

  // With reverse complements: a gap in the linearized prefix.
  // Original strings AC (linearized) and GA (in superstring GATT).
  {
    DedupStringArray array;
    array.reverse_complements = true;
    array.total_length = 8;
    array.linearized_length = 6;
    array.data_array_length = 14;
    array.offsets = sdsl::int_vector<0>(5, 0, 4);
    array.offsets[0] = 0; array.offsets[1] = 2; array.offsets[2] = 6; array.offsets[3] = 12; array.offsets[4] = 14;
    array.lengths = sdsl::int_vector<0>(4, 2, 2);
    std::string data = "ACGTNNGATTAATC";
    array.data = std::vector<char>(data.begin(), data.end());
    this->expect_invalid(array, "Gap in the linearized prefix");
  }
}

//------------------------------------------------------------------------------

} // namespace
