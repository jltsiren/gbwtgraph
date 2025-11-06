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
  std::vector<std::pair<std::string, std::string>> test_vectors =
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
