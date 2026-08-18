#include <gtest/gtest.h>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>

#include "shared.h"

using namespace gbwtgraph;

namespace
{

// Subgraph constructor is tested in test_gfa.cpp.
// Chunking and merging GBZs is tested in test_algorithms.cpp.

//------------------------------------------------------------------------------

class GBZSerialization : public ::testing::Test
{
public:
  std::unique_ptr<GBZ> create_gbz()
  {
    NaiveGraph source = build_naive_graph(false);
    return std::make_unique<GBZ>(build_gbwt_index(), source);
  }

  void check_gbz(const GBZ& gbz, const GBZ& truth, bool check_tags = true) const
  {
    // GBZ
    ASSERT_EQ(gbz.header, truth.header) << "GBZ: Invalid header";
    if(check_tags)
    {
      ASSERT_EQ(gbz.tags, truth.tags) << "GBZ: Invalid tags";
    }

    // GBWT
    ASSERT_EQ(gbz.index.size(), truth.index.size()) << "GBWT: Invalid size";
    ASSERT_EQ(gbz.index.sequences(), truth.index.sequences()) << "GBWT: Invalid number of sequences";
    ASSERT_EQ(gbz.index.sigma(), truth.index.sigma()) << "GBWT: Invalid alphabet size";
    ASSERT_EQ(gbz.index.effective(), truth.index.effective()) << "GBWT: Invalid effective alphabet size";
    ASSERT_EQ(gbz.index.samples(), truth.index.samples()) << "GBWT: Invalid number of samples";

    // Graph
    ASSERT_EQ(gbz.graph.header, truth.graph.header) << "Graph: Invalid header";
    ASSERT_EQ(gbz.graph.sequences, truth.graph.sequences) << "Graph: Invalid sequences";
    ASSERT_EQ(gbz.graph.real_nodes, truth.graph.real_nodes) << "Graph: Invalid real nodes";
    ASSERT_EQ(gbz.graph.segments, truth.graph.segments) << "Graph: Invalid segments";
    ASSERT_EQ(gbz.graph.node_to_segment, truth.graph.node_to_segment) << "Graph: Invalid node-to-segment mapping";
  }

  void simple_sds_serialize_v1(const GBZ& graph, const std::string& filename)
  {
    std::ofstream out(filename, std::ios_base::binary);
    if(!out)
    {
      throw sdsl::simple_sds::CannotOpenFile(filename, true);
    }
    out.exceptions(std::ios::failbit | std::ios::badbit);
    graph.simple_sds_serialize_v1(out);
    out.close();
  }
};

TEST_F(GBZSerialization, Empty)
{
  GBZ empty;
  size_t expected_size = empty.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
  std::string filename = gbwt::TempFile::getName("gbz");
  sdsl::simple_sds::serialize_to(empty, filename);

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  size_t bytes = gbwt::fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid file size";
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, empty);

  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, EmptyV1)
{
  GBZ empty;
  std::string filename = gbwt::TempFile::getName("gbz");
  this->simple_sds_serialize_v1(empty, filename);

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, empty);

  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, NonEmpty)
{
  std::unique_ptr<GBZ> original = this->create_gbz();
  size_t expected_size = original->simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
  std::string filename = gbwt::TempFile::getName("gbz");
  sdsl::simple_sds::serialize_to(*original, filename);

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  size_t bytes = gbwt::fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid file size";
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, *original);

  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, NonEmptyV1)
{
  std::unique_ptr<GBZ> original = this->create_gbz();
  std::string filename = gbwt::TempFile::getName("gbz");
  this->simple_sds_serialize_v1(*original, filename);

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, *original);

  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, ExternalObjects)
{
  // Serialization into separate GBWT / GBWTGraph files does not preserve tags.
  // We therefore do not check file sizes.
  std::unique_ptr<GBZ> original = this->create_gbz();
  std::string filename = gbwt::TempFile::getName("gbz");
  std::ofstream out(filename, std::ios_base::binary);
  GBZ::simple_sds_serialize(original->index, original->graph, out);
  out.close();

  GBZ duplicate;
  std::ifstream in(filename, std::ios_base::binary);
  duplicate.simple_sds_load(in);
  in.close();
  this->check_gbz(duplicate, *original, false);

  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, CopyAndSerialize)
{
  std::string filename = gbwt::TempFile::getName("gbz");
  {
    std::unique_ptr<GBZ> original = this->create_gbz();
    GBZ copied = *original; original.reset();
    sdsl::simple_sds::serialize_to(copied, filename);
  }
  {
    std::unique_ptr<GBZ> truth = this->create_gbz();
    GBZ loaded; sdsl::simple_sds::load_from(loaded, filename);
    this->check_gbz(loaded, *truth);
  }
  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, MoveAndSerialize)
{
  std::string filename = gbwt::TempFile::getName("gbz");
  {
    std::unique_ptr<GBZ> original = this->create_gbz();
    GBZ moved = std::move(*original); original.reset();
    sdsl::simple_sds::serialize_to(moved, filename);
  }
  {
    std::unique_ptr<GBZ> truth = this->create_gbz();
    GBZ loaded; sdsl::simple_sds::load_from(loaded, filename);
    this->check_gbz(loaded, *truth);
  }
  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, SwapAndSerialize)
{
  std::string filename = gbwt::TempFile::getName("gbz");
  {
    std::unique_ptr<GBZ> original = this->create_gbz();
    GBZ swapped;
    swapped.swap(*original); original.reset();
    sdsl::simple_sds::serialize_to(swapped, filename);
  }
  {
    std::unique_ptr<GBZ> truth = this->create_gbz();
    GBZ loaded; sdsl::simple_sds::load_from(loaded, filename);
    this->check_gbz(loaded, *truth);
  }
  gbwt::TempFile::remove(filename);
}

TEST_F(GBZSerialization, LoadTags)
{
  std::unique_ptr<GBZ> original = this->create_gbz();
  std::string filename = gbwt::TempFile::getName("gbz");
  sdsl::simple_sds::serialize_to(*original, filename);

  gbwt::Tags loaded_tags = GBZ::simple_sds_load_tags(filename);
  ASSERT_EQ(loaded_tags, original->tags) << "Invalid tags loaded from file";

  gbwt::TempFile::remove(filename);
}

//------------------------------------------------------------------------------

class GBZFunctionality : public ::testing::Test
{
public:
  GBZ build_gbz(const std::string& graph_name)
  {
    auto parse = gfa_to_gbwt(graph_name);
    return GBZ(parse.first, parse.second);
  }

  void set_reference_samples(GBZ& gbz, const sample_name_set& samples, size_t expected, const std::string& test)
  {
    size_t present = gbz.set_reference_samples(samples);
    ASSERT_EQ(present, expected) << test << ": Unexpected number of sample names present in the graph";
  }

  void check_named_paths(const GBZ& gbz, const sample_name_set& true_samples, size_t expected_paths, const std::string& test)
  {
    ASSERT_EQ(gbz.named_paths(), expected_paths) << test << ": Invalid number of named paths";

    const sample_name_set& samples = gbz.get_reference_samples();
    ASSERT_EQ(samples.size(), true_samples.size()) << test << ": Invalid number of reference samples";
    for(const std::string& sample : true_samples)
    {
      ASSERT_TRUE(samples.find(sample) != samples.end()) << test << ": Missing reference sample " << sample;
    }
  }

  // The test GFA files have their graph names computed with the reference pggname implementation.
  // Tests in test_gfa.cpp ensure that our implementation computes the same names.
  void check_graph_name(GBZ& gbz, const GraphName* parent, bool missing_name, const std::string& test_case)
  {
    EXPECT_EQ(gbz.pggname().empty(), missing_name) << "Unexpected graph name presence for " << test_case;
    EXPECT_EQ(gbz.translation_target(), "") << "Translation target should not be set for " << test_case;
    GraphName old = gbz.graph_name();
    EXPECT_EQ(gbz.pggname(), old.name()) << "Graph name mismatch for " << test_case;

    gbz.compute_pggname(parent);
    GraphName recomputed = gbz.graph_name();
    EXPECT_TRUE(recomputed.has_name()) << "Recomputed graph name missing for " << test_case;
    if(!missing_name)
    {
      EXPECT_EQ(recomputed, old) << "Graph name changed after recomputation for " << test_case;
    }
  }
};

TEST_F(GBZFunctionality, ReferenceSamples)
{
  GBZ gbz = this->build_gbz("gfas/components_ref.gfa");
  sample_name_set samples { "ref" };
  sample_name_set true_samples = samples;
  this->check_named_paths(gbz, true_samples, 2, "Initial graph");

  samples.erase("ref");
  samples.insert("sample");
  true_samples = samples;
  this->set_reference_samples(gbz, samples, 1, "New sample");
  this->check_named_paths(gbz, true_samples, 4, "New sample");

  samples.insert("ref");
  true_samples = samples;
  this->set_reference_samples(gbz, samples, 2, "Both samples");
  this->check_named_paths(gbz, true_samples, 6, "Both samples");

  samples.erase("sample");
  true_samples = samples;
  samples.insert("missing");
  this->set_reference_samples(gbz, samples, 1, "Invalid sample");
  this->check_named_paths(gbz, true_samples, 2, "Invalid sample");
}

TEST_F(GBZFunctionality, GraphNames)
{
  // Additional tests with parent graph relationships are in test_gfa.cpp.

  // Constructor from gfa_to_gbwt() output.
  {
    std::unique_ptr<gbwt::GBWT> index = std::make_unique<gbwt::GBWT>(build_gbwt_index());
    std::unique_ptr<NaiveGraph> graph = std::make_unique<NaiveGraph>(std::move(build_naive_graph(false)));
    GBZ gbz(index, graph);
    this->check_graph_name(gbz, nullptr, false, "gfa_to_gbwt() output");
  }

  // Constructor from GBWT and NaiveGraph.
  {
    gbwt::GBWT index = build_gbwt_index();
    NaiveGraph graph = build_naive_graph(false);
    GBZ gbz(index, graph);
    this->check_graph_name(gbz, nullptr, false, "GBWT and NaiveGraph");
  }

  // Subgraph construction.
  {
    gbwt::GBWT index = build_gbwt_index();
    NaiveGraph graph = build_naive_graph(false);
    GBZ supergraph(index, graph);
    GraphName parent = supergraph.graph_name();

    std::vector<gbwt::vector_type> paths { alt_path };
    gbwt::GBWT sub_index = build_gbwt(paths);
    GBZ gbz(std::move(sub_index), supergraph);
    this->check_graph_name(gbz, &parent, false, "subgraph construction");
  }

  // Constructor from GBWT and HandleGraph; graph name must be set manually.
  {
    gbwt::GBWT parent_index = build_gbwt_index();
    NaiveGraph graph = build_naive_graph(false);
    GBWTGraph parent_graph(parent_index, graph);

    gbwt::GBWT index = parent_index;
    GBZ gbz(std::move(index), parent_graph, nullptr);
    this->check_graph_name(gbz, nullptr, true, "GBWT and HandleGraph");
  }
}

//------------------------------------------------------------------------------

#ifdef GBWTGRAPH_ENABLE_SHARED_MEMORY

// These use two independent bi::managed_shared_memory handles to the same
// named segment to stand in for two separate processes: nothing here is
// shared except the segment name, so a real cross-process attach has to
// work the same way this does.
//
// This exercises GBWTGraph's construction path (via GBZ) actually
// publishing real per-node sequence data under a discoverable name in the
// segment. gbwt::StringArray itself is already covered by gbwt's own
// StringArraySharedMemoryTest (see gbwt/tests/test_support.cpp).
class GBZSharedMemoryTest : public ::testing::Test
{
public:
  std::string segment_name;

  void SetUp() override
  {
    this->segment_name = "gbwtgraph_test_shared_memory_" + std::to_string(::getpid());
    bi::shared_memory_object::remove(this->segment_name.c_str());
  }

  void TearDown() override
  {
    bi::shared_memory_object::remove(this->segment_name.c_str());
  }
};

/*
  A segment is only meaningful for the allocator that puts characters in one,
  and the constructors that build sequences on the heap are only meaningful for
  the allocators that do not. Which instantiation has which is part of the
  interface, so it is checked here rather than left to whatever the rest of the
  tests happen to name.

  A requires-expression only turns a failed requirement into `false` while
  substituting template arguments, so these are written against a parameter
  rather than naming the instantiations directly.
*/

template<typename GBZType>
concept BuildableFromSegment = requires(bi::managed_shared_memory* segment) { GBZType(segment); };
template<typename GBZType>
concept BuildableOnHeap = requires(const gbwt::GBWT& index, const NaiveGraph& graph) { GBZType(index, graph); };
template<typename GBZType>
concept MergeableFromSubgraphs = requires(std::vector<GBZ>&& subgraphs) { GBZType(std::move(subgraphs)); };

static_assert(BuildableFromSegment<CoreGBZ<SharedMemCharAllocatorType>>,
              "A shared-memory GBZ should be constructible from a segment");
static_assert(!BuildableFromSegment<GBZ>,
              "A plain GBZ should not be constructible from a segment");

static_assert(BuildableOnHeap<GBZ>,
              "A plain GBZ should be constructible from an index and a graph");
static_assert(!BuildableOnHeap<CoreGBZ<SharedMemCharAllocatorType>>,
              "A shared-memory GBZ should not be constructible without being given a segment");

// Subgraphs are always plain (see CoreGBWTGraph::subgraph()), so merging them
// can only produce a plain GBZ.
static_assert(MergeableFromSubgraphs<GBZ>,
              "A plain GBZ should be constructible by merging subgraphs");
static_assert(!MergeableFromSubgraphs<CoreGBZ<SharedMemCharAllocatorType>>,
              "A shared-memory GBZ should not be constructible by merging subgraphs");

TEST_F(GBZSharedMemoryTest, SequencesAttachFromIndependentHandle)
{
  // An independently built, plain GBZ is the ground truth for what the
  // node sequences should be; build_naive_graph()/build_gbwt_index() are
  // assumed deterministic across calls (other tests in this file rely on
  // the same assumption).
  GBZ truth(build_gbwt_index(), build_naive_graph(false));

  bi::managed_shared_memory writer_segment(bi::create_only, this->segment_name.c_str(), 16 * 1024 * 1024);
  CoreGBZ<SharedMemCharAllocatorType> writer(build_gbwt_index(), build_naive_graph(false), &writer_segment);
  ASSERT_EQ(writer.graph.sequences.size(), truth.graph.sequences.size()) << "Writer GBZ has the wrong number of sequences";
  for(size_t i = 0; i < truth.graph.sequences.size(); i++)
  {
    EXPECT_EQ(writer.graph.sequences.str(i), truth.graph.sequences.str(i)) << "Writer GBZ has the wrong sequence at offset " << i;
  }

  // A second, independent handle to the same segment stands in for a
  // second process. GBZ has no constructor that only attaches, so this names
  // the published sequences directly, the way a reader would have to.
  bi::managed_shared_memory reader_segment(bi::open_only, this->segment_name.c_str());
  gbwt::StringArray<SharedMemCharAllocatorType> attached(&reader_segment, "sequences");
  ASSERT_EQ(attached.size(), truth.graph.sequences.size()) << "Attached sequences have the wrong size";
  for(size_t i = 0; i < truth.graph.sequences.size(); i++)
  {
    EXPECT_EQ(attached.str(i), truth.graph.sequences.str(i)) << "Attached sequence at offset " << i << " does not match the source graph";
  }

  // Write-through visibility between two handles to the same segment is
  // exercised directly on gbwt::StringArray by StringArraySharedMemoryTest
  // (see gbwt/tests/test_support.cpp), so it is not repeated here.
}

// A segment nothing has published sequences into yet gives an empty array,
// which is what a GBZ built against it starts from.
TEST_F(GBZSharedMemoryTest, MissingSequencesAreEmpty)
{
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 1024 * 1024);
  gbwt::StringArray<SharedMemCharAllocatorType> sequences(&segment, "sequences");
  EXPECT_TRUE(sequences.empty()) << "Sequences nothing has published should be empty";
}

// Loading a serialized GBZ into a CoreGBZ whose segment another process has
// already published node sequences into must end up with those sequences
// rather than a second copy of them.
TEST_F(GBZSharedMemoryTest, LoadThenAttach)
{
  GBZ truth(build_gbwt_index(), build_naive_graph(false));

  bi::managed_shared_memory writer_segment(bi::create_only, this->segment_name.c_str(), 32 * 1024 * 1024);
  CoreGBZ<SharedMemCharAllocatorType> writer(build_gbwt_index(), build_naive_graph(false), &writer_segment);

  std::string filename = gbwt::TempFile::getName("gbz");
  sdsl::simple_sds::serialize_to(writer, filename);

  // A second, independent handle to the same segment stands in for a
  // second process, loading the same serialized bytes.
  bi::managed_shared_memory reader_segment(bi::open_only, this->segment_name.c_str());
  CoreGBZ<SharedMemCharAllocatorType> reader(&reader_segment);
  {
    std::ifstream in(filename, std::ios_base::binary);
    reader.simple_sds_load(in);
  }

  ASSERT_EQ(reader.graph.sequences.size(), truth.graph.sequences.size()) << "Reader has the wrong number of sequences";
  for(size_t i = 0; i < truth.graph.sequences.size(); i++)
  {
    EXPECT_EQ(reader.graph.sequences.str(i), truth.graph.sequences.str(i)) << "Reader has the wrong sequence at offset " << i;
  }

  gbwt::TempFile::remove(filename);
}

#endif // GBWTGRAPH_ENABLE_SHARED_MEMORY

//------------------------------------------------------------------------------

} // namespace
