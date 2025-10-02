// search_state_key.hpp
#pragma once

#include <gbwt/gbwt.h>
#include <utility>
#include <cstddef>
#include <functional>

namespace gbwtgraph {

// Auxiliary structures and functions for hasing search states and paths in the GBWT.
namespace detail {

// Hashable struct for search states in the GBWT.
// It stores a node and a range, representing the specific search state in the GBWT.
struct SearchStateKey {
    gbwt::node_type node;
    std::pair<size_t, size_t> range;

    inline bool operator==(const SearchStateKey& other) const {
        return node == other.node && range == other.range;
    }
};

// Which fields are we using for hashing the GBWT paths, suppose the path
// follows the structure: SAMPLE:CONTIG:PHASE:COUNT
// - SampleOnly: Only sample is used for hashing.
// - SamplePhase: Sample and phase are used for hashing.
// - Full: Sample, phase, and contig are used for hashing.
enum class HashMode {
    SampleOnly, 
    SamplePhase,
    Full,
};

// Key wrapper for hashing and equality based on the current mode
struct PathKey {
    gbwt::PathName path;
    HashMode mode;

    inline bool operator==(const PathKey& other) const {
        if (mode != other.mode) return false;
        if (path.sample != other.path.sample) return false;
        if (mode == HashMode::SampleOnly) return true;
        if (path.phase != other.path.phase) return false;
        if (mode == HashMode::SamplePhase) return true;
        return path.contig == other.path.contig;
    }
};

// Hash function for PathKey, which uses the mode to determine which fields to hash.
// The hashing function is based on boost's hash_combine pattern.
struct PathKeyHasher {
    inline size_t operator()(const PathKey& key) const {
        size_t h = std::hash<size_t>{}(key.path.sample);
        if (key.mode == HashMode::SampleOnly) return h;
        h ^= std::hash<size_t>{}(key.path.phase) + 0x9e3779b9 + (h << 6) + (h >> 2);
        if (key.mode == HashMode::SamplePhase) return h;
        h ^= std::hash<size_t>{}(key.path.contig) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

// PathIDMap is a class that maps paths to unique IDs based on the current HashMode.
class PathIDMap {
public:
    // Keep the default constructor
    PathIDMap() = default;

    // We also use the default assignment operators, but we get those automatically.

    inline explicit PathIDMap(HashMode mode) : mode(mode) {}

    inline explicit PathIDMap(const gbwt::Metadata& metadata) {
        for (HashMode try_mode : { HashMode::Full, HashMode::SamplePhase, HashMode::SampleOnly }) {
            if (build_map(metadata, try_mode)) {
                mode = try_mode;
                return;
            }
        }
        std::cerr << "[Minimizer Index] PathIDMap: Too many distinct path combinations (>64); collapsing haplotype info to 0\n";
        mode = HashMode::SampleOnly;
        collapse_all = true;
        map.clear(); 
    }

    inline uint8_t id(const gbwt::PathName& path) const {
        if (collapse_all) return 0;
        PathKey key{path, mode};
        auto it = map.find(key);
        if (it != map.end()) return it->second;
        std::cerr << "PathIDMap: Path not found!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    inline HashMode current_mode() const {
        return mode;
    }

    inline size_t size() const {
        return map.size();
    }

private:
    using Map = std::unordered_map<PathKey, uint8_t, PathKeyHasher>;

    HashMode mode = HashMode::Full;
    bool collapse_all = false;
    Map map;

    inline bool build_map(const gbwt::Metadata& metadata, HashMode try_mode) {
        map.clear();
        for (size_t i = 0; i < metadata.paths(); ++i) {
            PathKey key{metadata.path(i), try_mode};
            if (map.count(key)) continue;
            if (map.size() >= 64) return false;
            map[key] = static_cast<uint8_t>(map.size());
        }
        return true;
    }
};

} // namespace detail
} // namespace gbwtgraph

// Provide a hash specialization for SearchStateKey 
namespace std {
template<>
struct hash<gbwtgraph::detail::SearchStateKey> {
    inline std::size_t operator()(const gbwtgraph::detail::SearchStateKey& key) const {
        std::size_t seed = std::hash<gbwt::node_type>{}(key.node);
        seed ^= std::hash<size_t>{}(key.range.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<size_t>{}(key.range.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};
} // namespace std
