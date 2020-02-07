#ifndef GBWTGRAPH_MINIMIZER_H
#define GBWTGRAPH_MINIMIZER_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include <gbwt/utils.h>
#include <sdsl/int_vector.hpp>

#include <gbwtgraph/io.h>
#include <gbwtgraph/utils.h>

/*
  minimizer.h: Minimizer index.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

struct MinimizerHeader
{
  std::uint32_t tag, version;
  std::uint64_t k, w;
  std::uint64_t keys, capacity, max_keys;
  std::uint64_t values;
  std::uint64_t unique;
  std::uint64_t flags;

  constexpr static std::uint32_t TAG = 0x31513151;
  constexpr static std::uint32_t VERSION = Version::MINIMIZER_VERSION;

  constexpr static std::uint64_t FLAG_MASK       = 0x00FF;
  constexpr static std::uint64_t FLAG_KEY_MASK   = 0x00FF;
  constexpr static size_t        FLAG_KEY_OFFSET = 0;

  // Flags for old compatible versions.
  constexpr static std::uint32_t FIRST_VERSION   = 4;
  constexpr static std::uint64_t FIRST_FLAG_MASK = 0x0000;

  MinimizerHeader();
  MinimizerHeader(size_t kmer_length, size_t window_length, size_t initial_capacity, double max_load_factor, size_t key_bits);
  void sanitize(size_t kmer_max_length);
  bool check() const;

  void update_version(size_t key_bits);

  // Boolean flags.
  void set(std::uint64_t flag) { this->flags |= flag; }
  void unset(std::uint64_t flag) { this->flags &= ~flag; }
  bool get_flag(std::uint64_t flag) const { return (this->flags & flag); }

  // Integer flags.
  void set_int(std::uint64_t mask, size_t offset, size_t value);
  size_t get_int(std::uint64_t mask, size_t offset) const;

  size_t key_bits() const;

  bool operator==(const MinimizerHeader& another) const;
  bool operator!=(const MinimizerHeader& another) const { return !(this->operator==(another)); }
};

//------------------------------------------------------------------------------

/*
  A kmer encoded using 2 bits/character in a 64-bit integer. The encoding is only
  defined if all characters in the kmer are valid.
*/

struct Key128;

struct Key64
{
public:
  // Internal representation.
  typedef std::uint64_t key_type;
  key_type key;

  // Empty key.
  constexpr Key64() : key(EMPTY_KEY) {}

  // No key.
  constexpr static Key64 no_key() { return Key64(NO_KEY); }

  // Implicit conversion for testing.
  constexpr Key64(key_type value) : key(value) {}

  // Comparisons.
  bool operator==(Key64 another) const { return (this->key == another.key); }
  bool operator!=(Key64 another) const { return (this->key != another.key); }
  bool operator<(Key64 another) const { return (this->key < another.key); }

  // Hash of the key.
  size_t hash() const { return gbwt::wang_hash_64(this->key); }

  // Move the kmer forward, with c as the next character. Update the key, assuming that
  // it encodes the kmer in forward orientation.
  void forward(size_t k, unsigned char c, size_t& valid_chars)
  {
    key_type packed = CHAR_TO_PACK[c];
    if(packed > PACK_MASK) { this->key = EMPTY_KEY; valid_chars = 0; }
    else
    {
      this->key = ((this->key << PACK_WIDTH) | packed) & KMER_MASK[k];
      valid_chars++;
    }
  }

  // Move the kmer forward, with c as the next character. Update the key, assuming that
  // it encodes the kmer in reverse orientation.
  void reverse(size_t k, unsigned char c)
  {
    key_type packed = CHAR_TO_PACK[c];
    if(packed > PACK_MASK) { this->key = EMPTY_KEY; }
    else
    {
      packed ^= PACK_MASK; // The complement of the base.
      this->key = (packed << ((k - 1) * PACK_WIDTH)) | (this->key >> PACK_WIDTH);
    }
  }

  /// Encode a string of size k to a key.
  static Key64 encode(const std::string& sequence);

  /// Decode the key back to a string, given the kmer size used.
  std::string decode(size_t k) const;

  // Required numeric constants.
  constexpr static std::size_t KEY_BITS = sizeof(key_type) * gbwt::BYTE_BITS;
  constexpr static std::size_t KMER_LENGTH = 21;
  constexpr static std::size_t WINDOW_LENGTH = 11;
  constexpr static std::size_t KMER_MAX_LENGTH = 31;

private:
  // Specific key values.
  constexpr static key_type EMPTY_KEY = 0;
  constexpr static key_type NO_KEY = std::numeric_limits<key_type>::max();

  // Constants for the encoding between std::string and the key.
  constexpr static size_t   PACK_WIDTH = 2;
  constexpr static key_type PACK_MASK  = 0x3;

  // Arrays for the encoding between std::string and the key.
  const static std::vector<unsigned char> CHAR_TO_PACK;
  const static std::vector<char>          PACK_TO_CHAR;
  const static std::vector<key_type>      KMER_MASK;

  friend Key128;
};

// Required for printing keys.
std::ostream& operator<<(std::ostream& out, Key64 value);

//------------------------------------------------------------------------------

/*
  A kmer encoded using 2 bits/character in a pair of 64-bit integers. The encoding is
  only defined if all characters in the kmer are valid.
*/

struct Key128
{
public:
  // Internal representation.
  typedef std::uint64_t key_type;
  constexpr static std::size_t FIELD_BITS = sizeof(key_type) * gbwt::BYTE_BITS;
  key_type high, low;

  // Empty key.
  constexpr Key128() : high(EMPTY_KEY), low(EMPTY_KEY) {}

  // No key.
  constexpr static Key128 no_key() { return Key128(NO_KEY, NO_KEY); }

  // Implicit conversion for testing.
  constexpr Key128(key_type key) : high(EMPTY_KEY), low(key) {}

  // For testing.
  constexpr Key128(key_type high, key_type low) : high(high), low(low) {}

  // Comparisons.
  bool operator==(Key128 another) const { return (this->high == another.high && this->low == another.low); }
  bool operator!=(Key128 another) const { return (this->high != another.high || this->low != another.low); }
  bool operator<(Key128 another) const
  {
    if(this->high != another.high) { return (this->high < another.high); }
    return (this->low < another.low);
  }

  // Hash of the key. Essentially boost::hash_combine.
  size_t hash() const
  {
    size_t result = gbwt::wang_hash_64(this->high);
    result ^= gbwt::wang_hash_64(this->low) + 0x9e3779b9 + (result << 6) + (result >> 2);
    return result;
  }

  // Move the kmer forward, with c as the next character. Update the key, assuming that
  // it encodes the kmer in forward orientation.
  void forward(size_t k, unsigned char c, size_t& valid_chars)
  {
    key_type packed = CHAR_TO_PACK[c];
    if(packed > PACK_MASK) { this->high = EMPTY_KEY; this->low = EMPTY_KEY; valid_chars = 0; }
    else
    {
      this->high = ((this->high << PACK_WIDTH) | (this->low >> PACK_OVERFLOW)) & HIGH_MASK[k];
      this->low = ((this->low << PACK_WIDTH) | packed) & LOW_MASK[k];
      valid_chars++;
    }
  }

  // Move the kmer forward, with c as the next character. Update the key, assuming that
  // it encodes the kmer in reverse orientation.
  void reverse(size_t k, unsigned char c)
  {
    key_type packed = CHAR_TO_PACK[c];
    if(packed > PACK_MASK) { this->high = EMPTY_KEY; this->low = EMPTY_KEY; }
    else
    {
      packed ^= PACK_MASK; // The complement of the base.
      if(k > FIELD_CHARS)
      {
        this->low = ((this->high & PACK_MASK) << PACK_OVERFLOW) | (this->low >> PACK_WIDTH);
        this->high = (packed << ((k - FIELD_CHARS - 1) * PACK_WIDTH)) | (this->high >> PACK_WIDTH);
      }
      else // The entire kmer is in the lower part of the key.
      {
        this->low = (packed << ((k - 1) * PACK_WIDTH)) | (this->low >> PACK_WIDTH);
      }
    }
  }

  /// Encode a string of size k to a key.
  static Key128 encode(const std::string& sequence);

  /// Decode the key back to a string, given the kmer size used.
  std::string decode(size_t k) const;

  // Required numeric constants.
  constexpr static std::size_t KEY_BITS = 2 * FIELD_BITS;
  constexpr static std::size_t KMER_LENGTH = 39;
  constexpr static std::size_t WINDOW_LENGTH = 15;
  constexpr static std::size_t KMER_MAX_LENGTH = 63;

private:
  // Specific key values.
  constexpr static key_type EMPTY_KEY = 0;
  constexpr static key_type NO_KEY = std::numeric_limits<key_type>::max();

  // Constants for the encoding between std::string and the key.
  constexpr static size_t   PACK_WIDTH    = 2;
  constexpr static size_t   PACK_OVERFLOW = FIELD_BITS - PACK_WIDTH;
  constexpr static size_t   FIELD_CHARS   = FIELD_BITS / PACK_WIDTH;
  constexpr static key_type PACK_MASK     = 0x3;

  // Arrays for the encoding between std::string and the key.
  const static std::vector<unsigned char> CHAR_TO_PACK;
  const static std::vector<char>          PACK_TO_CHAR;
  const static std::vector<key_type>      HIGH_MASK, LOW_MASK;
};

// Required for printing keys.
std::ostream& operator<<(std::ostream& out, Key128 value);

//------------------------------------------------------------------------------

/*
  A class that implements the minimizer index as a hash table mapping kmers to sets of pos_t.
  The hash table uses quadratic probing with power-of-two size.
  We encode kmers using 2 bits/character and hash the encoding. A minimizer is the kmer with
  the smallest hash in a window of w consecutive kmers and their reverse complements.

  Index versions (this should be in the wiki):

    1  The initial version.

    2  Minimizer selection is based on hashes instead of lexicographic order. A sequence and
       its reverse complement have the same minimizers, reducing index size by 50%. Not
       compatible with version 1.

    3  Construction-time hit cap is no longer used. Compatible with version 2.

    4  Use SDSL data structures where appropriate. Not compatible with version 3.

    5  An option to use 128-bit keys. Compatible with version 4.
*/

template<class KeyType = Key64>
class MinimizerIndex
{
public:
  typedef KeyType key_type;
  typedef std::uint64_t code_type;
  typedef std::uint32_t offset_type;

  // Public constants.
  constexpr static size_t    INITIAL_CAPACITY = 1024;
  constexpr static double    MAX_LOAD_FACTOR  = 0.77;
  constexpr static code_type NO_VALUE         = 0;

  // Serialize the hash table in blocks of this many cells.
  constexpr static size_t    BLOCK_SIZE       = 4 * gbwt::MEGABYTE;

  const static std::string EXTENSION; // ".min"

  union value_type
  {
    code_type               value;
    std::vector<code_type>* pointer;
  };

  typedef std::pair<key_type, value_type> cell_type;

  // TODO: This should be constexpr in C++14.
  static cell_type empty_cell() { return cell_type(key_type::no_key(), { NO_VALUE }); }

  /*
    The sequence offset of a minimizer is the base that corresponds to the start of the
    minimizer. For minimizers in forward orientation, the offset is the first base
    covered by the minimizer. For reverse orientation, the offset is the last base
    covered by the minimizer.
  */
  struct minimizer_type
  {
    key_type    key;        // Encoded minimizer.
    size_t      hash;       // Hash of the minimizer.
    offset_type offset;     // Sequence offset.
    bool        is_reverse; // The minimizer is the reverse complement of the kmer.

    // Is the minimizer empty?
    bool empty() const { return (this->key == key_type::no_key()); }

    // Sort by (offset, !is_reverse). When the offsets are equal, a reverse complement
    // minimizer is earlier in the sequence than a forward minimizer.
    bool operator<(const minimizer_type& another) const
    {
      return ((this->offset < another.offset) ||
              (this->offset == another.offset && this->is_reverse > another.is_reverse));
    }

    bool operator==(const minimizer_type& another) const
    {
      return (this->key == another.key && this->offset == another.offset && this->is_reverse == another.is_reverse);
    }

    bool operator!=(const minimizer_type& another) const
    {
      return !(this->operator==(another));
    }
  };

//------------------------------------------------------------------------------

  MinimizerIndex() :
    header(KeyType::KMER_LENGTH, KeyType::WINDOW_LENGTH, INITIAL_CAPACITY, MAX_LOAD_FACTOR, KeyType::KEY_BITS),
    hash_table(this->header.capacity, empty_cell()),
    is_pointer(this->header.capacity, 0)
  {
    this->header.sanitize(KeyType::KMER_MAX_LENGTH);
  }

  MinimizerIndex(size_t kmer_length, size_t window_length) :
    header(kmer_length, window_length, INITIAL_CAPACITY, MAX_LOAD_FACTOR, KeyType::KEY_BITS),
    hash_table(this->header.capacity, empty_cell()),
    is_pointer(this->header.capacity, 0)
  {
    this->header.sanitize(KeyType::KMER_MAX_LENGTH);
  }

  MinimizerIndex(const MinimizerIndex& source)
  {
    this->copy(source);
  }

  MinimizerIndex(MinimizerIndex&& source)
  {
    *this = std::move(source);
  }

  ~MinimizerIndex()
  {
    this->clear();
  }

  void swap(MinimizerIndex& another)
  {
    if(&another == this) { return; }
    std::swap(this->header, another.header);
    this->hash_table.swap(another.hash_table);
    this->is_pointer.swap(another.is_pointer);
  }

  MinimizerIndex& operator=(const MinimizerIndex& source)
  {
    if(&source != this) { this->copy(source); }
    return *this;
  }

  MinimizerIndex& operator=(MinimizerIndex&& source)
  {
    if(&source != this)
    {
      this->header = std::move(source.header);
      this->hash_table = std::move(source.hash_table);
      this->is_pointer = std::move(source.is_pointer);
    }
    return *this;
  }

  // Serialize the index to the ostream. Returns the number of bytes written and
  // true if the serialization was successful.
  std::pair<size_t, bool> serialize(std::ostream& out) const
  {
    size_t bytes = 0;
    bool ok = true;

    bytes += io::serialize(out, this->header, ok);
    bytes += io::serialize_hash_table(out, this->hash_table, this->is_pointer, NO_VALUE, ok);
    bytes += this->is_pointer.serialize(out); // We do not know if this was successful.

    // Serialize the occurrence lists.
    for(size_t i = 0; i < this->capacity(); i++)
    {
      if(this->is_pointer[i]) { bytes += io::serialize_vector(out, *(this->hash_table[i].second.pointer), ok); }
    }

    if(!ok)
    {
      std::cerr << "MinimizerIndex::serialize(): Serialization failed" << std::endl;
    }

    return std::make_pair(bytes, ok);
  }

  // Load the index from the istream and return true if successful.
  bool deserialize(std::istream& in)
  {
    bool ok = true;

    // Load and check the header.
    ok &= io::load(in, this->header);
    if(!(this->header.check()))
    {
      std::cerr << "MinimizerIndex::deserialize(): Invalid or old index file" << std::endl;
      std::cerr << "MinimizerIndex::deserialize(): Index version is " << this->header.version << "; expected " << MinimizerHeader::VERSION << std::endl;
      return false;
    }
    // TODO: This depends on index version.
    if(this->header.version == MinimizerHeader::FIRST_VERSION)
    {
      if(KeyType::KEY_BITS != Key64::KEY_BITS)
      {
        std::cerr << "MinimizerIndex::deserialize(): Cannot load version " << this->header.version << " into a " << KeyType::KEY_BITS << "-bit index" << std::endl;
        return false;
      }
    }
    else if(this->header.key_bits() != KeyType::KEY_BITS)
    {
      std::cerr << "MinimizerIndex::deserialize(): Expected " << KeyType::KEY_BITS << "-bit keys, got " << this->header.key_bits() << "-bit keys" << std::endl;
      return false;
    }
    this->header.update_version(KeyType::KEY_BITS);

    // Load the hash table.
    if(ok) { ok &= io::load_vector(in, this->hash_table); }
    if(ok) { this->is_pointer.load(in); }

    // Load the occurrence lists.
    if(ok)
    {
      for(size_t i = 0; i < this->capacity(); i++)
      {
        if(this->is_pointer[i])
        {
          this->hash_table[i].second.pointer = new std::vector<code_type>();
          ok &= io::load_vector(in, *(this->hash_table[i].second.pointer));
        }
      }
    }

    if(!ok)
    {
      std::cerr << "MinimizerIndex::deserialize(): Index loading failed" << std::endl;
    }

    return ok;
  }

  // For testing.
  bool operator==(const MinimizerIndex& another) const
  {
    if(this->header != another.header || this->is_pointer != another.is_pointer) { return false; }

    for(size_t i = 0; i < this->capacity(); i++)
    {
      cell_type a = this->hash_table[i], b = another.hash_table[i];
      if(a.first != b.first) { return false; }
      if(this->is_pointer[i])
      {
        if(*(a.second.pointer) != *(b.second.pointer)) { return false; }
      }
      else
      {
        if(a.second.value != b.second.value) { return false; }
      }
    }

    return true;
  }

  bool operator!=(const MinimizerIndex& another) const { return !(this->operator==(another)); }

//------------------------------------------------------------------------------

  /*
    Returns the minimizer in the window specified by the iterators. If no minimizer
    exists (e.g. because all kmers contain invalid characters), the return value is
    an empty minimizer. If there are multiple occurrences of the minimizer, return
    the leftmost one.
  */
  minimizer_type minimizer(std::string::const_iterator begin, std::string::const_iterator end) const
  {
    minimizer_type result { key_type::no_key(), static_cast<size_t>(0), static_cast<offset_type>(0), false };
    if(static_cast<size_t>(end - begin) < this->k()) { return result; }

    size_t valid_chars = 0;
    bool found = false;
    key_type forward_key, reverse_key;
    for(std::string::const_iterator iter = begin; iter != end; ++iter)
    {
      forward_key.forward(this->k(), *iter, valid_chars);
      reverse_key.reverse(this->k(), *iter);
      if(valid_chars >= this->k())
      {
        size_t forward_hash = forward_key.hash(), reverse_hash = reverse_key.hash();
        size_t hash = std::min(forward_hash, reverse_hash);
        if(!found || hash < result.hash)
        {
          offset_type pos = (iter - begin) + 1 - this->k();
          if(reverse_hash < forward_hash) { result = { reverse_key, reverse_hash, pos, true }; }
          else                            { result = { forward_key, forward_hash, pos, false }; }
          found = true;
        }
      }
    }

    // It was more convenient to use the first offset of the kmer, regardless of the orientation.
    // If the minimizer is a reverse complement, we must return the last offset instead.
    if(result.is_reverse) { result.offset += this->k() - 1; }

    return result;
  }

  /*
    A circular buffer of size 2^i for all minimizer candidates. The candidates are sorted
    by both key and sequence offset. The candidate at offset j is removed when we reach
    offset j + w. Candidates at the tail are removed when we advance with a smaller key.
  */
  struct CircularBuffer
  {
    std::vector<minimizer_type> buffer;
    size_t head, tail;
    size_t w;

    constexpr static size_t BUFFER_SIZE = 16;

    CircularBuffer(size_t capacity) :
      buffer(),
      head(0), tail(0), w(capacity)
    {
      size_t buffer_size = BUFFER_SIZE;
      while(buffer_size < this->w) { buffer_size *= 2; }
      this->buffer.resize(buffer_size);
    }

    bool empty() const { return (this->head >= this->tail); }

    size_t begin() const { return this->head; }
    size_t end() const { return this->tail; }
    minimizer_type& at(size_t i) { return this->buffer[i & (this->buffer.size() - 1)]; }
    minimizer_type& front() { return this->at(this->begin()); }
    minimizer_type& back() { return this->at(this->end() - 1); }

    // Advance to the next offset (pos) with a valid kmer.
    void advance(offset_type pos, key_type forward_key, key_type reverse_key)
    {
      if(!(this->empty()) && this->front().offset + this->w <= pos) { this->head++; }
      size_t forward_hash = forward_key.hash(), reverse_hash = reverse_key.hash();
      size_t hash = std::min(forward_hash, reverse_hash);
      while(!(this->empty()) && this->back().hash > hash) { this->tail--; }
      this->tail++;
      if(reverse_hash < forward_hash) { this->back() = { reverse_key, reverse_hash, pos, true }; }
      else                            { this->back() = { forward_key, forward_hash, pos, false }; }
    }

    // Advance to the next offset (pos) without a valid kmer.
    void advance(offset_type pos)
    {
      if(!(this->empty()) && this->front().offset + this->w <= pos) { this->head++; }
    }
  };

  /*
    Returns all minimizers in the string specified by the iterators. The return
    value is a vector of minimizers sorted by their offsets. If there are multiple
    occurrences of one or more minimizer keys with the same hash in a window,
    return all of them.
  */
  std::vector<minimizer_type> minimizers(std::string::const_iterator begin, std::string::const_iterator end) const
  {
    std::vector<minimizer_type> result;
    size_t window_length = this->k() + this->w() - 1, total_length = end - begin;
    if(total_length < window_length) { return result; }

    // Find the minimizers.
    CircularBuffer buffer(this->w());
    size_t valid_chars = 0, start_pos = 0;
    size_t next_read_offset = 0;  // The first read offset that may contain a new minimizer.
    key_type forward_key, reverse_key;
    std::string::const_iterator iter = begin;
    while(iter != end)
    {
      forward_key.forward(this->k(), *iter, valid_chars);
      reverse_key.reverse(this->k(), *iter);
      if(valid_chars >= this->k()) { buffer.advance(start_pos, forward_key, reverse_key); }
      else                         { buffer.advance(start_pos); }
      ++iter;
      if(static_cast<size_t>(iter - begin) >= this->k()) { start_pos++; }
      // We have a full window with a minimizer.
      if(static_cast<size_t>(iter - begin) >= window_length && !buffer.empty())
      {
        // Insert the candidates if:
        // 1) this is the first minimizer we encounter;
        // 2) the last reported minimizer had the same key (we may have new occurrences); or
        // 3) the first candidate is located after the last reported minimizer.
        if(result.empty() || result.back().hash == buffer.front().hash || result.back().offset < buffer.front().offset)
        {
          // Insert all new occurrences of the minimizer in the window.
          for(size_t i = buffer.begin(); i < buffer.end() && buffer.at(i).hash == buffer.front().hash; i++)
          {
            if(buffer.at(i).offset >= next_read_offset)
            {
              result.emplace_back(buffer.at(i));
              next_read_offset = buffer.at(i).offset + 1;
            }
          }
        }
      }
    }

    // It was more convenient to use the first offset of the kmer, regardless of the orientation.
    // If the minimizer is a reverse complement, we must return the last offset instead.
    for(minimizer_type& minimizer : result)
    {
      if(minimizer.is_reverse) { minimizer.offset += this->k() - 1; }
    }
    std::sort(result.begin(), result.end());

    return result;
  }

  /*
    Returns all minimizers in the string. The return value is a vector of
    minimizers sorted by their offsets.
  */
  std::vector<minimizer_type> minimizers(const std::string& str) const
  {
    return this->minimizers(str.begin(), str.end());
  }
  
  /*
    Returns all minimizers in the string specified by the iterators, together
    with the weight of how many windows they arise from. The return value is a
    vector of pairs of minimizers and window counts sorted by their offsets. If
    there are multiple occurrences of one or more minimizer keys with the same
    hash in a window, they are all returned, but the window's weight is all
    assigned to an arbitrary minimizer that it contains.
  */
  std::vector<std::pair<minimizer_type, size_t>> weighted_minimizers(std::string::const_iterator begin, std::string::const_iterator end) const
  {
    std::vector<std::pair<minimizer_type, size_t>> result;
    size_t window_length = this->k() + this->w() - 1, total_length = end - begin;
    if(total_length < window_length) { return result; }
    
    // Find the minimizers.
    CircularBuffer buffer(this->w());
    size_t valid_chars = 0, start_pos = 0;
    size_t next_read_offset = 0;  // The first read offset that may contain a new minimizer.
    key_type forward_key, reverse_key;
    std::string::const_iterator iter = begin;
    while(iter != end)
    {
      // Get the forward and reverse strand minimizer candidates
      forward_key.forward(this->k(), *iter, valid_chars);
      reverse_key.reverse(this->k(), *iter);
      // If they don't have any Ns or anything in them, throw them into the sliding window tracked by buffer.
      // Otherwise just slide it along.
      if(valid_chars >= this->k()) { buffer.advance(start_pos, forward_key, reverse_key); }
      else                         { buffer.advance(start_pos); }
      ++iter;
      // If we have passed at least k characters, we must advance the starting position of the next kmer.
      if(static_cast<size_t>(iter - begin) >= this->k()) { start_pos++; }
      // We have a full window with a minimizer.
      if(static_cast<size_t>(iter - begin) >= window_length && !buffer.empty())
      {
        // Insert the candidates if:
        // 1) this is the first minimizer we encounter;
        // 2) the last reported minimizer had the same hash (we may have new occurrences); or
        // 3) the first candidate is located after the last reported minimizer.
        if(result.empty() || result.back().first.hash == buffer.front().hash || result.back().first.offset < buffer.front().offset)
        {
          // Insert all new occurrences of the minimizer in the window.
          for(size_t i = buffer.begin(); i < buffer.end() && buffer.at(i).hash == buffer.front().hash; i++)
          {
            if(buffer.at(i).offset >= next_read_offset)
            {
              result.emplace_back(buffer.at(i), 0);
              next_read_offset = buffer.at(i).offset + 1;
            }
          }
        }
        
        // Assign the window's weight to an arbitrary minimizer that occured in it.
        // Whatever is last in result right now will work.
        result.back().second++;
      }
    }

    // It was more convenient to use the first offset of the kmer, regardless of the orientation.
    // If the minimizer is a reverse complement, we must return the last offset instead.
    for(auto& weighted_minimizer : result)
    {
      if(weighted_minimizer.first.is_reverse) { weighted_minimizer.first.offset += this->k() - 1; }
    }
    std::sort(result.begin(), result.end());

    return result;
  }
  
  /*
    Returns all minimizers in the string. The return value is a vector of
    minimizers and window counts sorted by their offsets.
  */
  std::vector<std::pair<minimizer_type, size_t>> weighted_minimizers(const std::string& str) const
  {
    return this->weighted_minimizers(str.begin(), str.end());
  }
  
  /*
    Returns all minimizers in the string specified by the iterators, together
    with the start and length of the run of windows they arise from. The return
    value is a vector of tuples of minimizers, starts, and lengths, sorted by
    minimizer offset.
  */
  std::vector<std::tuple<minimizer_type, size_t, size_t>> minimizer_regions(std::string::const_iterator begin, std::string::const_iterator end) const
  {
    std::vector<std::tuple<minimizer_type, size_t, size_t>> result;
    size_t window_length = this->k() + this->w() - 1, total_length = end - begin;
    if(total_length < window_length) { return result; }
    
    // Find the minimizers.
    CircularBuffer buffer(this->w());
    // Note that start_pos isn't meaningfully the start of the window we are
    // looking at.
    size_t valid_chars = 0, start_pos = 0;
    size_t next_read_offset = 0;  // The first read offset that may contain a new minimizer.
    // All results before this are finished and have their lengths filled in.
    // All results after are current winning minimizers of the current window.
    size_t finished_through = 0; 
    key_type forward_key, reverse_key;
    std::string::const_iterator iter = begin;
    while(iter != end)
    {
    
      
    
      // Get the forward and reverse strand minimizer candidates
      forward_key.forward(this->k(), *iter, valid_chars);
      reverse_key.reverse(this->k(), *iter);
      // If they don't have any Ns or anything in them, throw them into the sliding window tracked by buffer.
      // Otherwise just slide it along.
      if(valid_chars >= this->k()) { buffer.advance(start_pos, forward_key, reverse_key); }
      else                         { buffer.advance(start_pos); }
      ++iter;
      if(static_cast<size_t>(iter - begin) >= this->k()) { start_pos++; }
      
      // We have a full window.
      if(static_cast<size_t>(iter - begin) >= window_length)
      {
        // Work out where the window we are minimizing in began
        size_t window_start = static_cast<size_t>(iter - begin) - window_length;
        
        // Work out the past-the-end index of the window we have just finished (not the current window)
        size_t prev_past_end_pos = window_start + window_length - 1;
      
        std::cerr << "At window start " << window_start << " (current kmer start " << (start_pos - 1)
          << ", prev past-end " << prev_past_end_pos << ")" << std::endl;
      
        // Finish off end positions for results that weren't replaced but are going out of range
        while(finished_through < result.size() &&
          std::get<0>(result[finished_through]).offset < window_start)
        {
          // Compute region length based on it stopping at the previous step
          std::get<2>(result[finished_through]) = prev_past_end_pos - std::get<1>(result[finished_through]);
          std::cerr << "Minimizer " << std::get<0>(result[finished_through]).key.decode(this->k())
            << " finished after covering " << std::get<2>(result[finished_through])
            << " bases due to out of range offset " << std::get<0>(result[finished_through]).offset
            << " < window start " << window_start << std::endl;
          finished_through++;
        }
      
        // Our full window has a minimizer in it
        if (!buffer.empty())
        {
        
          // Insert the candidates if:
          // 1) this is the first minimizer we encounter;
          // 2) the last reported minimizer had the same hash (we may have new occurrences); or
          // 3) the first candidate is located after the last reported minimizer.
          if(result.empty() ||
            std::get<0>(result.back()).hash == buffer.front().hash ||
            std::get<0>(result.back()).offset < buffer.front().offset)
          {
            // Insert all new occurrences of the minimizer in the window.
            for(size_t i = buffer.begin(); i < buffer.end() && buffer.at(i).hash == buffer.front().hash; i++)
            {
              if(buffer.at(i).offset >= next_read_offset)
              {
                // Insert the minimizer instance, with its region starting
                // where the window covered by the buffer starts.
                result.emplace_back(buffer.at(i), window_start, 0);
                // There can only ever really be one minimizer at a given start
                // position. So look for the next one 1 base to the right.
                next_read_offset = buffer.at(i).offset + 1;
                
                std::cerr << "Minimizer " << std::get<0>(result.back()).key.decode(this->k())
                  << " added at offset " <<  std::get<0>(result.back()).offset << std::endl;
              }
            }
            
            // If new minimizers beat out old ones, finish off the old ones.
            while(!result.empty() &&
              finished_through < result.size() &&
              std::get<0>(result.back()).hash != std::get<0>(result[finished_through]).hash)
            {
              // The window before the one we are looking at was the last one for this minimizer.
              std::get<2>(result[finished_through]) = prev_past_end_pos - std::get<1>(result[finished_through]);
              std::cerr << "Minimizer " << std::get<0>(result[finished_through]).key.decode(this->k())
                << " finished due to replacement after covering " << std::get<2>(result[finished_through]) << " bases" << std::endl;
              finished_through++;
            }
          }
        }
      }
    }

    // Now close off the minimizers left active when we hit the end of the string.
    while(finished_through < result.size())
    {
      // The region length is from the region start to the string end
      std::get<2>(result[finished_through]) = total_length - std::get<1>(result[finished_through]);
      std::cerr << "Minimizer " << std::get<0>(result[finished_through]).key.decode(this->k())
        << " finished after covering " << std::get<2>(result[finished_through])
        << " bases due to end of string" << std::endl;
      finished_through++;
    }

    // It was more convenient to use the first offset of the kmer, regardless of the orientation.
    // If the minimizer is a reverse complement, we must return the last offset instead.
    for(auto& record : result)
    {
      if(std::get<0>(record).is_reverse) { std::get<0>(record).offset += this->k() - 1; }
    }
    std::sort(result.begin(), result.end());

    return result;
  }
  
  /*
    Returns all minimizers in the string. The return value is a vector of
    minimizers, region starts, and region lengths sorted by their offsets.
  */
  std::vector<std::tuple<minimizer_type, size_t, size_t>> minimizer_regions(const std::string& str) const
  {
    return this->minimizer_regions(str.begin(), str.end());
  }

//------------------------------------------------------------------------------

  /*
    Inserts the position into the index, using minimizer.key as the key and
    minimizer.hash as its hash. Does not insert empty minimizers or positions.
    The offset of the position will be truncated to fit in OFFSET_BITS bits.
    Use minimizer() or minimizers() to get the minimizer and valid_offset() to check
    if the offset fits in the available space.
    The position should match the orientation of the minimizer: a path label
    starting from the position should have the minimizer as its prefix.
  */
  void insert(const minimizer_type& minimizer, const pos_t& pos)
  {
    if(minimizer.empty() || is_empty(pos)) { return; }

    size_t offset = this->find_offset(minimizer.key, minimizer.hash);
    code_type code = encode(pos);
    if(this->hash_table[offset].first == key_type::no_key())
    {
      this->insert(minimizer.key, code, offset);
    }
    else if(this->hash_table[offset].first == minimizer.key)
    {
      this->append(code, offset);
    }
  }

  /*
    Returns the sorted set of occurrences of the minimizer.
    Use minimizer() or minimizers() to get the minimizer.
    If the minimizer is in reverse orientation, use reverse_base_pos() to reverse
    the reported occurrences.
  */
  std::vector<pos_t> find(const minimizer_type& minimizer) const
  {
    std::vector<pos_t> result;
    if(minimizer.empty()) { return result; }

    size_t offset = this->find_offset(minimizer.key, minimizer.hash);
    cell_type cell = this->hash_table[offset];
    if(cell.first == minimizer.key)
    {
      if(this->is_pointer[offset])
      {
        result.reserve(cell.second.pointer->size());
        for(code_type pos : *(cell.second.pointer)) { result.emplace_back(decode(pos)); }
      }
      else { result.emplace_back(decode(cell.second.value)); }
    }

    return result;
  }

  /*
    Returns the occurrence count of the minimizer.
    Use minimizer() or minimizers() to get the minimizer.
  */
  size_t count(const minimizer_type& minimizer) const
  {
    if(minimizer.empty()) { return 0; }

    size_t offset = this->find_offset(minimizer.key, minimizer.hash);
    if(this->hash_table[offset].first == minimizer.key)
    {
      return (this->is_pointer[offset] ? this->hash_table[offset].second.pointer->size() : 1);
    }

    return 0;
  }

  /*
    Returns the occurrence count of the minimizer and a pointer to the internal
    representation of the occurrences (which are in sorted order). The pointer may be
    invalidated if new positions are inserted into the index.
    Use minimizer() or minimizers() to get the minimizer and decode() to decode the
    occurrences.
  */
  std::pair<size_t, const code_type*> count_and_find(const minimizer_type& minimizer) const
  {
    std::pair<size_t, const code_type*> result(0, nullptr);
    if(minimizer.empty()) { return result; }

    size_t offset = this->find_offset(minimizer.key, minimizer.hash);
    if(this->hash_table[offset].first == minimizer.key)
    {
      const cell_type& cell = this->hash_table[offset];
      if(this->is_pointer[offset])
      {
        result.first = cell.second.pointer->size();
        result.second = cell.second.pointer->data();
      }
      else
      {
        result.first = 1; result.second = &(cell.second.value);
      }
    }
    return result;
  }

//------------------------------------------------------------------------------

  // Length of the kmers in the index.
  size_t k() const { return this->header.k; }

  // Window length for the minimizers.
  size_t w() const { return this->header.w; }

  // Number of keys in the index.
  size_t size() const { return this->header.keys; }

  // Is the index empty.
  bool empty() const { return (this->size() == 0); }

  // Number of values (minimizer occurrences) in the index.
  size_t values() const { return this->header.values; }

  // Size of the hash table.
  size_t capacity() const { return this->header.capacity; }

  // Actual capacity of the hash table. Exceeding it will initiate rehashing.
  size_t max_keys() const { return this->header.max_keys; }

  // Current load factor of the hash table.
  double load_factor() const { return static_cast<double>(this->size()) / static_cast<double>(this->capacity()); }

  // Number of minimizers with a single occurrence.
  size_t unique_keys() const { return this->header.unique; }

//------------------------------------------------------------------------------

private:
  MinimizerHeader        header;
  std::vector<cell_type> hash_table;
  sdsl::bit_vector       is_pointer;

//------------------------------------------------------------------------------

public:

  // Constants for the encoding between pos_t and code_type.
  constexpr static size_t    OFFSET_BITS = 10;
  constexpr static size_t    ID_OFFSET   = OFFSET_BITS + 1;
  constexpr static code_type REV_MASK    = static_cast<code_type>(1) << OFFSET_BITS;
  constexpr static code_type OFF_MASK    = REV_MASK - 1;

  // Is the offset small enough to fit in the low-order bits of the encoding?
  static bool valid_offset(const pos_t& pos) { return (offset(pos) <= OFF_MASK); }

  // Encode pos_t as code_type.
  static code_type encode(const pos_t& pos)
  {
    return (static_cast<code_type>(id(pos)) << ID_OFFSET) |
           (static_cast<code_type>(is_rev(pos)) << OFFSET_BITS) |
           (static_cast<code_type>(offset(pos)) & OFF_MASK);
  }

  // Decode code_type as pos_t.
  static pos_t decode(code_type pos) { return make_pos_t(pos >> ID_OFFSET, pos & REV_MASK, pos & OFF_MASK); }

//------------------------------------------------------------------------------

private:
  void copy(const MinimizerIndex& source)
  {
    this->clear();
    this->header = source.header;
    this->hash_table = source.hash_table;
    this->is_pointer = source.is_pointer;
  }

  // Delete all pointers in the hash table.
  void clear()
  {
    for(size_t i = 0; i < this->hash_table.size(); i++)
    {
      if(this->is_pointer[i])
      {
        delete this->hash_table[i].second.pointer;
        this->hash_table[i].second.value = NO_VALUE;
        this->is_pointer[i] = false;
      }
    }
  }

  // Find the hash table offset for the key with the given hash value.
  size_t find_offset(key_type key, size_t hash) const
  {
    size_t offset = hash & (this->capacity() - 1);
    for(size_t attempt = 0; attempt < this->capacity(); attempt++)
    {
      if(this->hash_table[offset].first == key_type::no_key() || this->hash_table[offset].first == key) { return offset; }

      // Quadratic probing with triangular numbers.
      offset = (offset + attempt + 1) & (this->capacity() - 1);
    }

    // This should not happen.
    std::cerr << "MinimizerIndex::find_offset(): Cannot find the offset for key " << key << std::endl;
    return 0;
  }

  // Insert (key, pos) to hash_table[offset], which is assumed to be empty.
  // Rehashing may be necessary.
  void insert(key_type key, code_type pos, size_t offset)
  {
    this->hash_table[offset].first = key;
    this->hash_table[offset].second.value = pos;
    this->header.keys++;
    this->header.values++;
    this->header.unique++;

    if(this->size() > this->max_keys()) { this->rehash(); }
  }

  // Add pos to the list of occurrences of key at hash_table[offset].
  void append(code_type pos, size_t offset)
  {
    if(this->contains(offset, pos)) { return; }

    if(this->is_pointer[offset])
    {
      std::vector<code_type>* occs = this->hash_table[offset].second.pointer;
      occs->push_back(pos);
      size_t offset = occs->size() - 1;
      while(offset > 0 && occs->at(offset - 1) > occs->at(offset))
      {
        std::swap(occs->at(offset - 1), occs->at(offset));
        offset--;
      }
    }
    else
    {
      std::vector<code_type>* occs = new std::vector<code_type>(2);
      occs->at(0) = this->hash_table[offset].second.value;
      occs->at(1) = pos;
      if(occs->at(0) > occs->at(1)) { std::swap(occs->at(0), occs->at(1)); }
      this->hash_table[offset].second.pointer = occs;
      this->is_pointer[offset] = true;
      this->header.unique--;
    }
    this->header.values++;
  }

  // Does the list of occurrences at hash_table[offset] contain pos?
  bool contains(size_t offset, code_type pos) const
  {
    if(this->is_pointer[offset])
    {
      const std::vector<code_type>* occs = this->hash_table[offset].second.pointer;
      return std::binary_search(occs->begin(), occs->end(), pos);
    }
    else
    {
      return (this->hash_table[offset].second.value == pos);
    }
  }

  // Double the size of the hash table.
  void rehash()
  {
    // Reinitialize with a larger hash table.
    std::vector<cell_type> old_hash_table(2 * this->capacity(), empty_cell());
    sdsl::bit_vector old_is_pointer(2 * this->capacity(), 0);
    this->hash_table.swap(old_hash_table);
    this->is_pointer.swap(old_is_pointer);
    this->header.capacity = this->hash_table.size();
    this->header.max_keys = this->capacity() * MAX_LOAD_FACTOR;

    // Move the keys to the new hash table.
    for(size_t i = 0; i < old_hash_table.size(); i++)
    {
      key_type key = old_hash_table[i].first;
      if(key == key_type::no_key()) { continue; }

      size_t offset = this->find_offset(key, key.hash());
      this->hash_table[offset] = old_hash_table[i];
      this->is_pointer[offset] = old_is_pointer[i];
    }
  }
};

template<class KeyType>
std::ostream&
operator<<(std::ostream& out, const typename MinimizerIndex<KeyType>::minimizer_type& minimizer)
{
  out << "(" << minimizer.key << ", " << (minimizer.is_reverse ? "-" : "+") << minimizer.offset << ")";
  return out;
}

//------------------------------------------------------------------------------

// Choose the default index type.
//typedef MinimizerIndex<Key64> DefaultMinimizerIndex;
typedef MinimizerIndex<Key128> DefaultMinimizerIndex;

//------------------------------------------------------------------------------

// Numerical template class constants.

template<class KeyType> constexpr size_t MinimizerIndex<KeyType>::INITIAL_CAPACITY;
template<class KeyType> constexpr double MinimizerIndex<KeyType>::MAX_LOAD_FACTOR;
template<class KeyType> constexpr typename MinimizerIndex<KeyType>::code_type MinimizerIndex<KeyType>::NO_VALUE;

template<class KeyType> constexpr size_t MinimizerIndex<KeyType>::OFFSET_BITS;
template<class KeyType> constexpr size_t MinimizerIndex<KeyType>::ID_OFFSET;
template<class KeyType> constexpr typename MinimizerIndex<KeyType>::code_type MinimizerIndex<KeyType>::REV_MASK;
template<class KeyType> constexpr typename MinimizerIndex<KeyType>::code_type MinimizerIndex<KeyType>::OFF_MASK;

// Other template class variables.

template<class KeyType> const std::string MinimizerIndex<KeyType>::EXTENSION = ".min";

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_MINIMIZER_H
