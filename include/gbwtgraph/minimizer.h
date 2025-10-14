#ifndef GBWTGRAPH_MINIMIZER_H
#define GBWTGRAPH_MINIMIZER_H

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <utility>
#include <type_traits>
#include <vector>

#include <gbwt/utils.h>

#include <gbwtgraph/io.h>
#include <gbwtgraph/utils.h>

/*
  minimizer.h: Kmer indexes and minimizer indexes.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  A graph position encoded as a 64-bit integer. Uses the lowest 10 bits for
  node offset, limiting node length to 1024 bp. The next bit is used for node
  orientation and the remaining bits for node id.

  NOTE: This encoding must be consistent with the top-level MAX_NODE_LENGTH.
*/
struct Position
{
  typedef std::uint64_t code_type;
  code_type value;

  // Constants that define the encoding.
  constexpr static size_t    OFFSET_BITS = 10;
  constexpr static size_t    ID_OFFSET   = OFFSET_BITS + 1;
  constexpr static code_type REV_MASK    = static_cast<code_type>(1) << OFFSET_BITS;
  constexpr static code_type OFF_MASK    = REV_MASK - 1;

  // As node id 0 is not allowed, we can use 0 as an empty position.
  constexpr static code_type NO_POSITION = 0;

  // Default constructor for technical purposes. Creates an empty position.
  constexpr Position() : value(NO_POSITION) {}

  // Returns an empty position.
  constexpr static Position no_pos() { return Position(); }

  // This reinterprets the encoded value as a Position.
  explicit constexpr Position(code_type value) : value(value) {}

  // Writes the position to the given array pointer.
  void write(code_type* array_ptr) const { *array_ptr = this->value; }

  // Converts pos_t to Position.
  explicit Position(pos_t pos)
  {
    this->value =
      (static_cast<code_type>(gbwtgraph::id(pos)) << ID_OFFSET) |
      (static_cast<code_type>(gbwtgraph::is_rev(pos)) << OFFSET_BITS) |
      (static_cast<code_type>(gbwtgraph::offset(pos)) & OFF_MASK);
  }

  // Basic comparisons.
  bool operator==(Position another) const { return (this->value == another.value); }
  bool operator!=(Position another) const { return (this->value != another.value); }
  bool operator<(Position another) const { return (this->value < another.value); }
  bool operator>(Position another) const { return (this->value > another.value); }

  // Is the offset small enough to fit in the low-order bits of the encoding?
  static bool valid_offset(const pos_t& pos) { return (gbwtgraph::offset(pos) <= OFF_MASK); }

  // Returns the node identifier.
  nid_t id() const { return this->value >> ID_OFFSET; }

  // Returns the orientation.
  bool is_rev() const { return this->value & REV_MASK; }

  // Returns the node offset.
  size_t offset() const { return this->value & OFF_MASK; }

  // Decodes the position as pos_t.
  pos_t decode() const { return make_pos_t(this->id(), this->is_rev(), this->offset()); }
};

//------------------------------------------------------------------------------

/*
  Kmer encoding using 2 bits per base in 64-bit integers. If two kmers have the
  same length, comparing the integers or tuples of integers gives the same
  result as lexicographic comparison of the kmers.

  Also provides some type definitions shared between the kmer and minimizer
  indexes.
*/
struct KmerEncoding {
  // The type used for storing encoded keys and values.
  typedef std::uint64_t code_type;

  // A value consisting of a position and a pointer to the payload,
  // which is assumed to have the same size as the index payload size.
  typedef std::pair<Position, const code_type*> value_type;

  // A pointer to a list of encoded values and the number of values.
  // Use get_value() in the index to access and decode individual values.
  typedef std::pair<const code_type*, size_t> multi_value_type;

  // Conversion from characters to packed characters.
  const static std::vector<unsigned char> CHAR_TO_PACK;

  // Conversion from packed characters to upper-case characters.
  const static std::vector<char> PACK_TO_CHAR;

  // Masks for the bits used in the encoding. If kmer length is `k`,
  // `LOW_MASK[k]` marks the bits that are in use in the least significant word
  // of the encoding. `HIGH_MASK[k]` marks the bits that are in use in the
  // penultimate word of the encoding.
  const static std::vector<code_type> LOW_MASK;
  const static std::vector<code_type> HIGH_MASK;

  // Constants used in the encoding.
  constexpr static size_t    FIELD_BITS = sizeof(code_type) * gbwt::BYTE_BITS;
  constexpr static size_t    PACK_WIDTH    = 2;
  constexpr static size_t    PACK_OVERFLOW = FIELD_BITS - PACK_WIDTH;
  constexpr static size_t    FIELD_CHARS   = FIELD_BITS / PACK_WIDTH;
  constexpr static code_type PACK_MASK     = 0x3;
};

//------------------------------------------------------------------------------

/*
  A kmer encoded using 2 bits/character in a 64-bit integer. The encoding is only
  defined if all characters in the kmer are valid. The highest bit indicates
  whether the corresponding value in the hash table is a position or a pointer.
*/
struct Key64
{
public:
  // Internal representation.
  typedef KmerEncoding::code_type code_type;
  typedef code_type value_type;
  code_type key;

  // Empty key.
  constexpr Key64() : key(EMPTY_KEY) {}

  // Reads the key starting from the given pointer.
  static Key64 read(const code_type* ptr) { return Key64(*ptr); }

  // Writes the key to the given pointer.
  void write(code_type* ptr) const { *ptr = this->key; }

  // No key, with 63 set bits.
  constexpr static Key64 no_key() { return Key64(NO_KEY & KEY_MASK); }

  // For testing.
  explicit constexpr Key64(code_type key) : key(key) {}

  // Get a representation of the actual key.
  value_type get_key() const { return (this->key & KEY_MASK); }

  // Comparisons.
  bool operator==(Key64 another) const { return (this->get_key() == another.get_key()); }
  bool operator!=(Key64 another) const { return (this->get_key() != another.get_key()); }
  bool operator<(Key64 another) const { return (this->get_key() < another.get_key()); }

  // Is the corresponding value a pointer?
  bool is_pointer() const { return (this->key & IS_POINTER); }
  void clear_pointer() { this->key &= ~IS_POINTER; }
  void set_pointer() { this->key |= IS_POINTER; }

  // Hash of the key. Do not call this directly, as the index may use a derived
  // hash function.
  size_t hash() const { return wang_hash_64(this->key & KEY_MASK); }

  // Move the kmer forward, with c as the next character. Update the key, assuming that
  // it encodes the kmer in forward orientation.
  void forward(size_t k, unsigned char c, size_t& valid_chars)
  {
    code_type packed = KmerEncoding::CHAR_TO_PACK[c];
    if(packed > KmerEncoding::PACK_MASK) { this->key = EMPTY_KEY; valid_chars = 0; }
    else
    {
      this->key = ((this->key << KmerEncoding::PACK_WIDTH) | packed) & KmerEncoding::LOW_MASK[k];
      valid_chars++;
    }
  }

  // Move the kmer forward, with c as the next character. Update the key, assuming that
  // it encodes the kmer in reverse orientation.
  void reverse(size_t k, unsigned char c)
  {
    code_type packed = KmerEncoding::CHAR_TO_PACK[c];
    if(packed > KmerEncoding::PACK_MASK) { this->key = EMPTY_KEY; }
    else
    {
      packed ^= KmerEncoding::PACK_MASK; // The complement of the base.
      this->key = (packed << ((k - 1) * KmerEncoding::PACK_WIDTH)) | (this->key >> KmerEncoding::PACK_WIDTH);
    }
  }

  /// Encode a string of size k to a key.
  static Key64 encode(const std::string& sequence);

  /// Decode the key back to a string, given the kmer size used.
  std::string decode(size_t k) const;

  /// Returns the packed character value for kmer[i].
  code_type access_raw(size_t k, size_t i) const
  {
    size_t offset = (k - 1 - i) * KmerEncoding::PACK_WIDTH;
    return (this->key >> offset) & KmerEncoding::PACK_MASK;
  }

  /// Returns kmer[i].
  char access(size_t k, size_t i) const
  {
    return KmerEncoding::PACK_TO_CHAR[this->access_raw(k, i)];
  }

  /// Returns the reverse complement of the kmer.
  Key64 reverse_complement(size_t k) const;

  // Required numeric constants.
  constexpr static std::size_t KEY_BITS = KmerEncoding::FIELD_BITS;
  constexpr static std::size_t KMER_LENGTH = 29;
  constexpr static std::size_t WINDOW_LENGTH = 11;
  constexpr static std::size_t SMER_LENGTH = KMER_LENGTH - WINDOW_LENGTH;
  constexpr static std::size_t KMER_MAX_LENGTH = 31;

private:
  // Specific key values. Note that the highest bit is not a part of the key.
  constexpr static code_type EMPTY_KEY = 0;
  constexpr static code_type NO_KEY = std::numeric_limits<code_type>::max();
  constexpr static code_type KEY_MASK = NO_KEY >> 1;
  constexpr static code_type IS_POINTER = NO_KEY ^ KEY_MASK; // High bit.
};

// Required for printing keys.
std::ostream& operator<<(std::ostream& out, Key64 value);

//------------------------------------------------------------------------------

/*
  A kmer encoded using 2 bits/character in a pair of 64-bit integers. The encoding is
  only defined if all characters in the kmer are valid. The highest bit indicates
  whether the corresponding value in the hash table is a position or a pointer.
*/
struct Key128
{
public:
  // Internal representation.
  typedef KmerEncoding::code_type code_type;
  typedef std::pair<code_type, code_type> value_type;
  code_type high, low;

  // Empty key.
  constexpr Key128() : high(EMPTY_KEY), low(EMPTY_KEY) {}

  // Reads the key starting from the given pointer.
  static Key128 read(const code_type* ptr) { return Key128(ptr[0], ptr[1]); }

  // Writes the key to the given pointer.
  void write(code_type* ptr) const { ptr[0] = this->high; ptr[1] = this->low; }

  // No key, with 127 set bits.
  constexpr static Key128 no_key() { return Key128(NO_KEY & KEY_MASK, NO_KEY); }

  // For testing.
  explicit constexpr Key128(code_type key) : high(EMPTY_KEY), low(key) {}

  // For testing.
  constexpr Key128(code_type high, code_type low) : high(high), low(low) {}

  // Get a representation of the actual key.
  value_type get_key() const { return value_type(this->high & KEY_MASK, this->low); }

  // Comparisons.
  bool operator==(Key128 another) const { return (this->get_key() == another.get_key()); }
  bool operator!=(Key128 another) const { return (this->get_key() != another.get_key()); }
  bool operator<(Key128 another) const { return (this->get_key() < another.get_key()); }

  // Is the corresponding value a pointer?
  bool is_pointer() const { return (this->high & IS_POINTER); }
  void clear_pointer() { this->high &= ~IS_POINTER; }
  void set_pointer() { this->high |= IS_POINTER; }

  // Hash of the key. Essentially boost::hash_combine. Do not call this directly,
  // as the index may use a derived hash function.
  size_t hash() const
  {
    size_t result = wang_hash_64(this->high & KEY_MASK);
    result ^= wang_hash_64(this->low) + 0x9e3779b9 + (result << 6) + (result >> 2);
    return result;
  }

  // Move the kmer forward, with c as the next character. Update the key, assuming that
  // it encodes the kmer in forward orientation.
  void forward(size_t k, unsigned char c, size_t& valid_chars)
  {
    code_type packed = KmerEncoding::CHAR_TO_PACK[c];
    if(packed > KmerEncoding::PACK_MASK) { this->high = EMPTY_KEY; this->low = EMPTY_KEY; valid_chars = 0; }
    else
    {
      this->high = ((this->high << KmerEncoding::PACK_WIDTH) | (this->low >> KmerEncoding::PACK_OVERFLOW)) & KmerEncoding::HIGH_MASK[k];
      this->low = ((this->low << KmerEncoding::PACK_WIDTH) | packed) & KmerEncoding::LOW_MASK[k];
      valid_chars++;
    }
  }

  // Move the kmer forward, with c as the next character. Update the key, assuming that
  // it encodes the kmer in reverse orientation.
  void reverse(size_t k, unsigned char c)
  {
    code_type packed = KmerEncoding::CHAR_TO_PACK[c];
    if(packed > KmerEncoding::PACK_MASK) { this->high = EMPTY_KEY; this->low = EMPTY_KEY; }
    else
    {
      packed ^= KmerEncoding::PACK_MASK; // The complement of the base.
      if(k > KmerEncoding::FIELD_CHARS)
      {
        this->low = ((this->high & KmerEncoding::PACK_MASK) << KmerEncoding::PACK_OVERFLOW) | (this->low >> KmerEncoding::PACK_WIDTH);
        this->high = (packed << ((k - KmerEncoding::FIELD_CHARS - 1) * KmerEncoding::PACK_WIDTH)) | (this->high >> KmerEncoding::PACK_WIDTH);
      }
      else // The entire kmer is in the lower part of the key.
      {
        this->low = (packed << ((k - 1) * KmerEncoding::PACK_WIDTH)) | (this->low >> KmerEncoding::PACK_WIDTH);
      }
    }
  }

  /// Encode a string of size k to a key.
  static Key128 encode(const std::string& sequence);

  /// Decode the key back to a string, given the kmer size used.
  std::string decode(size_t k) const;

  /// Returns the packed character value for kmer[i].
  code_type access_raw(size_t k, size_t i) const
  {
    size_t offset = (k - 1 - i) * KmerEncoding::PACK_WIDTH;
    code_type shifted;
    if(offset >= KmerEncoding::FIELD_BITS) { shifted = this->high >> (offset - KmerEncoding::FIELD_BITS); }
    else { shifted = this->low >> offset; }
    return shifted & KmerEncoding::PACK_MASK;
  }

  /// Returns kmer[i].
  char access(size_t k, size_t i) const
  {
    return KmerEncoding::PACK_TO_CHAR[this->access_raw(k, i)];
  }

  /// Returns the reverse complement of the kmer.
  Key128 reverse_complement(size_t k) const;

  // Required numeric constants.
  constexpr static std::size_t KEY_BITS = 2 * KmerEncoding::FIELD_BITS;
  constexpr static std::size_t KMER_LENGTH = 39;
  constexpr static std::size_t WINDOW_LENGTH = 15;
  constexpr static std::size_t SMER_LENGTH = KMER_LENGTH - WINDOW_LENGTH;
  constexpr static std::size_t KMER_MAX_LENGTH = 63;

private:
  // Specific key values. Note that the highest bit is not a part of the key.
  constexpr static code_type EMPTY_KEY = 0;
  constexpr static code_type NO_KEY = std::numeric_limits<code_type>::max();
  constexpr static code_type KEY_MASK = NO_KEY >> 1;
  constexpr static code_type IS_POINTER = NO_KEY ^ KEY_MASK; // High bit.
};

// Required for printing keys.
std::ostream& operator<<(std::ostream& out, Key128 value);

//------------------------------------------------------------------------------

// Sequence offset in kmer extraction.
typedef std::uint32_t offset_type;

/*
  The sequence offset of a kmer is the base that corresponds to the start of
  the kmer. If the kmer is a reverse complement, that means the last base of
  the substring covered by the kmer.

  Values stored in a kmer index are graph positions for the start of the kmer
  used as the key. The kmer always extends forward from the position. When a
  kmer is in reverse orientation in the query sequence, the graph positions
  from the index must be flipped with reverse_base_pos(). The last base of the
  forward orientation of the query sequence covered by the kmer then matches
  the flipped graph position.
*/
template<class KeyType>
struct Kmer
{
  typedef KeyType key_type;

  key_type    key;        // Encoded kmer.
  size_t      hash;       // Hash of the kmer.
  offset_type offset;     // Sequence offset.
  bool        is_reverse; // The kmer is a reverse complement.

  // Is the kmer empty?
  bool empty() const { return (this->key == key_type::no_key()); }

  // Sort by (offset, !is_reverse). When the offsets are equal, a reverse complement
  // kmer is earlier in the sequence than a forward kmer.
  bool operator<(const Kmer& another) const
  {
    return ((this->offset < another.offset) ||
            (this->offset == another.offset && this->is_reverse > another.is_reverse));
  }

  bool operator==(const Kmer& another) const
  {
    return (this->key == another.key && this->offset == another.offset && this->is_reverse == another.is_reverse);
  }

  bool operator!=(const Kmer& another) const
  {
    return !(this->operator==(another));
  }
};

template<class KeyType>
std::ostream&
operator<<(std::ostream& out, const Kmer<KeyType>& kmer)
{
  out << "(" << kmer.key << ", " << (kmer.is_reverse ? "-" : "+") << kmer.offset << ")";
  return out;
}

//------------------------------------------------------------------------------

// Forward declarations for friend declarations.

class MinimizerHeader;

template<class KeyType>
class MinimizerIndex;

/*
  This is a general-purpose kmer index parameterized with a key type. The key
  can be either Key64 (64 bits, up to 31 bp) or Key128 (128 bits, up to 63 bp).

  Values consist of a graph position (Position) and payload. Payload size
  (in 64-bit words) is defined when the index is constructed and can be 0.

  The kmers are stored in a hash table with 2^n cells that uses quadratic
  probing. If a kmer is associated with a single value, it is stored directly
  within the hash table. When there are multiple values, the hash table stores
  a pointer to a vector of values. The values are sorted by position.

  Limitations:

  * The index does not support serialization.

  * Any key not equal to KeyType::no_key() can be inserted into the index.
    There are no guarantees that all kmers have the same length.
*/
template<class KeyType>
class KmerIndex
{
public:
  typedef KmerEncoding::code_type code_type;
  typedef KeyType key_type;

  // A cell consists of cell_size() words that encode a key and an
  // encoded value using value_size() words. The encoded value may
  // contain either a pointer to a vector of values or a single value.
  // A value consists of a graph position (Position) and 0 to
  // MAX_PAYLOAD_SIZE words of payload. For any given index, you can
  // get payload size using payload_size().
  typedef std::pair<key_type&, code_type*> cell_type;
  typedef std::pair<key_type, const code_type*> const_cell_type;

  // A value consisting of a position and a pointer to the payload,
  // which is assumed to have the same size as the index payload size.
  typedef KmerEncoding::value_type value_type;

  // A pointer to a list of encoded values and the number of values.
  // Use get_value() to access and decode individual values.
  typedef KmerEncoding::multi_value_type multi_value_type;

  // Public constants.
  constexpr static size_t INITIAL_CAPACITY = 1024;
  constexpr static double MAX_LOAD_FACTOR = 0.77;
  constexpr static size_t MAX_PAYLOAD_SIZE = 15; // Must be consistent with the number of bits reserved in MinimizerHeader.

  constexpr static size_t KEY_SIZE = sizeof(key_type) / sizeof(code_type);
  constexpr static size_t POS_SIZE = sizeof(Position) / sizeof(code_type);

  // Returns the size of the smallest hash table that can hold this many keys.
  static size_t minimum_size(size_t keys)
  {
    size_t bound = std::ceil(keys / MAX_LOAD_FACTOR);
    size_t capacity = INITIAL_CAPACITY;
    while(capacity < bound) { capacity *= 2; }
    return capacity;
  }

//------------------------------------------------------------------------------

  /*
    Constructors, destructors, and object handling.
  */

  // Throws std::runtime_error if payload size exceeds MAX_PAYLOAD_SIZE.
  explicit KmerIndex(size_t payload_size) :
    KmerIndex(INITIAL_CAPACITY, payload_size)
  {
  }

  // Hash table size must be a power of 2 and at least INITIAL_CAPACITY.
  // Otherwise it will revert to INITIAL_CAPACITY.
  // Throws std::runtime_error if payload size exceeds MAX_PAYLOAD_SIZE.
  KmerIndex(size_t cell_count, size_t payload_size) :
    keys(0), values(0), unique(0),
    payload(payload_size),
    downweight(0), frequent_kmers()
  {
    if(sdsl::bits::cnt(cell_count) != 1)
    {
      std::cerr << "KmerIndex::KmerIndex(): Cell count (" << cell_count << ") must be a power of 2; reverting to " << INITIAL_CAPACITY << std::endl;
      cell_count = INITIAL_CAPACITY;
    }
    if(cell_count < INITIAL_CAPACITY)
    {
      std::cerr << "KmerIndex::KmerIndex(): Cell count (" << cell_count << ") is too small; reverting to " << INITIAL_CAPACITY << std::endl;
      cell_count = INITIAL_CAPACITY;
    }
    this->max_keys = cell_count * MAX_LOAD_FACTOR;

    if(payload_size > MAX_PAYLOAD_SIZE)
    {
      std::string msg = "KmerIndex::KmerIndex(): Payload size (" + std::to_string(payload_size) + ") exceeds maximum (" + std::to_string(MAX_PAYLOAD_SIZE) + ")";
      throw std::runtime_error(msg);
    }
    this->hash_table = empty_hash_table(cell_count, this->cell_size());
  }

  KmerIndex(const KmerIndex& source)
  {
    this->copy(source);
  }

  KmerIndex(KmerIndex&& source)
  {
    *this = std::move(source);
  }

  ~KmerIndex()
  {
    this->clear();
  }

  void swap(KmerIndex& another)
  {
    if(&another == this) { return; }
    std::swap(this->keys, another.keys);
    std::swap(this->max_keys, another.max_keys);
    std::swap(this->values, another.values);
    std::swap(this->unique, another.unique);
    std::swap(this->payload, another.payload);
    this->hash_table.swap(another.hash_table);
    std::swap(this->downweight, another.downweight);
    this->frequent_kmers.swap(another.frequent_kmers);
  }

  KmerIndex& operator=(const KmerIndex& source)
  {
    if(&source != this) { this->copy(source); }
    return *this;
  }

  KmerIndex& operator=(KmerIndex&& source)
  {
    if(&source != this)
    {
      this->keys = std::move(source.keys);
      this->max_keys = std::move(source.max_keys);
      this->values = std::move(source.values);
      this->unique = std::move(source.unique);
      this->payload = std::move(source.payload);
      this->hash_table = std::move(source.hash_table);
      this->downweight = std::move(source.downweight);
      this->frequent_kmers = std::move(source.frequent_kmers);
    }
    return *this;
  }

  // For testing.
  bool operator==(const KmerIndex& another) const
  {
    if(this->keys != another.keys) { return false; }
    if(this->max_keys != another.max_keys) { return false; }
    if(this->values != another.values) { return false; }
    if(this->unique != another.unique) { return false; }
    if(this->payload != another.payload) { return false; }

    if(this->hash_table.size() != another.hash_table.size()) { return false; }
    for(size_t array_offset = 0; array_offset < this->hash_table.size(); array_offset += this->cell_size())
    {
      const_cell_type a = this->const_cell(array_offset);
      const_cell_type b = another.const_cell(array_offset);
      if(a.first != b.first) { return false; } // Key comparison.
      if(a.first.is_pointer() != b.first.is_pointer()) { return false; }
      multi_value_type a_values = this->get_values(a);
      multi_value_type b_values = another.get_values(b);
      if(a_values.second != b_values.second) { return false; }
      for(size_t j = 0; j < a_values.second * this->value_size(); j++)
      {
        if(a_values.first[j] != b_values.first[j]) { return false; }
      }
    }

    if(this->downweight != another.downweight) { return false; }
    if(this->frequent_kmers != another.frequent_kmers) { return false; }

    return true;
  }

  bool operator!=(const KmerIndex& another) const { return !(this->operator==(another)); }

//------------------------------------------------------------------------------

  /*
    Statistics and iteration.
  */

  // Number of keys in the index.
  size_t size() const { return this->keys; }

  // Is the index empty?
  bool empty() const { return (this->size() == 0); }

  // Number of values (kmer occurrences) in the index.
  size_t number_of_values() const { return this->values; }

  // Size of the hash table in cells.
  size_t cell_count() const { return this->hash_table.size() / this->cell_size(); }

  // Number of keys that can fit into the hash table. Exceeding it will initiate rehashing.
  size_t capacity() const { return this->max_keys; }

  // Current load factor of the hash table.
  double load_factor() const { return static_cast<double>(this->size()) / static_cast<double>(this->cell_count()); }

  // Number of kmers with a single occurrence.
  size_t unique_keys() const { return this->unique; }

  // Returns the size of a payload in words.
  size_t payload_size() const { return this->payload; }

  // Returns the size of a value in words.
  size_t value_size() const { return POS_SIZE + this->payload_size(); }

  // Returns the size of a cell in words.
  size_t cell_size() const { return KEY_SIZE + POS_SIZE + this->payload_size(); }

  // Call `callback` for every non-empty hash table cell.
  // The callback gets the key and encoded values.
  void for_each_kmer(const std::function<void(key_type, multi_value_type)>& callback) const
  {
    for(size_t array_offset = 0; array_offset < this->hash_table.size(); array_offset += this->cell_size())
    {
      const_cell_type cell = this->const_cell(array_offset);
      if(cell.first == key_type::no_key()) { continue; }
      multi_value_type values = this->get_values(cell);
      callback(cell.first, values);
    }
  }

  // Call `callback` for every non-empty hash table cell.
  // The callback gets the key and encoded values.
  // If callback returns false, then stop iterating.
  // Returns false if the iteration stopped early, true otherwise.
  bool for_each_kmer(const std::function<bool(key_type, multi_value_type)>& callback) const
  {
    for(size_t array_offset = 0; array_offset < this->hash_table.size(); array_offset += this->cell_size())
    {
      const_cell_type cell = this->const_cell(array_offset);
      if(cell.first == key_type::no_key()) { continue; }
      multi_value_type values = this->get_values(cell);
      if(!callback(cell.first, values)) { return false; }
    }
    return true;
  }

//------------------------------------------------------------------------------

  /*
    Operations.
  */

  /*
    Inserts (key, (position, payload)) into the index if it is not already there.
    Does not insert keys equal to key_type::no_key() or positions equal to
    Position::no_pos().
  */
  void insert(key_type key, value_type value)
  {
    this->insert(key, value, this->hash(key));
  }

  /*
    Inserts (key, (position, payload)) into the index if it is not already there
    using the given hash value. Does not insert keys equal to key_type::no_key()
    or positions equal to Position::no_pos().
  */
  void insert(key_type key, value_type value, size_t hash)
  {
    if(key == key_type::no_key() || value.first == Position::no_pos()) { return; }

    size_t array_offset = this->find_offset(key, hash);
    cell_type cell = this->cell(array_offset);
    if(cell.first == key_type::no_key())
    {
      this->insert_new(key, value, cell);
    }
    else if(cell.first == key)
    {
      this->append(value, cell);
    }
  }

  // Returns the occurrence count for the kmer.
  size_t count(key_type key) const
  {
    return this->count(key, this->hash(key));
  }

  // Returns the occurrence count for the kmer with the given hash value.
  size_t count(key_type key, size_t hash) const
  {
    if(key == key_type::no_key()) { return 0; }
    size_t array_offset = this->find_offset(key, hash);
    const_cell_type cell = this->const_cell(array_offset);
    if(cell.first == key)
    {
      multi_value_type values = this->get_values(cell);
      return values.second;
    }
    return 0;
  }

  /*
    Returns the sorted list of the occurrences for the given kmer and the
    number of occurrences. Any insertions into the index may invalidate the
    returned pointer.
  */
  multi_value_type find(key_type key) const
  {
    return this->find(key, this->hash(key));
  }

  /*
    Returns the sorted list of the occurrences for the kmer with the given hash
    value and the number of occurrences. Any insertions into the index may
    invalidate the returned pointer.
  */
  multi_value_type find(key_type key, size_t hash) const
  {
    if(key == key_type::no_key()) { return multi_value_type(nullptr, 0); }

    size_t array_offset = this->find_offset(key, hash);
    const_cell_type cell = this->const_cell(array_offset);
    if(cell.first == key)
    {
      return this->get_values(cell);
    }
    return multi_value_type(nullptr, 0);
  }

  // FIXME: add to tests
  /*
    Decodes the given value from the given list of values. Returns an empty
    value with no position and an empty payload if the index is out of bounds.
  */
  value_type get_value(multi_value_type values, size_t index) const
  {
    if(index >= values.second) { return std::make_pair(Position::no_pos(), nullptr); }
    const code_type* ptr = values.first + index * this->value_size();
    Position pos(*ptr);
    const code_type* payload = ptr + POS_SIZE;
    return std::make_pair(pos, payload);
  }

//------------------------------------------------------------------------------

  /*
    Internal implementation.
  */

private:
  size_t keys, max_keys;
  size_t values, unique;
  size_t payload; // Payload size in words.
  std::vector<code_type> hash_table;

  // Downweight hashes for frequent kmers. Note that frequent_kmers is a hash set similar to hash_table.
  size_t downweight;
  std::vector<key_type> frequent_kmers;

  // Needed for serialization.
  friend class MinimizerHeader;
  friend class MinimizerIndex<KeyType>;

  size_t serialize_hash_table(std::ostream& out) const
  {
    size_t words_in_block = io::BLOCK_SIZE * this->cell_size();
    size_t cell_size = this->cell_size();

    size_t bytes = 0;
    bytes += io::serialize_size(out, this->hash_table);

    // Data in blocks of BLOCK_SIZE cells. Replace pointers with empty values
    // to ensure that the serialization is deterministic.
    for(size_t array_offset = 0; array_offset < this->hash_table.size(); array_offset += words_in_block)
    {
      size_t current_block_size = std::min(this->hash_table.size() - array_offset, words_in_block);
      size_t byte_size = current_block_size * sizeof(code_type);
      std::vector<code_type> buffer(this->hash_table.begin() + array_offset, this->hash_table.begin() + array_offset + current_block_size);
      for(size_t buffer_offset = 0; buffer_offset < buffer.size(); buffer_offset += cell_size)
      {
        cell_type cell = hash_table_cell(buffer, buffer_offset);
        if(cell.first.is_pointer()) { this->clear_value(cell, false); }
      }
      out.write(reinterpret_cast<const char*>(buffer.data()), byte_size);
      bytes += byte_size;
    }

    return bytes;
  }

  static std::vector<code_type> empty_hash_table(size_t capacity, size_t cell_size)
  {
    std::vector<code_type> result = std::vector<code_type>(capacity * cell_size, 0);
    key_type no_key = key_type::no_key();
    Position no_pos = Position::no_pos();
    for(size_t array_offset = 0; array_offset < result.size(); array_offset += cell_size)
    {
      no_key.write(&result[array_offset]);
      no_pos.write(&result[array_offset + KEY_SIZE]);
    }
    return result;
  }

  void copy(const KmerIndex& source)
  {
    this->clear();
    this->keys = source.keys;
    this->max_keys = source.max_keys;
    this->values = source.values;
    this->unique = source.unique;
    this->payload = source.payload;

    // First make a shallow copy and then copy the occurrence lists.
    this->hash_table = source.hash_table;
    for(size_t array_offset = 0; array_offset < this->hash_table.size(); array_offset += this->cell_size())
    {
      cell_type cell = this->cell(array_offset);
      std::vector<code_type>* old = get_vector_pointer(cell);
      if(old != nullptr)
      {
        this->write_vector_pointer(new std::vector<code_type>(*old), cell);
      }
    }

    this->downweight = source.downweight;
    this->frequent_kmers = source.frequent_kmers;
  }

  // Clears the value stored in the given cell.
  // In the destructive mode, this also deletes any pointer stored in the cell
  // and clears the pointer bit in the key.
  void clear_value(cell_type cell, bool destructive) const
  {
    if(destructive && cell.first.is_pointer())
    {
      std::vector<code_type>* ptr = get_vector_pointer(cell);
      delete ptr;
      cell.first.clear_pointer();
    }
    Position::no_pos().write(cell.second);
    for(size_t value_offset = 1; value_offset < this->value_size(); value_offset++) { cell.second[value_offset] = 0; }
  }

  // Delete all pointers in the hash table and replace them with empty values.
  void clear()
  {
    for(size_t array_offset = 0; array_offset < this->hash_table.size(); array_offset += this->cell_size())
    {
      cell_type cell = this->cell(array_offset);
      this->clear_value(cell, true);
    }
  }

  // Returns a mutable view of the cell starting at the given array offset
  // in the given array.
  static cell_type hash_table_cell(std::vector<code_type>& array, size_t array_offset)
  {
    key_type* key = reinterpret_cast<key_type*>(&array[array_offset]);
    code_type* values = &array[array_offset + KEY_SIZE];
    return cell_type(*key, values);
  }

  // Returns a mutable view of the cell starting at the given array offset.
  cell_type cell(size_t array_offset)
  {
    return hash_table_cell(this->hash_table, array_offset);
  }

  // Returns an immutable view of the cell starting at the given array offset
  // in the given array.
  static const_cell_type hash_table_const_cell(const std::vector<code_type>& array, size_t array_offset)
  {
    key_type key = key_type::read(&array[array_offset]);
    const code_type* values = &array[array_offset + KEY_SIZE];
    return const_cell_type(key, values);
  }

  // Returns an immutable view of the cell starting at the given array offset.
  const_cell_type const_cell(size_t array_offset) const
  {
    return hash_table_const_cell(this->hash_table, array_offset);
  }

  // Converts a mutable cell to an immutable cell.
  static const_cell_type const_cell(cell_type cell)
  {
    return const_cell_type(cell.first, cell.second);
  }

  // Returns the values stored in given cell. The cell is assumed to be non-empty.
  multi_value_type get_values(const_cell_type cell) const
  {
    if(cell.first.is_pointer())
    {
      const std::vector<code_type>* ptr = reinterpret_cast<const std::vector<code_type>*>(cell.second[0]);
      return multi_value_type(ptr->data(), ptr->size() / this->value_size());
    }
    else { return multi_value_type(cell.second, 1); }
  }

  // Returns the vector storing the values for the hash table cell starting at the
  // given offset of the underlying array. This returns nullptr if the cell is empty
  // or if it contains a single value.
  static std::vector<code_type>* get_vector_pointer(cell_type cell)
  {
    if(!cell.first.is_pointer()) { return nullptr; }
    return reinterpret_cast<std::vector<code_type>*>(cell.second[0]);
  }

  // Returns the vector storing the values for the hash table cell starting at the
  // given offset of the underlying array. This returns nullptr if the cell is empty
  // or if it contains a single value.
  static const std::vector<code_type>* get_vector_pointer(const_cell_type cell)
  {
    if(!cell.first.is_pointer()) { return nullptr; }
    return reinterpret_cast<const std::vector<code_type>*>(cell.second[0]);
  }

  // Find the array offset for the key with the given hash value.
  // The cell may either be empty or contain the key.
  size_t find_offset(key_type key, size_t hash) const
  {
    size_t cell_count = this->cell_count();
    size_t cell_offset = hash & (cell_count - 1);
    size_t array_offset = cell_offset * this->cell_size();
    for(size_t attempt = 0; attempt < cell_count; attempt++)
    {
      const_cell_type cell = this->const_cell(array_offset);
      if(cell.first == key_type::no_key() || cell.first == key) { return array_offset; }

      // Quadratic probing with triangular numbers.
      cell_offset = (cell_offset + attempt + 1) & (cell_count - 1);
      array_offset = cell_offset * this->cell_size();
    }

    // This should not happen.
    std::cerr << "KmerIndex::find_offset(): Cannot find the offset for key " << key << std::endl;
    assert(false);
    return 0;
  }

  // Encodes the given value to the given array pointer.
  void encode_value(value_type value, code_type* array_ptr) const
  {
    value.first.write(array_ptr);
    for(size_t i = 0; i < this->payload_size(); i++) { array_ptr[i + POS_SIZE] = value.second[i]; }
  }

  // Stores the given vector pointer in the given cell. Assumes that
  // the cell does not already contain a pointer.
  void write_vector_pointer(std::vector<code_type>* ptr, cell_type cell) const
  {
    static_assert(sizeof(std::vector<code_type>*) == sizeof(code_type), "Pointer size does not match code_type size.");
    cell.second[0] = reinterpret_cast<code_type>(ptr);
    for(size_t i = 1; i < this->value_size(); i++) { cell.second[i] = 0; }
  }

  // Swaps the values starting at the given pointers, assuming that the values do not overlap.
  void swap_values(code_type* a, code_type* b) const
  {
    for(size_t i = 0; i < this->value_size(); i++) { std::swap(a[i], b[i]); }
  }

  // Encodes the given value, appends it to the given vector, and bubbles it down to its correct position.
  void insert_value(value_type value, std::vector<code_type>& vec) const
  {
    size_t value_size = this->value_size();
    assert(vec.size() % value_size == 0);
    size_t curr_offset = vec.size();
    vec.resize(vec.size() + value_size, 0);
    this->encode_value(value, vec.data() + curr_offset);

    while(curr_offset >= value_size)
    {
      size_t prev_offset = curr_offset - value_size;
      Position candidate(vec[prev_offset]);
      if(candidate > value.first)
      {
        this->swap_values(vec.data() + prev_offset, vec.data() + curr_offset);
      } else { break; }
      curr_offset = prev_offset;
    }
  }

  // Insert (key, value) the given cell, which is assumed to be empty.
  // Rehashing may be necessary, which invalidates the cell.
  void insert_new(key_type key, value_type value, cell_type cell)
  {
    cell.first = key;
    this->encode_value(value, cell.second);
    this->keys++;
    this->values++;
    this->unique++;

    if(this->size() > this->capacity()) { this->rehash(); }
  }

  // Inserts the given value into the given cell, which is assumed to
  // already contain a key.
  void append(value_type value, cell_type cell)
  {
    if(this->contains(const_cell(cell), value.first)) { return; }

    std::vector<code_type>* ptr = get_vector_pointer(cell);
    if(ptr != nullptr)
    {
      // We already have multiple occurrences.
      this->insert_value(value, *ptr);
    }
    else
    {
      // Now we have two occurrences of this kmer.
      std::vector<code_type>* occs = new std::vector<code_type>(cell.second, cell.second + this->value_size());
      this->insert_value(value, *occs);
      this->write_vector_pointer(occs, cell);
      cell.first.set_pointer();
      this->unique--;
    }
    this->values++;
  }

  // Returns true if the given cell contains the given graph position.
  bool contains(const_cell_type cell, Position pos) const
  {
    multi_value_type values = this->get_values(cell);
    size_t low = 0, high = values.second;
    while(low < high)
    {
      size_t mid = low + (high - low) / 2;
      size_t array_offset = mid * this->value_size();
      Position candidate(values.first[array_offset]);
      if(candidate < pos) { low = mid + 1; }
      else if(candidate > pos) { high = mid; }
      else { return true; }
    }
    return false;
  }

  // Get the hash value for the key. If downweighting is not in use, or
  // the kmer is not frequent, the hash value reported by the key is used
  // directly. Otherwise we downweight the hash value by the specified
  // number of iterations.
  size_t hash(key_type key) const
  {
    if(this->frequent_kmers.empty()) { return key.hash(); }

    size_t h = key.hash();
    size_t offset = h & (this->frequent_kmers.size() - 1);
    for(size_t attempt = 0; attempt < this->frequent_kmers.size(); attempt++)
    {
      if(this->frequent_kmers[offset] == key_type::no_key()) { break; }
      if(this->frequent_kmers[offset] == key)
      {
        h = std::numeric_limits<size_t>::max() - h;
        for(size_t i = 0; i < this->downweight; i++)
        {
          h = h >> 32; h = h * h;
        }
        return std::numeric_limits<size_t>::max() - h;
      }
      // Quadratic probing with triangular numbers.
      offset = (offset + attempt + 1) & (this->frequent_kmers.size() - 1);
    }

    return h;
  }

  // Double the size of the hash table.
  void rehash()
  {
    // Reinitialize with a larger hash table.
    size_t cell_size = this->cell_size();
    std::vector<code_type> old_hash_table = empty_hash_table(this->cell_count() * 2, cell_size);
    this->hash_table.swap(old_hash_table);
    this->max_keys = this->cell_count() * MAX_LOAD_FACTOR;

    // Move the keys to the new hash table.
    for(size_t old_array_offset = 0; old_array_offset < old_hash_table.size(); old_array_offset += cell_size)
    {
      cell_type source = hash_table_cell(old_hash_table, old_array_offset);
      if(source.first == key_type::no_key()) { continue; }

      size_t offset = this->find_offset(source.first, this->hash(source.first));
      cell_type target = this->cell(offset);
      assert(target.first == key_type::no_key());
      target.first = source.first;
      this->swap_values(target.second, source.second);
    }
  }

  // Inserts a key into the frequent kmers hash table.
  void insert_frequent(key_type key)
  {
    size_t h = key.hash();
    size_t offset = h & (this->frequent_kmers.size() - 1);
    for(size_t attempt = 0; attempt < this->frequent_kmers.size(); attempt++)
    {
      if(this->frequent_kmers[offset] == key) { return; }
      if(this->frequent_kmers[offset] == key_type::no_key())
      {
        this->frequent_kmers[offset] = key;
        return;
      }
      // Quadratic probing with triangular numbers.
      offset = (offset + attempt + 1) & (this->frequent_kmers.size() - 1);
    }

    // This should not happen.
    std::cerr << "KmerIndex::insert_frequent(): Cannot find the offset for key " << key << std::endl;
    assert(false);
  }
};

//------------------------------------------------------------------------------

struct MinimizerHeader
{
  std::uint32_t tag, version;
  std::uint64_t k, w_or_s;    // Minimizer parameters.
  std::uint64_t keys;         // Number of keys in the hash table.
  std::uint64_t unused;       // Currently unused.
  std::uint64_t capacity;     // Number of keys that can fit in the hash table without initiating rehashing.
  std::uint64_t values;       // Number of values in the index.
  std::uint64_t unique;       // Number of keys with a single value.
  std::uint64_t flags;

  constexpr static std::uint32_t TAG = 0x31513151;
  constexpr static std::uint32_t VERSION = Version::MINIMIZER_VERSION;

  constexpr static std::uint64_t FLAG_MASK           = 0xFFFF;
  constexpr static std::uint64_t FLAG_KEY_MASK       = 0x00FF;
  constexpr static size_t        FLAG_KEY_OFFSET     = 0;
  constexpr static std::uint64_t FLAG_SYNCMERS       = 0x0100;
  constexpr static std::uint64_t FLAG_WEIGHT_MASK    = 0x0E00;
  constexpr static size_t        FLAG_WEIGHT_OFFSET  = 9;
  constexpr static std::uint64_t FLAG_PAYLOAD_MASK   = 0xF000;
  constexpr static size_t        FLAG_PAYLOAD_OFFSET = 12;

  MinimizerHeader();
  MinimizerHeader(size_t kmer_length, size_t window_length, size_t key_bits, size_t payload_size);
  void sanitize(size_t kmer_max_length);

  // Throws `sdsl::simple_sds::InvalidData` if the header is invalid.
  void check() const;

  void update_version();

  // Boolean flags.
  void set(std::uint64_t flag) { this->flags |= flag; }
  void unset(std::uint64_t flag) { this->flags &= ~flag; }
  bool get_flag(std::uint64_t flag) const { return (this->flags & flag); }

  // Integer flags.
  void set_int(std::uint64_t mask, size_t offset, size_t value);
  size_t get_int(std::uint64_t mask, size_t offset) const;

  size_t key_bits() const;

  // If weighted minimizers are in use, deprioritize frequent kmers by
  // this many iterations.
  size_t downweight() const;

  // Returns payload size in 64-bit words.
  size_t payload_size() const;

  // These do not compare the statistics fields, because the actual values are
  // now stored within the kmer index.
  bool operator==(const MinimizerHeader& another) const;
  bool operator!=(const MinimizerHeader& another) const { return !(this->operator==(another)); }

  // Read the statistics from the kmer index.
  template<class KeyType>
  void set_statistics(const KmerIndex<KeyType>& index)
  {
    this->keys = index.size();
    this->unused = 0;
    this->capacity = index.capacity();
    this->values = index.number_of_values();
    this->unique = index.unique_keys();
    this->set_int(FLAG_WEIGHT_MASK, FLAG_WEIGHT_OFFSET, index.downweight);
    this->set_int(FLAG_PAYLOAD_MASK, FLAG_PAYLOAD_OFFSET, index.payload_size());
  }

  // Write the statistics into the kmer index.
  template<class KeyType>
  void fill_statistics(KmerIndex<KeyType>& index) const
  {
    index.keys = this->keys;
    index.max_keys = this->capacity;
    index.values = this->values;
    index.unique = this->unique;
    index.downweight = this->downweight();
    index.payload = this->payload_size();
  }
};

//------------------------------------------------------------------------------

// FIXME: operations for tags
/*
  A class that implements the minimizer index as a hash table mapping kmers to sets of pos_t.
  For each stored position, we also store 64 bits of payload for external purposes.
  The hash table uses quadratic probing with power-of-two size.
  We encode kmers using 2 bits/character and hash the encoding. A minimizer is the kmer with
  the smallest hash in a window of w consecutive kmers and their reverse complements.

  The index can also use closed syncmers:

    Edgar: Syncmers are more sensitive than minimizers for selecting conserved k-mers in
    biological sequences. PeerJ 9:e10805, 2021.

  A closed syncmer is a kmer where one of the occurrences of the smer (or its reverse
  complement) with the smallest hash is the first or the last. Because we want to select
  the same kmers in both orientations, we use the kmer or its reverse complement,
  depending on which has the smaller hash.

  Minimizers and closed syncmers should have roughly the same seed density when w = k - s.

  There is also an option to use weighted minimizers:

    Jain, Rhie, Zhang, Chu, Walenz, Koren, and Philippy: Weighted minimizer sampling improves
    long read mapping. Bioinformatics, 2020.

  Normally, a minimizer is the kmer with the smallest hash value. With weighted minimizers,
  we want to discourage selecting some (frequent) kmers as minimizers. Let x be the hash
  value interpreted as a number in [0.0, 1.0]. If the kmer is listed as one to be avoided,
  we downweight it by replacing its hash value with 1.0 - (1.0 - x)^(2^iterations), where
  iterations is a parameter between 1 and 7 (default 3).

  Index versions (this should be in the wiki):

    1  The initial version.

    2  Minimizer selection is based on hashes instead of lexicographic order. A sequence and
       its reverse complement have the same minimizers, reducing index size by 50%. Not
       compatible with version 1.

    3  Construction-time hit cap is no longer used. Compatible with version 2.

    4  Use SDSL data structures where appropriate. Not compatible with version 3.

    5  An option to use 128-bit keys. Compatible with version 4.

    6  Store the position/pointer bit in the key instead of a separate bit vector. Store
       64 bits of payload for each position. Not compatible with version 5.

    7  Option to use closed syncmers instead of minimizers. Compatible with version 6.

    8  Payload is now 128 bits per position. Not compatible with earlier versions.

    9  Option to provide a set of frequent kmers that should be avoided as minimizers.
       The capacity field in the header is no longer in use. Compatible with version 8.

    10 Replace the old distance index payload with zipcodes.
       Not compatible with earlier versions.

    11 The index is no longer parameterized by the value type. Instead, the value consists
       of a position and 0-15 words of payload. Payload size is specified when creating
       the index. There is also a tag (key-value) structure that can be used for storing
       the type of payload. Not compatible with earlier versions.
*/

template<class KeyType>
class MinimizerIndex
{
public:
  typedef KeyType key_type;
  typedef typename KeyType::code_type code_type;
  typedef typename KmerIndex<key_type>::cell_type cell_type;
  typedef typename KmerIndex<key_type>::const_cell_type const_cell_type;
  typedef typename KmerIndex<key_type>::value_type value_type;
  typedef typename KmerIndex<key_type>::multi_value_type multi_value_type;
  typedef Kmer<key_type> minimizer_type;

  const static std::string EXTENSION; // ".min"

//------------------------------------------------------------------------------

  explicit MinimizerIndex(size_t payload_size, bool use_syncmers = false) :
    header(KeyType::KMER_LENGTH, (use_syncmers ? KeyType::SMER_LENGTH : KeyType::WINDOW_LENGTH), KeyType::KEY_BITS, payload_size),
    tags(), index(payload_size)
  {
    if(use_syncmers) { this->header.set(MinimizerHeader::FLAG_SYNCMERS); }
    this->header.sanitize(KeyType::KMER_MAX_LENGTH);
  }

  MinimizerIndex(size_t kmer_length, size_t window_or_smer_length, size_t payload_size, bool use_syncmers = false) :
    header(kmer_length, window_or_smer_length, KeyType::KEY_BITS, payload_size),
    tags(), index(payload_size)
  {
    if(use_syncmers) { this->header.set(MinimizerHeader::FLAG_SYNCMERS); }
    this->header.sanitize(KeyType::KMER_MAX_LENGTH);
  }

  MinimizerIndex(const MinimizerIndex& source) :
    index(source.payload_size())
  {
    this->copy(source);
  }

  MinimizerIndex(MinimizerIndex&& source) :
    index(source.payload_size())
  {
    *this = std::move(source);
  }

  void swap(MinimizerIndex& another)
  {
    if(&another == this) { return; }
    std::swap(this->header, another.header);
    this->tags.swap(another.tags);
    this->index.swap(another.index);
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
      this->tags = std::move(source.tags);
      this->index = std::move(source.index);
    }
    return *this;
  }

  // Serializes the index to the ostream. Returns the number of bytes written.
  size_t serialize(std::ostream& out) const
  {
    size_t bytes = 0;

    // Get statistics from the index and serialize the header.
    MinimizerHeader copy = this->header;
    copy.set_statistics(this->index);
    bytes += io::serialize(out, copy);

    // Serialize the tags in the SDSL format.
    bytes += this->tags.serialize(out);

    // Serialize the hash table and the occurrence lists.
    bytes += this->index.serialize_hash_table(out);

    // Serialize the occurrence lists.
    for(size_t array_offset = 0; array_offset < this->index.hash_table.size(); array_offset += this->index.cell_size())
    {
      const_cell_type cell = this->index.const_cell(array_offset);
      const std::vector<code_type>* ptr = KmerIndex<key_type>::get_vector_pointer(cell);
      if(ptr != nullptr)
      {
        bytes += io::serialize_vector(out, *ptr);
      }
    }

    // Serialize the frequent kmers.
    bytes += io::serialize_vector(out, this->index.frequent_kmers);

    return bytes;
  }

  // Load the index from the istream.
  // Throws sdsl::simple_sds::InvalidData if sanity checks fail.
  void deserialize(std::istream& in)
  {
    // Get rid of the existing pointers in the hash table.
    this->index.clear();

    // Load and check the header and fill in the statistics in the kmer index.
    io::load(in, this->header);
    this->header.check();
    if(this->header.key_bits() != KeyType::KEY_BITS)
    {
      std::string msg = "MinimizerIndex::deserialize(): Expected " + std::to_string(KeyType::KEY_BITS) + "-bit keys, got " + std::to_string(this->header.key_bits()) + "-bit keys";
      throw sdsl::simple_sds::InvalidData(msg);
    }
    this->header.update_version();
    this->header.fill_statistics(this->index);

    // Load the tags in the SDSL format.
    this->tags.load(in);

    // Load the hash table and the occurrence lists.
    io::load_vector(in, this->index.hash_table);
    size_t cell_size = this->index.cell_size();
    for(size_t array_offset = 0; array_offset < this->index.hash_table.size(); array_offset += cell_size)
    {
      cell_type cell = this->index.cell(array_offset);
      if(cell.first.is_pointer())
      {
        std::vector<code_type>* ptr = new std::vector<code_type>();
        io::load_vector(in, *ptr);
        this->index.write_vector_pointer(ptr, cell);
      }
    }

    // Load the frequent kmers.
    io::load_vector(in, this->index.frequent_kmers);
  }

  // For testing.
  bool operator==(const MinimizerIndex& another) const
  {
    if(this->header != another.header) { return false; }
    if(this->tags != another.tags) { return false; }
    return (this->index == another.index);
  }

  bool operator!=(const MinimizerIndex& another) const { return !(this->operator==(another)); }

//------------------------------------------------------------------------------

  /*
    A circular buffer of size 2^i for all minimizer candidates. The candidates are sorted
    by both key and sequence offset. The candidate at offset j is removed when we reach
    offset j + w. Candidates at the tail are removed when we advance with a smaller hash.
  */
  struct CircularBuffer
  {
    const KmerIndex<key_type>& parent;
    std::vector<minimizer_type> buffer;
    size_t head, tail;
    size_t w;

    constexpr static size_t BUFFER_SIZE = 16;

    CircularBuffer(const KmerIndex<key_type>& parent, size_t capacity) :
      parent(parent), buffer(),
      head(0), tail(0), w(capacity)
    {
      size_t buffer_size = BUFFER_SIZE;
      while(buffer_size < this->w) { buffer_size *= 2; }
      this->buffer.resize(buffer_size);
    }

    bool empty() const { return (this->head >= this->tail); }
    size_t size() const { return this->tail - this->head; }

    size_t begin() const { return this->head; }
    size_t end() const { return this->tail; }
    minimizer_type& at(size_t i) { return this->buffer[i & (this->buffer.size() - 1)]; }
    minimizer_type& front() { return this->at(this->begin()); }
    minimizer_type& back() { return this->at(this->end() - 1); }

    // Advance to the next offset (pos) with a valid kmer.
    void advance(offset_type pos, key_type forward_key, key_type reverse_key)
    {
      if(!(this->empty()) && this->front().offset + this->w <= pos) { this->head++; }
      size_t forward_hash = this->parent.hash(forward_key), reverse_hash = this->parent.hash(reverse_key);
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

//------------------------------------------------------------------------------

  /*
    Finding minimizers.
  */

  /*
    Returns all minimizers in the string specified by the iterators. The return
    value is a vector of minimizers sorted by their offsets. If there are multiple
    occurrences of one or more minimizer keys with the same hash in a window,
    return all of them.

    Calls syncmers() if the index uses closed syncmers.
  */
  std::vector<minimizer_type> minimizers(std::string::const_iterator begin, std::string::const_iterator end) const
  {
    if(this->uses_syncmers()) { return this->syncmers(begin, end); }
    std::vector<minimizer_type> result;
    size_t window_length = this->window_bp(), total_length = end - begin;
    if(total_length < window_length) { return result; }

    // Find the minimizers.
    CircularBuffer buffer(this->index, this->w());
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
        // 2) the last reported minimizer had the same hash (we may have new occurrences); or
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

    Calls syncmers() if the index uses closed syncmers.
  */
  std::vector<minimizer_type> minimizers(const std::string& str) const
  {
    return this->minimizers(str.begin(), str.end());
  }
  
  /*
    Returns all minimizers in the string specified by the iterators, together
    with the start and length of the run of windows they arise from. The return
    value is a vector of tuples of minimizers, starts, and lengths, sorted by
    minimizer offset.

    Calls syncmers() if the index uses closed syncmers but leaves the start
    and length fields empty.
  */
  std::vector<std::tuple<minimizer_type, size_t, size_t>> minimizer_regions(std::string::const_iterator begin, std::string::const_iterator end) const
  {
    std::vector<std::tuple<minimizer_type, size_t, size_t>> result;
    if(this->uses_syncmers())
    {
      std::vector<minimizer_type> res = this->syncmers(begin, end);
      result.reserve(res.size());
      for(const minimizer_type& m : res) { result.emplace_back(m, 0, 0); }
      return result;
    }
    size_t window_length = this->window_bp(), total_length = end - begin;
    if(total_length < window_length) { return result; }
    
    // Find the minimizers.
    CircularBuffer buffer(this->index, this->w());
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
      
        // Finish off end positions for results that weren't replaced but are going out of range
        while(finished_through < result.size() &&
          std::get<0>(result[finished_through]).offset < window_start)
        {
          // Compute region length based on it stopping at the previous step
          std::get<2>(result[finished_through]) = prev_past_end_pos - std::get<1>(result[finished_through]);
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
              }
            }
            
            // If new minimizers beat out old ones, finish off the old ones.
            while(!result.empty() &&
              finished_through < result.size() &&
              std::get<0>(result.back()).hash != std::get<0>(result[finished_through]).hash)
            {
              // The window before the one we are looking at was the last one for this minimizer.
              std::get<2>(result[finished_through]) = prev_past_end_pos - std::get<1>(result[finished_through]);
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

    Calls syncmers() if the index uses closed syncmers but leaves the start
    and length fields empty.
  */
  std::vector<std::tuple<minimizer_type, size_t, size_t>> minimizer_regions(const std::string& str) const
  {
    return this->minimizer_regions(str.begin(), str.end());
  }

//------------------------------------------------------------------------------

  /*
    Finding syncmers.
  */

  /*
    Returns all closed syncmers in the string specified by the iterators. The return
    value is a vector of minimizers sorted by their offsets.

    Calls minimizers() if the index uses minimizers.
  */
  std::vector<minimizer_type> syncmers(std::string::const_iterator begin, std::string::const_iterator end) const
  {
    if(!(this->uses_syncmers())) { return this->minimizers(begin, end); }
    std::vector<minimizer_type> result;
    size_t total_length = end - begin;
    if(total_length < this->k()) { return result; }

    // Find the closed syncmers.
    CircularBuffer buffer(this->index, this->k() + 1 - this->s());
    size_t processed_chars = 0, dummy_valid_chars = 0, valid_chars = 0, smer_start = 0;
    key_type forward_smer, reverse_smer, forward_kmer, reverse_kmer;
    std::string::const_iterator iter = begin;
    while(iter != end)
    {
      forward_smer.forward(this->s(), *iter, valid_chars);
      reverse_smer.reverse(this->s(), *iter);
      forward_kmer.forward(this->k(), *iter, dummy_valid_chars);
      reverse_kmer.reverse(this->k(), *iter);
      if(valid_chars >= this->s()) { buffer.advance(smer_start, forward_smer, reverse_smer); }
      else                         { buffer.advance(smer_start); }
      ++iter; processed_chars++;
      if(processed_chars >= this->s()) { smer_start++; }
      // We have a full kmer with a closed syncmer.
      if(valid_chars >= this->k())
      {
        // Insert the kmer if the first or the last smer is among the smallest, i.e. if
        // 1) the first smer in the buffer is at the start of the kmer;
        // 2) the first smer in the buffer is at the end of the kmer; or
        // 3) the last smer in the buffer is at the end of the kmer and shares shares
        //    the hash with the first smer.
        if(buffer.front().offset == processed_chars - this->k() ||
          buffer.front().offset == smer_start - 1 ||
          (buffer.back().offset == smer_start - 1 && buffer.back().hash == buffer.front().hash))
        {
          size_t forward_hash = this->index.hash(forward_kmer), reverse_hash = this->index.hash(reverse_kmer);
          offset_type pos = processed_chars - this->k();
          if(reverse_hash < forward_hash)
          {
            result.push_back({ reverse_kmer, reverse_hash, pos, true });
          }
          else
          {
            result.push_back({ forward_kmer, forward_hash, pos, false });
          }
        }
      }
    }

    // It was more convenient to use the first offset of the kmer, regardless of the orientation.
    // If the closed syncmer is a reverse complement, we must return the last offset instead.
    for(minimizer_type& minimizer : result)
    {
      if(minimizer.is_reverse) { minimizer.offset += this->k() - 1; }
    }
    std::sort(result.begin(), result.end());

    return result;
  }

  /*
    Returns all closed syncmers in the string. The return value is a vector of minimizers sorted
    by their offsets.

    Calls minimizers() if the index uses minimizers.
  */
  std::vector<minimizer_type> syncmers(const std::string& str) const
  {
    return this->syncmers(str.begin(), str.end());
  }

//------------------------------------------------------------------------------

  /*
    Operations.
  */

  /*
    Adds a set of frequent kmers that should be avoided (in both orientations)
    when using weighted minimizers. The index must be empty. Clears the
    frequent kmers if the set of keys is empty or the number of iterations is 0.

    Throws std::runtime_error on failure.
  */
  void add_frequent_kmers(const std::vector<key_type>& kmers, size_t iterations)
  {
    if(!this->empty())
    {
      throw std::runtime_error("MinimizerIndex::add_frequent_kmers(): The index is not empty");
    }
    if(this->uses_syncmers())
    {
      throw std::runtime_error("MinimizerIndex::add_frequent_kmers(): Cannot use frequent kmers with syncmers");
    }
    size_t max_iterations = MinimizerHeader::FLAG_WEIGHT_MASK >> MinimizerHeader::FLAG_WEIGHT_OFFSET;
    if(iterations > max_iterations)
    {
      std::string msg = "MinimizerIndex::add_frequent_kmers(): Number of iterations (" + std::to_string(iterations) + ") must be at most " + std::to_string(max_iterations);
      throw std::runtime_error(msg);
    }
    if(kmers.empty() || iterations == 0)
    {
      this->index.frequent_kmers.clear();
      this->index.downweight = 0;
      this->header.set_int(MinimizerHeader::FLAG_WEIGHT_MASK, MinimizerHeader::FLAG_WEIGHT_OFFSET, 0);
    }

    // We oversize the hash table a bit to make it more likely to find an empty
    // cell when querying with a non-frequent kmer.
    size_t capacity = KmerIndex<key_type>::minimum_size(4 * kmers.size());
    this->index.frequent_kmers = std::vector<key_type>(capacity, key_type::no_key());
    for(key_type key : kmers)
    {
      this->index.insert_frequent(key);
      this->index.insert_frequent(key.reverse_complement(this->k()));
    }
    this->index.downweight = iterations;
    this->header.set_int(MinimizerHeader::FLAG_WEIGHT_MASK, MinimizerHeader::FLAG_WEIGHT_OFFSET, iterations);
  }

  /*
    Inserts the value into the index, using minimizer.key as the key and
    minimizer.hash as its hash. Does not insert keys equal to key_type::no_key()
    or positions equal to Position::no_pos(). Does not update the payload if the
    value has already been inserted with the same key.
    Use minimizer() or minimizers() to get the minimizer.
  */
  void insert(const minimizer_type& minimizer, value_type value)
  {
    this->index.insert(minimizer.key, value, minimizer.hash);
  }

  /*
    Returns the occurrence count of the minimizer.
    Use minimizer() or minimizers() to get the minimizer.
  */
  size_t count(const minimizer_type& minimizer) const
  {
    return this->index.count(minimizer.key, minimizer.hash);
  }

  /*
    Returns a pointer to the internal representation of the occurrences of the
    given minimizer and the occurrence count. The occurrences are in sorted
    order. The pointer may be invalidated if new positions are inserted into the
    index. Use minimizer() or minimizers() to get the minimizer and
    position.decode() to decode the positions within the occurrences.
    If the minimizer is in reverse orientation, use reverse_base_pos() to reverse
    the reported occurrences.
  */
  multi_value_type find(const minimizer_type& minimizer) const
  {
    return this->index.find(minimizer.key, minimizer.hash);
  }

  // FIXME: add to tests
  /*
    Decodes the given value from the given list of values. Returns an empty
    value with no position and an empty payload if the index is out of bounds.
  */
  value_type get_value(multi_value_type values, size_t index) const
  {
    return this->index.get_value(values, index);
  }

//------------------------------------------------------------------------------

  /*
    Statistics.
  */

  // Length of the kmers in the index.
  size_t k() const { return this->header.k; }

  // Window length for the minimizers. Does not make sense when using closed syncmers.
  size_t w() const { return this->header.w_or_s; }

  // Length of the smers when using closed syncmers. Does not make sense when using minimizers.
  size_t s() const { return this->header.w_or_s; }

  // Does the index use closed syncmers instead of minimizers.
  bool uses_syncmers() const { return this->header.get_flag(MinimizerHeader::FLAG_SYNCMERS); }

  // Does the index use weighted minimizers.
  bool uses_weighted_minimizers() const { return !(this->index.frequent_kmers.empty()); }

  // Window length in bp. We are guaranteed to have at least one kmer from the window if
  // all characters within it are valid.
  size_t window_bp() const
  {
    return (this->uses_syncmers() ? 2 * this->k() - this->s() - 1 : this->k() + this->w() - 1);
  }

  // Number of keys in the index.
  size_t size() const { return this->index.size(); }

  // Is the index empty.
  bool empty() const { return (this->size() == 0); }

  // Number of values (minimizer occurrences) in the index.
  size_t number_of_values() const { return this->index.number_of_values(); }

  // Size of the hash table.
  size_t cell_count() const { return this->index.cell_count(); }

  // Number of keys that can fit into the hash table. Exceeding it will initiate rehashing.
  size_t capacity() const { return this->index.capacity(); }

  // Current load factor of the hash table.
  double load_factor() const { return this->index.load_factor(); }

  // Number of minimizers with a single occurrence.
  size_t unique_keys() const { return this->index.unique_keys(); }

  // Returns the size of a payload in words.
  size_t payload_size() const { return this->index.payload_size(); }

  // Returns the size of a value in words.
  size_t value_size() const { return KmerIndex<key_type>::POS_SIZE + this->payload_size(); }

  // Returns the size of a cell in words.
  size_t cell_size() const { return KmerIndex<key_type>::KEY_SIZE + KmerIndex<key_type>::POS_SIZE + this->payload_size(); }

  // Call `callback` for every non-empty hash table cell in the index.
  // The callback gets the key and encoded values.
  // If callback returns false, then stop iterating.
  // Returns false if the iteration stopped early, true otherwise.
  bool for_each_minimizer(const std::function<bool(key_type, multi_value_type)>& callback) const 
  { 
      return this->index.for_each_kmer(callback);
  }

//------------------------------------------------------------------------------

  /*
    Internal implementation.
  */

private:
  MinimizerHeader    header;
  gbwt::Tags         tags;
  KmerIndex<KeyType> index;

  void copy(const MinimizerIndex& source)
  {
    this->header = source.header;
    this->tags = source.tags;
    this->index = source.index;
  }
};

//------------------------------------------------------------------------------

/*
  Decode the subset of minimizer hits and their payloads in the given subgraph induced
  by node identifiers.
  This version should only be used when the number of hits is small.
  If the minimizer is in reverse orientation, use reverse_base_pos() to reverse
  the reported occurrences.
*/
template<class KeyType>
void hits_in_subgraph
(
  const MinimizerIndex<KeyType>& index,
  KmerEncoding::multi_value_type hits,
  const std::unordered_set<nid_t>& subgraph,
  const std::function<void(KmerEncoding::value_type)>& report_hit
);

/*
  Decode the subset of minimizer hits and their payloads in the given subgraph induced
  by node identifiers. The set of node ids must be in sorted order.
  This version uses exponential search on the larger list, so it should be efficient
  regardless of the size of the subgraph and the number of hits.
  If the minimizer is in reverse orientation, use reverse_base_pos() to reverse
  the reported occurrences.
*/
template<class KeyType>
void hits_in_subgraph
(
  const MinimizerIndex<KeyType>& index,
  KmerEncoding::multi_value_type hits,
  const std::vector<nid_t>& subgraph,
  const std::function<void(KmerEncoding::value_type)>& report_hit
);

//------------------------------------------------------------------------------

// Choose the default index type.
using DefaultMinimizerIndex = MinimizerIndex<Key64>;

//------------------------------------------------------------------------------

// Numerical template class constants.

template<class KeyType>
constexpr size_t KmerIndex<KeyType>::INITIAL_CAPACITY;

template<class KeyType>
constexpr double KmerIndex<KeyType>::MAX_LOAD_FACTOR;

template<class KeyType>
constexpr size_t KmerIndex<KeyType>::MAX_PAYLOAD_SIZE;

template<class KeyType>
constexpr size_t KmerIndex<KeyType>::KEY_SIZE;

template<class KeyType>
constexpr size_t KmerIndex<KeyType>::POS_SIZE;

// Other template class variables.

template<class KeyType>
const std::string MinimizerIndex<KeyType>::EXTENSION = ".min";

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_MINIMIZER_H
