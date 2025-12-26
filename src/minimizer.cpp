#include <gbwtgraph/minimizer.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// KmerEncoding: Numerical class constants.
constexpr size_t KmerEncoding::FIELD_BITS;
constexpr size_t KmerEncoding::PACK_WIDTH;
constexpr size_t KmerEncoding::PACK_OVERFLOW;
constexpr size_t KmerEncoding::FIELD_CHARS;
constexpr KmerEncoding::code_type KmerEncoding::PACK_MASK;

// KmerEncoding: Other class variables.

const std::vector<unsigned char> KmerEncoding::CHAR_TO_PACK =
{
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

const std::vector<char> KmerEncoding::PACK_TO_CHAR = { 'A', 'C', 'G', 'T' };

const std::vector<KmerEncoding::code_type> KmerEncoding::LOW_MASK =
{
  // k = 0
  0x0000000000000000ull,

  // k = 1 to 32 in the low part of the key.
  0x0000000000000003ull,
  0x000000000000000Full,
  0x000000000000003Full,
  0x00000000000000FFull,
  0x00000000000003FFull,
  0x0000000000000FFFull,
  0x0000000000003FFFull,
  0x000000000000FFFFull,
  0x000000000003FFFFull,
  0x00000000000FFFFFull,
  0x00000000003FFFFFull,
  0x0000000000FFFFFFull,
  0x0000000003FFFFFFull,
  0x000000000FFFFFFFull,
  0x000000003FFFFFFFull,
  0x00000000FFFFFFFFull,
  0x00000003FFFFFFFFull,
  0x0000000FFFFFFFFFull,
  0x0000003FFFFFFFFFull,
  0x000000FFFFFFFFFFull,
  0x000003FFFFFFFFFFull,
  0x00000FFFFFFFFFFFull,
  0x00003FFFFFFFFFFFull,
  0x0000FFFFFFFFFFFFull,
  0x0003FFFFFFFFFFFFull,
  0x000FFFFFFFFFFFFFull,
  0x003FFFFFFFFFFFFFull,
  0x00FFFFFFFFFFFFFFull,
  0x03FFFFFFFFFFFFFFull,
  0x0FFFFFFFFFFFFFFFull,
  0x3FFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,

  // k = 33 to 63 in the high part of the key
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull
};

const std::vector<KmerEncoding::code_type> KmerEncoding::HIGH_MASK =
{
  // k = 0
  0x0000000000000000ull,

  // k = 1 to 32 in the low part of the key.
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,
  0x0000000000000000ull,

  // k = 33 to 63 in the high part of the key
  0x0000000000000003ull,
  0x000000000000000Full,
  0x000000000000003Full,
  0x00000000000000FFull,
  0x00000000000003FFull,
  0x0000000000000FFFull,
  0x0000000000003FFFull,
  0x000000000000FFFFull,
  0x000000000003FFFFull,
  0x00000000000FFFFFull,
  0x00000000003FFFFFull,
  0x0000000000FFFFFFull,
  0x0000000003FFFFFFull,
  0x000000000FFFFFFFull,
  0x000000003FFFFFFFull,
  0x00000000FFFFFFFFull,
  0x00000003FFFFFFFFull,
  0x0000000FFFFFFFFFull,
  0x0000003FFFFFFFFFull,
  0x000000FFFFFFFFFFull,
  0x000003FFFFFFFFFFull,
  0x00000FFFFFFFFFFFull,
  0x00003FFFFFFFFFFFull,
  0x0000FFFFFFFFFFFFull,
  0x0003FFFFFFFFFFFFull,
  0x000FFFFFFFFFFFFFull,
  0x003FFFFFFFFFFFFFull,
  0x00FFFFFFFFFFFFFFull,
  0x03FFFFFFFFFFFFFFull,
  0x0FFFFFFFFFFFFFFFull,
  0x3FFFFFFFFFFFFFFFull
};

//------------------------------------------------------------------------------

// MinimizerHeader: Numerical class constants.

constexpr std::uint32_t MinimizerHeader::TAG;
constexpr std::uint32_t MinimizerHeader::VERSION;

constexpr std::uint64_t MinimizerHeader::FLAG_MASK;
constexpr std::uint64_t MinimizerHeader::FLAG_KEY_MASK;
constexpr size_t MinimizerHeader::FLAG_KEY_OFFSET;
constexpr std::uint64_t MinimizerHeader::FLAG_WEIGHT_MASK;
constexpr size_t MinimizerHeader::FLAG_WEIGHT_OFFSET;
constexpr std::uint64_t MinimizerHeader::FLAG_SYNCMERS;
constexpr std::uint64_t MinimizerHeader::FLAG_PAYLOAD_MASK;
constexpr size_t MinimizerHeader::FLAG_PAYLOAD_OFFSET;

//------------------------------------------------------------------------------

// Position: Numerical class constants.

constexpr size_t Position::OFFSET_BITS;
constexpr size_t Position::ID_OFFSET;
constexpr Position::code_type Position::REV_MASK;
constexpr Position::code_type Position::OFF_MASK;
constexpr Position::code_type Position::NO_POSITION;

//------------------------------------------------------------------------------

// Key64: Numerical class constants.

constexpr std::size_t Key64::KEY_BITS;
constexpr std::size_t Key64::KMER_LENGTH;
constexpr std::size_t Key64::WINDOW_LENGTH;
constexpr std::size_t Key64::SMER_LENGTH;
constexpr std::size_t Key64::KMER_MAX_LENGTH;

constexpr Key64::code_type Key64::EMPTY_KEY;
constexpr Key64::code_type Key64::NO_KEY;
constexpr Key64::code_type Key64::KEY_MASK;
constexpr Key64::code_type Key64::IS_POINTER;

//------------------------------------------------------------------------------

// Key128: Numerical class constants.

constexpr std::size_t Key128::KEY_BITS;
constexpr std::size_t Key128::KMER_LENGTH;
constexpr std::size_t Key128::WINDOW_LENGTH;
constexpr std::size_t Key128::SMER_LENGTH;
constexpr std::size_t Key128::KMER_MAX_LENGTH;

constexpr Key128::code_type Key128::EMPTY_KEY;
constexpr Key128::code_type Key128::NO_KEY;
constexpr Key128::code_type Key128::KEY_MASK;
constexpr Key128::code_type Key128::IS_POINTER;

//------------------------------------------------------------------------------

MinimizerHeader::MinimizerHeader() :
  tag(TAG), version(VERSION),
  k(0), w_or_s(0),
  keys(0), unused(0), capacity(0),
  values(0), unique(0),
  flags(0)
{
}

MinimizerHeader::MinimizerHeader(size_t kmer_length, size_t window_length, size_t key_bits, size_t payload_size) :
  tag(TAG), version(VERSION),
  k(kmer_length), w_or_s(window_length),
  keys(0), unused(0), capacity(0),
  values(0), unique(0),
  flags(0)
{
    this->set_int(FLAG_KEY_MASK, FLAG_KEY_OFFSET, key_bits);
    this->set_int(FLAG_PAYLOAD_MASK, FLAG_PAYLOAD_OFFSET, payload_size);
}

void
MinimizerHeader::sanitize(size_t kmer_max_length)
{
  if(this->k > kmer_max_length)
  {
    std::cerr << "MinimizerHeader::sanitize(): Adjusting k from " << this->k << " to " << kmer_max_length << std::endl;
    this->k = kmer_max_length;
  }
  if(this->k <= 1)
  {
    std::cerr << "MinimizerHeader::sanitize(): Adjusting k from " << this->k << " to " << 2 << std::endl;
    this->k = 2;
  }

  if(this->get_flag(FLAG_SYNCMERS))
  {
    if(this->w_or_s == 0)
    {
      std::cerr << "MinimizerHeader::sanitize(): Adjusting s from " << this->w_or_s << " to " << 1 << std::endl;
      this->w_or_s = 1;
    }
    if(this->w_or_s >= this->k)
    {
      std::cerr << "MinimizerHeader::sanitize(): Adjusting s from " << this->w_or_s << " to " << (this->k - 1) << std::endl;
      this->w_or_s = this->k - 1;
    }
    if(this->downweight() > 0)
    {
      std::cerr << "MinimizerHeader::sanizize(): Weights cannot be used with syncmers" << std::endl;
      this->set_int(FLAG_WEIGHT_MASK, FLAG_WEIGHT_OFFSET, 0);
    }
  }
  else
  {
    if(this->w_or_s == 0)
    {
      std::cerr << "MinimizerHeader::sanitize(): Adjusting w from " << this->w_or_s << " to " << 1 << std::endl;
      this->w_or_s = 1;
    }
  }
}

void
MinimizerHeader::check() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("MinimizerHeader: Invalid tag");
  }

  if(this->version < VERSION || this->version > VERSION)
  {
    std::string msg = "MinimizerHeader: Expected version " + std::to_string(VERSION) + ", got version " + std::to_string(this->version);
    throw sdsl::simple_sds::InvalidData(msg);
  }

  std::uint64_t mask = (FLAG_MASK);
  if((this->flags & mask) != this->flags)
  {
    throw sdsl::simple_sds::InvalidData("MinimizerHeader: Invalid flags");
  }
}

void
MinimizerHeader::update_version()
{
  this->version = VERSION;
}

void
MinimizerHeader::set_int(std::uint64_t mask, size_t offset, size_t value)
{
  this->unset(mask);
  this->set((value << offset) & mask);
}

size_t
MinimizerHeader::get_int(std::uint64_t mask, size_t offset) const
{
  return (this->flags & mask) >> offset;
}

size_t
MinimizerHeader::key_bits() const
{
  return this->get_int(FLAG_KEY_MASK, FLAG_KEY_OFFSET);
}

size_t
MinimizerHeader::downweight() const
{
  return this->get_int(FLAG_WEIGHT_MASK, FLAG_WEIGHT_OFFSET);
}

size_t
MinimizerHeader::payload_size() const
{
  return this->get_int(FLAG_PAYLOAD_MASK, FLAG_PAYLOAD_OFFSET);
}

bool
MinimizerHeader::operator==(const MinimizerHeader& another) const
{
  return (this->tag == another.tag && this->version == another.version &&
          this->k == another.k && this->w_or_s == another.w_or_s &&
          this->flags == another.flags);
}

//------------------------------------------------------------------------------

Key64
Key64::encode(const std::string& sequence)
{
  code_type packed = 0;
  for(auto c : sequence)
  {
    auto packed_char = KmerEncoding::CHAR_TO_PACK[static_cast<std::uint8_t>(c)];
    if(packed_char > KmerEncoding::PACK_MASK)
    {
      throw std::runtime_error("Key64::encode(): Cannot encode character '" + std::to_string(c) + "'");
    }
    packed = (packed << KmerEncoding::PACK_WIDTH) | packed_char;
  }
  return Key64(packed);
}

std::string
Key64::decode(size_t k) const
{
  std::string result; result.reserve(k);
  for(size_t i = 0; i < k; i++)
  {
    result.push_back(KmerEncoding::PACK_TO_CHAR[(this->key >> ((k - i - 1) * KmerEncoding::PACK_WIDTH)) & KmerEncoding::PACK_MASK]);
  }
  return result;
}

Key64
Key64::reverse_complement(size_t k) const
{
  value_type source = ~(this->get_key()); // Complement of the kmer, plus weird high-order bits.
  code_type result = 0;
  for(size_t i = 0; i < k; i++)
  {
    result = (result << KmerEncoding::PACK_WIDTH) | (source & KmerEncoding::PACK_MASK);
    source >>= KmerEncoding::PACK_WIDTH;
  }
  return Key64(result);
}

std::ostream&
operator<<(std::ostream& out, Key64 value)
{
  out << value.key;
  return out;
}

Key128
Key128::encode(const std::string& sequence)
{
  size_t low_limit = (sequence.size() > KmerEncoding::FIELD_CHARS ? KmerEncoding::FIELD_CHARS : sequence.size());
  
  code_type packed_high = 0;
  code_type packed_low = 0;
  
  for(size_t i = 0; i < sequence.size(); i++)
  {
    auto c = sequence[i];
    auto packed_char = KmerEncoding::CHAR_TO_PACK[static_cast<std::uint8_t>(c)];
    if(packed_char > KmerEncoding::PACK_MASK)
    {
      throw std::runtime_error("Key128::encode(): Cannot encode character '" + std::to_string(c) + "'");
    }
    
    code_type& pack_to = (i < sequence.size() - low_limit) ? packed_high : packed_low;
    
    pack_to = (pack_to << KmerEncoding::PACK_WIDTH) | packed_char;
  }
  
  return Key128(packed_high, packed_low);
}

std::string
Key128::decode(size_t k) const
{
  std::string result; result.reserve(k);
  size_t low_limit = (k > KmerEncoding::FIELD_CHARS ? KmerEncoding::FIELD_CHARS : k);
  for(size_t i = KmerEncoding::FIELD_CHARS; i < k; i++)
  {
    result.push_back(KmerEncoding::PACK_TO_CHAR[(this->high >> ((k - i - 1) * KmerEncoding::PACK_WIDTH)) & KmerEncoding::PACK_MASK]);
  }
  for(size_t i = 0; i < low_limit; i++)
  {
    result.push_back(KmerEncoding::PACK_TO_CHAR[(this->low >> ((low_limit - i - 1) * KmerEncoding::PACK_WIDTH)) & KmerEncoding::PACK_MASK]);
  }
  return result;
}

Key128
Key128::reverse_complement(size_t k) const
{
  constexpr size_t HIGH_SHIFT = KmerEncoding::FIELD_BITS - KmerEncoding::PACK_WIDTH;

  // Get the complement of the kmer, plus weird high-order bits.
  value_type source = this->get_key();
  source.first = ~(source.first); source.second = ~(source.second);

  code_type high = 0, low = 0;
  for(size_t i = 0; i < k; i++)
  {
    high = (high << KmerEncoding::PACK_WIDTH) | ((low >> HIGH_SHIFT) & KmerEncoding::PACK_MASK);
    low = (low << KmerEncoding::PACK_WIDTH) | (source.second & KmerEncoding::PACK_MASK);
    source.second = ((source.first & KmerEncoding::PACK_MASK) << HIGH_SHIFT) | (source.second >> KmerEncoding::PACK_WIDTH);
    source.first >>= KmerEncoding::PACK_WIDTH;
  }

  return Key128(high, low);
}

std::ostream&
operator<<(std::ostream& out, Key128 value)
{
  out << "(" << value.high << ", " << value.low << ")";
  return out;
}

//------------------------------------------------------------------------------

template<class KeyType>
void hits_in_subgraph
(
  const MinimizerIndex<KeyType>& index,
  KmerEncoding::multi_value_type hits,
  const std::unordered_set<nid_t>& subgraph,
  const std::function<void(KmerEncoding::value_type)>& report_hit
)
{
  for(size_t i = 0; i < hits.second; i++)
  {
    auto value = index.get_value(hits, i);
    auto iter = subgraph.find(value.first.id());
    if(iter != subgraph.end()) { report_hit(value); }
  }
}

/*
  Exponential search that returns the first offset with get_value(offset) >= target.
  We assume that start < limit and get_value(start) < target.
  Returns limit if get_value(offset) < target for all offset < limit.
*/
size_t
exponential_search(size_t start, size_t limit, nid_t target, const std::function<nid_t(size_t)>& get_value)
{
  // Exponential search: low is too early.
  size_t step = 1;
  size_t low = start, candidate = start + step;
  while(candidate < limit && get_value(candidate) < target)
  {
    step *= 2;
    low = candidate; candidate += step;
  }

  // Binary search: low + 1 is the first candidate while candidate is the last.
  low++;
  size_t count = std::min(limit, candidate + 1) - low;
  while(count > 0)
  {
    step = count / 2;
    candidate = low + step;
    if(get_value(candidate) < target) { low = candidate + 1; count -= step + 1; }
    else { count = step; }
  }
  return low;
}

template<class KeyType>
void hits_in_subgraph
(
  const MinimizerIndex<KeyType>& index,
  KmerEncoding::multi_value_type hits,
  const std::vector<nid_t>& subgraph,
  const std::function<void(KmerEncoding::value_type)>& report_hit
)
{
  size_t hit_offset = 0, subgraph_offset = 0;
  while(hit_offset < hits.second && subgraph_offset < subgraph.size())
  {
    auto value = index.get_value(hits, hit_offset);
    nid_t node = value.first.id();
    if(node < subgraph[subgraph_offset])
    {
      hit_offset = exponential_search(hit_offset, hits.second, subgraph[subgraph_offset], [&](size_t offset) -> nid_t
      {
        return index.get_value(hits, offset).first.id();
      });
    }
    else if(node > subgraph[subgraph_offset])
    {
      subgraph_offset = exponential_search(subgraph_offset, subgraph.size(), node, [&](size_t offset) -> nid_t
      {
        return subgraph[offset];
      });
    }
    else
    {
      report_hit(value);
      hit_offset++;
    }
  }
}

// Instantiate templates, as we did not define the functions in the header.

template void hits_in_subgraph<Key64>
(
  const MinimizerIndex<Key64>& index,
  KmerEncoding::multi_value_type hits,
  const std::unordered_set<nid_t>& subgraph,
  const std::function<void(KmerEncoding::value_type)>& report_hit
);

template void hits_in_subgraph<Key128>
(
  const MinimizerIndex<Key128>& index,
  KmerEncoding::multi_value_type hits,
  const std::unordered_set<nid_t>& subgraph,
  const std::function<void(KmerEncoding::value_type)>& report_hit
);

template void hits_in_subgraph<Key64>
(
  const MinimizerIndex<Key64>& index,
  KmerEncoding::multi_value_type hits,
  const std::vector<nid_t>& subgraph,
  const std::function<void(KmerEncoding::value_type)>& report_hit
);

template void hits_in_subgraph<Key128>
(
  const MinimizerIndex<Key128>& index,
  KmerEncoding::multi_value_type hits,
  const std::vector<nid_t>& subgraph,
  const std::function<void(KmerEncoding::value_type)>& report_hit
);

//------------------------------------------------------------------------------

} // namespace gbwtgraph
