#include <gbwtgraph/minimizer.h>

namespace gbwtgraph
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr std::uint32_t MinimizerHeader::TAG;
constexpr std::uint32_t MinimizerHeader::VERSION;
constexpr std::uint64_t MinimizerHeader::FLAG_MASK;
constexpr std::uint64_t MinimizerHeader::FLAG_KEY_MASK;
constexpr size_t MinimizerHeader::FLAG_KEY_OFFSET;
constexpr std::uint32_t MinimizerHeader::FIRST_VERSION;
constexpr std::uint64_t MinimizerHeader::FIRST_FLAG_MASK;

constexpr std::size_t Key64::KEY_BITS;
constexpr std::size_t Key64::KMER_LENGTH;
constexpr std::size_t Key64::WINDOW_LENGTH;
constexpr std::size_t Key64::KMER_MAX_LENGTH;

constexpr Key64::key_type Key64::EMPTY_KEY;
constexpr Key64::key_type Key64::NO_KEY;

constexpr size_t Key64::PACK_WIDTH;
constexpr Key64::key_type Key64::PACK_MASK;

//------------------------------------------------------------------------------

// Other class variables.

const std::vector<unsigned char> Key64::CHAR_TO_PACK =
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

const std::vector<char> Key64::PACK_TO_CHAR = { 'A', 'C', 'G', 'T' };

const std::vector<Key64::key_type> Key64::KMER_MASK =
{
  0x0000000000000000ull,
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

MinimizerHeader::MinimizerHeader() :
  tag(TAG), version(VERSION),
  k(0), w(0),
  keys(0), capacity(0), max_keys(0),
  values(0), unique(0),
  flags(0)
{
}

MinimizerHeader::MinimizerHeader(size_t kmer_length, size_t window_length, size_t initial_capacity, double max_load_factor, size_t key_bits) :
  tag(TAG), version(VERSION),
  k(kmer_length), w(window_length),
  keys(0), capacity(initial_capacity), max_keys(initial_capacity * max_load_factor),
  values(0), unique(0),
  flags(0)
{
  this->set_int(FLAG_KEY_MASK, FLAG_KEY_OFFSET, key_bits);
}

void
MinimizerHeader::sanitize(size_t kmer_max_length)
{
  if(this->k > kmer_max_length)
  {
    std::cerr << "MinimizerHeader::sanitize(): Adjusting k from " << this->k << " to " << kmer_max_length << std::endl;
    this->k = kmer_max_length;
  }
  if(this->k == 0)
  {
    std::cerr << "MinimizerHeader::sanitize(): Adjusting k from " << this->k << " to " << 1 << std::endl;
    this->k = 1;
  }

  if(this->w == 0)
  {
    std::cerr << "MinimizerHeader::sanitize(): Adjusting w from " << this->w << " to " << 1 << std::endl;
    this->w = 1;
  }
}

bool
MinimizerHeader::check() const
{
  if(this->TAG != TAG) { return false; }
  switch(this->version)
  {
  case VERSION:
    return ((this->flags & FLAG_MASK) == this->flags);
  case FIRST_VERSION:
    return ((this->flags & FIRST_FLAG_MASK) == this->flags);
  default:
    return false;
  }
}

void
MinimizerHeader::update_version(size_t key_bits)
{
  this->version = VERSION;
  this->set_int(FLAG_KEY_MASK, FLAG_KEY_OFFSET, key_bits);
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
  switch(this->version)
  {
  case VERSION:
    return this->get_int(FLAG_KEY_MASK, FLAG_KEY_OFFSET);
  case FIRST_VERSION:
    return Key64::KEY_BITS;
  default:
    return 0;
  }
}

bool
MinimizerHeader::operator==(const MinimizerHeader& another) const
{
  return (this->tag == another.tag && this->version == another.version &&
          this->k == another.k && this->w == another.w &&
          this->keys == another.keys && this->capacity == another.capacity && this->max_keys == another.max_keys &&
          this->values == another.values && this->unique == another.unique &&
          this->flags == another.flags);
}

//------------------------------------------------------------------------------

} // namespace gbwtgraph
