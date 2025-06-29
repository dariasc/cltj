#pragma once
#include "types.h"
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>

struct bdz_internal_bit_vector {
  typedef usize size_type;

  size_type m_n, m_arg_cnt = 0;
  sdsl::int_vector<64> m_v;

  bdz_internal_bit_vector() = default;

  bdz_internal_bit_vector(sdsl::int_vector<> &G)
      : m_n(G.size()), m_v(2 * ((m_n + 63) / 64)) {
    for (usize i = 0; i < G.size(); i++) {
      usize j = i >> 6;
      u64 mask = u64(1) << (i & 0x3f);
      bool hi = (G[i] & 0b10) != 0;
      bool lo = (G[i] & 0b01) != 0;

      if (hi)
        m_v[2 * j + 1] |= mask;
      else
        m_v[2 * j + 1] &= ~mask;

      if (lo)
        m_v[2 * j] |= mask;
      else
        m_v[2 * j] &= ~mask;

      if (!(lo && hi)) {
        m_arg_cnt++;
      }
    }
  }

  u8 G(usize i) const {
    usize j = i >> 6;
    u64 mask = u64(1) << (i & 0x3f);
    return 2 * ((m_v[2 * j + 1] & mask) != 0) + ((m_v[2 * j] & mask) != 0);
  }

  u64 B(usize j) const { return ~(m_v[2 * j] & m_v[2 * j + 1]); }

  bool operator[](usize i) const {
    usize j = i >> 6;
    u64 mask = u64(1) << (i & 0x3f);
    return (B(j) & mask) != 0;
  }

  bool empty() const { return m_n == 0; }

  size_type size() const { return m_n; }

  size_type bit_size() const { return m_n; }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
    size_type written_bytes = 0;
    written_bytes += sdsl::write_member(m_n, out);
    written_bytes += sdsl::write_member(m_arg_cnt, out);
    sdsl::structure_tree_node * child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    written_bytes += m_v.serialize(out, child, "double_bitvector");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  void load(std::istream &in) {
    sdsl::read_member(m_n, in);
    sdsl::read_member(m_arg_cnt, in);
    m_v.load(in);
  }
};

template <uint8_t t_b = 1, uint8_t t_pat_len = 1>
class select_support_bdz_internal {
  static_assert(t_b == 1u,
                "select_support_bdz_internal: bit pattern must be `1`");
  static_assert(t_pat_len == 1u,
                "select_support_bdz_internal: bit pattern length must be 1");

public:
  typedef usize size_type;
  typedef bdz_internal_bit_vector bit_vector_type;
  typedef sdsl::select_support_trait<t_b, t_pat_len> trait_type;

  enum { bit_pat = t_b };
  enum { bit_pat_len = (uint8_t)1 };

private:
  bit_vector_type const *m_v;
  u32 m_logn = 0,  // \f$ log(size) \f$
      m_logn2 = 0, // \f$ log^2(size) \f$
      m_logn4 = 0; // \f$ log^4(size) \f$
  // entry i of m_superblock equals the answer to select_1(B,i*4096)
  sdsl::int_vector<0> m_superblock;
  sdsl::int_vector<0> *m_longsuperblock = nullptr;
  sdsl::int_vector<0> *m_miniblock = nullptr;
  size_type m_arg_cnt = 0;

public:
  select_support_bdz_internal &operator=(select_support_bdz_internal &&ss) {
    if (this != &ss) {
      m_logn = ss.m_logn;                        // copy log n
      m_logn2 = ss.m_logn2;                      // copy (logn)^2
      m_logn4 = ss.m_logn4;                      // copy (logn)^4
      m_superblock = std::move(ss.m_superblock); // move long superblock
      m_arg_cnt = ss.m_arg_cnt;                  // copy count of 1-bits
      m_v = ss.m_v; // copy pointer to the supported bit vector

      delete[] m_longsuperblock;
      m_longsuperblock = ss.m_longsuperblock;
      ss.m_longsuperblock = nullptr;

      delete[] m_miniblock;
      m_miniblock = ss.m_miniblock;
      ss.m_miniblock = nullptr;
    }
    return *this;
  }

  explicit select_support_bdz_internal(bit_vector_type const *v = nullptr) {
    set_vector(v);
    init_data();

    if (m_v == nullptr)
      return;
    // Count the number of arguments in the bit vector
    m_arg_cnt = v->m_arg_cnt;

    const size_type SUPER_BLOCK_SIZE = 64 * 64;

    if (m_arg_cnt ==
        0) // if there are no arguments in the vector we are done...
      return;

    //    size_type sb = (m_arg_cnt+63+SUPER_BLOCK_SIZE-1)/SUPER_BLOCK_SIZE; //
    //    number of superblocks, add 63 as the last block could contain 63
    //    uninitialized bits
    size_type sb = (m_arg_cnt + SUPER_BLOCK_SIZE - 1) /
                   SUPER_BLOCK_SIZE; // number of superblocks
    delete[] m_miniblock;
    m_miniblock = new sdsl::int_vector<0>[sb];

    m_superblock = sdsl::int_vector<0>(
        sb, 0, m_logn); // TODO: hier koennte man logn noch optimieren...s

    size_type arg_position[SUPER_BLOCK_SIZE];
    size_type l = 0;
    // uint64_t const * data = v->data();
    uint64_t carry_new = 0;
    size_type last_k64 = 1, sb_cnt = 0;
    for (size_type i = 0, cnt_old = 0, cnt_new = 0, last_k64_sum = 1;
         i < (((v->bit_size() + 63) >> 6) << 6); i += 64, ++l) {
      cnt_new += trait_type::args_in_the_word(m_v->B(l), carry_new);
      cnt_new = std::min(cnt_new,
                         m_arg_cnt); // For (0, 1), we may find nonexistent args
                                     // in the padding after the bitvector.
      if (cnt_new >= last_k64_sum) {
        arg_position[last_k64 - 1] =
            i + trait_type::ith_arg_pos_in_the_word(
                    m_v->B(l), last_k64_sum - cnt_old, carry_new);
        last_k64 += 64;
        last_k64_sum += 64;

        if (last_k64 == SUPER_BLOCK_SIZE + 1) {
          m_superblock[sb_cnt] = arg_position[0];
          size_type pos_of_last_arg_in_the_block = arg_position[last_k64 - 65];

          for (size_type ii = arg_position[last_k64 - 65] + 1,
                         j = last_k64 - 65;
               ii < v->size() and j < SUPER_BLOCK_SIZE; ++ii)
            if ((*v)[ii] == 0b1) {
              pos_of_last_arg_in_the_block = ii;
              ++j;
            }
          size_type pos_diff = pos_of_last_arg_in_the_block - arg_position[0];
          if (pos_diff > m_logn4) { // long block
            if (m_longsuperblock == nullptr)
              m_longsuperblock =
                  new sdsl::int_vector<0>[sb + 1]; // create longsuperblock
            // GEANDERT am 2010-07-17 +1 nach pos_of_last_arg..
            m_longsuperblock[sb_cnt] = sdsl::int_vector<0>(
                SUPER_BLOCK_SIZE, 0,
                sdsl::bits::hi(pos_of_last_arg_in_the_block) + 1);
            for (size_type j = arg_position[0], k = 0;
                 k < SUPER_BLOCK_SIZE and j <= pos_of_last_arg_in_the_block;
                 ++j)
              if ((*v)[j] == 0b1) {
                if (k >= SUPER_BLOCK_SIZE) {
                  for (size_type ii = 0; ii < SUPER_BLOCK_SIZE; ++ii) {
                    std::cout << "(" << ii << ","
                              << m_longsuperblock[sb_cnt][ii] << ") ";
                  }
                  std::cout << std::endl;
                  std::cout << "k=" << k
                            << " SUPER_BLOCK_SIZE=" << SUPER_BLOCK_SIZE
                            << std::endl;
                  std::cout << "pos_of_last_arg_in_the_block"
                            << pos_of_last_arg_in_the_block << std::endl;
                  std::cout.flush();
                }
                m_longsuperblock[sb_cnt][k++] = j;
              }
          } else {
            m_miniblock[sb_cnt] =
                sdsl::int_vector<0>(64, 0, sdsl::bits::hi(pos_diff) + 1);
            for (size_type j = 0; j < SUPER_BLOCK_SIZE; j += 64) {
              m_miniblock[sb_cnt][j / 64] = arg_position[j] - arg_position[0];
            }
          }
          ++sb_cnt;
          last_k64 = 1;
        }
      }
      cnt_old = cnt_new;
    }
    // handle last block: append long superblock
    if (last_k64 > 1) {
      if (m_longsuperblock == nullptr)
        m_longsuperblock =
            new sdsl::int_vector<0>[sb + 1]; // create longsuperblock
      m_longsuperblock[sb_cnt] = sdsl::int_vector<0>(
          SUPER_BLOCK_SIZE, 0, sdsl::bits::hi(v->size() - 1) + 1);
      for (size_type i = arg_position[0], k = 0; i < v->size(); ++i) {
        if ((*v)[i] == 0b1) {
          m_longsuperblock[sb_cnt][k++] = i;
        }
      }
      ++sb_cnt;
    }
  }

  ~select_support_bdz_internal() {
    delete[] m_longsuperblock;
    delete[] m_miniblock;
  }

  size_type select(size_type i) const {
    assert(i > 0 and i <= m_arg_cnt);

    i = i - 1;
    size_type sb_idx = i >> 12;   // i/4096
    size_type offset = i & 0xFFF; // i%4096
    if (m_longsuperblock != nullptr and !m_longsuperblock[sb_idx].empty()) {
      return m_longsuperblock[sb_idx][offset];
    } else {
      if ((offset & 0x3F) == 0) {
        assert(sb_idx < m_superblock.size());
        assert((offset >> 6) < m_miniblock[sb_idx].size());
        return m_superblock[sb_idx] + m_miniblock[sb_idx][offset >> 6 /*/64*/];
      } else {
        i = i - (sb_idx << 12) - ((offset >> 6) << 6);
        // now i > 0 and i <= 64
        assert(i > 0);
        size_type pos =
            m_superblock[sb_idx] + m_miniblock[sb_idx][offset >> 6] + 1;

        // now pos is the position from where we search for the ith argument
        size_type word_pos = pos >> 6;
        size_type word_off = pos & 0x3F;
        size_type l = word_pos;
        // uint64_t const * data = m_v->data() + word_pos;
        // uint64_t carry = trait_type::init_carry(data, word_pos);
        uint64_t carry = 0;
        size_type args =
            trait_type::args_in_the_first_word(m_v->B(l), word_off, carry);

        if (args >= i) {
          return (word_pos << 6) + trait_type::ith_arg_pos_in_the_first_word(
                                       m_v->B(l), i, word_off, carry);
        }
        word_pos += 1;
        size_type sum_args = args;
        // carry = trait_type::get_carry(*data);
        uint64_t old_carry = carry;
        args = trait_type::args_in_the_word(m_v->B(++l), carry);
        while (sum_args + args < i) {
          sum_args += args;
          // assert(data + 1 < m_v->data() + ((m_v->bit_size() + 63) >> 6));
          // old_carry = carry;
          args = trait_type::args_in_the_word(m_v->B(++l), carry);
          word_pos += 1;
        }
        return (word_pos << 6) + trait_type::ith_arg_pos_in_the_word(
                                     m_v->B(l), i - sum_args, old_carry);
      }
    }
  }

  size_type operator()(size_type i) const { return select(i); }

  size_type size() const { return m_v->size(); }

  void set_vector(bit_vector_type const *v = nullptr) { m_v = v; }

  void init_data() {
    m_arg_cnt = 0;
    if (nullptr == m_v) {
      m_logn = m_logn2 = m_logn4 = 0;
    } else {
      m_logn = sdsl::bits::hi(((m_v->bit_size() + 63) >> 6) << 6) +
               1; // TODO maybe it's better here to take a max(...,12)
      m_logn2 = m_logn * m_logn;
      m_logn4 = m_logn2 * m_logn2;
    }
    delete[] m_longsuperblock;
    m_longsuperblock = nullptr;
    delete[] m_miniblock;
    m_miniblock = nullptr;
  }

  void load(std::istream &in, bit_vector_type const *v = nullptr) {
    set_vector(v);
    init_data();
    // read the number of 1-bits in the supported bit_vector
    in.read((char *)&m_arg_cnt, sizeof(size_type) / sizeof(char));
    size_type sb = (m_arg_cnt + 4095) >> 12;

    if (m_arg_cnt) {         // if there exists 1-bits to be supported
      m_superblock.load(in); // load superblocks

      delete[] m_miniblock;
      m_miniblock = nullptr;
      delete[] m_longsuperblock;
      m_longsuperblock = nullptr;

      sdsl::bit_vector mini_or_long; // Helper vector: mini or long block?
      mini_or_long.load(in);         // Load the helper vector
      m_miniblock =
          new sdsl::int_vector<0>[sb]; // Create miniblock int_vector<0>
      if (!mini_or_long.empty())
        m_longsuperblock =
            new sdsl::int_vector<0>[sb]; // Create longsuperblock int_vector<0>

      for (size_type i = 0; i < sb; ++i)
        if (!mini_or_long.empty() and not mini_or_long[i]) {
          m_longsuperblock[i].load(in);
        } else {
          m_miniblock[i].load(in);
        }
    }
  }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                      std::string name = "") const {
    sdsl::structure_tree_node *child =
        sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;
    // write the number of 1-bits in the supported bit_vector
    out.write((char *)&m_arg_cnt, sizeof(size_type) / sizeof(char));
    written_bytes = sizeof(size_type) / sizeof(char);
    // number of superblocks in the data structure
    size_type sb = (m_arg_cnt + 4095) >> 12;

    if (m_arg_cnt) { // if there exists 1-bits to be supported
      written_bytes += m_superblock.serialize(
          out, child, "superblock"); // serialize superblocks
      sdsl::bit_vector mini_or_long; // Helper vector: mini or long block?
      if (m_longsuperblock != nullptr) {
        mini_or_long.resize(
            sb); // resize indicator bit_vector to the number of superblocks
        for (size_type i = 0; i < sb; ++i)
          mini_or_long[i] = !m_miniblock[i].empty();
      }
      written_bytes += mini_or_long.serialize(out, child, "mini_or_long");
      size_type written_bytes_long = 0;
      size_type written_bytes_mini = 0;
      for (size_type i = 0; i < sb; ++i)
        if (!mini_or_long.empty() and !mini_or_long[i]) {
          written_bytes_long += m_longsuperblock[i].serialize(out);
        } else {
          written_bytes_mini += m_miniblock[i].serialize(out);
        }
      written_bytes += written_bytes_long;
      written_bytes += written_bytes_mini;
      sdsl::structure_tree_node *child_long = sdsl::structure_tree::add_child(
          child, "longsuperblock", sdsl::util::class_name(m_longsuperblock));
      sdsl::structure_tree::add_size(child_long, written_bytes_long);
      sdsl::structure_tree_node *child_mini = sdsl::structure_tree::add_child(
          child, "minisuperblock", sdsl::util::class_name(m_miniblock));
      sdsl::structure_tree::add_size(child_mini, written_bytes_mini);
    }
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  void swap(select_support_bdz_internal& ss) {
    if (this != &ss) {
      std::swap(m_logn, ss.m_logn);
      std::swap(m_logn2, ss.m_logn2);
      std::swap(m_logn4, ss.m_logn4);
      m_superblock.swap(ss.m_superblock);
      std::swap(m_longsuperblock, ss.m_longsuperblock);
      std::swap(m_miniblock, ss.m_miniblock);
      std::swap(m_arg_cnt, ss.m_arg_cnt);
    }
  }
};

template <uint8_t t_b = 1, uint8_t t_pat_len = 1>
class rank_support_bdz_internal {
private:
  static_assert(t_b == 1u,
                "rank_support_bdz_internal: bit pattern must be `1`");
  static_assert(t_pat_len == 1u or t_pat_len == 2u,
                "rank_support_bdz_internal: bit pattern length must be 1 or 2");

public:
  typedef usize size_type;
  typedef bdz_internal_bit_vector bit_vector_type;
  typedef sdsl::rank_support_trait<t_b, t_pat_len> trait_type;
  enum { bit_pat = t_b };
  enum { bit_pat_len = t_pat_len };

private:
  // basic block for interleaved storage of superblockrank and blockrank
  bit_vector_type const *m_v;
  sdsl::int_vector<64> m_basic_block;

public:
  explicit rank_support_bdz_internal(bit_vector_type const *v = nullptr) {
    set_vector(v);
    if (v == nullptr) {
      return;
    } else if (v->empty()) {
      m_basic_block =
          sdsl::int_vector<64>(2, 0); // resize structure for basic_blocks
      return;
    }
    size_type basic_block_size = (((v->bit_size() + 63) >> 9) + 1) << 1;
    m_basic_block.resize(basic_block_size); // resize structure for basic_blocks
    if (m_basic_block.empty())
      return;
    size_type l = 0;
    // uint64_t const * data = m_v->data();
    size_type i, j = 0;
    m_basic_block[0] = m_basic_block[1] = 0;

    uint64_t carry = trait_type::init_carry();
    uint64_t sum = trait_type::args_in_the_word(v->B(l), carry);
    uint64_t second_level_cnt = 0;
    for (i = 1; i < ((m_v->bit_size() + 63) >> 6); ++i) {
      if (!(i & 0x7)) { // if i%8==0
        j += 2;
        m_basic_block[j - 1] = second_level_cnt;
        m_basic_block[j] = m_basic_block[j - 2] + sum;
        second_level_cnt = sum = 0;
      } else {
        second_level_cnt |=
            sum << (63 - 9 * (i & 0x7)); //  54, 45, 36, 27, 18, 9, 0
      }
      sum += trait_type::args_in_the_word(v->B(++l), carry);
    }
    if (i & 0x7) { // if i%8 != 0
      second_level_cnt |= sum << (63 - 9 * (i & 0x7));
      m_basic_block[j + 1] = second_level_cnt;
    } else { // if i%8 == 0
      j += 2;
      m_basic_block[j - 1] = second_level_cnt;
      m_basic_block[j] = m_basic_block[j - 2] + sum;
      m_basic_block[j + 1] = 0;
    }
  }

  rank_support_bdz_internal(rank_support_bdz_internal const &) = default;
  rank_support_bdz_internal(rank_support_bdz_internal &&) = default;
  rank_support_bdz_internal &
  operator=(rank_support_bdz_internal const &) = default;
  rank_support_bdz_internal &operator=(rank_support_bdz_internal &&) = default;

  size_type rank(size_type idx) const {
    assert(m_v != nullptr);
    assert(idx <= m_v->size());
    uint64_t const *p = m_basic_block.data() +
                        ((idx >> 8) & 0xFFFFFFFFFFFFFFFEULL); // (idx/512)*2
    if (idx & 0x3F) {                                         // if (idx%64)!=0
      u64 data = m_v->B(idx >> 6);
      return *p + ((*(p + 1) >> (63 - 9 * ((idx & 0x1FF) >> 6))) & 0x1FF) +
             trait_type::word_rank(&data, idx & 0x3F);
    } else {
      return *p + ((*(p + 1) >> (63 - 9 * ((idx & 0x1FF) >> 6))) & 0x1FF);
    }
  }

  inline size_type operator()(size_type idx) const { return rank(idx); }

  size_type size() const { return m_v->size(); }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                      std::string name = "") const {
    size_type written_bytes = 0;
    sdsl::structure_tree_node *child =
        sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    written_bytes += m_basic_block.serialize(out, child, "cumulative_counts");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  void load(std::istream &in, bit_vector_type const *v = nullptr) {
    set_vector(v);
    m_basic_block.load(in);
  }

  void set_vector(bit_vector_type const *v = nullptr) { m_v = v; }

  void swap(rank_support_bdz_internal& rs) {
    if (this != &rs) { // if rs and _this_ are not the same object
      m_basic_block.swap(rs.m_basic_block);
    }
  }
};
