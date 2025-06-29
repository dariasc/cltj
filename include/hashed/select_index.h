#pragma once
#include "types.h"
#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>

template <
  auto string_access, 
  usize k, 
  usize s,
  bool full_sigma = 0,
  typename x_type = sdsl::bit_vector,
  typename b_type = sdsl::bit_vector
>
struct select_index {
  static_assert((k & (k - 1)) == 0, "k must be power of 2");

  using T = typename member_function_traits<string_access>::class_type;
  T *v;
  usize n, sigma;
  // char freq
  x_type X;
  typename x_type::rank_0_type rank0_X;
  typename x_type::select_1_type select1_X;
  // inverse permutation
  sdsl::int_vector<> S;
  b_type B; 
  typename b_type::rank_1_type rank1_B;
  // rank D prime
  sdsl::int_vector<> R;
  using f_type = sdsl::sd_vector<>;
  f_type F;
  typename f_type::rank_1_type rank1_F;

  select_index() = default;

  select_index(T* v, usize n, usize sigma) : v(v), n(n), sigma(sigma) {
    R = sdsl::int_vector<>(n, 0, sdsl::bits::hi(k));
    sdsl::int_vector<> freq(sigma);
    std::vector<std::vector<usize>> F_pos(sigma);
    usize F_ones = 0;
    for (usize i = 0; i < n; i++) {
      usize c = access(i); 
      assert(c < sigma);
      freq[c]++;
      R[i] = freq[c] & (k-1);
      if ((freq[c] & (k-1)) == 0) {
        F_ones++;
        F_pos[c].emplace_back(i);
      }
    }
    sdsl::sd_vector_builder F_builder(sigma*n, F_ones);
    for (usize c = 0; c < sigma; c++) {
      for (auto i : F_pos[c]) {
        F_builder.set(c*n + i);
      }
    }
    F = f_type(F_builder);
    sdsl::util::init_support(rank1_F, &F);

    usize X_len = n+sigma+1;
    if constexpr (full_sigma) {
      X_len = n+1;
    }
    sdsl::bit_vector X_tmp(X_len);
    usize w = 0;
    for (usize c = 0; c < sigma; c++) {
      X_tmp[w++] = 1;
      if constexpr (full_sigma) {
        while (--freq[c]) {
          X_tmp[w++] = 0;
        }
      } else {
        while (freq[c]--) {
          X_tmp[w++] = 0;
        }
      }
    }
    X_tmp[w++] = 1;
    X = x_type(X_tmp);
    sdsl::util::init_support(rank0_X, &X);
    sdsl::util::init_support(select1_X, &X);

    sdsl::bit_vector B_tmp(n), V(n);
    std::vector<std::array<usize, 2>> P;
    for (usize i = 0; i < n; i++) {
      if (V[i] == 0) {
        V[i] = 1;
        usize l = 0;
        usize last = i;
        usize j = i;
        usize next = pi_inv(j);
        while (next != i) {
          j = next;
          next = pi_inv(j);
          V[j] = 1;
          l++;
          if (l % s == 0) {
            P.push_back({j, last});
            last = j;
            B_tmp[j] = 1;
          }
        }
        if (l >= s) {
          P.push_back({i, last});
          B_tmp[i] = 1;
        }
      }
    }
    B = b_type(B_tmp);
    sdsl::util::init_support(rank1_B, &B);
    S = sdsl::int_vector<>(rank1_B(n), 0, sdsl::bits::hi(n-1) + 1); 
    for (auto [i, last] : P) {
      S[rank1_B(i)] = last;
    }
  }

  void operator=(const select_index&) = delete;

  select_index& operator=(select_index&& other) {
    v = other.v;
    n = other.n;
    sigma = other.sigma;
    X = std::move(other.X);
    rank0_X = std::move(other.rank0_X);
    rank0_X.set_vector(&X);
    select1_X = std::move(other.select1_X);
    select1_X.set_vector(&X);
    B = std::move(other.B);
    rank1_B = std::move(other.rank1_B);
    rank1_B.set_vector(&B);
    S = std::move(other.S);
    R = std::move(other.R);
    F = std::move(other.F);
    rank1_F = std::move(other.rank1_F);
    rank1_F.set_vector(&F);
    return *this;
  }

  usize access(usize i) {
    return ((*v).*string_access)(i);
  }

  usize select(usize c, usize r) {
    usize l = count_less_than_c(c);
    return S[l + r - 1];
    // return pi(l + r - 1);
  }

  usize count_less_than_c(usize c) {
    if constexpr (full_sigma) {
      return rank0_X(select1_X(c+1)) + c;
    }
    return rank0_X(select1_X(c+1));
  }

  usize count_c(usize c) {
    return count_less_than_c(c+1) - count_less_than_c(c);
  }

  usize rank_D_prime(usize i, usize c) {
    usize upper = rank1_F(n*c + i + 1) - rank1_F(n*c);
    return k*upper + R[i];
  }

  usize pi_inv(usize i) {
    usize c = access(i);
    return count_less_than_c(c) + rank_D_prime(i, c) - 1;
  }

  usize pi(usize i) {
    bool done = 0;
    usize j = i;
    usize next = pi_inv(j);
    while (next != i) {
      if (B[j] && !done) {
        done = 1;
        j = S[rank1_B(j)];
      } else {
        j = next;
      }
      next = pi_inv(j);
    }
    return j;
  }

  usize serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    usize written_bytes = 0;
    written_bytes += X.serialize(out, child, "X");
    written_bytes += rank0_X.serialize(out, child, "rank0_X");
    written_bytes += select1_X.serialize(out, child, "select1_X");
    written_bytes += R.serialize(out, child, "R");
    written_bytes += F.serialize(out, child, "F");
    written_bytes += rank1_F.serialize(out, child, "rank1_F");
    written_bytes += S.serialize(out, child, "S");
    written_bytes += B.serialize(out, child, "B");
    written_bytes += rank1_B.serialize(out, child, "rank1_B");
    written_bytes += sdsl::write_member(n, out);
    written_bytes += sdsl::write_member(sigma, out);
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  void load(std::istream &in) {
    X.load(in);
    rank0_X.load(in);
    rank0_X.set_vector(&X);
    select1_X.load(in);
    select1_X.set_vector(&X);
    R.load(in);
    F.load(in);
    rank1_F.load(in);
    rank1_F.set_vector(&F);
    S.load(in);
    B.load(in);
    rank1_B.load(in);
    rank1_B.set_vector(&B);
    sdsl::read_member(n, in);
    sdsl::read_member(sigma, in);
  }
};
