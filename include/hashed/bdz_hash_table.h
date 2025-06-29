#pragma once
#include "types.h"
#include "util.h"
#include "bdz_internal_vector.h"
#include <vector>
#include <sdsl/int_vector.hpp>

template <bool quotient = 1>
struct bdz_hash_table {
  usize n, m;
  using E = std::array<u64, 3>;
  E r, d, a, a_inv;
  bdz_internal_bit_vector W;
  rank_support_bdz_internal<1> rank1_B;
  select_support_bdz_internal<1> select1_B;
  sdsl::int_vector<> Q;

  bdz_hash_table() = default;

  bdz_hash_table(std::vector<u64> &keys) {
    i32 skip = 0;
    while (!build(keys, skip)) {
      skip += 3;
      std::cerr << "bdz_hash_table::build failed" << std::endl;
    }
  }

  void operator=(const bdz_hash_table&) = delete;

  bdz_hash_table& operator=(bdz_hash_table&& other) {
    n = other.n;
    m = other.m;
    r = other.r;
    d = other.d;
    a = other.a;
    a_inv = other.a_inv;
    W = std::move(other.W);
    rank1_B = std::move(other.rank1_B);
    rank1_B.set_vector(&W);
    select1_B = std::move(other.select1_B);
    select1_B.set_vector(&W);
    Q = std::move(other.Q);
    return *this;
  }

  bool build(std::vector<u64> &keys, i32 skip) {
    n = keys.size();
    m = (5*n+3)/4;
    usize rk = m/3;
    while (skip-- > 0) { // skip some primes if they didn't work last time
      while (!is_prime(rk) || !all_coprime(keys, rk)) {
        rk++;
      }
      rk++;
    }
    d[0] = 0;
    for (usize k = 0; k < 3; k++) {
      while (!is_prime(rk) || !all_coprime(keys, rk)) {
        rk++;
      }
      r[k] = rk;
      rk++;
      if (k > 0) {
        d[k] = r[k-1] + d[k-1];
      }
    }
    m = r[0] + r[1] + r[2];
    // std::cout << r[0] << " " << r[1] << " " << r[2] << " " << m << " " << n << std::endl;

    std::random_device dev;
    std::mt19937_64 rng(dev());
    for (usize k = 0; k < 3; k++) {
      a[k] = rng() % r[k];
      a_inv[k] = mod_pow(a[k], r[k] - 2, r[k]);
    }

    usize logn = sdsl::bits::hi(n-1) + 1;
    usize logm = sdsl::bits::hi(m-1) + 1;
    std::vector<E> vals(n);
    sdsl::int_vector<> N(m, 0, logn);
    // std::vector<sdsl::int_vector<>> L(m, sdsl::int_vector<>(0, 0, logn));
    std::vector<std::vector<u32>> L(m);
    for (usize i = 0; i < n; i++) {
      E v;
      for (usize k = 0; k < 3; k++) {
        v[k] = h(keys[i], k);
        N[v[k]]++;
        L[v[k]].push_back(i);
      }
      vals[i] = v;
    }

    // initialize "priority_queue"
    // sdsl::int_vector<> q(0, 0, logm);
    std::vector<u32> q;
    q.reserve(m);
    for (usize j = 0; j < m; j++) {
      if (N[j] == 1) {
        q.push_back(j);
      }
    }

    // topological sort?
    sdsl::bit_vector M(n);
    // sdsl::int_vector<> S(0, 0, logn);
    std::vector<u32> S;
    S.reserve(n);
    while (!q.empty()){
      auto j = q.back();
      q.pop_back();
      if (N[j] == 0) {
        continue;
      }
      assert(N[j] == 1);
      usize I = n;
      for (usize i : L[j]) {
        if (!M[i]) {
          I = i;
          break;
        }
      }
      assert(I != n);
      auto v = vals[I];
      for (usize k = 0; k < 3; k++) {
        N[v[k]]--;
        if (N[v[k]] == 1) {
          q.push_back(v[k]);
        }
      }
      S.push_back(I);
      M[I] = 1;
    }

    // check for cycles
    for (usize j = 0; j < m; j++) {
      if (N[j] > 0) {
        return 0;
      }
    }

    // build G
    sdsl::int_vector<> G(m, 3, 2);
    sdsl::bit_vector V(m);
    while (!S.empty()) {
      auto i = S.back();
      S.pop_back();
      auto v = vals[i];
      usize K = 3;
      for (usize k = 0; k < 3; k++) {
        if (V[v[k]] == 0) {
          K = k;
          break;
        }
      }
      assert(K != 3);
      G[v[K]] = (K - G[v[0]] - G[v[1]] - G[v[2]] + 9) % 3;
      for (usize k = 0; k < 3; k++) {
        V[v[k]] = 1;
      }
    }
    W = bdz_internal_bit_vector(G);
    sdsl::util::init_support(rank1_B, &W);
    if constexpr (quotient) {
      sdsl::util::init_support(select1_B, &W);
    }
    Q = sdsl::int_vector<>(n);
    for (auto key : keys) {
      if (quotient) {
        Q[rank1_B(h(key))] = key / r[index_of(key)];
      } else {
        Q[rank1_B(h(key))] = key;
      }
    }
    sdsl::util::bit_compress(Q);
    return 1;
  }

  u8 G(usize j) {
    return W.G(j);
  }

  usize h(u64 key, usize k) {
    return d[k] + mod_mul(key % r[k], a[k], r[k]);
  }

  usize index_of(u64 key) {
    return (G(h(key, 0)) + G(h(key, 1)) + G(h(key, 2))) % 3;
  }

  usize h(u64 key) {
    return h(key, index_of(key));
  }

  u64 access(usize i) {
    if constexpr (quotient) {
      usize j = select1_B(i + 1);
      usize k = 2;
      if (j < d[1]) {
        k = 0;
      } else if (j < d[2]){
        k = 1;
      }
      u64 lo = j - d[k];
      return Q[i] * r[k] + mod_mul(lo, a_inv[k], r[k]);
    } else {
      return Q[i];
    }
  }

  usize index_of_key(u64 key) { // assumes key exists
    usize j = h(key);
    assert(W[j]);
    return rank1_B(j);
  }

  bool contains(u64 key) {
    usize j = h(key);
    if (!W[j]) {
      return 0;
    }
    usize i = rank1_B(j);
    if constexpr (quotient) {
      usize k = 2;
      if (j < d[1]) {
        k = 0;
      } else if (j < d[2]){
        k = 1;
      }
      return Q[i] == key / r[k];
    } else {
      return Q[i] == key;
    }
  }

  usize serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    usize written_bytes = 0;
    written_bytes += W.serialize(out, child, "W");
    written_bytes += rank1_B.serialize(out, child, "rank1_B");
    written_bytes += select1_B.serialize(out, child, "select1_B");
    written_bytes += Q.serialize(out, child, "Q");
    written_bytes += sdsl::write_member(n, out);
    written_bytes += sdsl::write_member(m, out);
    written_bytes += sdsl::write_member(r, out);
    written_bytes += sdsl::write_member(d, out);
    written_bytes += sdsl::write_member(a, out);
    written_bytes += sdsl::write_member(a_inv, out);
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  void load(std::istream &in) {
    W.load(in);
    rank1_B.load(in);
    rank1_B.set_vector(&W);
    select1_B.load(in);
    select1_B.set_vector(&W);
    Q.load(in);
    sdsl::read_member(n, in);
    sdsl::read_member(m, in);
    sdsl::read_member(r, in);
    sdsl::read_member(d, in);
    sdsl::read_member(a, in);
    sdsl::read_member(a_inv, in);
  }
};
/*
#include "types.h"
#include "util.h"
#include "bdz_internal_vector.h"
#include <vector>
#include <sdsl/int_vector.hpp>

#define assert(x) if (!(x)) { exit(1); }

template <bool quotient = 1>
struct bdz_hash_table {
  usize n, m;
  using E = std::array<u64, 3>;
  E r, d, a, a_inv;
  bdz_internal_bit_vector W;
  rank_support_bdz_internal<1> rank1_B;
  select_support_bdz_internal<1> select1_B;
  sdsl::int_vector<> Q;

  bdz_hash_table() = default;

  explicit bdz_hash_table(std::vector<u64> &keys) {
    i32 skip = 0;
    while (!build(keys, skip)) {
      skip += 3;
      std::cerr << "bdz_hash_table::build failed" << std::endl;
    }
  }

  void operator=(const bdz_hash_table&) = delete;

  bdz_hash_table& operator=(bdz_hash_table&& other) {
    n = other.n;
    m = other.m;
    r = other.r;
    d = other.d;
    a = other.a;
    a_inv = other.a_inv;
    W = std::move(other.W);
    rank1_B = std::move(other.rank1_B);
    rank1_B.set_vector(&W);
    select1_B = std::move(other.select1_B);
    select1_B.set_vector(&W);
    Q = std::move(other.Q);
    return *this;
  }

  bool build(std::vector<u64> &keys, i32 skip) {
    n = keys.size();
    m = (5*n-3)/4;
    usize rk = m/3;
    while (skip-- >= 0) { // skip some primes if they didn't work last time
      while (!is_prime(rk) || !all_coprime(keys, rk)) {
        rk++;
      }
      rk++;
    }
    d[0] = 0;
    for (usize k = 0; k < 3; k++) {
      while (!is_prime(rk) || !all_coprime(keys, rk)) {
        rk++;
      }
      r[k] = rk;
      rk++;
      if (k > 0) {
        d[k] = r[k-1] + d[k-1];
      }
    }
    m = r[0] + r[1] + r[2];
    std::cout << r[0] << " " << r[1] << " " << r[2] << " " << m << " " << n << std::endl;

    std::random_device dev;
    std::mt19937_64 rng(dev());
    for (usize k = 0; k < 3; k++) {
      a[k] = rng() % r[k];
      a_inv[k] = mod_pow(a[k], r[k] - 2, r[k]);
    }

    usize logn = sdsl::bits::hi(n-1) + 1;
    usize logm = sdsl::bits::hi(m-1) + 1;
    std::vector<E> vals(n);
    sdsl::int_vector<> N(m, 0, logn);
    std::vector<std::vector<usize>> L(m);
    for (usize i = 0; i < n; i++) {
      E v; 
      for (usize k = 0; k < 3; k++) {
        v[k] = h(keys[i], k);
        ++N[v[k]];
        L[v[k]].push_back(i);
      }
      vals[i] = v;
    }

    // initialize "priority_queue"
    std::vector<usize> q;
    q.reserve(m);
    for (usize j = 0; j < m; j++) {
      if (N[j] == 1) {
        q.push_back(j);
      }
    }

    // topological sort?
    sdsl::bit_vector M(n);
    std::vector<usize> S;
    S.reserve(n);
    while (!q.empty()){
      auto j = q.back();
      q.pop_back();
      if (N[j] == 0) {
        continue;
      }
      assert(N[j] == 1);
      usize I = n;
      for (usize i : L[j]) {
        if (!M[i]) {
          I = i;
          break;
        }
      }
      assert(I != n);
      auto v = vals[I];
      for (usize k = 0; k < 3; k++) {
        --N[v[k]];
        if (N[v[k]] == 1) {
          q.push_back(v[k]);
        }
      }
      S.push_back(I);
      M[I] = 1;
    }

    // check for cycles
    bool fail = false;
    for (usize j = 0; j < m; j++) {
      if (N[j] > 0) {
        std::cout << "fail at " << j << std::endl;
        for (auto x : L[j]) {
          if (!M[x]) {
            std::cout << keys[x] << " " << h(keys[x], 0) << " " << h(keys[x], 1) << " " << h(keys[x], 2) << std::endl;
          }
        }
        fail = true;
      }
    }
    if (fail) {
      return 0;
    }

    // build G
    sdsl::int_vector<> G(m, 3, 2);
    sdsl::bit_vector V(m);
    while (!S.empty()) {
      auto i = S.back();
      S.pop_back();
      auto v = vals[i];
      usize K = 3;
      for (usize k = 0; k < 3; k++) {
        if (V[v[k]] == 0) {
          K = k;
          break;
        }
      }
      assert(K != 3);
      G[v[K]] = (K - G[v[0]] - G[v[1]] - G[v[2]] + 9) % 3;
      for (usize k = 0; k < 3; k++) {
        V[v[k]] = 1;
      }
    }
    W = bdz_internal_bit_vector(G);
    sdsl::util::init_support(rank1_B, &W);
    if constexpr (quotient) {
      sdsl::util::init_support(select1_B, &W);
    }
    Q = sdsl::int_vector<>(n);
    for (auto key : keys) {
      if (quotient) {
        Q[rank1_B(h(key))] = key / r[index_of(key)];
      } else {
        Q[rank1_B(h(key))] = key;
      }
    }
    sdsl::util::bit_compress(Q);
    return 1;
  }

  u8 G(usize j) {
    return W.G(j);
  }

  usize h(u64 key, usize k) {
    return d[k] + mod_mul(key % r[k], a[k], r[k]);
  }

  usize index_of(u64 key) {
    return (G(h(key, 0)) + G(h(key, 1)) + G(h(key, 2))) % 3;
  }

  usize h(u64 key) {
    return h(key, index_of(key));
  }

  u64 access(usize i) {
    if constexpr (quotient) {
      usize j = select1_B(i + 1);
      usize k = 2;
      if (j < d[1]) {
        k = 0;
      } else if (j < d[2]){
        k = 1;
      }
      u64 lo = j - d[k];
      return Q[i] * r[k] + mod_mul(lo, a_inv[k], r[k]);
    } else {
      return Q[i];
    }
  }

  usize index_of_key(u64 key) { // assumes key exists
    usize j = h(key);
    assert(W[j]);
    return rank1_B(j);
  }

  bool contains(u64 key) {
    usize j = h(key);
    if (!W[j]) {
      return 0;
    }
    usize i = rank1_B(j);
    if constexpr (quotient) {
      usize k = 2;
      if (j < d[1]) {
        k = 0;
      } else if (j < d[2]){
        k = 1;
      }
      return Q[i] == key / r[k];
    } else {
      return Q[i] == key;
    }
  }

  usize serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    usize written_bytes = 0;
    written_bytes += W.serialize(out, child, "W");
    written_bytes += rank1_B.serialize(out, child, "rank1_B");
    written_bytes += select1_B.serialize(out, child, "select1_B");
    written_bytes += Q.serialize(out, child, "Q");
    written_bytes += sdsl::write_member(n, out);
    written_bytes += sdsl::write_member(m, out);
    written_bytes += sdsl::write_member(r, out);
    written_bytes += sdsl::write_member(d, out);
    written_bytes += sdsl::write_member(a, out);
    written_bytes += sdsl::write_member(a_inv, out);
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  void load(std::istream &in) {
    W.load(in);
    rank1_B.load(in);
    rank1_B.set_vector(&W);
    select1_B.load(in);
    select1_B.set_vector(&W);
    Q.load(in);
    sdsl::read_member(n, in);
    sdsl::read_member(m, in);
    sdsl::read_member(r, in);
    sdsl::read_member(d, in);
    sdsl::read_member(a, in);
    sdsl::read_member(a_inv, in);
  }
};
*/