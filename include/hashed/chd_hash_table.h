#include "types.h"
#include <sdsl/int_vector.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/coder_elias_gamma.hpp>
#include "util.h"

template <bool quotient = 1> 
struct chd_hash_table {
  usize n, m, r;
  u64 a;
  sdsl::int_vector<> Q, K, K_inv;
  sdsl::vlc_vector<sdsl::coder::elias_gamma<>> J, J_inv;
  sdsl::bit_vector B;
  sdsl::rank_support_v<1> rank1_B;
  sdsl::select_support_mcl<1> select1_B;

  chd_hash_table() = default;

  chd_hash_table(std::vector<u64> &keys, const u32 lambda = 2, const double alpha = 0.5) {
    while (!build(keys, lambda, alpha)) {
      std::cerr << "chd_hash_table::build failed" << std::endl;
    }
  }

  void operator=(const chd_hash_table&) = delete;

  chd_hash_table& operator=(chd_hash_table&& other) {
    n = other.n;
    m = other.m;
    r = other.r;
    a = other.a;
    Q = std::move(other.Q);
    K = std::move(other.K);
    K_inv = std::move(other.K_inv);
    J = std::move(other.J);
    J_inv = std::move(other.J_inv);
    B = std::move(other.B);
    rank1_B = std::move(other.rank1_B);
    rank1_B.set_vector(&B);
    select1_B = std::move(other.select1_B);
    select1_B.set_vector(&B);
    return *this;
  }

  bool build(std::vector<u64> &keys, const u32 lambda, const double alpha) {
    n = keys.size();
    m = usize(n / alpha);
    while (!is_prime(m) || !all_coprime(keys, m)) {
      m++;
    }
    r = (n + lambda - 1) / lambda;
    while (!is_prime(r)) {
      r++;
    }

    std::random_device dev;
    std::mt19937_64 rng(dev());
    a = (rng()) % r;

    usize max_bucket_size = 0;
    std::vector<sdsl::int_vector<>> buckets(r);
    for (auto key : keys) {
      buckets[g(key)].push_back(key);
      max_bucket_size = std::max(max_bucket_size, buckets[g(key)].size());
    }

    std::vector<sdsl::int_vector<>> I(max_bucket_size+1);
    for (usize i = 0; i < r; i++) {
      I[max_bucket_size - buckets[i].size()].push_back(i);
    }

    sdsl::int_vector<> S, J_tmp(r);
    B = sdsl::bit_vector(m);
    K = sdsl::int_vector<>();
    K_inv = sdsl::int_vector<>();
    for (auto &v : I) {
      for (auto i : v) {
        bool done = 0;
        for (usize j = 0; j < m; j++) {
          if (j == K.size()) {
            K.push_back(rng() % (m-1) + 1);
            K_inv.push_back(mod_pow(K.back(), m-2, m));
          }
          auto k = K[j];
          for (auto key : buckets[i]) {
            if (B[h(key, k)]) {
              while (!S.empty()) {
                B[S.back()] = 0;
                S.pop_back();
              }
              break;
            }
            B[h(key, k)] = 1;
            S.push_back(h(key, k));
          }
          if (S.size() == buckets[i].size()) {
            S.clear();
            done = 1;
            J_tmp[i] = j;
            break;
          }
        }
        if (!done) {
          return 0;
        }
      }
    }
    J = decltype(J)(J_tmp);
    rank1_B = sdsl::rank_support_v<1>(&B);
    select1_B = sdsl::select_support_mcl<1>(&B);

    Q = sdsl::int_vector<>(n);
    sdsl::int_vector<> J_inv_tmp(n);
    for (auto key : keys) {
      if (quotient) {
        Q[rank1_B(hash(key))] = key / m;
        J_inv_tmp[rank1_B(hash(key))] = J[g(key)];
      } else {
        Q[rank1_B(hash(key))] = key;
      }
    }
    if (quotient) {
      J_inv = decltype(J_inv)(J_inv_tmp);
    }
    sdsl::util::bit_compress(Q);
    return 1;
  }

  usize g(u64 key) {
    return ((key % r) * a) % r;
  }

  usize h(u64 key, u64 k) {
    return ((key % m) * k) % m;
  }

  usize hash(u64 key) {
    return h(key, K[J[g(key)]]);
  }

  u64 access(usize i) {
    if (quotient) {
      return Q[i] * m + (select1_B(i+1) * K_inv[J_inv[i]]) % m;
    } else {
      return Q[i];
    }
  }

  bool contains(u64 key) {
    usize j = hash(key); 
    if (!B[j]) {
      return 0;
    }
    usize i = rank1_B(j);
    if (quotient) {
      return key == Q[i] * m + key % m;
    } else {
      return key == Q[i];
    }
  }

  usize size_in_bytes() {
    return sdsl::size_in_bytes(J) 
      + sdsl::size_in_bytes(J_inv) 
      + sdsl::size_in_bytes(K)
      + sdsl::size_in_bytes(K_inv)
      + sdsl::size_in_bytes(Q)
      + sdsl::size_in_bytes(B)
      + sdsl::size_in_bytes(rank1_B)
      + sdsl::size_in_bytes(select1_B);
  }

  void decomp() {
    usize total = size_in_bytes();
    std::cout << "bits per elem: " << (total*8)/double(n) << '\n';
    std::cout << "J: " << sdsl::size_in_bytes(J)*8/double(n) << std::endl;
    std::cout << "J_inv: " << sdsl::size_in_bytes(J_inv)*8/double(n) << std::endl;
    std::cout << "K: " << sdsl::size_in_bytes(K)*8/double(n) << std::endl;
    std::cout << "K_inv: " << sdsl::size_in_bytes(K_inv)*8/double(n) << std::endl;
    std::cout << "Q: " << sdsl::size_in_bytes(Q)*8/double(n) << std::endl;
    std::cout << "B: " << sdsl::size_in_bytes(B)*8/double(n) << std::endl;
    std::cout << "B rank1: " << sdsl::size_in_bytes(rank1_B)*8/double(n) << std::endl;
    std::cout << "B select1: " << sdsl::size_in_bytes(select1_B)*8/double(n) << std::endl;
  }
};
