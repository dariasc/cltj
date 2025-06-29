#pragma once
#include <cltj_config.hpp>
#include <array>

#include "triple_pattern.hpp"
#include "tuple_index.h"

namespace cltj {

  struct hashed_trie {
    typedef uint64_t size_type;
    typedef uint64_t value_type;

    size_type sigma_x = 0;
    pair_table *pair;
    std::array<uint8_t, 2> pair_perm;
    triple_table *triple;
    std::array<uint8_t, 3> triple_perm;

    hashed_trie() : pair(nullptr), triple(nullptr) {
    }

    hashed_trie(pair_table *pair, decltype(pair_perm) pair_perm, triple_table *triple,
                     decltype(triple_perm) triple_perm) : pair(pair), pair_perm(pair_perm), triple(triple),
                                                          triple_perm(triple_perm) {
      sigma_x = pair_perm[0] == 0 ? pair->sigma_x : pair->sigma_y;
    }

    bool exists(size_type len, const array<value_type, 3> &value) const {
      if (len == 0) {
        return true;
      } else if (len <= 1) {
        return value[0]-1 < sigma_x;
      } else if (len == 2) {
        return pair->contains(value[pair_perm[0]]-1, value[pair_perm[1]]-1);
      } else if (len == 3) {
        return triple->contains(value[triple_perm[0]]-1, value[triple_perm[1]]-1, value[triple_perm[2]]-1);
      }
    }

    size_type children(size_type len, const array<value_type, 3> &value) const {
      if (len == 0) {
        return sigma_x;
      } else if (len == 1) {
        return pair->children(value[0]-1, pair_perm[0]);
      } else if (len == 2) {
        return triple->children(value[pair_perm[0]]-1, value[pair_perm[1]]-1, triple_perm[2]);
      } else if (len == 3) {
        return 0;
      }
    }

    value_type seek(size_type len, const array<value_type, 3> &value, size_type i) const {
      if (len == 0) {
        return i;
      } else if (len == 1) {
        return pair->seek(value[0]-1, i, pair_perm[0])+1;
      } else if (len == 2) {
        return triple->seek(value[pair_perm[0]]-1, value[pair_perm[1]]-1, i, triple_perm[2])+1;
      }
    }
  };

  template <class trie_t>
  struct cltj_hashed_index {
    typedef uint64_t size_type;
    typedef trie_t trie_type;

    std::array<pair_table, 3> pairs;
    triple_table spo;

    std::array<trie_type, 6> tries;

    explicit cltj_hashed_index() {}

    explicit cltj_hashed_index(std::vector<cltj::spo_triple> &v) {
      u32 sigma_s = 0, sigma_p = 0, sigma_o = 0;
      auto v_cpy = v;
      for (auto &[s, p, o]: v_cpy) {
        s--;
        p--;
        o--;
        sigma_s = std::max(sigma_s, s + 1);
        sigma_p = std::max(sigma_p, p + 1);
        sigma_o = std::max(sigma_o, o + 1);
      }

      pairs[0] = pair_table(v_cpy, 0, 1, sigma_s, sigma_p);
      cout << "pair 0 done" << endl;
      pairs[1] = pair_table(v_cpy, 0, 2, sigma_s, sigma_o);
      cout << "pair 1 done" << endl;
      pairs[2] = pair_table(v_cpy, 1, 2, sigma_p, sigma_o);
      cout << "pair 2 done" << endl;
      spo = triple_table(v_cpy, &pairs[0], &pairs[2], &pairs[1], sigma_p);
      cout << "triple done" << endl;
      initialize_tries();
    }

    void initialize_tries() {
      for (size_type i = 0; i < 6; ++i) {
        std::array<uint8_t, 2> pair_perm = {0, 1};
        if (i / 3) {
          pair_perm = {1, 0};
        }
        std::array triple_perm = {cltj::spo_orders[i][0], cltj::spo_orders[i][1], cltj::spo_orders[i][2]};
        tries[i] = trie_type(&pairs[i % 3], pair_perm, &spo, triple_perm);
      }
    }

    trie_type* get_trie(size_type i) const {
      return &tries[i];
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
      sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
      size_type written_bytes = 0;
      for(const auto & pair : pairs){
        written_bytes += pair.serialize(out, child, "pairs");
      }
      written_bytes += spo.serialize(out, child, "spo");
      sdsl::structure_tree::add_size(child, written_bytes);
      return written_bytes;
    }

    void load(std::istream &in) {
      for(auto & pair : pairs){
        pair.load(in);
      }
      spo.load(in);
      spo.xy = &pairs[0];
      spo.xz = &pairs[1];
      spo.yz = &pairs[2];
      initialize_tries();
    }
  };

  typedef cltj::cltj_hashed_index<cltj::hashed_trie> hashed_ltj;
}

