#pragma once
#include "types.h"
#include <sdsl/int_vector.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/select_support_mcl.hpp>
#include "bdz_hash_table.h"
#include "select_index.h"
#include <cltj_config.hpp>

struct pair_table {
  static const usize k = 16;
  static const usize s = 1;
  typedef bdz_hash_table<1> hash_table;
  typedef uint64_t size_type;

  usize n, sigma_x = 0, sigma_y = 0;
  hash_table table;

  u64 access_x(usize i) {
    return access(i, 0);
  }

  u64 access_y(usize i) {
    return access(i, 1);
  }

  select_index<&pair_table::access_x, k, s> index_x;
  select_index<&pair_table::access_y, k, s> index_y;

  pair_table() = default;

  pair_table(std::vector<cltj::spo_triple> &v, usize i, usize j, usize sigma_x, usize sigma_y) : sigma_x(sigma_x), sigma_y(sigma_y) {
    std::vector<u64> t;
    t.reserve(n);
    for (auto arr : v) {
      t.emplace_back(arr[i] + sigma_x*arr[j]);
    }
    std::sort(t.begin(), t.end());
    auto last = unique(t.begin(), t.end());
    t.erase(last, t.end());
    n = t.size();
    table = hash_table(t);
    t.clear();
    t.shrink_to_fit();
    index_x = decltype(index_x)(this, n, sigma_x);
    index_y = decltype(index_y)(this, n, sigma_y);
  }

  pair_table& operator=(pair_table&& other) noexcept {
    if (this != &other) {
      n = other.n;
      sigma_x = other.sigma_x;
      sigma_y = other.sigma_y;
      table = std::move(other.table);
      index_x = std::move(other.index_x);
      index_x.v = this;
      index_y = std::move(other.index_y);
      index_y.v = this;
    }
    return *this;
  }

  u64 access(usize i, usize k) {
    assert(0 <= k && k < 2);
    if (k == 0) {
      return table.access(i) % sigma_x;
    }
    return table.access(i) / sigma_x;
  }

  bool contains(u64 x, u64 y) {
    return table.contains(x + y*sigma_x);
  }

  usize index_of(u64 x, u64 y) {
    return table.index_of_key(x + y*sigma_x);
  }

  usize children(u64 c, usize k) {
    assert(0 <= k && k < 2);
    if (k == 0) {
      return index_x.count_c(c);
    }
    return index_y.count_c(c);
  }

  usize select(u64 c, usize i, usize k) {
    assert(0 <= k && k < 2);
    if (k == 0) {
      return index_x.select(c, i);
    }
    return index_y.select(c, i);
  }

  usize seek(u64 c, usize i, usize k) {
    assert(0 <= k && k < 2);
    if (k == 0) {
      return access_y(select(c, i, k));
    }
    return access_x(select(c, i, k));
  }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += table.serialize(out, child, "table");
    written_bytes += index_x.serialize(out, child, "index_y");
    written_bytes += index_y.serialize(out, child, "index_x");
    written_bytes += sdsl::write_member(n, out);
    written_bytes += sdsl::write_member(sigma_x, out);
    written_bytes += sdsl::write_member(sigma_y, out);
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  void load(std::istream &in) {
    table.load(in);
    index_x.load(in);
    index_x.v = this;
    index_y.load(in);
    index_y.v = this;
    sdsl::read_member(n, in);
    sdsl::read_member(sigma_x, in);
    sdsl::read_member(sigma_y, in);
  }
};

struct triple_table {
  static const usize k = 16;
  static const usize s = 1;
  typedef bdz_hash_table<1> hash_table;
  typedef uint64_t size_type;

  usize n, sigma_y = 0;

  pair_table *xy, *yz, *xz;
  hash_table table; // (y, xz_i)

  u64 access_yz(usize i) {
    return access(i, 0);
  }

  u64 access_xz(usize i) {
    return access(i, 1);
  }

  u64 access_xy(usize i) {
    return access(i, 2);
  }

  select_index<&triple_table::access_xy, k, s> index_xy;
  select_index<&triple_table::access_yz, k, s> index_yz;
  select_index<&triple_table::access_xz, k, s> index_xz;

  triple_table() = default;

  triple_table(std::vector<cltj::spo_triple> &v, pair_table *xy, pair_table *yz, pair_table *xz, usize sigma_y) : xy(xy), yz(yz), xz(xz), sigma_y(sigma_y) {
    std::vector<u64> t;
    t.reserve(n);
    for (auto [x, y, z] : v) {
      t.emplace_back(sigma_y * xz->index_of(x, z) + y);
    }
    std::sort(t.begin(), t.end());
    auto last = unique(t.begin(), t.end());
    t.erase(last, t.end());
    n = t.size();
    table = hash_table(t);
    t.clear();
    t.shrink_to_fit();

    index_xy = decltype(index_xy)(this, n, xy->n);
    index_yz = decltype(index_yz)(this, n, yz->n);
    index_xz = decltype(index_xz)(this, n, xz->n);
  }

  triple_table& operator=(triple_table&& other) noexcept {
    if (this != &other) {
      n = other.n;
      sigma_y = other.sigma_y;
      xy = other.xy;
      yz = other.yz;
      xz = other.xz;
      table = std::move(other.table);
      index_xy = std::move(other.index_xy);
      index_xy.v = this;
      index_xz = std::move(other.index_xz);
      index_xz.v = this;
      index_yz = std::move(other.index_yz);
      index_yz.v = this;
    }
    return *this;
  }

  u64 access(usize i, usize k) {
    assert(0 <= k && k < 3);
    auto res = table.access(i);
    usize xz_i = res / sigma_y;
    if (k == 1) {
      return xz_i;
    }
    u64 y = res % sigma_y;
    if (k == 0) {
      return yz->index_of(y, xz->access_y(xz_i));
    }
    return xy->index_of(xz->access_x(xz_i), y);
  }

  bool contains(u64 x, u64 y, u64 z) {
    return xz->contains(x, z) && table.contains(y + sigma_y * xz->index_of(x, z));
  }

  usize children(u64 a, u64 b, usize k) {
    assert(0 <= k && k < 3);
    if (k == 0) {
      return index_yz.count_c(yz->index_of(a, b));
    } else if (k == 1) {
      return index_xz.count_c(xz->index_of(a, b));
    }
    return index_xy.count_c(xy->index_of(a, b));
  }

  usize select(u64 a, u64 b, usize i, usize k) {
    assert(0 <= k && k < 3);
    if (k == 0) {
      return index_yz.select(yz->index_of(a, b), i);
    } else if (k == 1) {
      return index_xz.select(xz->index_of(a, b), i);
    } else {
      return index_xy.select(xy->index_of(a, b), i);
    }
  }

  usize seek(u64 a, u64 b, usize i, usize k) {
    assert(0 <= k && k < 3);
    auto val = table.access(select(a, b, i, k));
    if (k == 0) {
      return xz->access_x(val / sigma_y);
    } else if (k == 1) {
      return val % sigma_y;
    } else {
      return xz->access_y(val / sigma_y);
    }
  }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += table.serialize(out, child, "table");
    written_bytes += index_yz.serialize(out, child, "index_yz");
    written_bytes += index_xz.serialize(out, child, "index_xz");
    written_bytes += index_xy.serialize(out, child, "index_xy");
    written_bytes += sdsl::write_member(n, out);
    written_bytes += sdsl::write_member(sigma_y, out);
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  void load(std::istream &in) {
    table.load(in);
    index_yz.load(in);
    index_yz.v = this;
    index_xz.load(in);
    index_xz.v = this;
    index_xy.load(in);
    index_xy.v = this;
    sdsl::read_member(n, in);
    sdsl::read_member(sigma_y, in);
  }
};
