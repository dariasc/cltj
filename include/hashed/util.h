#pragma once
#include "types.h"
#include <vector>

u64 mod_mul(u64 a, u64 b, u64 M) {
	i64 ret = a * b - M * u64(1.L / M * a * b);
	return ret + M * (ret < 0) - M * (ret >= (i64)M);
}

u64 mod_pow(u64 b, u64 e, u64 M) {
	u64 ans = 1;
	for (; e; b = mod_mul(b, b, M), e /= 2)
		if (e & 1) ans = mod_mul(ans, b, M);
	return ans;
}

bool is_prime(u64 n) {
	if (n < 2 || n % 6 % 4 != 1) return (n | 1) == 3;
	u64 A[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022},
	    s = __builtin_ctzll(n-1), d = n >> s;
	for (u64 a : A) {   // ^ count trailing zeroes
		u64 p = mod_pow(a%n, d, n), i = s;
		while (p != 1 && p != n - 1 && a % n && i--)
			p = mod_mul(p, p, n);
		if (p != n-1 && i != s) return 0;
	}
	return 1;
}

bool all_coprime(std::vector<u64> &keys, usize m) {
  for (auto x : keys) {
    if (x != 0 && x % m == 0) {
      return 0;
    }
  }
  return 1;
}

template <auto>
struct member_function_traits;

template <typename Class, typename Ret, typename... Args, Ret (Class::*Func)(Args...)>
struct member_function_traits<Func> {
    using class_type = Class;
};
