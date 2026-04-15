#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <tuple>
#include <vector>

template <typename T>
T powmod(T base, T exp, T mod) {
  T result = 1;
  T b = base % mod;

  // Binary exponentiation
  while (exp > 0) {
    if (exp & 1) {
      result = result * b % mod;
    }
    b = b * b % mod;
    exp >>= 1;
  }
  return result;
}

template <typename T>
T pow6(T x) {
  T x2 = x * x;
  T x4 = x2 * x2;
  return x4 * x2;
}

template <typename T>
T integer_sixth_root(T n) {
  if (n == 0) {
    return 0;
  }

  // Quickly find upper bound
  T high = 2;
  while (pow6(high) <= n) {
    high <<= 1;
  }

  // Binary search for the floor
  T low = 1;
  T ans = 1;
  while (low <= high) {
    T mid = low + (high - low) / 2;
    if (pow6(mid) <= n) {
      ans = mid;
      low = mid + 1;
    } else {
      high = mid - 1;
    }
  }
  return ans;
}

template <typename T, T Mod>
class ModularFilter {
 public:
  ModularFilter() {
    std::array<bool, Mod> is_sixth_power{};
    for (T a = 0; a < Mod; ++a) {
      is_sixth_power[powmod(a, T{6}, Mod)] = true;
    }

    std::vector<T> res;
    for (T a = 0; a < Mod; ++a) {
      if (is_sixth_power[a]) {
        res.push_back(a);
      }
    }

    for (T a : res) {
      for (T b : res) {
        for (T c : res) {
          reachable_[(a + b + c) % Mod] = true;
        }
      }
    }
  }

  [[nodiscard]] bool is_impossible(T n) const { return !reachable_[n % Mod]; }

 private:
  std::array<bool, Mod> reachable_{};
};

template <typename T, T Mod>
class DiophantineSolver {
 public:
  DiophantineSolver() { precompute(); }

  [[nodiscard]] size_t get_precomputed_size() const {
    size_t count = 0;
    for (const auto& v : sols_) {
      count += !v.empty();
    }
    return count;
  }

  void solve(T a_max) {
    // Check that a1^6 + a2^6 will not overflow
    auto lim = static_cast<double>(std::numeric_limits<T>::max());
    assert(a_max < std::pow(0.5 * lim, 1.0 / 6) &&
           "a_max is too large and will cause overflow");

    // Maximum value of v_div that can occur:
    // v <= a1^6 + a2^6 <= 2 * a_max^6
    T max_v_div = (pow6(a_max) * 2) / Mod;
    precompute_pair_sums(max_v_div);

#pragma omp parallel for schedule(guided)
    for (T a1 = 1; a1 <= a_max; ++a1) {
      T a16 = pow6(a1);
      for (T a2 = 1; a2 <= a1; ++a2) {
        // Prune for only non-primitive solutions
        if ((a1 % 2 == 0 && a2 % 2 == 0) || (a1 % 3 == 0 && a2 % 3 == 0)) {
          continue;
        }

        T a26 = pow6(a2);

        // For b1 such that b1^6 <= a1^6 + a2^6
        T b_max = integer_sixth_root(a16 + a26);
        for (T b1 = 1; b1 <= b_max; ++b1) {
          T b16 = pow6(b1);
          T t = a16 + a26 - b16;

          if (t % Mod == 0) {
            continue;
          }

          // Number theory tells there are exactly 6 solutions for each key
          const auto& b2_candidates = sols_[t % Mod];
          for (auto b2 : b2_candidates) {
            if (b2 > b_max) {
              break;
            }

            T b26 = pow6(b2);

            // Enforce that v is positive
            if (b26 >= t) {
              break;  // b2_candidates is sorted; all further b2 values are
                      // larger
            }

            T v = t - b26;
            T v_div = v / Mod;

            // Eliminate as much as possible using modular restrictions
            if (filter16_.is_impossible(v_div) ||
                filter9_.is_impossible(v_div) ||
                filter13_.is_impossible(v_div) ||
                filter7_.is_impossible(v_div) ||
                filter31_.is_impossible(v_div) ||
                filter19_.is_impossible(v_div)) {
              continue;
            }

            // Try to express v_div as a sum of three sixth powers
            auto cs = try_decompose(v_div);
            if (cs.has_value()) {
              // If we can, we've found a solution!
              auto [c1, c2, c3] = cs.value();
              T b3 = 7 * c1;
              T b4 = 7 * c2;
              T b5 = 7 * c3;

              // Mark this as a critical section to avoid garbled output
#pragma omp critical
              report_solution(a1, a2, b1, b2, b3, b4, b5);
            }
          }
        }
      }
    }
  }

 private:
  // This can't be unordered_map because T may not be hashable, e.g. __uint128_t
  // TODO Better to use array of arrays? The inner vector should contain exactly
  // 6 elements provided that the index into the outer vector is 1 mod 7.
  std::array<std::vector<T>, Mod> sols_;

  // Precompute solutions to x^6 = t mod 7^6 for 0 < t < 7^6 and t = 1 mod 7
  void precompute() {
    for (T x = 0; x < Mod; ++x) {
      const T t = powmod<T>(x, 6, Mod);
      if (t != 0) {
        // Ascending by construction
        sols_[t].push_back(x);
      }
    }
  }

  // TODO Can we somehow combine these for efficiency?
  ModularFilter<T, 16> filter16_;
  ModularFilter<T, 9> filter9_;
  ModularFilter<T, 13> filter13_;
  ModularFilter<T, 7> filter7_;
  ModularFilter<T, 31> filter31_;
  ModularFilter<T, 19> filter19_;

  // TODO Is it better to move decomposition logic (pair sum computations and
  // try_decomposte) into its own class like ModularFilter?
  struct PairSum {
    T sum;
    T c1;
    T c2;
  };

  std::vector<PairSum> pair_sums_;

  void precompute_pair_sums(T max_n) {
    pair_sums_.clear();

    T c_max = integer_sixth_root(max_n);
    for (T c1 = 1; c1 <= c_max; ++c1) {
      T c16 = pow6(c1);
      for (T c2 = 1; c2 <= c1; ++c2) {
        T sum = c16 + pow6(c2);
        if (sum <= max_n) {
          pair_sums_.push_back({sum, c1, c2});
        }
      }
    }

    // Sort by the sum for binary search
    std::sort(pair_sums_.begin(), pair_sums_.end(),
              [](const PairSum& a, const PairSum& b) { return a.sum < b.sum; });
  }

  // Try to express n as a sum of three sixth powers
  std::optional<std::tuple<T, T, T>> try_decompose(T n) {
    T c3_max = integer_sixth_root(n);
    for (T c3 = 1; c3 <= c3_max; ++c3) {
      T c36 = pow6(c3);
      T target = n - c36;

      // Binary search for all pair sums equal to target
      auto lower = std::lower_bound(
          pair_sums_.begin(), pair_sums_.end(), target,
          [](const PairSum& ps, const T& value) { return ps.sum < value; });

      for (auto it = lower; it != pair_sums_.end() && it->sum == target; ++it) {
        T c1 = it->c1;
        T c2 = it->c2;

        // Enforce ordering to avoid duplicate permutations
        if (c3 <= c2) {
          return std::tuple{c1, c2, c3};
        }
      }
    }
    return std::nullopt;
  }

  void report_solution(T a1, T a2, T b1, T b2, T b3, T b4, T b5) {
    std::cout << static_cast<uint64_t>(a1) << "^6 + "
              << static_cast<uint64_t>(a2)
              << "^6 = " << static_cast<uint64_t>(b1) << "^6 + "
              << static_cast<uint64_t>(b2) << "^6 + "
              << static_cast<uint64_t>(b3) << "^6 + "
              << static_cast<uint64_t>(b4) << "^6 + "
              << static_cast<uint64_t>(b5) << "^6\n";
  }
};

int main() {
  auto solver = DiophantineSolver<__uint128_t, 117649>();  // 7^6
  assert(solver.get_precomputed_size() == 16807);          // 7^5

  auto start = std::chrono::high_resolution_clock::now();
  solver.solve(1200);
  auto stop = std::chrono::high_resolution_clock::now();

  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Finished search in " << duration.count() << " ms!\n";
  return 0;
}
