#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// ============================= Helper Functions ==============================

template <typename T>
constexpr T powmod(T base, T exp, T mod) {
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

// TODO If we precompute sixth powers, we can delete this function and just use
// an std::upper_bound instead
template <typename T>
size_t integer_sixth_root(T n) {
  if (n == 0) {
    return 0;
  }

  // Quickly find upper bound
  size_t high = 2;
  while (pow6<T>(high) <= n) {
    high <<= 1;
  }

  // Binary search for the floor
  size_t low = 1;
  size_t ans = 1;
  while (low <= high) {
    size_t mid = low + ((high - low) / 2);
    if (pow6<T>(mid) <= n) {
      ans = mid;
      low = mid + 1;
    } else {
      high = mid - 1;
    }
  }
  return ans;
}

// ============================ Power Residue Table ============================

template <size_t Mod>
class PowerResidueTable {
 public:
  PowerResidueTable() {
    for (size_t x = 0; x < Mod; ++x) {
      const size_t t = powmod(x, 6UL, Mod);
      if (t != 0) {
        // Ascending by construction
        residues_[t].push_back(x);
      }
    }
  }

  std::span<const size_t> operator[](size_t t) const { return residues_[t]; }

 private:
  struct SolutionSet {
   public:
    void push_back(size_t v) { values_[count_++] = v; }

    [[nodiscard]] const size_t* begin() const { return values_.data(); }
    [[nodiscard]] const size_t* end() const { return values_.data() + count_; }

   private:
    std::array<size_t, 6> values_{};
    uint8_t count_ = 0;
  };

  // TODO This only stores data for indices which are 1 mod 7. Are we using too
  // much memory?
  std::array<SolutionSet, Mod> residues_;
};

// ============================= Modular Filtering =============================

// TODO Can we generalize to sums of n kth powers? Here, n=3 and k=6
template <size_t Mod>
constexpr uint64_t compute_reachable_mask_impl() {
  static_assert(Mod <= 63, "Mod must be <= 63");

  // TODO Can we eliminate containers entirely and only use bitmasks?
  std::array<bool, Mod> is_sixth_power{};
  for (size_t a = 0; a < Mod; ++a) {
    is_sixth_power[powmod(a, 6UL, Mod)] = true;
  }

  uint64_t mask = 0;
  for (size_t a = 0; a < Mod; ++a) {
    if (!is_sixth_power[a]) {
      continue;
    }
    for (size_t b = 0; b < Mod; ++b) {
      if (!is_sixth_power[b]) {
        continue;
      }
      for (size_t c = 0; c < Mod; ++c) {
        if (!is_sixth_power[c]) {
          continue;
        }
        mask |= 1UL << ((a + b + c) % Mod);
      }
    }
  }
  return mask;
}

template <typename T, size_t Mod>
class ModularFilter {
 public:
  [[nodiscard]] bool is_impossible(T n) const {
    return !((kReachable >> (n % Mod)) & 1);
  }

 private:
  static constexpr uint64_t kReachable = compute_reachable_mask_impl<Mod>();
};

// ============================= Power Decomposer ==============================

template <typename T>
class PowerDecomposer {
 public:
  explicit PowerDecomposer(T max_n) {
    // Precompute all x^6 + y^6 under a limit
    size_t c_max = integer_sixth_root(max_n);
    for (size_t x = 1; x <= c_max; ++x) {
      T x6 = pow6<T>(x);
      for (size_t y = 1; y <= x; ++y) {
        T sum = x6 + pow6<T>(y);
        if (sum <= max_n) {
          pair_sum_map_[sum].push_back({x, y});
        }
      }
    }
  }

  // Try to express y as a sum of three sixth powers x1^6 + x2^6 + x3^6
  std::optional<std::tuple<size_t, size_t, size_t>> try_decompose(T y) const {
    // Deep pruning for 7, 8, and 9
    if ((y % 7 == 0 && y % 117649 != 0) || (y % 8 == 0 && y % 64 != 0) ||
        (y % 9 == 0 && y % 729 != 0)) {
      return std::nullopt;
    }

    size_t x3_max = integer_sixth_root(y);
    for (size_t x3 = 1; x3 <= x3_max; ++x3) {
      T x36 = pow6<T>(x3);
      T target = y - x36;

      auto it = pair_sum_map_.find(target);
      if (it != pair_sum_map_.end()) {
        for (auto [x1, x2] : it->second) {
          // Enforce ordering to avoid duplicate permutations
          if (x3 <= x2) {
            return std::tuple{x1, x2, x3};
          }
        }
      }
    }
    return std::nullopt;
  }

 private:
  // TODO AI is insisting that this would be better as a vector with equal_range
  std::unordered_map<T, std::vector<std::pair<size_t, size_t>>> pair_sum_map_;
};

// ============================ Diophantine Solver =============================

// TODO Does this make sense as a class or should it just be its own free
// function?
template <typename T, size_t Mod>
class DiophantineSolver {
 public:
  void solve(size_t a_max) {
    // Check that a1^6 + a2^6 will not overflow
    auto lim = static_cast<double>(std::numeric_limits<T>::max());
    assert(a_max < std::pow(0.5 * lim, 1.0 / 6) &&
           "a_max is too large and will cause overflow");

    // Maximum value of v_div that can occur:
    // v <= a1^6 + a2^6 <= 2 * a_max^6
    T max_v_div = (pow6<T>(a_max) * 2) / Mod;
    PowerDecomposer<T> decomposer(max_v_div);

#pragma omp parallel for schedule(guided)
    for (size_t a1 = 1; a1 <= a_max; ++a1) {
      T a16 = pow6<T>(a1);
      for (size_t a2 = 1; a2 <= a1; ++a2) {
        // a1 and a2 both being divisible by either 2 or 3 imply the solution is
        // non-primitive
        if ((a1 % 2 == 0 && a2 % 2 == 0) || (a1 % 3 == 0 && a2 % 3 == 0)) {
          continue;
        }

        T a26 = pow6<T>(a2);

        // For b1 such that b1^6 <= a1^6 + a2^6
        size_t b_max = integer_sixth_root(a16 + a26);
        for (size_t b1 = 1; b1 <= b_max; ++b1) {
          T b16 = pow6<T>(b1);
          T t = a16 + a26 - b16;

          // Number theory tells there are either 0 or 6 solutions for each key
          const auto& b2_candidates = power_residues_[t % Mod];
          // TODO We are missing some candidates, e.g. when a1 = 1117, a2 =
          // 770, and b1 = 1092, we never get b2 = 861.
          for (auto b2 : b2_candidates) {
            if (b2 > b1) {
              break;  // Enforce b1 >= b2; candidates are ascending
            }

            T b26 = pow6<T>(b2);
            if (b26 >= t) {
              break;  // Enforce v > 0; candidates are ascending
            }

            T v = t - b26;
            T v_div = v / Mod;

            // Eliminate as much as possible using modular filters, manually
            // ordered from most-pruning to least-pruning
            if (filter16_.is_impossible(v_div) ||
                filter9_.is_impossible(v_div) ||
                filter13_.is_impossible(v_div) ||
                filter7_.is_impossible(v_div) ||
                filter31_.is_impossible(v_div) ||
                filter19_.is_impossible(v_div)) {
              continue;
            }

            // Try to express v_div as a sum of three sixth powers
            auto cs = decomposer.try_decompose(v_div);
            if (cs.has_value()) {
              // If we can, we've found a solution!
              auto [c1, c2, c3] = cs.value();
              size_t b3 = 7 * c1;
              size_t b4 = 7 * c2;
              size_t b5 = 7 * c3;

              std::ostringstream oss;
              oss << a1 << "^6 + " << a2 << "^6 = " << b1 << "^6 + " << b2
                  << "^6 + " << b3 << "^6 + " << b4 << "^6 + " << b5 << "^6\n";

              // Lock while writing to avoid garbled output
#pragma omp critical
              std::cout << oss.str();
            }
          }
        }
      }
    }
  }

 private:
  // Precomputed solutions to x^6 = t mod 7^6 for 0 < t < 7^6 and t = 1 mod 7
  PowerResidueTable<Mod> power_residues_;

  ModularFilter<T, 16> filter16_;
  ModularFilter<T, 9> filter9_;
  ModularFilter<T, 13> filter13_;
  ModularFilter<T, 7> filter7_;
  ModularFilter<T, 31> filter31_;
  ModularFilter<T, 19> filter19_;
};

namespace std {
template <>
struct hash<__uint128_t> {
  size_t operator()(__uint128_t x) const {
    auto hi = static_cast<uint64_t>(x >> 64);
    auto lo = static_cast<uint64_t>(x);
    return std::hash<uint64_t>{}(hi ^ lo);
  }
};
}  // namespace std

int main(int /*argc*/, char* argv[]) {
  // Read in the search limit from command line arguments
  const size_t a_max = std::stoull(argv[1]);

  DiophantineSolver<__uint128_t, 117649> solver;  // 7^6

  auto start = std::chrono::high_resolution_clock::now();
  solver.solve(a_max);
  auto stop = std::chrono::high_resolution_clock::now();

  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Finished search in " << duration.count() << " ms!\n";
  return 0;
}
