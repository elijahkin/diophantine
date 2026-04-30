#include <array>
#include <cassert>
#include <chrono>
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

// Maps each t which is 1 mod 7 to the six solutions of x^6 = t mod 7^6
// TODO Too many magic numbers inside here...
template <size_t Mod>
class PowerResidueTable {
 public:
  PowerResidueTable() {
    for (size_t x = 0; x < Mod; ++x) {
      const size_t t = powmod(x, 6UL, Mod);
      if (t != 0) {
        // Ascending by construction
        residues_[t / 7].push_back(x);
      }
    }
  }

  std::span<const size_t> operator[](size_t t) const {
    if (t % 7 != 1) {
      return {};
    }
    return residues_[t / 7];
  }

 private:
  // A simple data structure for emulating the behavior of std::vector when a
  // fixed max size is known at compile time
  template <typename T, size_t N>
  struct SmallVector {
   public:
    void push_back(T v) { values_[count_++] = v; }

    [[nodiscard]] const T* begin() const { return values_.data(); }
    [[nodiscard]] const T* end() const { return values_.data() + count_; }

   private:
    std::array<T, N> values_{};
    size_t count_ = 0;
  };

  std::array<SmallVector<size_t, 6>, Mod / 7> residues_;
};

// ============================= Modular Filtering =============================

// TODO Can we used __uint128_t instead for larger moduli?
template <size_t N, size_t K>
struct ImpossibleSumPowers {
  // TODO Can we compute this function or make it more clear/simpler?
  template <size_t Mod>
  static constexpr uint64_t add_residue_sets(uint64_t lhs, uint64_t rhs) {
    constexpr uint64_t kMask = (1UL << Mod) - 1;
    uint64_t result = 0;
    while (rhs) {
      const int i = __builtin_ctzll(rhs);
      result |= ((lhs << i) | (lhs >> (Mod - i))) & kMask;
      rhs &= rhs - 1;
    }
    return result;
  }

  template <size_t Mod>
  static constexpr uint64_t get_mask() {
    static_assert(Mod <= 63, "Mod must be <= 63 for uint64_t mask");

    uint64_t kth_powers = 0;
    for (size_t a = 0; a < Mod; ++a) {
      kth_powers |= 1UL << powmod(a, K, Mod);
    }

    // Start from {0}; after i iterations, `reachable` holds all residues
    // expressible as a sum of exactly i Kth powers
    uint64_t reachable = 1UL;  // {0}
    for (size_t i = 0; i < N; ++i) {
      reachable = add_residue_sets<Mod>(reachable, kth_powers);
    }

    return ~reachable;
  }
};

template <typename Policy, size_t... Moduli>
struct ModularFilter {
  // Uses a fold expression to apply all filters in the provided order,
  // execution short-circuits as soon as any filter returns true
  template <typename T>
  static constexpr bool includes(T n) {
    return (check<Moduli>(n) || ...);
  }

 private:
  template <size_t Mod, typename T>
  static constexpr bool check(T n) {
    constexpr auto kMask = Policy::template get_mask<Mod>();
    return (kMask >> (n % Mod)) & 1;
  }
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

  // Try to express y as a sum of five sixth powers
  std::optional<std::tuple<size_t, size_t, size_t, size_t, size_t>>try_decompose(T y) const {
    // Deep pruning for 7, 8, and 9
    // TODO Generalize deep pruning for arbitrary sixth-power sums
    if ((y % 7 == 0 && y % 117649 != 0) || (y % 8 == 0 && y % 64 != 0) ||
        (y % 9 == 0 && y % 729 != 0)) {
      return std::nullopt;
    }

    // Outermost loop: x5 (smallest value)
    size_t x5_max = integer_sixth_root(y);
    for (size_t x5 = 1; x5 <= x5_max; ++x5) {
      T rem5 = y - pow6<T>(x5);

      // Middle loop: x4 >= x5
      size_t x4_max = integer_sixth_root(rem5);
      for (size_t x4 = x5; x4 <= x4_max; ++x4) {
        T rem4 = rem5 - pow6<T>(x4);

        // Inner loop: x3 >= x4
        size_t x3_max = integer_sixth_root(rem4);
        for (size_t x3 = x4; x3 <= x3_max; ++x3) {
          T target = rem4 - pow6<T>(x3);

          if (ModularFilter<ImpossibleSumPowers<2, 6>, 13, 19, 31, 37, 43,
                            61>::includes(target)) {
            continue;
          }

          // Look up target = x1^6 + x2^6 in the pair map
          auto it = pair_sum_map_.find(target);
          if (it != pair_sum_map_.end()) {
            for (auto [x1, x2] : it->second) {
              // Enforce ordering: x2 >= x3
              if (x3 <= x2) {
                return std::tuple{x1, x2, x3, x4, x5};
              }
            }
          }
        }
      }
    }
    return std::nullopt;
  }

 private:
  // Faster than std::vector for a_max=3000, binary search kills performance
  std::unordered_map<T, std::vector<std::pair<size_t, size_t>>> pair_sum_map_;
};

// ============================ Diophantine Solver =============================

template <typename T, size_t Mod>
void solve_diophantine(size_t a_max) {
  // Check that a1^6 will not overflow
  auto lim = static_cast<double>(std::numeric_limits<T>::max());
  assert(a_max < integer_sixth_root(lim) &&
         "a_max is too large and will cause overflow");

  // Precompute solutions to x^6 = t mod 7^6 for 0 < t < 7^6 and t = 1 mod 7
  PowerResidueTable<Mod> power_residues;

  // Maximum value of v_div that can occur:
  // v <= a1^6 <= a_max^6
  T max_v_div = pow6<T>(a_max) / Mod;
  PowerDecomposer<T> decomposer(max_v_div);

#pragma omp parallel for schedule(guided)
  for (size_t a1 = 1; a1 <= a_max; ++a1) {
    // Skip non-primitive solutions.
    if (a1 % 7 == 0) continue;

    T a16 = pow6<T>(a1);

    // Since our congruence is now a1^6 \equiv b1^6 (mod 7^6), 
    // we query the table directly using a16.
    const auto& b1_candidates = power_residues[a16 % Mod];
    for (auto b1 : b1_candidates) {
      if (b1 >= a1) {
        break;  // Enforce a1 > b1 to avoid negative v 
      }

      T b16 = pow6<T>(b1);
      T v_div = (a16 - b16) / Mod;

      // Eliminate as much as possible using modular filters
      if (ModularFilter<ImpossibleSumPowers<5, 6>, 13, 19, 31, 37, 43,
                        61>::includes(v_div)) {
        continue;
      }

      // Try to express v_div as a sum of five sixth powers
      auto cs = decomposer.try_decompose(v_div);
      if (cs.has_value()) {
        // If we can, we've found a solution!
        auto [c1, c2, c3, c4, c5] = cs.value();
        size_t b2 = 7 * c1;
        size_t b3 = 7 * c2;
        size_t b4 = 7 * c3;
        size_t b5 = 7 * c4;
        size_t b6 = 7 * c5;

        std::ostringstream oss;
        oss << a1 << "^6 = " << b1 << "^6 + " << b2 << "^6 + " << b3
            << "^6 + " << b4 << "^6 + " << b5 << "^6 + " << b6 << "^6\n";

        // Lock while writing to avoid garbled output
#pragma omp critical
        std::cout << oss.str();
      }
    }
  }
}

// =================================== Usage ===================================

// libstdc++ (Linux) does not provide hash<__uint128_t>, so we define it there.
// libc++ (macOS) already provides it.
#ifdef __GLIBCXX__
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
#endif

int main(int /*argc*/, char* argv[]) {
  // Read in the search limit from command line arguments
  const size_t a_max = std::stoull(argv[1]);

  const size_t mod = 117649;  // 7^6

  auto start = std::chrono::high_resolution_clock::now();
  solve_diophantine<__uint128_t, mod>(a_max);
  auto stop = std::chrono::high_resolution_clock::now();

  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Finished search in " << duration.count() << " ms!\n";
  return 0;
}
