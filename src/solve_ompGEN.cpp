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
        residues_[t].push_back(x);
      }
    }
  }

  std::span<const size_t> operator[](size_t t) const {
    return residues_[t];
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

  std::array<SmallVector<size_t, 6>, Mod> residues_;
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

// Generalized decomposer: tries to express a value as a sum of N sixth powers.
template <size_t N, typename T>
class PowerDecomposer {
  static_assert(N >= 2, "Need at least 2 terms to decompose");

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

  // Try to express y as a sum of N sixth powers
  std::optional<std::array<size_t, N>> try_decompose(T y) const {
    // Deep pruning for 7, 8, and 9
    // TODO Generalize deep pruning for arbitrary sixth-power sums
    if ((y % 7 == 0 && y % 117649 != 0) || (y % 8 == 0 && y % 64 != 0) ||
        (y % 9 == 0 && y % 729 != 0)) {
      return std::nullopt;
    }

    std::array<size_t, N> result{};
    size_t start_max = integer_sixth_root(y);
    if (scan<0>(result, start_max, y)) {
      return result;
    }
    return std::nullopt;
  }

 private:
  // Replacing the nested for loops for picking the sixth-power values and adding them to the result array.
  // Depth:       which scan level we're at (0 = outermost, N-2 = base case)
  // result:      the array being filled in
  // upper_bound: ordering constraint from previous level (current x <= this)
  // remainder:      remaining value after subtracting all previously chosen powers
  template <size_t Depth>
  __attribute__((always_inline)) inline bool scan(std::array<size_t, N>& result, size_t upper_bound, T remainder) const {

    // Base case
    if constexpr (Depth == N - 2) {
      // Now we only have 2 terms to decompose (x1^6 + x2^6 = remainder)
      auto it = pair_sum_map_.find(remainder);
      if (it != pair_sum_map_.end()) {
        for (auto [x1, x2] : it->second) {
          // Enforce ordering 
          if (upper_bound <= x2) {
            result[0] = x1;
            result[1] = x2;
            return true;
          }
        }
      }
      return false;
    } else {
      // Recursive case: pick x for slot (N - 1 - Depth)
      constexpr size_t slot = N - 1 - Depth;
      size_t x_max = integer_sixth_root(remainder);
      if (x_max > upper_bound) {
        x_max = upper_bound;
      }

      for (size_t x = 1; x <= x_max; ++x) {
        result[slot] = x;
        if (scan<Depth + 1>(result, x, remainder - pow6<T>(x))) {
          return true;
        }
      }
      return false;
    }
  }

  // Faster than std::vector for a_max=3000, binary search kills performance
  std::unordered_map<T, std::vector<std::pair<size_t, size_t>>> pair_sum_map_;
};

// ============================= Diophantine Solver ==============================
// Wnat to generalize for other diophanite equations expressed as sums of 6th powers
// We have two functions to handle the general recursive looping logic (scan_lhs and scan_rhs)

// Solves equations of the form:
//   a1^6 + ... + a_{NL}^6 = b1^6 + ... + b_{NR}^6
template <typename T, size_t NL, size_t NR, size_t Mod>
class DiophantineSolver {
  static_assert(NR > NL, "RHS must have more terms than LHS w/o loss of generality");
  static constexpr size_t N_FREE = NR - NL;

  const PowerResidueTable<Mod>& power_residues_;
  const PowerDecomposer<N_FREE, T>& decomposer_;

 public:
  DiophantineSolver(const PowerResidueTable<Mod>& pr,
                    const PowerDecomposer<N_FREE, T>& dec)
      : power_residues_(pr), decomposer_(dec) {}

  void solve(size_t a_max) const {
#pragma omp parallel for schedule(guided)
    for (size_t a1 = 1; a1 <= a_max; ++a1) {
      // NL=1: skip a1 divisible by 7 (non-primitive)
      // TODO CAM Generalize primitivity check for NL > 2
      if constexpr (NL == 1) {
        if (a1 % 7 == 0) continue;
      }

      std::array<size_t, NL> a_vals{};
      a_vals[0] = a1;
      T lhs_sum = pow6<T>(a1);

      // Fill remaining LHS terms (if NL > 1), then scan RHS
      scan_lhs<1>(a_vals, a1, lhs_sum);
    }
  }

 private:
  // Recursive LHS scan to find a_i values (Note: depth 0 is handled by parallel for loop, hence we start with 1)
  template <size_t Depth>
  __attribute__((always_inline)) inline void scan_lhs(std::array<size_t, NL>& a_vals, size_t upper_bound, T lhs_sum) const {
    if constexpr (Depth == NL) {
      // NL=2: skip when both a-values share a factor of 2 or 3 (non-primitive)
      // TODO CAM Generalize primitiv check for NL > 2
      if constexpr (NL == 2) {
        if ((a_vals[0] % 2 == 0 && a_vals[1] % 2 == 0) ||
            (a_vals[0] % 3 == 0 && a_vals[1] % 3 == 0)) {
          return;
        }
      }

      // Start scanning brute-force RHS terms.
      std::array<size_t, NR> b_vals{};
      size_t b1_max = integer_sixth_root(lhs_sum);
      scan_rhs<0>(a_vals, b_vals, lhs_sum, b1_max);
    } else {
      for (size_t a = 1; a <= upper_bound; ++a) {
        a_vals[Depth] = a;
        scan_lhs<Depth + 1>(a_vals, a, lhs_sum + pow6<T>(a));
      }
    }
  }

  // Recursive RHS scan ot find b_i values: (Note: base case is NL-1 because the last b term is handled by the residue table)
  template <size_t Depth>
  __attribute__((always_inline)) inline void scan_rhs(const std::array<size_t, NL>& a_vals, std::array<size_t, NR>& b_vals, T remainder, size_t upper_bound) const {
    if constexpr (Depth == NL - 1) {
      // Base case: do residue table lookup + decomposition
      residue_and_decompose(a_vals, b_vals, remainder);
    } else {
      for (size_t b = 1; b <= upper_bound; ++b) {
        T b_pow = pow6<T>(b);
        if (b_pow >= remainder) break;

        b_vals[Depth] = b;
        scan_rhs<Depth + 1>(a_vals, b_vals, remainder - b_pow, b);
      }
    }
  }

  // Uses residue table for b_vals[NL-1], then calls decomposer for the rest.
  __attribute__((always_inline)) inline void residue_and_decompose(const std::array<size_t, NL>& a_vals, std::array<size_t, NR>& b_vals, T remainder) const {
    constexpr size_t res_slot = NL - 1;

    const auto& candidates = power_residues_[remainder % Mod];
    for (auto b_res : candidates) {
      // Enforce ordering to avoid duplicates
      if constexpr (NL == 1) {
        if (b_res >= a_vals[0]) break;
      } else {
        if (b_res > b_vals[res_slot - 1]) break;
      }

      T b_res_pow = pow6<T>(b_res);
      if (b_res_pow >= remainder) break;

      T v_div = (remainder - b_res_pow) / Mod;

      if (ModularFilter<ImpossibleSumPowers<N_FREE, 6>, 13, 19, 27, 31, 32,
                        49>::includes(v_div)) {
        continue;
      }

      // Try to express the remaining sums of 6th pwoers
      auto cs = decomposer_.try_decompose(v_div);
      if (cs.has_value()) {
        const auto& free_vals = cs.value();

        b_vals[res_slot] = b_res;
        for (size_t i = 0; i < N_FREE; ++i) {
          b_vals[NL + i] = 7 * free_vals[i];
        }

        std::ostringstream oss;
        oss << a_vals[0] << "^6";
        for (size_t i = 1; i < NL; ++i) {
          oss << " + " << a_vals[i] << "^6";
        }
        oss << " = " << b_vals[0] << "^6";
        for (size_t i = 1; i < NR; ++i) {
          oss << " + " << b_vals[i] << "^6";
        }
        oss << "\n";

#pragma omp critical
        std::cout << oss.str();
      }
    }
  }
};

// Wrapper
template <typename T, size_t NL, size_t NR, size_t Mod>
void solve_diophantine(size_t a_max) {
  constexpr size_t N_FREE = NR - NL;

  auto lim = static_cast<double>(std::numeric_limits<T>::max());
  assert(a_max < integer_sixth_root(lim / NL) &&
         "a_max is too large and will cause overflow");

  PowerResidueTable<Mod> power_residues;
  T max_v_div = pow6<T>(a_max) * NL / Mod;
  PowerDecomposer<N_FREE, T> decomposer(max_v_div);

  DiophantineSolver<T, NL, NR, Mod> solver(power_residues, decomposer);
  solver.solve(a_max);
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

  solve_diophantine<__uint128_t, 2, 5, mod>(a_max);

  auto stop = std::chrono::high_resolution_clock::now();

  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Finished search in " << duration.count() << " ms!\n";
  return 0;
}
