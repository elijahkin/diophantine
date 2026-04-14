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

// Returns true if we're sure that n is impossible express as a sum of three
// sixth powers, and false otherwise
// TODO Could this function be automatically generated from a list of moduli?
template <typename T>
bool modularly_impossible(T n) {
  T r;

  r = n % 16;
  if (r == 4 || r == 5 || r == 6 || r == 7 || r == 8 || r == 12 || r == 13 ||
      r == 14 || r == 15) {
    return true;
  }

  r = n % 9;
  if (r == 4 || r == 5 || r == 6 || r == 7 || r == 8) {
    return true;
  }

  r = n % 13;
  if (r == 4 || r == 5 || r == 6 || r == 7 || r == 8 || r == 9) {
    return true;
  }

  r = n % 7;
  if (r == 4 || r == 5 || r == 6) {
    return true;
  }

  r = n % 31;
  if (r == 15 || r == 23 || r == 27 || r == 29 || r == 30) {
    return true;
  }

  r = n % 19;
  if (r == 5 || r == 16 || r == 17) {
    return true;
  }

  // Test is inconclusive
  return false;
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

  // Uses a float estimate then corrects by +/- 1 to handle rounding error
  T r = static_cast<T>(std::pow(static_cast<double>(n), 1.0 / 6.0));

  while (pow6(r + 1) <= n) {
    ++r;  // Step up if float rounded down
  }
  while (r > 0 && pow6(r) > n) {
    --r;  // Step down if float rounded up
  }
  return r;
}

template <typename T>
class DiophantineSolver {
 public:
  explicit DiophantineSolver(T mod) : mod_(mod) { precompute(); }

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

#pragma omp parallel for schedule(dynamic, 1)
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

          if (t % mod_ == 0) {
            continue;
          }

          // Number theory tells there are exactly 6 solutions for each key
          const auto& b2_candidates = sols_[t % mod_];
          for (auto b2 : b2_candidates) {
            if (b2 > b_max) {
              break;
            }

            T b26 = pow6(b2);

            // Enforce that v is positive
            if (t <= b26) {
              continue;
            }

            T v = t - b26;
            T v_div = v / mod_;

            // Eliminate as much as possible using modular restrictions
            if (modularly_impossible(v_div)) {
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
              report_solution(a1, a2, b1, b2, b3, b4, b5);
            }
          }
        }
      }
    }
  }

 private:
  T mod_;

  // This can't be unordered_map because T may not be hashable, e.g. __uint128_t
  // TODO Better to use a vector of arrays?
  std::vector<std::vector<T>> sols_;

  // Precompute solutions to x^6 = t mod 7^6 for 0 < t < 7^6 and t = 1 mod 7
  void precompute() {
    sols_.resize(mod_);
    for (T x = 0; x < mod_; ++x) {
      const T t = powmod<T>(x, 6, mod_);
      if (t != 0) {
        // Ascending by construction
        sols_[t].push_back(x);
      }
    }
  }

  // Try to express n as a sum of three sixth powers
  // TODO Faster to keep precomputed sorted array of all x^6 + y^6?
  std::optional<std::tuple<T, T, T>> try_decompose(T n) {
    T c1_max = integer_sixth_root(n);
    for (T c1 = 1; c1 <= c1_max; ++c1) {
      T c16 = pow6(c1);
      T rem1 = n - c16;

      // Constrain c2 <= c1 to avoid duplicate permutations
      T c2_max = std::min(integer_sixth_root(rem1), c1);
      for (T c2 = 1; c2 <= c2_max; ++c2) {
        T c26 = pow6(c2);
        T rem2 = rem1 - c26;

        T c3 = integer_sixth_root(rem2);

        // Check if remaining value is exactly a 6th root and c3 <= c2
        if (c3 > 0 && c3 <= c2 && pow6(c3) == rem2) {
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
  auto solver = DiophantineSolver<__uint128_t>(117649);  // 7^6
  assert(solver.get_precomputed_size() == 16807);        // 7^5

  auto start = std::chrono::high_resolution_clock::now();
  solver.solve(400);
  auto stop = std::chrono::high_resolution_clock::now();

  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Finished search in " << duration.count() << " ms!\n";
  return 0;
}
