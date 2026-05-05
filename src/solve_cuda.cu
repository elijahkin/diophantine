#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using u128 = unsigned __int128;

constexpr uint64_t kMod = 117649;  // 7^6
constexpr uint64_t kResidueKeys = kMod / 7;
constexpr uint32_t kResiduesPerKey = 6;
constexpr uint32_t kDefaultMaxSolutions = 65536;
constexpr int kDefaultThreadsPerBlock = 128;
constexpr int kMaxThreadsPerBlock = 1024;
constexpr uint64_t kMaxUint64 = 0xffffffffffffffffULL;

struct PairEntry {
  u128 sum;
  uint32_t x;
  uint32_t y;
};

struct Solution {
  uint64_t a1;
  uint64_t a2;
  uint64_t b1;
  uint64_t b2;
  uint64_t b3;
  uint64_t b4;
  uint64_t b5;
};

void cuda_check(cudaError_t err, const char* call, const char* file, int line) {
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string(file) + ":" + std::to_string(line) +
                             ": " + call + " failed: " +
                             cudaGetErrorString(err));
  }
}

#define CUDA_CHECK(call) cuda_check((call), #call, __FILE__, __LINE__)

__host__ __device__ constexpr uint64_t powmod(uint64_t base, uint64_t exp,
                                              uint64_t mod) {
  uint64_t result = 1;
  uint64_t b = base % mod;
  while (exp > 0) {
    if ((exp & 1) != 0) {
      result = (result * b) % mod;
    }
    b = (b * b) % mod;
    exp >>= 1;
  }
  return result;
}

__host__ __device__ constexpr u128 pow6(uint64_t x) {
  const u128 x2 = static_cast<u128>(x) * x;
  const u128 x4 = x2 * x2;
  return x4 * x2;
}

__host__ __device__ uint64_t integer_sixth_root(u128 n) {
  if (n == 0) {
    return 0;
  }

  uint64_t high = 2;
  while (pow6(high) <= n && high <= (kMaxUint64 / 2)) {
    high <<= 1;
  }

  uint64_t low = 1;
  uint64_t ans = 1;
  while (low <= high) {
    const uint64_t mid = low + ((high - low) / 2);
    if (pow6(mid) <= n) {
      ans = mid;
      low = mid + 1;
    } else {
      high = mid - 1;
    }
  }
  return ans;
}

template <size_t N, size_t K>
struct ImpossibleSumPowers {
  template <size_t Mod>
  __host__ __device__ static constexpr uint64_t add_residue_sets(uint64_t lhs,
                                                                 uint64_t rhs) {
    constexpr uint64_t kMask = (1ULL << Mod) - 1;
    uint64_t result = 0;
    for (size_t i = 0; i < Mod; ++i) {
      if (((rhs >> i) & 1ULL) != 0) {
        result |= ((lhs << i) | (lhs >> (Mod - i))) & kMask;
      }
    }
    return result;
  }

  template <size_t Mod>
  __host__ __device__ static constexpr uint64_t get_mask() {
    static_assert(Mod <= 63, "Mod must be <= 63 for uint64_t masks");

    uint64_t kth_powers = 0;
    for (size_t a = 0; a < Mod; ++a) {
      kth_powers |= 1ULL << powmod(a, K, Mod);
    }

    uint64_t reachable = 1ULL;
    for (size_t i = 0; i < N; ++i) {
      reachable = add_residue_sets<Mod>(reachable, kth_powers);
    }
    return ~reachable;
  }
};

template <typename Policy, size_t... Moduli>
struct ModularFilter {
  template <typename T>
  __host__ __device__ static bool includes(T n) {
    return (check<Moduli>(n) || ...);
  }

 private:
  template <size_t Mod, typename T>
  __host__ __device__ static bool check(T n) {
    constexpr uint64_t kMask = Policy::template get_mask<Mod>();
    return ((kMask >> static_cast<uint64_t>(n % Mod)) & 1ULL) != 0;
  }
};

__device__ size_t lower_bound_sum(const PairEntry* pairs, size_t pair_count,
                                  u128 target) {
  size_t low = 0;
  size_t high = pair_count;
  while (low < high) {
    const size_t mid = low + ((high - low) / 2);
    if (pairs[mid].sum < target) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }
  return low;
}

__device__ bool try_decompose(u128 y, const PairEntry* pairs, size_t pair_count,
                              uint32_t* c1, uint32_t* c2, uint32_t* c3) {
  if ((y % 7 == 0 && y % 117649 != 0) ||
      (y % 8 == 0 && y % 64 != 0) || (y % 9 == 0 && y % 729 != 0)) {
    return false;
  }

  const uint64_t x3_max = integer_sixth_root(y);
  for (uint64_t x3 = 1; x3 <= x3_max; ++x3) {
    const u128 target = y - pow6(x3);

    if (ModularFilter<ImpossibleSumPowers<2, 6>, 13, 19, 31, 37, 43,
                      61>::includes(target)) {
      continue;
    }

    const size_t first = lower_bound_sum(pairs, pair_count, target);
    for (size_t i = first; i < pair_count && pairs[i].sum == target; ++i) {
      if (x3 <= pairs[i].y) {
        *c1 = pairs[i].x;
        *c2 = pairs[i].y;
        *c3 = static_cast<uint32_t>(x3);
        return true;
      }
    }
  }
  return false;
}

__device__ uint64_t triangular(uint64_t n) { return (n * (n + 1)) / 2; }

__device__ void index_to_a_pair(uint64_t index, uint64_t a_max, uint64_t* a1,
                                uint64_t* a2) {
  uint64_t low = 1;
  uint64_t high = a_max;
  while (low < high) {
    const uint64_t mid = low + ((high - low) / 2);
    if (triangular(mid) > index) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }

  *a1 = low;
  *a2 = index - triangular(low - 1) + 1;
}

__global__ void solve_kernel(uint64_t a_max, const uint32_t* residues,
                             const uint8_t* residue_counts,
                             const PairEntry* pairs, size_t pair_count,
                             Solution* solutions,
                             unsigned int* solution_count,
                             unsigned int max_solutions) {
  const uint64_t start =
      static_cast<uint64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  const uint64_t stride = static_cast<uint64_t>(gridDim.x) * blockDim.x;
  const uint64_t total_pairs = (a_max * (a_max + 1)) / 2;

  for (uint64_t pair_index = start; pair_index < total_pairs;
       pair_index += stride) {
    uint64_t a1 = 0;
    uint64_t a2 = 0;
    index_to_a_pair(pair_index, a_max, &a1, &a2);

    if ((a1 % 2 == 0 && a2 % 2 == 0) ||
        (a1 % 3 == 0 && a2 % 3 == 0)) {
      continue;
    }

    const u128 a16 = pow6(a1);
    const u128 a26 = pow6(a2);
    const u128 lhs = a16 + a26;
    const uint64_t b_max = integer_sixth_root(lhs);

    for (uint64_t b1 = 1; b1 <= b_max; ++b1) {
      const u128 t = lhs - pow6(b1);
      const uint64_t t_mod = static_cast<uint64_t>(t % kMod);
      if (t_mod % 7 != 1) {
        continue;
      }

      const uint64_t key = t_mod / 7;
      const uint8_t count = residue_counts[key];
      for (uint8_t i = 0; i < count; ++i) {
        const uint64_t b2 = residues[(key * kResiduesPerKey) + i];
        if (b2 > b1) {
          break;
        }

        const u128 b26 = pow6(b2);
        if (b26 >= t) {
          break;
        }

        const u128 v_div = (t - b26) / kMod;
        if (ModularFilter<ImpossibleSumPowers<3, 6>, 13, 19, 27, 31, 32,
                          49>::includes(v_div)) {
          continue;
        }

        uint32_t c1 = 0;
        uint32_t c2 = 0;
        uint32_t c3 = 0;
        if (try_decompose(v_div, pairs, pair_count, &c1, &c2, &c3)) {
          const unsigned int idx = atomicAdd(solution_count, 1U);
          if (idx < max_solutions) {
            solutions[idx] = Solution{a1, a2, b1, b2, 7ULL * c1, 7ULL * c2,
                                      7ULL * c3};
          }
        }
      }
    }
  }
}

std::string to_string(u128 value) {
  if (value == 0) {
    return "0";
  }

  std::string result;
  while (value > 0) {
    const int digit = static_cast<int>(value % 10);
    result.push_back(static_cast<char>('0' + digit));
    value /= 10;
  }
  std::reverse(result.begin(), result.end());
  return result;
}

uint64_t parse_u64(const char* text, const char* name) {
  try {
    size_t pos = 0;
    const unsigned long long value = std::stoull(text, &pos);
    if (text[pos] != '\0') {
      throw std::invalid_argument("trailing characters");
    }
    return static_cast<uint64_t>(value);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("invalid ") + name + " '" + text +
                             "': " + e.what());
  }
}

void build_residue_table(std::vector<uint32_t>* residues,
                         std::vector<uint8_t>* counts) {
  residues->assign(kResidueKeys * kResiduesPerKey, 0);
  counts->assign(kResidueKeys, 0);

  for (uint64_t x = 0; x < kMod; ++x) {
    const uint64_t t = powmod(x, 6, kMod);
    if (t == 0) {
      continue;
    }

    const uint64_t key = t / 7;
    const uint8_t pos = (*counts)[key];
    if (pos >= kResiduesPerKey) {
      throw std::runtime_error("residue table overflow");
    }
    (*residues)[(key * kResiduesPerKey) + pos] = static_cast<uint32_t>(x);
    (*counts)[key] = static_cast<uint8_t>(pos + 1);
  }
}

std::vector<PairEntry> build_pair_entries(u128 max_n) {
  const uint64_t c_max = integer_sixth_root(max_n);
  if (c_max > std::numeric_limits<uint32_t>::max()) {
    throw std::runtime_error("pair table indices exceed uint32_t");
  }

  const u128 pair_count_128 =
      (static_cast<u128>(c_max) * static_cast<u128>(c_max + 1)) / 2;
  if (pair_count_128 > std::numeric_limits<size_t>::max()) {
    throw std::runtime_error("pair table is too large for this host");
  }

  std::vector<PairEntry> pairs;
  pairs.reserve(static_cast<size_t>(pair_count_128));
  for (uint64_t x = 1; x <= c_max; ++x) {
    const u128 x6 = pow6(x);
    for (uint64_t y = 1; y <= x; ++y) {
      const u128 sum = x6 + pow6(y);
      if (sum <= max_n) {
        pairs.push_back(PairEntry{sum, static_cast<uint32_t>(x),
                                  static_cast<uint32_t>(y)});
      }
    }
  }

  std::sort(pairs.begin(), pairs.end(),
            [](const PairEntry& lhs, const PairEntry& rhs) {
              if (lhs.sum != rhs.sum) {
                return lhs.sum < rhs.sum;
              }
              if (lhs.x != rhs.x) {
                return lhs.x < rhs.x;
              }
              return lhs.y < rhs.y;
            });
  return pairs;
}

void print_solution(const Solution& s) {
  std::cout << s.a1 << "^6 + " << s.a2 << "^6 = " << s.b1 << "^6 + " << s.b2
            << "^6 + " << s.b3 << "^6 + " << s.b4 << "^6 + " << s.b5
            << "^6\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    if (argc < 2 || argc > 4) {
      std::cerr << "Usage: " << argv[0]
                << " <a_max> [max_solutions] [threads_per_block]\n";
      return 2;
    }

    const uint64_t a_max = parse_u64(argv[1], "a_max");
    const uint64_t max_solutions_arg =
        argc >= 3 ? parse_u64(argv[2], "max_solutions") : kDefaultMaxSolutions;
    const uint64_t threads_per_block_arg =
        argc == 4 ? parse_u64(argv[3], "threads_per_block")
                  : kDefaultThreadsPerBlock;
    if (a_max == 0) {
      throw std::runtime_error("a_max must be positive");
    }
    if (a_max >= kMod) {
      throw std::runtime_error("a_max must be less than 7^6 (" +
                               std::to_string(kMod) +
                               ") for this residue-based search");
    }
    if (max_solutions_arg == 0 ||
        max_solutions_arg > std::numeric_limits<unsigned int>::max()) {
      throw std::runtime_error("max_solutions must fit in unsigned int");
    }
    const unsigned int max_solutions =
        static_cast<unsigned int>(max_solutions_arg);
    if (threads_per_block_arg == 0 ||
        threads_per_block_arg > static_cast<uint64_t>(kMaxThreadsPerBlock)) {
      throw std::runtime_error("threads_per_block must be in [1, 1024]");
    }
    const int threads_per_block = static_cast<int>(threads_per_block_arg);

    int device = 0;
    CUDA_CHECK(cudaGetDevice(&device));
    cudaDeviceProp props{};
    CUDA_CHECK(cudaGetDeviceProperties(&props, device));

    const auto start = std::chrono::high_resolution_clock::now();

    std::vector<uint32_t> h_residues;
    std::vector<uint8_t> h_residue_counts;
    build_residue_table(&h_residues, &h_residue_counts);

    const u128 max_v_div = (pow6(a_max) * 2) / kMod;
    std::vector<PairEntry> h_pairs = build_pair_entries(max_v_div);

    uint32_t* d_residues = nullptr;
    uint8_t* d_residue_counts = nullptr;
    PairEntry* d_pairs = nullptr;
    Solution* d_solutions = nullptr;
    unsigned int* d_solution_count = nullptr;

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_residues),
                          h_residues.size() * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_residue_counts),
                          h_residue_counts.size() * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_pairs),
                          h_pairs.size() * sizeof(PairEntry)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_solutions),
                          max_solutions * sizeof(Solution)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_solution_count),
                          sizeof(unsigned int)));

    CUDA_CHECK(cudaMemcpy(d_residues, h_residues.data(),
                          h_residues.size() * sizeof(uint32_t),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_residue_counts, h_residue_counts.data(),
                          h_residue_counts.size() * sizeof(uint8_t),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pairs, h_pairs.data(),
                          h_pairs.size() * sizeof(PairEntry),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_solution_count, 0, sizeof(unsigned int)));

    const uint64_t total_pairs = (a_max * (a_max + 1)) / 2;
    const uint64_t needed_blocks =
        (total_pairs + threads_per_block - 1) / threads_per_block;
    const int occupancy_blocks = std::max(1, props.multiProcessorCount * 16);
    const int blocks =
        static_cast<int>(std::min<uint64_t>(needed_blocks, occupancy_blocks));

    solve_kernel<<<blocks, threads_per_block>>>(
        a_max, d_residues, d_residue_counts, d_pairs, h_pairs.size(),
        d_solutions, d_solution_count, max_solutions);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    unsigned int solution_count = 0;
    CUDA_CHECK(cudaMemcpy(&solution_count, d_solution_count,
                          sizeof(unsigned int), cudaMemcpyDeviceToHost));

    const unsigned int stored_count = std::min(solution_count, max_solutions);
    std::vector<Solution> h_solutions(stored_count);
    if (stored_count > 0) {
      CUDA_CHECK(cudaMemcpy(h_solutions.data(), d_solutions,
                            stored_count * sizeof(Solution),
                            cudaMemcpyDeviceToHost));
    }

    CUDA_CHECK(cudaFree(d_residues));
    CUDA_CHECK(cudaFree(d_residue_counts));
    CUDA_CHECK(cudaFree(d_pairs));
    CUDA_CHECK(cudaFree(d_solutions));
    CUDA_CHECK(cudaFree(d_solution_count));

    std::sort(h_solutions.begin(), h_solutions.end(),
              [](const Solution& lhs, const Solution& rhs) {
                return std::tie(lhs.a1, lhs.a2, lhs.b1, lhs.b2, lhs.b3,
                                lhs.b4, lhs.b5) <
                       std::tie(rhs.a1, rhs.a2, rhs.b1, rhs.b2, rhs.b3,
                                rhs.b4, rhs.b5);
              });
    for (const Solution& solution : h_solutions) {
      print_solution(solution);
    }

    if (solution_count > max_solutions) {
      std::cerr << "Stored " << max_solutions << " of " << solution_count
                << " solutions; rerun with a larger max_solutions argument.\n";
    }

    const auto stop = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Finished CUDA search in " << duration.count() << " ms!\n";
    std::cerr << "CUDA device: " << props.name
              << ", blocks: " << blocks
              << ", threads/block: " << threads_per_block
              << ", total launched threads: "
              << (static_cast<uint64_t>(blocks) * threads_per_block)
              << ", pair entries: " << h_pairs.size()
              << ", max v/7^6: " << to_string(max_v_div) << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
