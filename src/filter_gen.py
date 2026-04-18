from itertools import product

# A quick and dirty script for determining which modular filters to use

# Return the set of remainders of sums of n kth powers mod m
def sum_powers_remainders(n, k, m):
  res = set(a**k % m for a in range(m))
  return set(sum(x) % m for x in product(res, repeat=n))

# All prime powers < 64
prime_powers = [
  2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43,
  47, 49, 53, 59, 61
  ]

# Look for moduli for which there are some unreachable remainders
restrictive_moduli = {}
for m in prime_powers:
  l = len(sum_powers_remainders(2, 6, m))
  if l < m:
    restrictive_moduli[m] = l / m

# Sort the moduli according to what eliminates the most
print(sorted(restrictive_moduli, key=restrictive_moduli.get))
