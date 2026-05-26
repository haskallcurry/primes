# Billion Primes

Small Rust experiments around fast prime counting and finding large nth primes.
The current target is the billionth prime:

```text
p_1,000,000,000 = 22,801,763,489
```

## Run

Find the default target:

```bash
cargo run --release
```

Find a chosen nth prime:

```bash
cargo run --release -- --target=1000000000
```

Count primes up to `x`:

```bash
cargo run --release -- --count=22801763489
```

## Prime Counters

The default counter is a Meissel-Lehmer prime-counting implementation:

```bash
cargo run --release -- --count=22801763489
```

There is also a Legendre counter matching the math used by the QueenJewels
inspiration submodule:

```bash
cargo run --release --features legendre -- --count=22801763489
```

The shared `PrimeCounter` trait lets the nth-prime search use either counter at
compile time.

## Current Approach

The nth-prime path estimates `li^-1(n)`, expands a small bracket around that
estimate, and then searches for the smallest `x` with:

```text
pi(x) >= n
```

For the billionth prime, the initial bracket is currently only about 228k
integers wide.

## References

- D. H. Lehmer, "On the Exact Number of Primes Less than a Given Limit", 1959.
- J. C. Lagarias, V. S. Miller, A. M. Odlyzko, "Computing pi(x): The
  Meissel-Lehmer Method", 1985.
- `inspiration/QueenJewels` is included as a git submodule for comparison.
