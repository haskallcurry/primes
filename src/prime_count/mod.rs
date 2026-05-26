#[cfg(not(feature = "legendre"))]
mod lehmer;

#[cfg(feature = "legendre")]
mod legendre;

#[cfg(feature = "legendre")]
pub use legendre::LegendreCounter as DefaultCounter;

#[cfg(not(feature = "legendre"))]
pub use lehmer::LehmerCounter as DefaultCounter;

pub trait PrimeCounter {
    fn pi(&mut self, x: u64) -> u64;
}

struct PrimeTable {
    primes: Vec<u64>,
    pi: Vec<u32>,
}

impl PrimeTable {
    fn new(limit: usize) -> Self {
        let limit = limit.max(2);
        let is_prime = sieve(limit);
        let mut primes = Vec::new();
        let mut pi = vec![0u32; limit + 1];
        let mut count = 0u32;

        for n in 2..=limit {
            if is_prime[n] {
                count += 1;
                primes.push(n as u64);
            }
            pi[n] = count;
        }

        Self { primes, pi }
    }
}

fn sieve(limit: usize) -> Vec<bool> {
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    for p in 2..=isqrt(limit as u64) as usize {
        if is_prime[p] {
            let mut multiple = p * p;
            while multiple <= limit {
                is_prime[multiple] = false;
                multiple += p;
            }
        }
    }

    is_prime
}

struct PhiTable {
    depth: usize,
    modulus: u64,
    values: Vec<Vec<u64>>,
}

impl PhiTable {
    const PRIMES: [usize; 6] = [2, 3, 5, 7, 11, 13];

    fn new() -> Self {
        let depth = Self::PRIMES.len();
        let modulus = Self::PRIMES.iter().product::<usize>();
        let mut values = vec![vec![0; modulus + 1]; depth + 1];

        for (n, value) in values[0].iter_mut().enumerate() {
            *value = n as u64;
        }

        for a in 1..=depth {
            let p = Self::PRIMES[a - 1];
            let (previous_rows, current_rows) = values.split_at_mut(a);
            let previous = &previous_rows[a - 1];

            for (n, value) in current_rows[0].iter_mut().enumerate() {
                *value = previous[n] - previous[n / p];
            }
        }

        Self {
            depth,
            modulus: modulus as u64,
            values,
        }
    }

    fn get(&self, x: u64, a: usize) -> u64 {
        let blocks = x / self.modulus;
        let remainder = (x % self.modulus) as usize;
        blocks * self.values[a][self.modulus as usize] + self.values[a][remainder]
    }
}

fn isqrt(n: u64) -> u64 {
    iroot(n, 2)
}

fn iroot(n: u64, degree: u32) -> u64 {
    let mut root = (n as f64).powf(1.0 / f64::from(degree)).floor() as u64;

    while pow_leq(root.saturating_add(1), degree, n) {
        root += 1;
    }

    while !pow_leq(root, degree, n) {
        root -= 1;
    }

    root
}

fn pow_leq(base: u64, degree: u32, limit: u64) -> bool {
    let mut value = 1u64;

    for _ in 0..degree {
        let Some(next) = value.checked_mul(base) else {
            return false;
        };

        if next > limit {
            return false;
        }

        value = next;
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_prime_count_table() {
        let table = PrimeTable::new(20);

        assert_eq!(table.primes, vec![2, 3, 5, 7, 11, 13, 17, 19]);
        assert_eq!(table.pi[0], 0);
        assert_eq!(table.pi[1], 0);
        assert_eq!(table.pi[2], 1);
        assert_eq!(table.pi[3], 2);
        assert_eq!(table.pi[4], 2);
        assert_eq!(table.pi[20], 8);
    }

    #[test]
    fn computes_roots_exactly() {
        assert_eq!(isqrt(0), 0);
        assert_eq!(isqrt(1), 1);
        assert_eq!(isqrt(2), 1);
        assert_eq!(isqrt(3), 1);
        assert_eq!(isqrt(4), 2);
        assert_eq!(isqrt(15), 3);
        assert_eq!(isqrt(16), 4);
        assert_eq!(isqrt(17), 4);
        assert_eq!(iroot(27, 3), 3);
        assert_eq!(iroot(28, 3), 3);
    }
}
